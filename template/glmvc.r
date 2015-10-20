##predefine_begin
setwd("H:/shengquanhu/projects/somaticmutation/TCGA_rsmc_positionInRead_RNA/TCGA-BH-A0B3-RNA-TP-NT")
inputfile<-"TCGA-BH-A0B3-RNA-TP-NT.bases"
outputfile<-"TCGA-BH-A0B3-RNA-TP-NT.tsv"
file<-"M_13710_A_A_T_2644_54_2032_387_1.5E-77"
errorrate<-0.01
pvalue<-0.05
isvalidation<-0
##predefine_end

library("brglm")
library("stringr")

data<-read.table(inputfile, header=T, stringsAsFactors=F)
pmax<-max(data$PositionInRead)
midx<-round(pmax/2)

isdebug<-0

files<-unique(data$Identity)
if(isdebug){
  filecount<-min(50, length(files))
}else{
  filecount<-length(files)
}

outcols<-c("chr", 
           "start", 
           "end", 
           "major_allele", 
           "minor_allele", 
           "ref_allele", 
           "normal_major_count", 
           "normal_minor_count", 
           "tumor_major_count", 
           "tumor_minor_count", 
           "fisher_group", 
           "fisher_normal", 
           "brglm_converged",
           "brglm_group", 
           "brglm_score", 
           "brglm_strand", 
           "brglm_position", 
           "brglm_group_fdr", 
           "filter", 
           "Identity") 
out<-matrix(NA, length(files), length(outcols))
colnames(out)<-outcols

for (i in 1:filecount) {
  file<-files[i]
  
  cat(paste0(i, "/" , filecount, ":", file, "\n"))
  
  parts<-unlist(strsplit(file,'_'))
  chr<-parts[1]
  locus<-parts[2]
  ref<-parts[3]
  major<-parts[4]
  minor<-parts[5]
  
  filedata<-data[data$Identity == file,]
  filedata$SAMPLE<-factor(filedata$SAMPLE,levels=c("NORMAL", "TUMOR"))
  filedata$Strand<-as.factor(filedata$Strand)
  filedata$Base<-factor(filedata$Base,levels=c(major,minor))
  filedata$Position<-as.factor(filedata$Position)
  
  majordata<-filedata[filedata$Base == major,]
  minordata<-filedata[filedata$Base == minor,]
  hasScore = (median(minordata$Score) - median(majordata$Score) < -2.0)
  
  tumordata<-filedata[filedata$SAMPLE=="TUMOR",]
  hasStrand <- length(unique(tumordata$Strand)) > 1
  hasPosition <- TRUE
  
  filter<-"PASS"
  
  tb<-table(filedata$Base, filedata$SAMPLE)
  pvalue.fisher.group<-fisher.test(tb)$p.value
  
  allnorm = sum(tb[,"NORMAL"])
  ntb<-data.frame(NORMAL=tb[,"NORMAL"], ERROR=round(c(allnorm * (1-errorrate), allnorm * errorrate)))
  if(ntb[minor,"NORMAL"] < ntb[minor,"ERROR"]){
    pvalue.fisher.normal<-1
  }else{
    pvalue.fisher.normal<-fisher.test(ntb)$p.value
  }
    
  pvalue.brglm.group<-""
  pvalue.brglm.strand<-""
  pvalue.brglm.position<-""
  pvalue.brglm.score<-""
  
  baseformula<-"Base ~ SAMPLE";
  res <- try(fit0<-brglm(baseformula, family=binomial, data=filedata))
  addfake<-FALSE
  failed<-FALSE
  if(class(res) == "try-error" || !fit0$converged){
    if(tb[minor, "NORMAL"] == 0){
      normaldata<-filedata[filedata$SAMPLE=="NORMAL",]
      
      fake<-minordata[1,]
      fake$SAMPLE<-"NORMAL"
      fake$Score<-round(mean(minordata$Score))
      tbstrand<-table(normaldata$Strand)
      fake$Strand<-ifelse(tbstrand["FORWARD"] < tbstrand["REVERSE"], "FORWARD", "REVERSE")
      fake$Position<-"MIDDLE"      
      fake$PositionInRead<-midx      
      filedata<-rbind(filedata, fake)
      print ("Add fake normal data.")
      addfake<-TRUE
    }
    
    if(tb[major, "TUMOR"] == 0){
      fake<-majordata[1,]
      fake$SAMPLE<-"TUMOR"
      fake$Score<-round(mean(majordata$Score))
      tbstrand<-table(tumordata$Strand)
      fake$Strand<-ifelse(tbstrand["FORWARD"] < tbstrand["REVERSE"], "FORWARD", "REVERSE")
      fake$Position<-"MIDDLE"      
      fake$PositionInRead<-midx     
      filedata<-rbind(filedata, fake)
      print ("Add fake tumor data.")
      addfake<-TRUE
    }
    
    if(addfake){
      res <- try(fit0<-brglm(baseformula, family=binomial, data=filedata))
      if(class(res) == "try-error"){
        failed<-TRUE
      }
    }else{
      failed<-TRUE
    }
  }
  
  if(!failed){
    if(fit0$converged && coef(summary(fit0))[2, 4] < pvalue){
      formulas<-c("Base ~ SAMPLE + Score + Strand + PositionInRead",
                  "Base ~ SAMPLE + Strand + PositionInRead",
                  "Base ~ SAMPLE + Score + PositionInRead",
                  "Base ~ SAMPLE + Score + Strand",
                  "Base ~ SAMPLE + PositionInRead",
                  "Base ~ SAMPLE + Strand",
                  "Base ~ SAMPLE + Score",
                  "Base ~ SAMPLE")
      
      if(!hasScore){
        formulas<-formulas[!str_detect(formulas, "Score")]
      }
      
      if(!hasStrand){
        formulas<-formulas[!str_detect(formulas, "Strand")]
      }
      
      if(!hasPosition){
        formulas<-formulas[!str_detect(formulas, "PositionInRead")]
      }
      
      for(curformula in formulas){
        cat("\t", curformula, "\n")
        
        if(curformula == baseformula){
          fit<-fit0
          break
        }
        
        res <- try(fit<-brglm(curformula, family=binomial, data=filedata))
        if(class(res) != "try-error" && fit$converged){
          break
        }
      }
    }else{
      cat("\t", baseformula, "\n")
      fit<-fit0
    }
  }
  
  pvalue.brglm.converged<-ifelse(failed, FALSE, fit$converged)
  
  if(pvalue.brglm.converged){
    fit.coef=coef(summary(fit))
    for(index in c(1:nrow(fit.coef))){
      if(rownames(fit.coef)[index] == "SAMPLETUMOR"){
        pvalue.brglm.group<-fit.coef[index, 4] 
        pvalue.brglm_logistiX_group<-pvalue.brglm.group
      }else if(rownames(fit.coef)[index] == "Score"){
        pvalue.brglm.score<-fit.coef[index, 4] 
      }else if(rownames(fit.coef)[index] == "StrandREVERSE"){
        pvalue.brglm.strand<-fit.coef[index, 4] 
      }else if(rownames(fit.coef)[index] == "PositionInRead"){
        pvalue.brglm.position<-fit.coef[index, 4] 
      }
    }
  }
  
  filter<-ifelse(failed, "GLM_FAILED", ifelse(pvalue.brglm.converged, "PASS", "NOT_CONVERGED"))
  
  v<-c(chr, 
       locus, 
       locus, 
       major, 
       minor, 
       ref,
       tb[major, "NORMAL"], 
       tb[minor, "NORMAL"],
       tb[major, "TUMOR"], 
       tb[minor, "TUMOR"],
       pvalue.fisher.group, 
       pvalue.fisher.normal,
       pvalue.brglm.converged,
       pvalue.brglm.group, 
       pvalue.brglm.score, 
       pvalue.brglm.strand, 
       pvalue.brglm.position, 
       NA,
       filter, file)
  
  out[i,]<-v
}

fout<-data.frame(out, stringsAsFactors=F)

fout$brglm_group<-as.numeric(fout$brglm_group)
passed<-fout$filter == "PASS"

passedout<-fout[passed,]
unpassedout<-fout[!passed,]

passedout$brglm_group_fdr<-p.adjust(passedout$brglm_group, method="fdr")

if(isvalidation){
  failed<-passedout$brglm_group >= pvalue
  passedout[failed,"filter"]<-"GLM_PVALUE"
}else{
  failed<-passedout$brglm_group_fdr >= pvalue
  passedout[failed,"filter"]<-"GLM_FDR"
}

gpassed<-passedout$filter == "PASS"
gpassedout<-passedout[gpassed,]
gunpassedout<-passedout[!gpassed,]

if(nrow(gpassedout) > 0){
  failed<-gpassedout$fisher_normal <= pvalue
  gpassedout[failed,"filter"]<-"NORMAL_FISHER"
}
  
filtered<-rbind(unpassedout, gunpassedout, gpassedout)
rownames(filtered)<-filtered$Identity

fout$brglm_group_fdr<-filtered[fout$Identity, "brglm_group_fdr" ]
fout$filter<-filtered[fout$Identity, "filter" ]

write.table(fout, file=outputfile, col.names=T, row.names=F, quote=F, sep="\t")
