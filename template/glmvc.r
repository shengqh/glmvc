##predefine_begin
setwd("H:/shengquanhu/projects/GuoYan/20151015_guoyan_glmvc/tcga_rna_nt_glmvc_tp0.1_tr5_score5_np0.02_g0.1/result/TCGA-BH-A0C0-RNA-TP-NT")
inputfile<-"TCGA-BH-A0C0-RNA-TP-NT.bases"
outputfile<-"TCGA-BH-A0C0-RNA-TP-NT.unfiltered.new.tsv"
file<-"9916873_A_A_G_17_0_0_12_1.9E-08"
errorrate<-0.01
pvalue<-0.05
israwpvalue<-0
checkscore<-1
min_median_score_diff<-5
use_zero_minor_allele_strategy<-0
zero_minor_allele_strategy_glm_pvalue<-0.2
##predefine_end

library("brglm")
library("stringr")

data<-read.table(inputfile, header=T, stringsAsFactors=F)
pmax<-max(data$PositionInRead)
midx<-round(pmax/2)

brconsts_for_not_converged<-0.5

isdebug<-0

files<-unique(data$Identity)
if(isdebug){
  filecount<-min(50, length(files))
}else{
  filecount<-length(files)
}

notconveraged<-NULL

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
  istart<-length(parts)
  if(is.na(as.numeric(parts[istart]))){
    istart<-istart-1
  }
  istart<-istart-9
  chr<-paste(parts[1:istart],collapse = '_')
  locus<-parts[istart+1]
  ref<-parts[istart+2]
  major<-parts[istart+3]
  minor<-parts[istart+4]
  
  filedata<-data[data$Identity == file,]
  filedata$SAMPLE<-factor(filedata$SAMPLE,levels=c("NORMAL", "TUMOR"))
  filedata$Strand<-as.factor(filedata$Strand)
  filedata$Base<-factor(filedata$Base,levels=c(major,minor))
  filedata$Position<-as.factor(filedata$Position)
  
  majordata<-filedata[filedata$Base == major,]
  minordata<-filedata[filedata$Base == minor,]
  hasScore = ifelse(checkscore, (median(majordata$Score) - median(minordata$Score) >= min_median_score_diff), 0)
  
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
    
    res <- try(fit<-brglm(curformula, family=binomial, data=filedata))
    if(class(res) != "try-error" && fit$converged){
      break
    }
    
    res <- try(fit<-brglm(curformula, family=binomial, data=filedata, br.consts=brconsts_for_not_converged))
    if(class(res) != "try-error" && fit$converged){
      break
    }
  }  
  
  if(class(res) == "try-error"){
    failed = TRUE
  }else{
    failed = FALSE
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
  }else{
    if(is.null(notconveraged)){
      notconveraged <-filedata
    }else{
      notconveraged<-rbind(notconveraged, filedata)
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


if(!is.null(notconveraged)){
  write.table(file=paste0(outputfile, ".not_coveraged"), notconveraged, row.names=FALSE, sep="\t", quote=FALSE)
}
fout<-data.frame(out, stringsAsFactors=F)

fout$brglm_group<-as.numeric(fout$brglm_group)
passed<-fout$filter == "PASS"

passedout<-fout[passed,]
unpassedout<-fout[!passed,]

passedout$brglm_group_fdr<-p.adjust(passedout$brglm_group, method="fdr")

if(israwpvalue){
  if(use_zero_minor_allele_strategy){
    failed<-(passedout$normal_minor_count == 0 & passedout$brglm_group > zero_minor_allele_strategy_glm_pvalue) | (passedout$normal_minor_count > 0 & passedout$brglm_group > pvalue)
  }else{
    failed<-passedout$brglm_group > pvalue
  }
  passedout[failed,"filter"]<-"GLM_PVALUE"
}else{
  if(use_zero_minor_allele_strategy){
    failed<-(passedout$normal_minor_count == 0 & passedout$brglm_group_fdr > zero_minor_allele_strategy_glm_pvalue) | (passedout$normal_minor_count > 0 & passedout$brglm_group_fdr > pvalue)
  }else{
    failed<-passedout$brglm_group_fdr > pvalue
  }
  passedout[failed,"filter"]<-"GLM_FDR"
}

gpassed<-passedout$filter == "PASS"
gpassedout<-passedout[gpassed,]
gunpassedout<-passedout[!gpassed,]

filtered<-rbind(unpassedout, gunpassedout, gpassedout)
rownames(filtered)<-filtered$Identity

fout$brglm_group_fdr<-filtered[fout$Identity, "brglm_group_fdr" ]
fout$filter<-filtered[fout$Identity, "filter" ]

write.table(fout, file=outputfile, col.names=T, row.names=F, quote=F, sep="\t")
