setwd("H:/shengquanhu/projects/somaticmutation/TCGA-A7-A0D9-TOPHAT2-SINGLE/candidates")
minscore<-20
pvalue<-0.01
filename<-"../TCGA-BH-A0H7-TOPHAT2-SINGLE.tsv"

isdebug<-FALSE

library("brglm") 
library("elrm")

fl<-list.files(, pattern=".wsm") 

print (paste0("Total file:",length(fl)))

if(isdebug){
  filecount<-min(50, length(fl))
}else{
  filecount<-length(fl)
}

out<-matrix(NA, filecount, 19) 

for (i in 1:filecount) { 
  file<-fl[i]
  
  #file<-"10_18969245_T.wsm"
  print(paste0(i,"/", filecount, ": ", file))
  
  #get chr and position from filename
  parts<-unlist(strsplit(file,'_'))
  chr<-parts[1]
  position<-parts[2]
  parts<-unlist(strsplit(parts[3],'[.]'))
  nucleotide<-parts[1]
  
  raw.data<-read.table(file, header=T)
  tb<-table(raw.data$Base)
  if(tb[1] < tb[2]){
    levels=names(tb)
  }else{
    levels=c(names(tb)[2], names(tb)[1])
  }
  
  fit<-brglm(Base ~ factor(SAMPLE) + Score , family=binomial, data=raw.data) 
  fit.coef=coef(summary(fit))
  pvalue.logistic.group<-fit.coef[2, 4] 
  pvalue.logistic.groupError<-fit.coef[2, 2] 
  if(2 == nrow(fit.coef)){
    pvalue.logistic.score<-NA
  }else{
    pvalue.logistic.score<-fit.coef[3, 4] 
  }
  pvalue.logistic.converged<-fit$converged
  
  if(pvalue.logistic.converged){
    pvalue.exactlog<-NA
    pvalue.final<-pvalue.logistic.group
  } else {
    #################################################  
    # exact logistic regression
    x <- xtabs(~Base + interaction(SAMPLE, Score), data = raw.data)
    n.score<-dim( table(raw.data$Score))
    lab.score<-as.numeric( names( table(raw.data$Score)))
    dd <- data.frame(SAMPLE = rep(1:0, n.score), 
                     Score = rep(lab.score, each = 2), 
                     y = x[1, ], ntrials = colSums(x) )
    updated.data<-dd[dd$ntrials!=0,]
    fit.elrm<-elrm(y/ntrials ~ SAMPLE+Score, interest = ~ SAMPLE, 
                   dataset = updated.data, iter=5000, burnIn = 3000)
    pvalue.exactlog<-fit.elrm$p.values
    pvalue.final<-pvalue.exactlog
    #################################################
  }
  
  tb2<-table( raw.data$Base, raw.data$SAMPLE)
  v<-c(chr, position, position, levels[2], levels[1], nucleotide,
       tb2[levels[1], "S1"], tb2[levels[2], "S1"], tb2[levels[1], "S2"], tb2[levels[2], "S2"],
       pvalue.logistic.group, pvalue.logistic.groupError, pvalue.logistic.score, pvalue.logistic.converged, pvalue.exactlog, pvalue.final,
       NA, NA, file)
  
  out[i, ]<-v 
} 

colnames(out)<-c("chr", "start", "end", "major_allele", "minor_allele", "ref", 
                 "S1_minor_allele", "S1_major_allele", "S2_minor_allele", "S2_major_allele", 
                 "logistic_group", "logistic_group_error", "logistic_score", "logistic_converged",
                 "exact_logistic_group", "final_group_pvalue", "final_group_pvalue_fdr", "passed", "filename") 

out[,"final_group_pvalue_fdr"]<-p.adjust(as.numeric(out[,"final_group_pvalue"]), method="BH")
out[,"passed"]<-apply(out, 1, FUN=function(x){
  if(is.na(x["final_group_pvalue_fdr"])){
    return (FALSE)
  }
  
  return (as.numeric(x["final_group_pvalue_fdr"]) <= pvalue)
})

passed<-out[as.logical(out[,"passed"]),]

passed<-passed[order(passed[,"chr"], as.numeric(passed[,"start"])),]

write.csv(passed, file=filename, col.names=T, row.names=F, quote=F)

