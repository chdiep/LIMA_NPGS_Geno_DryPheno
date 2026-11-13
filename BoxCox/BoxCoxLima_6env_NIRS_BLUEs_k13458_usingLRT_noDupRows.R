rm(list=ls())
# renamed this file from BoxCoxLimaDav2024_NIRS_HSW1repLargerSampleSize (1).R because it had been used for 6-env file

setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting") #CHD uncommented and added SpatialModelFitting

library(MASS)
library(lme4)

dat <- read.table("AllEnvsAllTraits_BLUEs_2026-09-06.txt",header=T)
meta_wExtraCols <- read.csv("2.4.26_predicted_5fCV_nopreteat_final6E_meta_updated_meta_added_withdb.csv",head=T)
meta_wDupRows <- meta_wExtraCols[,c("LDP","Kgroup")]
meta <- meta_wDupRows[!duplicated(meta_wDupRows),]
dat_withMeta <- merge(dat,meta,by="LDP")
dat_withMeta_13458_wExtraCols <- dat_withMeta[which(dat_withMeta$Kgroup %in% c("K1","K3","K4","K5","K8")),]
dat_withMeta_13458 = dat_withMeta_13458_wExtraCols[,c(1:61)]

best.lambdas = NULL #for recording lambdas used in boxcox transformations
transformed.vals.to.print.ALL.traits = NULL #to compile transformed values across all traits

for(col.idx in 2:ncol(dat_withMeta_13458)){
  
#Select trait and print to console
  this.trait = colnames(dat_withMeta_13458)[col.idx]
  print(this.trait)
  
  mean_this.trait = mean(dat_withMeta_13458[,col.idx],na.rm=T)
  sd_this.trait = sd(dat_withMeta_13458[,col.idx],na.rm=T)
  
  for(idx in 1:length(dat_withMeta_13458[,col.idx])){
    if(is.na(dat_withMeta_13458[idx,col.idx])){next}
    if(dat_withMeta_13458[idx,col.idx] > (mean_this.trait + 3*sd_this.trait)){
      dat_withMeta_13458[idx,col.idx] = NA 
      print('changed one high val to NA')
      print(dat_withMeta_13458$LDP[idx]) #JAA replaced X.Trait. with LDP
    }else if(dat_withMeta_13458[idx,col.idx] < (mean_this.trait - 3*sd_this.trait)){
      dat_withMeta_13458[idx,col.idx] = NA
      print('changed one low val to NA')
      print(dat_withMeta_13458$LDP[idx]) #JAA replaced X.Trait. with LDP
    }
  }
  
  ##Conduct Box-Cox transformation
#and print transformed values to file (save lambdas as well, to be printed after trait loop)
mod.for.boxcox = lm(dat_withMeta_13458[,col.idx] ~ 1)
boxcox.lambda.eval = 
  as.data.frame(boxcox(mod.for.boxcox,
                       lambda = seq(-2,2,by=0.5),plotit=T,interp=F)) #testing the convenient lambda values; CHD changed from -3/3 to -2/2, as -3 was resulting in values for which it was not clear all digits were being retained for mapping
convenient.lambda.likelihoods = 
  boxcox.lambda.eval[which(boxcox.lambda.eval$x %in% c(-2,-1,-0.5,0,0.5,1,2)),]

best.xval = 
  convenient.lambda.likelihoods[which.max(convenient.lambda.likelihoods$y),"x"]

if(best.xval == 0){
  dat_withMeta_13458.transformed = log(dat_withMeta_13458[,col.idx]) #CHD confirmed log is by default base e in R
} else {
  dat_withMeta_13458.transformed = (dat_withMeta_13458[,col.idx]^(best.xval)-1)/best.xval
}

transformed.vals.to.print <- data.frame(
  `taxa` = dat_withMeta_13458[["LDP"]],  # your genotype/ID column
  setNames(list(dat_withMeta_13458.transformed), paste0("TRANSF_predicted.value_", this.trait))
)

if(is.null(transformed.vals.to.print.ALL.traits)){
  transformed.vals.to.print.ALL.traits = transformed.vals.to.print
}else{
  transformed.vals.to.print.ALL.traits = cbind(transformed.vals.to.print.ALL.traits,transformed.vals.to.print[,2])
  colnames(transformed.vals.to.print.ALL.traits)[ncol(transformed.vals.to.print.ALL.traits)] = this.trait
}

best.lambdas = rbind(best.lambdas,c(this.trait, best.xval))

hist(as.numeric(dat_withMeta_13458[,col.idx]),main=this.trait,breaks=40)
hist(as.numeric(transformed.vals.to.print[,2]),main=this.trait,breaks=40)

}#end trait loop

write.table(transformed.vals.to.print.ALL.traits,
            paste("AllEnvsAllTraits_transformed.BLUEs_k13458_usingLRT_noDupRows.txt",sep=''),
            sep='\t',quote=F,row.names=F)

write.table(best.lambdas,"bestLambdas_usingLRT_noDupRows.txt",sep='\t',quote=F,row.names=F)