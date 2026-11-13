rm(list=ls())

library(rMVP)
library(BLINK)
library(qqman)
library(dplyr)
library(ggplot2)
library(rTASSEL)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
library(CMplot)
##Descriptive Agtrait file
setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting")
#2.2.26 updated trait data and geno data to exclude accessions with low qc/common bean 
#added trait data for color (white vs colored), removed green and combined as white because could not tell from video 
geno.obj = MVP.Data(fileVCF ="KHU056_kV_lima_bean.lisnaming.fixup_1034.accessions_no.filters_Filt360Descrip_NoHSW_NoSingTraitAcc_MAF0.05_PropHet0.8_SiteMinCount180.vcf", #"KHU056_kV_lima_bean.lisnaming.fixup.934acc_qcagtrait.vcf", "KHU056_kV_lima_bean.lisnaming.fixup.152filteredacc_descrip.vcf",#"
                    fileKin = TRUE, #computes kinship matrix
                    filePC = TRUE, #computes PC
                    filePhe = "rMVP_lima_descriptive_photo_color_filtered619_NoHSW_NoSingTraitAcc.txt",#"rMVP_lima_descriptive_photo_color_filtered152.txt",# # #"rMVP_lima_trait_Compiled_with_photo_filtered619.txt",#
                    sep.phe = "\t",
                    out="mvp_descrip_Obj360")
genotype1 <- attach.big.matrix("mvp_descrip_Obj360.geno.desc")
map1 = read.table("mvp_descrip_Obj360.geno.map",head=T)
Kinship1 <- attach.big.matrix("mvp_descrip_Obj360.kin.desc") #attach kinship object created above
Covariates <- bigmemory::as.matrix(attach.big.matrix("mvp_descrip_Obj360.pc.desc")) #attach pc object created above
Covariates_PC1 <- Covariates[,c(1:2)]
phenotype1<- read.table("mvp_descrip_Obj360.phe",head=TRUE)

## Load files/objects for rTASSEL
# Load geno file from the same VCF used for rMVP
tasGenoVCF <- readGenotypeTableFromPath(
  path = "./KHU056_kV_lima_bean.lisnaming.fixup_1034.accessions_no.filters_Filt360Descrip_NoHSW_NoSingTraitAcc_MAF0.05_PropHet0.8_SiteMinCount180.vcf" #"KHU056_kV_lima_bean.lisnaming.fixup.152filteredacc_descrip.vcf" #./KHU056_kV_lima_bean.lisnaming.fixup.619FILTEREDbyPCacc_fordescp.vcf"# #./KHU056_kV_lima_bean.lisnaming.fixup.934acc_qcagtrait.vcf" #"./KHU056_kV_lima_bean.lisnaming.fixup.1013acc_agtrait.vcf" #
)
tasGenoVCF 

dir.create(paste("GWAS_",Sys.Date(),"_DescripTraits360_AlsoFiltSiteMinCount",sep=''))
setwd(paste("GWAS_",Sys.Date(),"_DescripTraits360_AlsoFiltSiteMinCount",sep=''))

dir.create("rMVP")
dir.create("rTASSEL")
dir.create("BLINK")

#setwd("./GWAS_descriptive_rerun")

for(i in c(2:ncol(phenotype1))){
  #for(i in 2:)){ 
  phe_i1 <- phenotype1[, c(1, i)]
  trait.name = colnames(phenotype1)[i]
  #if(trait.name == "TRANSF_predicted.value_Moisture_percent.12.comps_CF23"){next}
  
  setwd("rMVP")
  
  mvp <- MVP(
    phe=phe_i1,
    geno=genotype1,
    map=map1,
    K=Kinship1,
    CV.MLM=Covariates_PC1,
    CV.FarmCPU=Covariates_PC1,
    maxLine=10000,
    method.bin="FaST-LMM", #had previously been using static as default setting
    method=c("MLM","FarmCPU"),
    permutation.threshold = TRUE,
    file.output= F #options if specifying a string: c("pmap", "pmap.signal", "plot", "log")
  )
  
  setwd("../")
  
  
  
  # Prepare the exact pheno object used for rMVP for use in rTASSEL, with PC1 appended
  phenoDF <- cbind(phe_i1,Covariates_PC1) #removed [,1] after Covariates_PC1 in order to use both PCs
  colnames(phenoDF)[1] <- "Taxon"
  colnames(phenoDF)[c((ncol(phenoDF)-1),ncol(phenoDF))] <- c("PC1","PC2") # added to ensure both PC cols are named
  head(phenoDF)
  
  # Read in that pheno object, and specify variable types (needed to ensure PCs are read and used as covariates)
  tasPhenoDF <- readPhenotypeFromDataFrame(
    phenotypeDF = phenoDF, taxaID = "Taxon", attributeTypes = c("data","covariate","covariate") #added second covariate mention here for PC2
  )
  tasPhenoDF
  
  # Join geno and pheno for use in rTASSEL
  tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = tasGenoVCF,
    phenoPathDFOrObj = tasPhenoDF
  )
  
  
  # Calculate kinship in rTASSEL (tried using rMVP-generated kinship matrix but format was not recognized)
  tasKin <- kinshipMatrix(tasObj = tasGenoPheno)
  
  # Conduct MLM via rTASSEL
  tasMLM <- assocModelFitter(
    tasObj = tasGenoPheno,
    formula = . ~ .,
    fitMarkers = TRUE,
    kinship = tasKin,
    fastAssociation = FALSE
  )
  tasMLM
  write.table(tasMLM@results$MLM_Stats,paste("rTASSEL/rTASSEL_MLM_Stats_",trait.name,sep=''),quote=F)
  
  ## Run BLINK with PC1 as a covariate
  myBlink=Blink(Y=phe_i1,GD=genotype1,GM=map1,CV = as.matrix(Covariates_PC1),maxLoop=10,time.cal=T,BIC.method="naive") #removed [,1] from after Covariates_PC1 in order to use both PCs
  write.table(myBlink$GWAS,paste("BLINK/BLINK_Stats_",trait.name,sep=''),quote=F)
  
  ## Calculate Bonferroni threshold and plot results from all four methods using this threshold (each to file, and combined plot within RStudio that also prints to file)
  Bonf_threshold = 0.05/ncol(genotype1)
  print(Bonf_threshold)
  
  
  #CM plot methods overlayed
  tasMLM.df <- tasMLM@results$MLM_Stats[, c(2,3,4,7)]
  colnames(tasMLM.df) <- c("SNP", "Chr", "Pos", "TASSEL_MLM")
  tasMLM.df$Pos<-as.numeric(tasMLM.df$Pos)
  
  rMVP_MLM.df <- cbind(
    mvp$map[, c(1:3)],
    mvp$mlm.results[, ncol(mvp$mlm.results)])
  colnames(rMVP_MLM.df) <- c("SNP", "Chr", "Pos", "MVP_MLM")
  rMVP_MLM.df$Pos<-as.numeric(rMVP_MLM.df$Pos)
  
  rMVP_FC.df <- cbind(
    mvp$map[, c(1:3)],
    mvp$farmcpu.results[, ncol(mvp$farmcpu.results)])
  colnames(rMVP_FC.df) <- c("SNP", "Chr", "Pos", "MVP_FarmCPU")
  rMVP_FC.df$Pos<-as.numeric(rMVP_FC.df$Pos)
  
  BLINK.df <- myBlink$GWAS[, c(1:3,6)]
  colnames(BLINK.df) <- c("SNP", "Chr", "Pos", "BLINK")
  BLINK.df$Pos<-as.numeric(BLINK.df$Pos)
  
  CMplot.df <- tasMLM.df %>%
    full_join(rMVP_MLM.df,  by = c("SNP","Chr","Pos")) %>%
    full_join(rMVP_FC.df,   by = c("SNP","Chr","Pos")) %>%
    full_join(BLINK.df,     by = c("SNP","Chr","Pos")) %>%
    arrange(Chr, as.numeric(Pos))
  
  write.csv(CMplot.df, paste("GWAS_hits_Bonf_",trait.name, ".csv"),row.names=F)
  
  CMplot(
    CMplot.df,
    plot.type = "m",
    multraits = TRUE,
    threshold = Bonf_threshold,
    amplify = FALSE,
    lab.cex = 1,
    dpi = 600,
    legend.ncol = 1,
    legend.pos = "left",
    file.output = T,
    file.name=paste("method_comp",trait.name,sep='')
  )
  
  CMplot(CMplot.df,plot.type="q",multraits=TRUE, threshold=Bonf_threshold,ylab.pos=2, signal.cex=1.2,signal.col="red",
         conf.int=F,box=FALSE,axis.cex=1,file="jpg",file.name=paste("qq_comp",trait.name,sep=''),dpi=600,file.output=TRUE,
         verbose=TRUE)
  
  
  png(paste("Results_all4Methods_",trait.name,".png",sep=''),width=10,height=8,units="in",res=300)
  #par(mfrow = c(4,1))
  par(mfrow = c(4,1),oma=c(1,3,1,1),mai=c(0.3,0.3,0.3,0.3))
  
  tasMLM.to.plot = tasMLM@results$MLM_Stats[,c(2,3,4,7)]
  CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
  setwd("rTASSEL/")
  CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("rTASSEL_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
  
  rMVP.to.plot.MLM = cbind(mvp$map[,c(1:3)], mvp$mlm.results[,ncol(mvp$mlm.results)])
  CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
  setwd("../rMVP/")
  CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
  
  rMVP.to.plot.farmcpu = cbind(mvp$map[,c(1:3)], mvp$farmcpu.results[,ncol(mvp$farmcpu.results)])
  CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
  CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_FarmCPU_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
  
  BLINK.to.plot = myBlink$GWAS[,c(1:3,6)]
  CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="BLINK",dpi=600)
  setwd("../BLINK/")
  CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("BLINK_",trait.name,sep=''),lab.cex=1,amplify=F,main="BLINK",dpi=600)
  dev.off()
  
  setwd("../")
}
