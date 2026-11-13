rm(list=ls())

setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting") #CHD uncommented and added SpatialModelFitting
#devtools::install_github("YaoZhou89/BLINK")
#install.packages("qqman")
#install.packages("pak") 
#pak::pak("maize-genetics/rTASSEL@v0.10.0")
#install.packages("CMplot")
library(rMVP)
library(BLINK)
library(qqman)
library(dplyr)
library(ggplot2)
library(devtools)
library(rTASSEL)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
library(CMplot)


#setwd("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/")

## Read in geno file and format it for rMVP; only needs to be run once
geno.obj = MVP.Data(fileVCF ="KHU056_kV_lima_bean.lisnaming.fixup.127acc_meso_ucd_fl_h.vcf", #CHD changed 9/20 from the 248 used for seed traits
                     fileKin = TRUE, #computes kinship matrix
                     filePC = TRUE, #computes PC
                     filePhe = "AllEnvsAllTraits_transformed.BLUEs_FT_usingLRT_k13458_noDupRows_forTASSEL.txt", #CHD added; changed to LRT 9/21
                     #filePhe = "transformed.vals_ALL_nirs6E_filtered_k13458_forCHiDO_alltraits_forTASSEL.txt", #"transformed.vals_ALL_nirs6E_filtered_k13458.txt",
                     sep.phe = "\t",
                     out="mvp_nirs_trait_Obj_1PC")

# Extract various objects from that MVP.Data object
genotype1 <- attach.big.matrix("mvp_nirs_trait_Obj_1PC.geno.desc")
map1 = read.table("mvp_nirs_trait_Obj_1PC.geno.map",head=T)
Kinship1 <- attach.big.matrix("mvp_nirs_trait_Obj_1PC.kin.desc") #attach kinship object created above
Covariates <- bigmemory::as.matrix(attach.big.matrix("mvp_nirs_trait_Obj_1PC.pc.desc")) #attach pc object created above
Covariates_PC1 <- Covariates[,c(1)]
phenotype1<- read.table("mvp_nirs_trait_Obj_1PC.phe",head=TRUE)

# Load geno file from the same VCF used for rMVP
tasGenoVCF <- readGenotypeTableFromPath(
  path = "./KHU056_kV_lima_bean.lisnaming.fixup.127acc_meso_ucd_fl_h.vcf" #CHD changed 9/20 from the 248 used for seed traits
)
tasGenoVCF 

dir.create(paste("GWAS_",Sys.Date(),"usingLRT_1PCForAllMethods",sep=''))
setwd(paste("GWAS_",Sys.Date(),"usingLRT_1PCForAllMethods",sep=''))

dir.create("rMVP")
dir.create("rTASSEL")
dir.create("BLINK")

#setwd("GWAS_nirs")

# Run rMVP for multiple columns in phenotype
for(i in 2:ncol(phenotype1)){
  
  phe_i1 <- phenotype1[, c(1, i)]
  trait.name = colnames(phenotype1)[i]

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
    file.output= TRUE #options if specifying a string: c("pmap", "pmap.signal", "plot", "log")
  )
  
  setwd("../")
  
  ## Load files/objects for rTASSEL

  
  # Prepare the exact pheno object used for rMVP for use in rTASSEL, with PCs appended (changed to PCs from PC1)
  phenoDF <- cbind(phe_i1,Covariates_PC1) #removed [,1] after Covariates_PC1 in order to use both PCs
  colnames(phenoDF)[1] <- "Taxon"
  colnames(phenoDF)[ncol(phenoDF)] <- c("PC1")
  head(phenoDF)
  
  # Read in that pheno object, and specify variable types (needed to ensure PCs are read and used as covariates)
  tasPhenoDF <- readPhenotypeFromDataFrame(
    phenotypeDF = phenoDF, taxaID = "Taxon", attributeTypes = c("data","covariate")
  )
  tasPhenoDF
  
  # Join geno and pheno for use in rTASSEL
  tasGenoPheno <- readGenotypePhenotype(
    genoPathOrObj = tasGenoVCF,
    phenoPathDFOrObj = tasPhenoDF
  )
  
  # Calculate kinship in rTASSEL (tried using rMVP-generated kinship matrix but format was not recognized)
  tasKin <- kinshipMatrix(tasObj = tasGenoPheno)
  
  # CHD adding try/catch loop around rTASSEL GWAS and write step as UCD25 GDD is throwing error; 'Error in validObject(.Object) : invalid class “AssociationResultsMLM” object: invalid object for slot "traits" in class "AssociationResultsMLM": got class "NULL", should be or extend class "character"'
  tryCatch({
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
  },error=function(x){print("TASSEL GWAS failed")}
  )

  ## Run BLINK with PCs as covariates
  myBlink=tryCatch(
    {Blink(Y=phe_i1,GD=genotype1,GM=map1,CV = as.matrix(Covariates_PC1),maxLoop=10,time.cal=T,BIC.method="naive")},
    error=function(e) {
      message(paste("BLINK failed for", trait.name,
                                     "with maxLoop=10; retrying maxLoop=1"))
      return(NULL)
    })
  if (is.null(myBlink)) {
    myBlink <- tryCatch(
      {Blink(Y =phe_i1,GD=genotype1,GM=map1,CV = as.matrix(Covariates_PC1), maxLoop = 1, time.cal = TRUE, BIC.method = "naive")},
      error = function(e) {message(paste("BLINK failed again for", trait.name,"- skipping trait"))
        return(NULL)}
    )}

  
  ## Write results only if valid
  if (!is.null(myBlink) &&
      !is.null(myBlink$GWAS) &&
      nrow(myBlink$GWAS) > 0) 
    {write.table(myBlink$GWAS, paste0("BLINK/BLINK_Stats_", trait.name), quote = FALSE,row.names = FALSE)}
  else {message(paste("No BLINK results for", trait.name))}

  
  #write.table(myBlink$GWAS,paste("BLINK/BLINK_Stats_",trait.name,sep=''),quote=F)
  
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
  
  write.csv(CMplot.df, paste("GWAS_hits_Bonf_",trait.name,".csv"),row.names=F)
  
  CMplot(CMplot.df,plot.type = "m",multraits = TRUE, threshold = Bonf_threshold, amplify = FALSE, lab.cex = 1, dpi = 600, legend.ncol = 1, legend.pos = "left", file.output = T,
    file.name=paste("rMVP_FarmCPU_",trait.name,sep=''))
  
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
  dev.off() #CHD added
  
  # QQ plot
  png(filename = paste("qqplot_",trait.name,".png",sep=''), width = 480, height = 480,
      pointsize = 12, bg = "white")
  
  qq(myBlink$GWAS$P.value)
  dev.off()
  
  setwd("../")
#} #CHD moved bracket below multiPlot portion of script below
## Did not run the script beyond this point on 9/7/26

#### CHD commented this entire large section 9/7/26 
# ##Descriptive Agtrait file
# setwd("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/")
# #2.2.26 updated trait data and geno data to exclude accessions with low qc/common bean 
# #added trait data for color (white vs colored), removed green and combined as white because could not tell from video 
# geno.obj = MVP.Data(fileVCF ="KHU056_kV_lima_bean.lisnaming.fixup.619FILTEREDbyPCacc_fordescp.vcf", #"KHU056_kV_lima_bean.lisnaming.fixup.934acc_qcagtrait.vcf", "KHU056_kV_lima_bean.lisnaming.fixup.152filteredacc_descrip.vcf",#"
#                     fileKin = TRUE, #computes kinship matrix
#                     filePC = TRUE, #computes PC
#                     filePhe = "rMVP_lima_descriptive_photo_color_filtered619.txt",#"rMVP_lima_descriptive_photo_color_filtered152.txt",# # #"rMVP_lima_trait_Compiled_with_photo_filtered619.txt",#
#                     sep.phe = "\t",
#                     out="mvp_agtrait_Obj")
# genotype1 <- attach.big.matrix("mvp_agtrait_Obj.geno.desc")
# map1 = read.table("mvp_agtrait_Obj.geno.map",head=T)
# Kinship1 <- attach.big.matrix("mvp_agtrait_Obj.kin.desc") #attach kinship object created above
# Covariates <- bigmemory::as.matrix(attach.big.matrix("mvp_agtrait_Obj.pc.desc")) #attach pc object created above
# Covariates_PC1 <- Covariates[,c(1:3)]
# phenotype1<- read.table("mvp_agtrait_Obj.phe",head=TRUE)
# 
# ## Load files/objects for rTASSEL
# # Load geno file from the same VCF used for rMVP
# tasGenoVCF <- readGenotypeTableFromPath(
#   path = "KHU056_kV_lima_bean.lisnaming.fixup.152filteredacc_descrip.vcf" #./KHU056_kV_lima_bean.lisnaming.fixup.619FILTEREDbyPCacc_fordescp.vcf"# #./KHU056_kV_lima_bean.lisnaming.fixup.934acc_qcagtrait.vcf" #"./KHU056_kV_lima_bean.lisnaming.fixup.1013acc_agtrait.vcf" #
# )
# tasGenoVCF 
# 
# dir.create(paste("GWAS_",Sys.Date(),sep=''))
# setwd(paste("GWAS_",Sys.Date(),sep=''))
# 
# dir.create("rMVP")
# dir.create("rTASSEL")
# dir.create("BLINK")
# 
# #setwd("./GWAS_descriptive_rerun")
# setwd(paste("GWAS_",Sys.Date(),sep=''))
# 
# for(i in 2:ncol(phenotype1)){
# #for(i in 2:)){ 
#   phe_i1 <- phenotype1[, c(1, i)]
#   trait.name = colnames(phenotype1)[i]
#   #if(trait.name == "TRANSF_predicted.value_Moisture_percent.12.comps_CF23"){next}
#   
#   setwd("rMVP")
#   
#   mvp <- MVP(
#     phe=phe_i1,
#     geno=genotype1,
#     map=map1,
#     K=Kinship1,
#     CV.MLM=Covariates_PC1,
#     CV.FarmCPU=Covariates_PC1,
#     maxLine=10000,
#     method.bin="FaST-LMM", #had previously been using static as default setting
#     method=c("MLM","FarmCPU"),
#     permutation.threshold = TRUE,
#     file.output= F #options if specifying a string: c("pmap", "pmap.signal", "plot", "log")
#   )
#   
#   setwd("../")
#   
# 
#   
#   # Prepare the exact pheno object used for rMVP for use in rTASSEL, with PC1 appended
#   phenoDF <- cbind(phe_i1,Covariates_PC1[,1])
#   colnames(phenoDF)[1] <- "Taxon"
#   colnames(phenoDF)[ncol(phenoDF)] <- "PC1"
#   head(phenoDF)
#   
#   # Read in that pheno object, and specify variable types (needed to ensure PC1 is read and used as a covariate)
#   tasPhenoDF <- readPhenotypeFromDataFrame(
#     phenotypeDF = phenoDF, taxaID = "Taxon", attributeTypes = c("data","covariate")
#   )
#   tasPhenoDF
#   
#   # Join geno and pheno for use in rTASSEL
#   tasGenoPheno <- readGenotypePhenotype(
#     genoPathOrObj = tasGenoVCF,
#     phenoPathDFOrObj = tasPhenoDF
#   )
#   
#   # Calculate kinship in rTASSEL (tried using rMVP-generated kinship matrix but format was not recognized)
#   tasKin <- kinshipMatrix(tasObj = tasGenoPheno)
#   
#   # Conduct MLM via rTASSEL
#   tasMLM <- assocModelFitter(
#     tasObj = tasGenoPheno,
#     formula = . ~ .,
#     fitMarkers = TRUE,
#     kinship = tasKin,
#     fastAssociation = FALSE
#   )
#   tasMLM
#   write.table(tasMLM@results$MLM_Stats,paste("rTASSEL/rTASSEL_MLM_Stats_",trait.name,sep=''),quote=F)
#   
#   ## Run BLINK with PC1 as a covariate
#   myBlink=Blink(Y=phe_i1,GD=genotype1,GM=map1,CV = as.matrix(Covariates_PC1[,1]),maxLoop=10,time.cal=T,BIC.method="naive")
#   write.table(myBlink$GWAS,paste("BLINK/BLINK_Stats_",trait.name,sep=''),quote=F)
#   
#   ## Calculate Bonferroni threshold and plot results from all four methods using this threshold (each to file, and combined plot within RStudio that also prints to file)
#   Bonf_threshold = 0.05/ncol(genotype1)
#   print(Bonf_threshold)
#   
#   
#   #CM plot methods overlayed
#   tasMLM.df <- tasMLM@results$MLM_Stats[, c(2,3,4,7)]
#   colnames(tasMLM.df) <- c("SNP", "Chr", "Pos", "TASSEL_MLM")
#   tasMLM.df$Pos<-as.numeric(tasMLM.df$Pos)
#   
#   rMVP_MLM.df <- cbind(
#     mvp$map[, c(1:3)],
#     mvp$mlm.results[, ncol(mvp$mlm.results)])
#   colnames(rMVP_MLM.df) <- c("SNP", "Chr", "Pos", "MVP_MLM")
#   rMVP_MLM.df$Pos<-as.numeric(rMVP_MLM.df$Pos)
#   
#   rMVP_FC.df <- cbind(
#     mvp$map[, c(1:3)],
#     mvp$farmcpu.results[, ncol(mvp$farmcpu.results)])
#   colnames(rMVP_FC.df) <- c("SNP", "Chr", "Pos", "MVP_FarmCPU")
#   rMVP_FC.df$Pos<-as.numeric(rMVP_FC.df$Pos)
#   
#   BLINK.df <- myBlink$GWAS[, c(1:3,6)]
#   colnames(BLINK.df) <- c("SNP", "Chr", "Pos", "BLINK")
#   BLINK.df$Pos<-as.numeric(BLINK.df$Pos)
#   
#   CMplot.df <- tasMLM.df %>%
#     full_join(rMVP_MLM.df,  by = c("SNP","Chr","Pos")) %>%
#     full_join(rMVP_FC.df,   by = c("SNP","Chr","Pos")) %>%
#     full_join(BLINK.df,     by = c("SNP","Chr","Pos")) %>%
#     arrange(Chr, as.numeric(Pos))
#   
#   write.csv(CMplot.df, paste("GWAS_hits_Bonf_",trait.name, ".csv"),row.names=F)
#   
#   CMplot(
#     CMplot.df,
#     plot.type = "m",
#     multraits = TRUE,
#     threshold = Bonf_threshold,
#     amplify = FALSE,
#     lab.cex = 1,
#     dpi = 600,
#     legend.ncol = 1,
#     legend.pos = "left",
#     file.output = T,
#     file.name=paste("method_comp",trait.name,sep='')
#   )
#   
#   CMplot(CMplot.df,plot.type="q",multraits=TRUE, threshold=Bonf_threshold,ylab.pos=2, signal.cex=1.2,signal.col="red",
#          conf.int=F,box=FALSE,axis.cex=1,file="jpg",file.name=paste("qq_comp",trait.name,sep=''),dpi=600,file.output=TRUE,
#          verbose=TRUE)
#   
#   
#   png(paste("Results_all4Methods_",trait.name,".png",sep=''),width=10,height=8,units="in",res=300)
#   #par(mfrow = c(4,1))
#   par(mfrow = c(4,1),oma=c(1,3,1,1),mai=c(0.3,0.3,0.3,0.3))
#   
#   tasMLM.to.plot = tasMLM@results$MLM_Stats[,c(2,3,4,7)]
#   CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
#   setwd("rTASSEL/")
#   CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("rTASSEL_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
#   
#   rMVP.to.plot.MLM = cbind(mvp$map[,c(1:3)], mvp$mlm.results[,ncol(mvp$mlm.results)])
#   CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
#   setwd("../rMVP/")
#   CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
#   
#   rMVP.to.plot.farmcpu = cbind(mvp$map[,c(1:3)], mvp$farmcpu.results[,ncol(mvp$farmcpu.results)])
#   CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
#   CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_FarmCPU_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
#   
#   BLINK.to.plot = myBlink$GWAS[,c(1:3,6)]
#   CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="BLINK",dpi=600)
#   setwd("../BLINK/")
#   CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("BLINK_",trait.name,sep=''),lab.cex=1,amplify=F,main="BLINK",dpi=600)
#   dev.off()
#   
#   setwd("../")
# }
#                                                              
# 
# 
# ##UCD flowering and height
# setwd("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/")
# 
# geno.obj = MVP.Data(fileVCF ="KHU056_kV_lima_bean.lisnaming.fixup.127acc_meso_ucd_fl_h.vcf" ,#"KHU056_kV_lima_bean.lisnaming.fixup.205acc_ucd_fl_h.vcf",
#                     fileKin = TRUE, #computes kinship matrix
#                     filePC = TRUE, #computes PC
#                     filePhe = "transformed.vals_ft_height_127_k13458_2.txt", #"trait_UCD_flower_height_205.txt", #untransformed, vine and bush together
#                     sep.phe = "\t",
#                     out="mvp_ucdagtrait_Obj")
# 
# #"transformed.vals_ALL_UCD_flowering_height.txt", #transformed, together
# #"transformed.vals_ALL_UCD_flowering_height_sep.txt", #separated vine and bush before transforming because bimodal
# 
# genotype1 <- attach.big.matrix("mvp_ucdagtrait_Obj.geno.desc")
# map1 = read.table("mvp_ucdagtrait_Obj.geno.map",head=T)
# Kinship1 <- attach.big.matrix("mvp_ucdagtrait_Obj.kin.desc") #attach kinship object created above
# Covariates <- bigmemory::as.matrix(attach.big.matrix("mvp_ucdagtrait_Obj.pc.desc")) #attach pc object created above
# Covariates_PC1 <- Covariates[,c(1:3)]
# phenotype1<- read.table("mvp_ucdagtrait_Obj.phe",head=TRUE)
# 
# ## Load files/objects for rTASSEL
# # Load geno file from the same VCF used for rMVP
# tasGenoVCF <- readGenotypeTableFromPath(
#   path = "./KHU056_kV_lima_bean.lisnaming.fixup.127acc_meso_ucd_fl_h.vcf" #"./KHU056_kV_lima_bean.lisnaming.fixup.205acc_ucd_fl_h.vcf"
# )
# tasGenoVCF 
# 
# setwd(paste("GWAS_",Sys.Date(),sep=''))
# setwd("GWAS_2025-12-26") #for the first run with comibined growth habits (failed for GDD UCD25 because tranformation was weird)
# setwd("GWAS_2025-12-30") #separated growth habits
# setwd("GWAS_UCD_flower_height_untransformed")
# for(i in 2:ncol(phenotype1)){
#   
#   phe_i1 <- phenotype1[, c(1, i)]
#   trait.name = colnames(phenotype1)[i]
#   #if(trait.name == "TRANSF_predicted.value_Moisture_percent.12.comps_CF23"){next}
#   
#   setwd("rMVP")
#   
#   mvp <- MVP(
#     phe=phe_i1,
#     geno=genotype1,
#     map=map1,
#     K=Kinship1,
#     CV.MLM=Covariates_PC1,
#     CV.FarmCPU=Covariates_PC1,
#     maxLine=10000,
#     method.bin="FaST-LMM", #had previously been using static as default setting
#     method=c("MLM","FarmCPU"),
#     permutation.threshold = TRUE,
#     file.output= TRUE #options if specifying a string: c("pmap", "pmap.signal", "plot", "log")
#   )
#   
#   setwd("../")
#   
#   
#   
#   # Prepare the exact pheno object used for rMVP for use in rTASSEL, with PC1 appended
#   phenoDF <- cbind(phe_i1,Covariates_PC1[,1])
#   colnames(phenoDF)[1] <- "Taxon"
#   colnames(phenoDF)[ncol(phenoDF)] <- "PC1"
#   head(phenoDF)
#   
#   # Read in that pheno object, and specify variable types (needed to ensure PC1 is read and used as a covariate)
#   tasPhenoDF <- readPhenotypeFromDataFrame(
#     phenotypeDF = phenoDF, taxaID = "Taxon", attributeTypes = c("data","covariate")
#   )
#   tasPhenoDF
#   
#   # Join geno and pheno for use in rTASSEL
#   tasGenoPheno <- readGenotypePhenotype(
#     genoPathOrObj = tasGenoVCF,
#     phenoPathDFOrObj = tasPhenoDF
#   )
#   
#   # Calculate kinship in rTASSEL (tried using rMVP-generated kinship matrix but format was not recognized)
#   tasKin <- kinshipMatrix(tasObj = tasGenoPheno)
#   
#   # Conduct MLM via rTASSEL
#   tasMLM <- assocModelFitter(
#     tasObj = tasGenoPheno,
#     formula = . ~ .,
#     fitMarkers = TRUE,
#     kinship = tasKin,
#     fastAssociation = FALSE
#   )
#   tasMLM
#   write.table(tasMLM@results$MLM_Stats,paste("rTASSEL/rTASSEL_MLM_Stats_",trait.name,sep=''),quote=F)
#   
#   ## Run BLINK with PC1 as a covariate
#   myBlink=Blink(Y=phe_i1,GD=genotype1,GM=map1,CV = as.matrix(Covariates_PC1[,1]),maxLoop=10,time.cal=T,BIC.method="naive")
#   write.table(myBlink$GWAS,paste("BLINK/BLINK_Stats_",trait.name,sep=''),quote=F)
#   
#   ## Calculate Bonferroni threshold and plot results from all four methods using this threshold (each to file, and combined plot within RStudio that also prints to file)
#   Bonf_threshold = 0.05/ncol(genotype1)
#   print(Bonf_threshold)
#   
#   png(paste("Results_all4Methods_",trait.name,".png",sep=''),width=10,height=8,units="in",res=300)
#   #par(mfrow = c(4,1))
#   par(mfrow = c(4,1),oma=c(1,3,1,1),mai=c(0.3,0.3,0.3,0.3))
#   
#   tasMLM.to.plot = tasMLM@results$MLM_Stats[,c(2,3,4,7)]
#   CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
#   setwd("rTASSEL/")
#   CMplot(tasMLM.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("rTASSEL_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rTASSEL MLM",dpi=600)
#   
#   rMVP.to.plot.MLM = cbind(mvp$map[,c(1:3)], mvp$mlm.results[,ncol(mvp$mlm.results)])
#   CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
#   setwd("../rMVP/")
#   CMplot(rMVP.to.plot.MLM, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_MLM_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP MLM",dpi=600)
#   
#   rMVP.to.plot.farmcpu = cbind(mvp$map[,c(1:3)], mvp$farmcpu.results[,ncol(mvp$farmcpu.results)])
#   CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
#   CMplot(rMVP.to.plot.farmcpu, plot.type="m",threshold=Bonf_threshold,file.name=paste("rMVP_FarmCPU_",trait.name,sep=''),lab.cex=1,amplify=F,main="rMVP FarmCPU",dpi=600)
#   
#   BLINK.to.plot = myBlink$GWAS[,c(1:3,6)]
#   CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.output = F,lab.cex=1,amplify=F,main="BLINK",dpi=600)
#   setwd("../BLINK/")
#   CMplot(BLINK.to.plot, plot.type="m",threshold=Bonf_threshold,file.name=paste("BLINK_",trait.name,sep=''),lab.cex=1,amplify=F,main="BLINK",dpi=600)
#   dev.off()
#   
#   setwd("../")
# }

#CMplot Mulitraits=TRUE plot models on top of each other 

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

png(
  paste("Results_all4Methods_", trait.name, ".png", sep = ""),
  width = 10, height = 8, units = "in", res = 300
)

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
  file.output = FALSE,
  file.name=paste("rMVP_FarmCPU_",trait.name,sep='')
)

dev.off()
}