# Jonathan Berlingeri; University of California, Davis; 12/14/2023
# edited by: Jaclyn Adaskaveg 3/2024, 4/1/25

##### Script Overview #####
# This script reads compiled NIRS spectra data from spectragryph and calculates
# kennard-stone algorithm for selecting samples for calibration based on spectra
# differences between samples

# Load packages and Data --------------------------------------------------
rm(list=ls())

# install.packages("prospectr") # has kenStone function
library(prospectr)
library(tidyverse)
setwd("~/Library/LIMA/Field_Trial")

#code from calibration script to pretreat data
dat_spectral_pretransp = read.csv("nir_spectra_CVARS_all_LIMA.csv", head=T, check.names = F)
dat_spectral_w.head = t(as.data.frame(dat_spectral_pretransp)) #includes wavelength header
dim(dat_spectral_w.head)
colnames(dat_spectral_w.head) = dat_spectral_w.head[1,] #sets wavelengths as column names
dat_spectral = dat_spectral_w.head[-1,] #removes first row that had contained wavelengths #update 4.1.25: and row that contained "Log(1/R)"
View(dat_spectral[1:15,1:15])
d <- dat_spectral
names <- rownames(dat_spectral)
rownames(d) <- NULL
dat_spectral_ID <- cbind(names,d)
#clean up data and filter
dat_spectral_cvars<-as.data.frame(dat_spectral_ID)%>%  filter(!grepl("Gem",names ) )
dat_spectral_cvars$names <- gsub("sample.no:.", "", dat_spectral_cvars$names)
write.csv(dat_spectral_cvars, file="nir_spectra_CVARS23_LIMAonly_formated.csv")


Spectra<-read.csv("nir_spectra_CVARS23_LIMAonly_formated_dups_removed.csv", head=T, check.names = F)
Spectra<-column_to_rownames(Spectra, "names")

#Spec<-read.csv("~/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Central_Ferry/spectra_CF_all.csv", , header=T, stringsAsFactors = T)

# Choose Samples With Kennard-Stone Algorithm ---------------------------

sel.m.80<-kenStone(Spectra, k = 80, pc = .99, metric = 'mahal')

#Selected from kenStone Mahalanobis distance
Results_Mah80<-rownames(sel.m.80$pc)[sel.m.80$model]
write.csv(Results_Mah80, "KenStone_results_2025_04_02_CVARS23.csv")

# Visualize data ----------------------------------------------------------

#Plot PCs
par(mfrow=c(2,2))


plot(sel.m.80$pc[, 1:2], xlab = "PC1", ylab = "PC2")
points(sel.m.80$pc[sel.m.80$model, 1:2], pch = 19, col = rgb(0,0,1, alpha = .5))

plot(sel.m.80$pc[, 2:3], xlab = "PC1", ylab = "PC2")
points(sel.m.80$pc[sel.m.80$model, 2:3], pch = 19, col = rgb(0,0,1, alpha = .5))

plot(sel$pc[, 1:2], xlab = "PC1", ylab = "PC2")
plot(sel$pc[, 2:3], xlab = "PC2", ylab = "PC3")
plot(sel$pc[, c(1,3)], xlab = "PC1", ylab = "PC3")
plot(sel$pc[, c(2,2)], xlab = "PC2", ylab = "PC2")#as a check
par(mfrow=c(1,1))
