#Jaclyn Adaskaveg #9/24/25
#Analysis of Descriptive Trait Data, Seed Trait Data, and Flowering Time Data 
library(ggplot2)
library(GGally)
library(ggbiplot)
library(tidyverse)
library(corrplot)
library(dplyr)
library(ggridges)

#Descriptive Trait Data

#pheno_descrip<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/GWAS_trait_data/rMVP_lima_descriptive_photo_color_qual934.txt",  header = T)
pheno_descrip<-read.table("~/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript Submission/Revisions/GWAS_2PCsForAllMethods/GWAS_2026-11-08_descrip_2PC_360_NoHSW_NoSingTraitAcc_alsoFiltSiteMinCount/rMVP_lima_descriptive_photo_color_filtered619_NoHSW_NoSingTraitAcc.txt", header=T)
pheno_descrip$PurpleFlower<-gsub("NaN","NA",pheno_descrip$PurpleFlower)
pheno_descrip$Vine<-gsub("NaN","NA",pheno_descrip$Vine)
pheno_descrip$Seed_white<-gsub("NaN","NA",pheno_descrip$Seed_white)
#pheno_descrip<-pheno_descrip%>%dplyr::select(-"GH_HSW")

pheno_descrip$PurpleFlower<-as.numeric(pheno_descrip$PurpleFlower)
pheno_descrip$Vine<-as.numeric(pheno_descrip$Vine)
pheno_descrip$Seed_white<-as.numeric(pheno_descrip$Seed_white)
pheno_descrip$photoperiod<-as.numeric(pheno_descrip$photoperiod)

#Correlations
#all accessions
cor.ag<-pheno_descrip%>%dplyr::select(2:5)%>%cor(use = "pairwise.complete.obs") #%>%na.omit()

cor.test(pheno_descrip$PurpleFlower, pheno_descrip$Seed_white,conf.level =.95)
cor.test(pheno_descrip$PurpleFlower, pheno_descrip$Vine,conf.level =.95)
cor.test(pheno_descrip$PurpleFlower, pheno_descrip$photoperiod,conf.level =.95)
cor.test(pheno_descrip$Vine, pheno_descrip$Seed_white,conf.level =.95)
cor.test(pheno_descrip$photoperiod, pheno_descrip$Seed_white,conf.level =.95)
cor.test(pheno_descrip$Vine, pheno_descrip$photoperiod,conf.level =.95)


#%>%na.omit()%>%cor.mtest(conf.level = 0.95) 

corrplot::corrplot(cor.ag, p.mat = cor.pval$p, method = 'number', 
                   type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', number.cex = 1.2, 
                   order = 'alphabet',diag=F,
                   #col = colorRampPalette(c("yellow", "red"))(200),
                   #tl.pos = "lt",
                   tl.cex = 1.5, cl.cex=1, tl.col = "black", cl.pos = "r")


#mesoamerican accessions
pheno_descrip_meso<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/GWAS_trait_data/rMVP_lima_descriptive_photo_color_filtered619.txt",  header = T)

pheno_descrip_meso$PurpleFlower<-gsub("NaN","NA",pheno_descrip_meso$PurpleFlower)
pheno_descrip_meso$Vine<-gsub("NaN","NA",pheno_descrip_meso$Vine)
pheno_descrip_meso$Seed_white<-gsub("NaN","NA",pheno_descrip_meso$Seed_white)
pheno_descrip_meso<-pheno_descrip_meso%>%dplyr::select(-"GH_HSW")

pheno_descrip_meso$PurpleFlower<-as.numeric(pheno_descrip_meso$PurpleFlower)
pheno_descrip_meso$Vine<-as.numeric(pheno_descrip_meso$Vine)
pheno_descrip_meso$Seed_white<-as.numeric(pheno_descrip_meso$Seed_white)
pheno_descrip_meso$photoperiod<-as.numeric(pheno_descrip_meso$photoperiod)

cor.meso<- pheno_descrip_meso%>%dplyr::select(2:6)%>%cor(use = "pairwise.complete.obs") #%>%na.omit()
cor.pval.meso <- pheno_descrip_meso%>%dplyr::select(2:6)%>%na.omit()%>%cor.mtest( conf.level = 0.95) 

corrplot::corrplot(cor.meso, p.mat = cor.pval.meso$p, method = 'number', 
                   type = 'lower', 
                   insig='blank',
                   #addCoef.col ='black', number.cex = 1.2, 
                   order = 'alphabet',diag=F,
                   #col = colorRampPalette(c("yellow", "red"))(200),
                   #tl.pos = "lt",
                   tl.cex = 1.5, cl.cex=1, tl.col = "black", cl.pos = "r")

#Seed Trait data across 6 environments
nirs_6E<-read.csv("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/nirs_multiE_data/NIRS_predicted_with_meta/2.4.26_predicted_5fCV_nopreteat_final6E_meta_updated_meta_added_withdb.csv", header=T, stringsAsFactors = T)

str(nirs_6E)

#sort into Population subgroups
nirs_6E$Pop = NA
nirs_6E[nirs_6E$Kgroup == "K1", "Pop"] = "Meso"
nirs_6E[nirs_6E$Kgroup == "K2", "Pop"] = "B"
nirs_6E[nirs_6E$Kgroup == "K3", "Pop"] = "Meso"
nirs_6E[nirs_6E$Kgroup == "K4", "Pop"] = "Meso"
nirs_6E[nirs_6E$Kgroup == "K5", "Pop"] = "Meso"
nirs_6E[nirs_6E$Kgroup == "K6", "Pop"] = "C"
nirs_6E[nirs_6E$Kgroup == "K7", "Pop"] = "C"
nirs_6E[nirs_6E$Kgroup == "K8", "Pop"] = "Meso"
nirs_6E[nirs_6E$Kgroup == "admixed", "Pop"] = "Admixed"

nirs_6E$Pop<-as.factor(nirs_6E$Pop)
nirs_6E$Plot<-as.factor(nirs_6E$Plot)
nirs_6E$Rep<-as.factor(nirs_6E$Rep)
nirs_6E$Block<-as.factor(nirs_6E$Block)
nirs_6E$Bed<-as.factor(nirs_6E$Bed)
nirs_6E$Tier<-as.factor(nirs_6E$Tier)
nirs_6E$Year<-as.factor(nirs_6E$Year)


#PCA of macronutrient data by dry basis for all accessions
nirs_6E.db.na<-nirs_6E%>% filter(Pop!= "C")%>%#filter(Pop %in% c("Andean","Meso"))%>% #filter(Kgroup %in% c("K1","K3","K4","K5","K8"))%>%
  dplyr::select(c("LDP","Field","Pop","Protein.db","Fat.db","Ash.db", "Starch.db","HSW"))%>%
  na.omit()
str(nirs_6E.db.na)
#PCA for all samples, colored by pop gorup
pca_res2 = prcomp(nirs_6E.db.na[4:7],scale=T) #removed HSW
ggbiplot::ggbiplot(pca_res2, choices= c(1,2), obs.scale=1, var.scale=1,
                   groups= nirs_6E.db.na$Pop, #nirs_6E.db.na$Field,
                   #labels = ,
                   point.size = 1,
                   ellipse=TRUE, labels.size=4,
                   varname.size = 6, varname.abbrev=F)+
  #geom_point(aes(shape = nirs_6E.db.na$Pop, color=nirs_6E.db.na$Field)) +
  #scale_shape_manual(values=c(17,16,15))+
  scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  #scale_color_brewer(palette = "Set1")+
  #scale_fill_brewer(palette = "Set1")+
  theme_classic()+ 
  labs(title = "")+
  guides(color=guide_legend("")) +
  theme(text = element_text(size=20))

#PCA of macronutrient data for Mesoamerican accessions by dry basis
nirs_6E.db.na_meso<-nirs_6E%>% filter(Pop== "Meso")%>%#filter(Pop %in% c("Andean","Meso"))%>% #filter(Kgroup %in% c("K1","K3","K4","K5","K8"))%>%
  dplyr::select(c("LDP","Field","Pop","Protein.db","Fat.db","Ash.db", "Starch.db","HSW"))%>%
  na.omit()
str(nirs_6E.db.na_meso)
#PCA for all samples, colored by field

pca_res2 = prcomp(nirs_6E.db.na_meso[4:7],scale=T)
ggbiplot::ggbiplot(pca_res2, choices= c(1,2), obs.scale=1, var.scale=1,
                  groups= nirs_6E.db.na_meso$Field, #nirs_6E.db.na$Field,
                  #labels = ,
                  point.size = 1,
                  ellipse=TRUE, labels.size=4,
                  varname.size = 6, varname.abbrev=F)+
  scale_colour_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  scale_fill_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  theme_classic()+ 
  labs(title = "")+
  guides(color=guide_legend("")) +
  theme(text = element_text(size=20))


#create a vector for each trait and year combo
#summarize reps across field
nirs_6E_sum<-nirs_6E %>%
  #filter(Source=="UCD")%>%
  group_by(Field,Kgroup, LDP)%>%
  summarise_if(is.numeric, mean, na.rm=T)%>%
  ungroup()

nirs_6E_sum$Pop = NA
nirs_6E_sum[nirs_6E_sum$Kgroup == "K1", "Pop"] = "Meso"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K2", "Pop"] = "Andean"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K3", "Pop"] = "Meso"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K4", "Pop"] = "Meso"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K5", "Pop"] = "Meso"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K6", "Pop"] = "C"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K7", "Pop"] = "C"
nirs_6E_sum[nirs_6E_sum$Kgroup == "K8", "Pop"] = "Meso"
nirs_6E_sum[nirs_6E_sum$Kgroup == "admixed", "Pop"] = "Admixed"

nirs_6E_sum.db.na<-nirs_6E_sum%>% filter(Pop!= "C")%>%#filter(Pop %in% c("Andean","Meso"))%>% #filter(Kgroup %in% c("K1","K3","K4","K5","K8"))%>%
  dplyr::select(c("LDP","Field","Pop","Protein.db","Fat.db","Ash.db", "Starch.db","HSW"))%>%
  na.omit()
str(nirs_6E_sum.db.na)

#PCA for all samples, colored by pop gorup
pca_res2 = prcomp(nirs_6E_sum.db.na[4:7],scale=T) #removed HSW
ggbiplot::ggbiplot(pca_res2, choices= c(1,2), obs.scale=1, var.scale=1,
                   groups= nirs_6E_sum.db.na$Pop, #nirs_6E.db.na$Field,
                   #labels = ,
                   point.size = 1,
                   ellipse=TRUE, labels.size=4,
                   varname.size = 6, varname.abbrev=F)+
  #geom_point(aes(shape = nirs_6E.db.na$Pop, color=nirs_6E.db.na$Field)) +
  #scale_shape_manual(values=c(17,16,15))+
  scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  #scale_color_brewer(palette = "Set1")+
  #scale_fill_brewer(palette = "Set1")+
  theme_classic()+ 
  labs(title = "")+
  guides(color=guide_legend("")) +
  theme(text = element_text(size=20))

nirs_6E.sum.db.na_meso<-nirs_6E_sum%>% filter(Pop== "Meso")%>%#filter(Pop %in% c("Andean","Meso"))%>% #filter(Kgroup %in% c("K1","K3","K4","K5","K8"))%>%
  dplyr::select(c("LDP","Field","Pop","Protein.db","Fat.db","Ash.db", "Starch.db","HSW"))%>%
  na.omit()
str(nirs_6E.sum.db.na_meso)

#PCA for all samples, colored by field

pca_res2 = prcomp(nirs_6E.sum.db.na_meso[4:7],scale=T)
ggbiplot::ggbiplot(pca_res2, choices= c(1,2), obs.scale=1, var.scale=1,
                   groups= nirs_6E_sum.db.na$Field, #nirs_6E.db.na$Field,
                   #labels = ,
                   point.size = 1,
                   ellipse=TRUE, labels.size=4,
                   varname.size = 6, varname.abbrev=F)+
  scale_colour_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  scale_fill_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  theme_classic()+ 
  labs(title = "")+
  guides(color=guide_legend("")) +
  theme(text = element_text(size=20))


nirs_6E$Plot<-as.factor(nirs_6E$Plot)
nirs_6E$Rep<-as.factor(nirs_6E$Rep)
nirs_6E$Block<-as.factor(nirs_6E$Block)
nirs_6E$Bed<-as.factor(nirs_6E$Bed)
nirs_6E$Tier<-as.factor(nirs_6E$Tier)
nirs_6E$Year<-as.factor(nirs_6E$Year)


nirs_6E_sum.db.long<-nirs_6E_sum%>%#select(1,2,10,16:19)%>%
  gather(key = Nutr, value = Percent, c(HSW,Fat.db,Protein.db,Ash.db,Starch.db))

labels<-c(
  Protein.db="Protein",
  Starch.db="Starch",
  Ash.db="Ash",
  Fat.db="Fat",
  HSW="HSW")

nirs_6E_sum.db.long%>% filter(Pop %in% c("Meso","Andean","Admixed"))%>% #%>%filter(Pop!="C")
ggplot(aes(x = Percent, y = Field, fill = Pop)) + #
  facet_grid(~factor(Nutr, levels = c("Protein.db", "Starch.db", "Ash.db","Fat.db", "HSW")),
             labeller=labeller(Nutr=labels),scales = "free_x") +
  geom_density_ridges( jittered_points = TRUE,rel_min_height = 0.01,
                       point_size = .3,alpha=0.6, point_alpha = .8,scale=1
                       )+
  theme_ridges() + 
  scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16, colour = "black"),
        axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(color = "black", linewidth =1.5, fill = "white"),
        strip.text = element_text(size = 16, color= "black"),
        legend.position = "bottom")


## Correlation plot for all and for enviro/trait combos

nirs_6E_sum %>% 
  dplyr::select(c("Protein.db","Starch.db", "Fat.db","Ash.db","HSW"))%>%
  ggpairs(aes(color=nirs_6E_sum$Field, alpha =0.5), columnLabels = c("Protein","Starch","Fat","Ash","HSW"))+ theme_bw()+
  scale_fill_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  scale_color_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  theme(strip.background = element_rect(color = "black", size=1.5, fill = "white"),
    strip.text = element_text(size = 18),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

meso<-nirs_6E_sum %>% 
  filter(Pop=="Meso")
meso %>% 
  dplyr::select(c("Protein.db","Starch.db","Fat.db", "Ash.db","HSW"))%>%
  ggpairs(aes(color=meso$Field, alpha =0.5), columnLabels = c("Protein","Starch","Fat","Ash","HSW"))+ theme_bw()+
  scale_fill_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  scale_color_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  theme(strip.background = element_rect(color = "black", size=1.5, fill = "white"),
    strip.text = element_text(size = 18),
    panel.grid.minor = element_blank())

#all accessions untransformed BLUES
df<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript\ Submission/Revisions/UntransformedBLUEs/AllEnvsAllTraits_BLUEs_2026-09-06.txt", header = T)
df<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript\ Submission/Revisions/UntransformedBLUEs/AllEnvsAllTraits_BLUEs_2026-09-08_usingBIC.txt", header = T)
#meso accessions untransformed BLUES
df<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript\ Submission/Revisions/UntransformedBLUEs/AllEnvsAllTraits_BLUEs_2026-09-06_usingLRT_k13458_noDupRows.txt", header = T)
df<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript\ Submission/Revisions/UntransformedBLUEs/AllEnvsAllTraits_BLUEs_2026-09-08_usingBIC_k13458_noDupRows.txt", header = T)

df<- df %>%
  select("LDP",contains("HSW"),contains(".db")) 

df_long <- df %>%
  pivot_longer(
    cols = -LDP,
    names_to = c("Environment", "Trait"),
    names_pattern = "^(.+?)_(HSW|Protein\\.db|Fat\\.db|Starch\\.db|Ash\\.db)$",
    values_to = "Value"
  ) %>%
  mutate(Trait = sub("\\.db$", "", Trait)) %>%
  pivot_wider(
    names_from = Trait,
    values_from = Value
  ) %>%
  select(LDP, Environment, HSW, Protein, Fat, Starch, Ash)



df_long %>% 
  dplyr::select(c("Protein","Starch","Fat","Ash","HSW"))%>%
  ggpairs(aes(color=df_long$Environment, alpha =0.5))+ theme_bw()+
  scale_fill_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  scale_color_manual(values = c("#33a02c","#b15928","#e31a1c","#fd8d3c","#6a3d9a","#1f78b4"))+
  theme(strip.background = element_rect(color = "black", size=1.5, fill = "white"),
        strip.text = element_text(size = 18),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


cor.seed<- df%>%dplyr::select(2:19)%>%cor(use = "pairwise.complete.obs") #%>%na.omit()
cor.pval <- df%>%dplyr::select(2:19)%>%na.omit()%>%cor.mtest( conf.level = 0.95) 

#Flowering Time Data

ft<-read.csv("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Flowering_Davis24_25_gxe.csv", header=T, stringsAsFactors = T)
str(ft)
ft$Pop = NA
ft[ft$Kgroup == "K1", "Pop"] = "Meso"
ft[ft$Kgroup == "K2", "Pop"] = "Andean"
ft[ft$Kgroup == "K3", "Pop"] = "Meso"
ft[ft$Kgroup == "K4", "Pop"] = "Meso"
ft[ft$Kgroup == "K5", "Pop"] = "Meso"
ft[ft$Kgroup == "K6", "Pop"] = "NA"
ft[ft$Kgroup == "K7", "Pop"] = "NA"
ft[ft$Kgroup == "K8", "Pop"] = "Meso"
ft[ft$Kgroup == "admixed", "Pop"] = "Admixed"

#raw flowering data histogram by subpop
ft%>%filter(Pop!="NA")%>%
  ggplot(aes(x = dpp, y = Field, fill = Pop)) +
  geom_density_ridges( jittered_points = TRUE,rel_min_height = 0.01,
                       point_size = .3,alpha=0.6, point_alpha = .8,scale=1)+
  theme_ridges() + 
  scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16, colour = "black"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(color = "black", linewidth =1.5, fill = "white"),
        strip.text = element_text(size = 16, color= "black"),
        legend.position = "bottom")

#raw flowering data histogram for meso by growth habit
ft%>%filter(Pop=="Meso")%>%
  ggplot(aes(x = dpp, y = Field, fill = Vine.Bush.1)) +
  #facet_wrap(~Nutr, scales = "free_x") +
  geom_density_ridges( jittered_points = TRUE,rel_min_height = 0.01,
                       point_size = .3,alpha=0.6, point_alpha = .8,scale=1)+
  theme_ridges() + 
  #scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  #scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16, colour = "black"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(color = "black", linewidth =1.5, fill = "white"),
        strip.text = element_text(size = 16, color= "black"),
        legend.position = "bottom")

#Raw averages 
sum_ft<-ft%>%
  #filter(Source=="UCD")%>%
  group_by(SCRI, Field,Pop, Vine.Bush.1)%>%
  summarise_if(is.numeric, mean, na.rm=T)%>%
  ungroup()

#raw average flowering data histogram by subpop
sum_ft%>%filter(Pop!="NA")%>%
  ggplot(aes(x = dpp, y = Field, fill = Pop)) +
  geom_density_ridges( jittered_points = TRUE,rel_min_height = 0.01,
                       point_size = .3,alpha=0.6, point_alpha = .8,scale=1)+
  theme_ridges() + 
  scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16, colour = "black"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(color = "black", linewidth =1.5, fill = "white"),
        strip.text = element_text(size = 16, color= "black"),
        legend.position = "bottom")

#raw flowering averages histogram for meso by growth habit
sum_ft%>%filter(Pop=="Meso")%>%
  ggplot(aes(x = dpp, y = Field, fill = Vine.Bush.1)) +
  #facet_wrap(~Nutr, scales = "free_x") +
  geom_density_ridges( jittered_points = TRUE,rel_min_height = 0.01,
                       point_size = .3,alpha=0.6, point_alpha = .8,scale=1)+
  theme_ridges() + 
  #scale_color_manual(values = c("#969696","#377eb8","#4daf4a"))+
  #scale_fill_manual(values = c("#969696","#377eb8","#4daf4a"))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16, colour = "black"),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        #axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        strip.background = element_rect(color = "black", linewidth =1.5, fill = "white"),
        strip.text = element_text(size = 16, color= "black"),
        legend.position = "bottom")

df<-read.table("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/Manuscript\ Submission/Revisions/UntransformedBLUEs/AllEnvsAllTraits_BLUEs_2026-09-20_FT_usingBIC.txt", header = T)

cor<- df%>%dplyr::select(2,4)%>%cor(use = "pairwise.complete.obs") #%>%na.omit()
cor.pval <- df%>%dplyr::select(2,4)%>%na.omit()%>%cor.mtest( conf.level = 0.95) 

