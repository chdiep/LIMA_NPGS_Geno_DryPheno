rm(list=ls())

library(asreml)
library(dplyr)
library(tidyr)

setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting")

#pheno_dat_noPlotID = read.table("transformed.vals_ft_height_127_k13458_2.txt",sep='\t',head=T) #was used for GWAS but doesn't have plot ID
#colnames(pheno_dat_noPlotID)[which(colnames(pheno_dat_noPlotID) == "TRANSF_predicted.value_dpp_UCD24")] = "dpp_UCD24"
#pheno_dat_noPlotID_long = pivot_longer(pheno_dat_noPlotID, cols=c(2:6),names_to = c("Trait","Field"),names_sep="_")
#meta_dat = read.csv("2.4.26_predicted_5fCV_nopreteat_final6E_meta_updated_meta_added_withdb.csv",head=T)
#names(meta_dat)
#meta_dat = meta_dat[,c(1:10)]

dat = read.csv("Flowering_Davis24_25_gxe.csv",head=T)
head(dat)
str(dat)

View(dat %>% group_by(Field) %>% summarise_all(list(min,max)))
colnames(dat)
dat[,c(2:15)] = lapply(dat[,c(2:15)],as.factor)

env_vect = unique(dat$Field)
trait_vect = colnames(dat)[c(16:17)]

#______________________________________________________________________________#
#Make grids merge with data----------
#______________________________________________________________________________#

# UCD24
grid.UCD24 = expand.grid(Field = 'UCD24', Bed = 57:78, Tier = 2:24)
dat.UCD24 = merge(dat, grid.UCD24, by = c('Field', 'Bed', 'Tier'), all.y = TRUE)

# UCD25
grid.UCD25 = expand.grid(Field = 'UCD25', Bed = 41:64, Tier = 2:19)
dat.UCD25 = merge(dat, grid.UCD25, by = c('Field', 'Bed', 'Tier'), all.y = TRUE)

# Now combine all datasets and reorder
dat.ALL = bind_rows(dat.UCD24,dat.UCD25) %>%
  arrange(Field, Bed, Tier)  # Final ordering

# Set data types
dat.ALL[,c(2:15)] = lapply(dat.ALL[,c(2:15)], as.factor)
dat.ALL[,c(16:ncol(dat.ALL))] = lapply(dat.ALL[,c(16:ncol(dat.ALL))], as.numeric)
str(dat.ALL)

Best.Res.Er.Str.LRT = NULL
Best.Res.Er.Str.BIC = NULL
LocYearTrait.vec = NULL
blues_toPrint_all = NULL

for(env in env_vect){
  print(env)
  dat_thisEnv = dat.ALL[which(dat.ALL$Field == env),]
  
  #__________________________________________________________________________#
  # Spatial Model Fitting for single environments ----
  #__________________________________________________________________________#
  
  # Define random structure
    rand_form <- as.formula("~ LDP + Rep + Bed + Tier")
    rand_form_BLUE <- as.formula("~ Rep + Bed + Tier")
  
  for(this.trait in trait_vect){
    
    LocYearTrait <- paste(env, this.trait, sep = "_")
    
    # Fit candidate models varying spatial error structures
    model.cand1 = asreml(
      fixed = as.formula(paste(this.trait, "~ 1")), 
      random = rand_form,
      residual = ~id(units), 
      na.action = na.method(y = 'include', x = 'include'), 
      data=dat_thisEnv, workspace='3.0gb', maxit = 50)
    
    model.cand2 = asreml(
      fixed = as.formula(paste(this.trait, "~ 1")), 
      random = rand_form,
      residual = ~id(Bed):ar1(Tier), 
      na.action=na.method(y = 'include', x = 'include'), 
      data=dat_thisEnv, workspace='3.0gb', maxit = 50)
    
    model.cand3 = asreml(
      fixed = as.formula(paste(this.trait, "~ 1")), 
      random = rand_form,
      residual = ~ar1(Bed):id(Tier), 
      na.action=na.method(y = 'include', x = 'include'), 
      data=dat_thisEnv, workspace='3.0gb', maxit = 50)
    
    model.cand4 = asreml(
      fixed = as.formula(paste(this.trait, "~ 1")), 
      random = rand_form, 
      residual = ~ar1(Bed):ar1(Tier), 
      na.action=na.method(y = 'include', x = 'include'), 
      data=dat_thisEnv, workspace='3.0gb', maxit = 50)
    
    #Combine models into list
    results.list = list(model.cand1, model.cand2, model.cand3, model.cand4)
    
    #Calculate BIC to choose best residual error structure
    ## Note: recommend using built in summary(mod)$bic functionality (AsREML-R 4) if using BIC to determine best residual structure, but note that that parameter seems to round;
    ## good to keep an eye on that rounding if models exhibit highly similar BIC values in case affects nuanced choice between models (or if lack of rounding helps break ties)
      BICs = sapply(results.list, function(mod) {
      log(mod$nedf) * length(mod$vparameters) - (2 * mod$loglik) #Changed to mod$vparameters (AsREML-R 4) from mod$gammas (AsREML-R 3)
    })
    
    #Check to find the best model with BIC
    best.index = which.min(BICs)
    best.mod = results.list[[best.index]] 
    
    #Likelihood ratio test to choose best residual error structure
    LRT4.2 = lrt(model.cand4, model.cand2)$`Pr(Chisq)` #step down bed
    LRT4.3 = lrt(model.cand4, model.cand3)$`Pr(Chisq)` #step down tier
    LRT2.1 = lrt(model.cand2, model.cand1)$`Pr(Chisq)` #step down both
    
    # Check model for param errors (B=boundary, F=fixed, C=constrained)
    has_boundary <- function(mod) {
      any(mod$vparameters.con %in% c("B", "F", "C"))
    }
    cand1_boundary <- has_boundary(model.cand1)
    cand2_boundary <- has_boundary(model.cand2)
    cand3_boundary <- has_boundary(model.cand3)
    cand4_boundary <- has_boundary(model.cand4)
    BoundSumm <- c(cand1_boundary, cand2_boundary, cand3_boundary, cand4_boundary)
    
    # Store results in a table
    model.result = data.frame(
      Trait = this.trait,
      Env = env,
      Model = c(1,2,3,4),
      Res_formula = c("~id(units)", 
                      "~(id(Bed):ar1(Tier))", 
                      "~ar1(Bed):id(Tier)", 
                      "~ar1(Bed):ar1(Tier)"),
      BIC = BICs,
      LRT_p = c(LRT2.1, LRT4.2, LRT4.3, NA),
      LRT_formula = c("LRT2.1", "LRT4.2", "LRT4.3","NA"),
      Boundary = BoundSumm
    )
    
    #Initialize or append
    if (!exists("model.results") || is.null(model.results)) {
      model.results <- model.result
    } else {
      model.results <- rbind(model.results, model.result)
    }
    
    # Set significance level
    alpha = 0.05
    
    # Start with the most complex model
    Lrt.best.mod = model.cand4
    Lrt.best.index = 4
    
    # Step 1: Compare Model 4 vs Model 2 (removing AR1 on Bed)
    if (LRT4.2 > alpha) {        #LRT tests if we are losing information  
      Lrt.best.mod = model.cand2 #by going to a simpler model (sig = yes) 
      Lrt.best.index = 2  
    }
    
    # Step 2: Compare Model 4 vs Model 3 (removing AR1 on Tier) 
    if (LRT4.3 > alpha && Lrt.best.index == 4) {  
      Lrt.best.mod = model.cand3  
      Lrt.best.index = 3  
    }
    
    # Step 3: Compare Model 2 vs Model 1 (removing AR1 on both) 
    if (Lrt.best.index == 2 && LRT2.1 > alpha) {  
      Lrt.best.mod = model.cand1  
      Lrt.best.index = 1  
    }
    
    #Re-assign best fit as scored by LRT to 'mod.to.fit' #changed from LRT 9/8/26, changed back 10/22 for clarity as this script uses LRT (as seen in next line)
    mod.to.fit = Lrt.best.mod #changed from best.mod 9/21/26
    
    #Check to see if LRT and BIC agree
    if(best.index == Lrt.best.index){
      cat(this.trait, env, ": Model selected by LRT is IDENTICAL to BIC")
    }else if(best.index != Lrt.best.index){
      cat(this.trait, env, ": Model selected by LRT is DIFFERENT from BIC")
    }
    
    #Save best residual error structure from BIC and LRT
    best.res.LRT = paste0(deparse(Lrt.best.mod$formulae$residual))
    best.res.BIC = paste0(deparse(best.mod$formulae$residual))
    #Grow/Append to vectors each iteration 
    Best.Res.Er.Str.LRT = c(Best.Res.Er.Str.LRT, best.res.LRT)
    Best.Res.Er.Str.BIC = c(Best.Res.Er.Str.BIC, best.res.BIC)
    LocYearTrait.vec = c(LocYearTrait.vec, LocYearTrait)
    
    if(is.null(rand_form_BLUE)){ # If excluding tier and bed from random effects
      
      blues.mod.to.fit <- asreml(
        fixed    = as.formula(paste(this.trait, "~ LDP")),
        residual = ~id(units), #changed from mod.to.fit$formulae$residual because singularity
        data     = dat_thisEnv,
        na.action = na.method(y = "include", x = "include")
      )
      
    } else {
      
      blues.mod.to.fit <- asreml(
        fixed    = as.formula(paste(this.trait, "~ LDP")),
        random   = rand_form_BLUE,
        residual = mod.to.fit$formulae$residual,
        data     = dat_thisEnv,
        na.action = na.method(y = "include", x = "include")
      )
    }
    
    # predict BLUEs
    blues <- predict(blues.mod.to.fit, classify = "LDP")
    blues_toPrint  <- data.frame(
      LDP   = as.character(blues$pvals$LDP),
      value = blues$pvals$predicted.value)
    colnames(blues_toPrint)[2] = LocYearTrait
    
    if (is.null(blues_toPrint_all)) {
      blues_toPrint_all <- blues_toPrint
      
    } else {
      blues_toPrint_all <- merge(blues_toPrint_all,
                                 blues_toPrint,
                                 by = "LDP", all = TRUE)
      
    }
  }
}

#Print model residual error structure testing performance metrics to file
LRTvsBIC_res = as.data.frame(
  cbind(LocYearTrait.vec, Best.Res.Er.Str.LRT, Best.Res.Er.Str.BIC))

model.results$LocYearTrait.vec <- paste(model.results$Env,
                                        model.results$Trait, sep = "_")

model.results.LRTvsBIC <- merge(model.results, LRTvsBIC_res, by = "LocYearTrait.vec", all=T)

write.table(model.results.LRTvsBIC,
            paste("Model_Performance_Comparison_BestResErrorStrucs_compiled",
                  Sys.Date(), "_FT_usingLRT.txt", sep = "_"),
            sep='\t',quote=F,row.names=F)

write.table(blues_toPrint_all,paste0("AllEnvsAllTraits_BLUEs_",
        Sys.Date(), "_FT_usingLRT.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
