rm(list=ls()) #added 10/7/26

### load libraries
library(statgenGxE)
library(ggplot2)
library(tidyr)
library(dplyr)
library(forcats)
library(ggpubr)
library(gridExtra)
library(statgenSTA) #added 10/7/26


### skip to 513

#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} methods
#' are available.
#'
#' @param fitMod A fitted variance components model.
#' @param modDat A data.frame containing the data used in fitting the model.
#'
#' @keywords internal
createVarComp <- function(fitMod,
                          modDat,
                          trait,
                          nestingFactor,
                          useLocYear,
                          useRegionLocYear,
                          fullRandVC,
                          aovFullFixedMod,
                          engine,
                          diagTabs) {
  varComp <- structure(list(fitMod = fitMod,
                            modDat = modDat,
                            trait = trait,
                            nestingFactor = nestingFactor,
                            useLocYear = useLocYear,
                            useRegionLocYear = useRegionLocYear,
                            fullRandVC = fullRandVC,
                            aovFullFixedMod = aovFullFixedMod,
                            engine = engine,
                            diagTabs = diagTabs),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @export
print.varComp <- function(x,
                          ...) {
  summary(x)
}

#' @export
summary.varComp <- function(object,
                            ...) {
  engine <- object$engine
  trait <- object$trait
  if (engine == "lme4") {
    ## Display model formula in text form.
    ## This might cut off the formula if it is longer than 500 character.
    ## This is however highly unlikely given the fixed structure of the models.
    fitModCall <- deparse(formula(getCall(object$fitMod)), width.cutoff = 500)
  } else if (engine == "asreml") {
    ## For asreml, rewrite the formula to match the display of formulas
    ## from lme4,
    ## No changes needed in the fixed part.
    fitModCallFixed <- deparse(object$fitMod$call$fixed, width.cutoff = 500)
    ## The terms in the random part need to be surrounded by (1 | ...).
    fitModCallRandTerms <- attr(x = terms(object$fitMod$call$random),
                                which = "term.labels")
    fitModCallRand <- paste0("(1 | ", fitModCallRandTerms, ")",
                             collapse = " + ")
    ## Concatenate fixed and random part to form the full formula.
    fitModCall <- paste(fitModCallFixed, "+", fitModCallRand)
  }
  ## Print source of variation as percentage.
  fullRandVC <- object$fullRandVC
  ## Prevent scientific notation in variance component.
  fullRandVC[["vcov"]] <- sprintf("%1.2f", fullRandVC[["vcov"]])
  if (engine == "asreml") {
    fullRandVC[["stdError"]] <- sprintf("%1.3f", fullRandVC[["stdError"]])
  }
  fullRandVC[["vcovPerc"]] <- sprintf("%1.2f %%", 100 * fullRandVC[["vcovPerc"]])
  colnames(fullRandVC) <- c("Component", if (engine == "asreml") "SE",
                            "% Variance expl.")
  aovFullFixedMod <- object$aovFullFixedMod
  ## Construct fully fixed and fully random model formula from ANOVA.
  modTerms <- rownames(aovFullFixedMod)[-nrow(aovFullFixedMod)]
  fixModCall <- paste(trait, "~", paste0(modTerms, collapse = " + "))
  randModCall <- paste(trait, "~",
                       paste0("(1 | ", modTerms, ")", collapse = " + "))
  ## Print ANOVA for fully fixed model with alternative header.
  attr(x = aovFullFixedMod, which = "heading") <-
    paste("Analysis of Variance Table for fully fixed model:\n",
          fixModCall, "\n")
  ## Print output.
  cat("Fitted model formula final mixed model\n\n")
  cat("", fitModCall, "\n\n")
  cat("Sources of variation for fully random model:\n")
  cat("", randModCall, "\n\n")
  print(fullRandVC)
  cat("\n")
  print(aovFullFixedMod)
}

plot.varComp <- function(x,
                         ...,
                         plotType = c("sd", "percVar"),
                         title = NULL,
                         output = TRUE) {
  plotType <- match.arg(plotType)
  chkChar(title, len = 1)
  ## Extract mu from the fitted model.
  mu <- round(mean(fitted(x$fitMod)))
  if (is.null(title)) {
    title <- paste0(ifelse(plotType == "sd", "Standard deviations",
                           "Percentage of variance explained"),
                    " (general mean = ", mu, ")")
  }
  ## The actual variable to plot depends on the plotType.
  plotVar <- if (plotType == "sd") "sd" else "vcovPerc"
  ## Extract var comps for random model and anova for fixed model.
  fullRandVC <- x$fullRandVC
  aovFullFixedMod <- x$aovFullFixedMod
  ## Degrees of freedom in the plot come from the anova for the fixed model.
  ## These need to be merged to the random model.
  ## lm (used for the anova) orders interaction terms alphabetically.
  ## To merge these interactions need to be split into their terms and
  ## matching can be done on the terms.
  fullRandVC[["vars"]] <- sapply(X = strsplit(x = rownames(fullRandVC),
                                              split = ":"),
                                 FUN = function(var) {
                                   paste0(sort(var), collapse = "_")
                                 })
  aovFullFixedMod[["vars"]] <- sapply(X = strsplit(x = rownames(aovFullFixedMod),
                                                   split = ":"),
                                      FUN = function(var) {
                                        paste0(sort(var), collapse = "_")
                                      })
  fullRandVC[["Df"]] <- aovFullFixedMod[["Df"]][match(aovFullFixedMod[["vars"]],
                                                      fullRandVC[["vars"]])]
  ## term will be used as y-axis label. It consists of term + df
  spaces <- sapply(X = 6 - nchar(fullRandVC[["Df"]]),
                   FUN = function(i) {
                     paste0(rep(" ", times = i), collapse = "")
                   })
  fullRandVC[["term"]] <- paste(rownames(fullRandVC), spaces,
                                fullRandVC[["Df"]])
  ## Revert levels term to get a nice ordering on the y-axis.
  fullRandVC[["term"]] <- factor(fullRandVC[["term"]],
                                 levels = rev(fullRandVC[["term"]]))
  ## Compute standard deviations.
  fullRandVC[["sd"]] <- sqrt(fullRandVC[["vcov"]])
  ## Compute position for annotation. 5e5 is found empirically.
  ## The idea is to annotate the y-axis just left of the x = 0.
  annoPosX <- -max(fullRandVC[[plotVar]]) / 5e5
  p <- ggplot2::ggplot(fullRandVC,
                       ggplot2::aes(x = .data[[plotVar]],
                                    y = .data[["term"]])) +
    ggplot2::geom_point(na.rm = TRUE, size = 2) +
    ## Add line from y-axis to points.
    ggplot2::geom_segment(ggplot2::aes(xend = .data[[plotVar]],
                                       yend = .data[["term"]]),
                          x = 0) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ## Set lower xlim to 0. This assures 0 is always displayed on the x-axis
    ## even if the lowest variance component is e.g. 1e-8.
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(color = "grey50"),
                   axis.line = ggplot2::element_line(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks.length.y = grid::unit(0, "mm"),
                   axis.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::annotation_custom(grid::textGrob("Source      df ", just = "right",
                                              gp = grid::gpar(size = 14)),
                               xmin = annoPosX, xmax = annoPosX,
                               ymin = Inf, ymax = Inf) +
    ggplot2::ggtitle(title)
  if (plotType == "sd") {
    p <- p + ggplot2::labs(x = "Square root of variance estimate")
  } else if (plotType == "percVar") {
    p <- p + ggplot2::labs(x = "Percentage of variance explained")
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

vc <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  fitMod <- varComp$fitMod
  ## Extract variance component and rename so rows/columns to assure
  ## matching outputs for lme4/asreml.
  if (varComp$engine == "lme4") {
    modTerms <- colnames(attr(x = terms(getCall(fitMod)$formula,
                                        keep.order = TRUE), which = "factors"))
    modTerms <- gsub(pattern = "1 | ", replacement = "", x = modTerms,
                     fixed = TRUE)
    varcomps <- as.data.frame(lme4::VarCorr(fitMod))
    rownames(varcomps) <- varcomps[["grp"]]
    rownames(varcomps)[nrow(varcomps)] <- "residuals"
    modTermsRand <- modTerms[modTerms %in% rownames(varcomps)]
    varcomps <- varcomps[c(modTermsRand, "residuals"), "vcov", drop = FALSE]
    colnames(varcomps) <- "Component"
  } else if (varComp$engine == "asreml") {
    modTerms <- colnames(attr(x = terms(fitMod$call$random, keep.order = TRUE),
                              which = "factors"))
    varcomps <- summary(fitMod)$varcomp
    rownames(varcomps)[nrow(varcomps)] <- "residuals"
    varcomps <- varcomps[c(modTerms, "residuals"), c("component", "std.error")]
    colnames(varcomps) <- c("Component", "SE")
  }
  return(varcomps)
}

#' Calculate heritability


herit <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Extract fitted model and model data.
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  ## Compute variance components.
  varcomps <- vc(varComp)
  ## Extract variance components for genotype and residual.
  sigmaG <- varcomps["genotype", "Component"]
  sigmaRes <- varcomps["residual", "Component"]
  ## Numerator is constructed by looping over all random model terms and
  ## Adding their share. It always includes sigmaG.
  numerator <- sigmaG
  ## Get the terms used in the random part of the model.
  modTerms <- rownames(varcomps)
  ## Extract all variables used in the random part of the model.
  ## The are needed for computing the contribution of the residual variance.
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]
    
  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random),
                             which = "factors"))[-1]
  }
  ## Get median number for times genotypes are tested within modVars.
  nLevModVars <- sapply(X = modVars, FUN = function(modVar) {
    median(rowSums(table(modDat[["genotype"]], modDat[[modVar]]) > 0))
  })
  for (term in modTerms[-c(1, length(modTerms))]) {
    ## Get variance for current term.
    sigmaTerm <- varcomps[term, "Component"]
    ## Get variables in current term, exclude genotype (always the first var).
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    ## Divide variance by product of #levels for all variables in current term.
    ## Add that to numerator.
    numerator <- numerator + sigmaTerm / prod(nLevModVars[termVars])
  }
  nReps <- median(table(modDat[["genotype"]], modDat[["trial"]]))
  if (length(modVars) > 0) {
    ## Contribution for residual variance is computed by dividing sigmaRes by
    ## product of #levels of all variables in random part of model and
    ## #replicates.
    numerator <- numerator + sigmaRes / prod(nLevModVars, nReps)
  } else {
    ## No other variables in random part.
    ## Just divide sigmaRes by #replicates.
    numerator <- numerator + sigmaRes / nReps
  }
  return(sigmaG / numerator)
}

#' Calculate the correlated response to selection
#'
#' Calculate the correlated response to selection (CRDR) based on the fitted
#' model. The CRDR is calculated as described by Atlin et al. E.g. for a model
#' with trials nested within scenarios, which has a random part that looks like
#' this: genotype + genotype:scenario + genotype:scenario:trial the CRDR is
#' calculated as:\cr\cr
#' \deqn{H1 = \sigma_G^2 / (\sigma_G^2 + \sigma_S^2 / s + \sigma_{ST}^2 / st +
#' \sigma_E^2 / str)}
#' \deqn{H2 = (\sigma_G^2 + \sigma_S^2) / (\sigma_G^2 + \sigma_S^2 +
#' \sigma_{ST}^2 / st + \sigma_E^2 / str)}
#' \deqn{CRDR = (\sigma_G^2 / (\sigma_G^2 + \sigma_S^2)) * sqrt(H1 / H2)}
#' In these formulas the \eqn{\sigma} terms stand for the standard deviations of
#' the respective model terms, and the lower case letters for the number of
#' levels for the respective model terms. So \eqn{\sigma_S} is the standard
#' deviation for the scenario term in the model and \eqn{s} is the number of
#' scenarios. \eqn{\sigma_E} corresponds to the residual standard deviation and
#' \eqn{r} to the number of replicates.
#'
#' @inheritParams herit
#'
#' @references Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000).
#' Selection response in subdivided target regions. Crop Science, 40(1), 7–13.
#' \doi{10.2135/cropsci2000.4017}
#'
#' @family Mixed model analysis
#'
#' @export
CRDR <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  if (is.null(varComp$nestingFactor) && isFALSE(varComp$useRegionLocYear)) {
    stop("CRDR can only be computed when a model is fitted with a nesting ",
         "structure or when regions are included in the model.\n")
  }
  H1 <- herit(varComp)
  ## Extract fitted model and model data.
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  ## Get factor for computing H2
  if (!is.null(varComp$nestingFactor)) {
    H2factor <- paste0("genotype:", varComp$nestingFactor)
  } else if (varComp$useRegionLocYear) {
    H2factor <- "genotype:region"
  }
  ## Compute variance components.
  varcomps <- vc(varComp)
  ## Extract variance components for genotype, H2factor and residual.
  sigmaG <- varcomps["genotype", "component"]
  sigmaH2 <- varcomps[H2factor, "component"]
  sigmaRes <- varcomps["residual", "component"]
  ## Numerator is constructed by looping over all random model terms and
  ## Adding their share. It always includes sigmaG and sigmaH2.
  numerator <- sigmaG + sigmaH2
  ## Get the terms used in the random part of the model.
  modTerms <- rownames(varcomps)
  ## Extract all variables used in the random part of the model.
  ## They are needed for computing the contribution of the residual variance.
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]
    
  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random),
                             which = "factors"))[-1]
  }
  ## Get median number for times genotypes are tested within modVars.
  nLevModVars <- sapply(X = modVars, FUN = function(modVar) {
    median(rowSums(table(modDat[["genotype"]], modDat[[modVar]]) > 0))
  })
  H2factorPos <- which(modTerms == H2factor)
  for (term in modTerms[-c(1, H2factorPos, length(modTerms))]) {
    ## Get variance for current term.
    sigmaTerm <- varcomps[term, "component"]
    ## Get variables in current term, exclude genotype (always the first var).
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    ## Divide variance by product of #levels for all variables in current term.
    ## Add that to numerator.
    numerator <- numerator + sigmaTerm / prod(nLevModVars[termVars])
  }
  nReps <- median(table(modDat[["genotype"]], modDat[["trial"]]))
  if (length(modVars) > 0) {
    ## Contribution for residual variance is computed by dividing sigmaRes by
    ## product of #levels of all variables in random part of model and
    ## #replicates.
    numerator <- numerator + sigmaRes / prod(nLevModVars, nReps)
  } else {
    ## No other variables in random part.
    ## Just divide sigmaRes by #replicates.
    numerator <- numerator + sigmaRes / nReps
  }
  H2 <- (sigmaG + sigmaH2) / numerator
  r <- sigmaG / (sigmaG + sigmaH2)
  return(r * sqrt(H1 / H2))
}

#' Compute different types of correlations.
#'
#' Compute three types of correlations for models fitted with a nesting factor.
#' \itemize{
#' \item correlation between scenarios or environment types:
#' \deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2)}
#' \item correlation between trials within scenarios or environment types:
#' \deqn{(\sigma_G^2 + \sigma_{GS}^2) / (\sigma_G^2 + \sigma_{GS}^2 +
#' \sigma_E^2)}
#' \item correlation between trials that belong to different
#' scenarios/environment types:
#' \deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2 + \sigma_E^2)}
#' }
#' In these formulas the \eqn{\sigma} terms stand for the standard deviations of
#' the respective model terms. So \eqn{\sigma_S} is the standard deviation for
#' the scenario term in the model, \eqn{\sigma_{GS}} for the standard deviation
#' of the genotype by scenario term and \eqn{\sigma_E} corresponds to the
#' residual standard deviation.
#'
#' @inheritParams herit
#'
#' @return A list with three correlations.
#'
#' @family Mixed model analysis
#'
#' @export
correlations <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Get factor for computing correlations.
  if (!is.null(varComp$nestingFactor)) {
    corFactor <- paste0("genotype:", varComp$nestingFactor)
  } else {
    stop("correlations can only be computed when a model is fitted with a ",
         "nesting structure.\n")
  }
  ## Compute variance components.
  varComps <- vc(varComp)
  varGeno <- varComps["genotype", "component"]
  varCorFactor <- varComps[corFactor, "component"]
  varRes <- varComps["residuals", "component"]
  ## Compute correlation between scenarios.
  rScen <- varGeno / (varGeno + varCorFactor)
  ## Compute correlation between trials within scenarios.
  rTrScen <- (varGeno + varCorFactor) / (varGeno + varCorFactor + varRes)
  ## Compute correlattion between trials that belong to different scenarios.
  rTrDiffScen <- varGeno / (varGeno + varCorFactor + varRes)
  return(list(rScen = rScen, rTrScen = rTrScen, rTrDiffScen = rTrDiffScen))
}

#' Get diagnostics for an object of class varComp
#'
#' Get the diagnostics for the model fitted. This will print a table of
#' combinations missing in the data. For each random factor in the model a
#' table is printed.
#'
#' @param varComp An object of class varComp.
#'
#' @return A list of tables is invisibly returned.
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Display diagnostics.
#' diagnostics(geVarComp)
#'
#' @family Mixed model analysis
#'
#' @export
diagnostics <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Get diagTabs from varComp.
  diagTabs <- varComp$diagTabs
  ## For each diagTab print either its content or display a message that
  ## it has no content.
  for (diagTab in diagTabs) {
    if (nrow(diagTab) > 0) {
      cat(nrow(diagTab), " missing combinations for ",
          paste(colnames(diagTab), collapse = " x "), ".\n", sep = "")
      print(diagTab, row.names = FALSE)
      cat("\n\n")
    } else {
      cat("No missing combinations for ",
          paste(colnames(diagTab), collapse = " x "), ".\n\n", sep = "")
    }
  }
  invisible(diagTabs)
}


###
setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting/")

#LIMA
#nirs_6E<-read.csv("/Users/jaadaskaveg/Library/CloudStorage/Box-Box/PostDoc/LIMA_USDA_SCRI/Lima-geno-pheno/2.4.26_predicted_5fCV_nopreteat_final6E_meta_updated.csv", header=T) #CHD commented
#nirs_6E<-read.csv("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting/2.4.26_predicted_5fCV_nopreteat_final6E_meta_updated_meta_added_withdb.csv",header=T) #CHD added
FT_2E<-read.csv("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting/Flowering_Davis24_25_gxe.csv",header=T,stringsAsFactors = T) #added stringsAsFactors = T on 10/14/26
str(FT_2E) # note experimental vars not factors, check that is what statgenGxE wants

FT_2E$Pop = NA
FT_2E[FT_2E$Kgroup == "K1", "Pop"] = "Meso"
FT_2E[FT_2E$Kgroup == "K2", "Pop"] = "B"
FT_2E[FT_2E$Kgroup == "K3", "Pop"] = "Meso"
FT_2E[FT_2E$Kgroup == "K4", "Pop"] = "Meso"
FT_2E[FT_2E$Kgroup == "K5", "Pop"] = "Meso"
FT_2E[FT_2E$Kgroup == "K6", "Pop"] = "C"
FT_2E[FT_2E$Kgroup == "K7", "Pop"] = "C"
FT_2E[FT_2E$Kgroup == "K8", "Pop"] = "Meso"
FT_2E[FT_2E$Kgroup == "admixed", "Pop"] = "Admixed"
FT_2E$Location = "Davis" #CHD added as this column was not present for FT since all in Davis

FT_all<-FT_2E%>%filter(Pop!="NA")
FT_meso<-FT_2E%>%filter(Pop=="Meso")
FT_and<-FT_2E%>%filter(Pop=="B")
FT_c<-FT_2E%>%filter(Pop=="C")
FT_adm<-FT_2E%>%filter(Pop=="Admixed")

#Next 3 lines were commented 10/11, not relevant for FT since only in Davis
#nirs_cvars<-nirs_meso%>%filter(Location=="Coachella")
#nirs_cf<-nirs_meso%>%filter(Location=="Central.Ferry")
#nirs_ucd<-nirs_meso%>%filter(Location=="Davis")

dat.TD <- createTD(data = FT_meso, #CHD changed from _all
                   genotype = "LDP", #CHD changed from Genotype
                   trial = "Field",
                   year = "Year",
                   loc = "Location",
                   repId = "Rep", #uncommented 10/7/26
                   subBlock = "Block",
                   rowCoord = "Tier", 
                   colCoord = "Bed")

TD = dat.TD

trials = names(TD)

traits <- names(FT_2E)[16] #dpp

TDTot <- Reduce(f = rbind, x = TD[trials])

TDTot[["wt"]] <- 1

# Commented next two lines 10/7/26 as they were throwing an error and are not used
# RepTrialTab <- table(TDTot[["Rep"]], TDTot[["trial"]])
# RepTrialTab <- all(RepTrialTab > 0)

 
TDTot <- droplevels(TDTot)
## Increase maximum number of iterations for asreml.
maxIter <- 200
  ## Asreml can't handle missing weights, so set them to 0 when missing.
TDTot[is.na(TDTot[["wt"]]), "wt"] <- 0
  ## Construct formula for fixed part - first as text.
  ## Trying to fit this in something 'smart' actually makes it unreadable.
  ## First create a vector with the separate terms.
  ## This avoids difficult constructions to get the +-es correct.
fixedTerms <- c("trial") #CHD changed from loc
  
## Check if the data contains replicates.
  repTab <- table(TDTot[c("genotype",
                          unlist(strsplit(x = tail(fixedTerms, 1),
                                          split = ":")))])
  hasReps <- any(repTab > 1)
  ## Random terms are genotype x fixedTerms.
  ## If there are no replicates or weights the final random term is the actual
  ## residual and therefore left out of the model.
  if (hasReps) { #} || useWt) {
    randTermIncl <- fixedTerms
  } else {
    randTermIncl <- fixedTerms[-length(fixedTerms)]
  }
  randTerms <- c("genotype","genotype:trial","trial:repId") #CHD changed from "genotype","loc:year","genotype:loc","loc:year:Rep","genotype:loc:year"
  ## First fit a model with all terms fixed to determine:
  ## - should all terms in the fixed part really be present.
  ## - Predict which terms in the random part of the model will probably
  ##   have a zero variance component.
  
  #used traits for lima (above)

  df <- data.frame(matrix(ncol = 0, nrow = 5)) #CHD changed from 7 to 5
  df_herit <- data.frame(matrix(ncol=0,nrow =1))
  
  
  #modified traits_of_interest to be traits in loop for(i in traits_of_interest) 
  for (i in traits){
    trait = i
    
    fullFixedTxt <- paste0("`", trait, "`~",
                           paste(c(randTerms), collapse = "+"))
    
    ## Fit the fully fixed model.
    fullFixedMod <- suppressWarnings(lm(formula(fullFixedTxt), data = TDTot))
    
    aovFullFixedMod <- anova(fullFixedMod)
    rownames(aovFullFixedMod)[nrow(aovFullFixedMod)] <- "residuals"
    ## Reorder terms to the original order in the model call.
    modTerms <- c(fixedTerms, randTerms)
    #modTerms <- c(randTerms)
    
    aovVars <- sapply(X = strsplit(x = rownames(aovFullFixedMod)[-nrow(aovFullFixedMod)],
                                   split = ":"), FUN = function(var) {
                                     paste0(sort(var), collapse = "_")
                                   })
    modVars <- sapply(X = strsplit(x = modTerms, split = ":"),
                      FUN = function(var) {
                        paste0(sort(var), collapse = "_")
                      })
    aovFullFixedMod <- aovFullFixedMod[c(match(modVars, aovVars),
                                         nrow(aovFullFixedMod)), ]
    rownames(aovFullFixedMod) <- c(modTerms, "residuals")
    ## Get all model terms as used by lm (might involve reordered terms).
    fullFixedLabs <- attr(x = terms(fullFixedMod), which = "term.labels")
    if (!all(fullFixedLabs %in% rownames(aovFullFixedMod))) {
      ## At least one terms missing from ANOVA.
      ## If this is a fixed term remove it from fixed.
      missTerms <- fullFixedLabs[!fullFixedLabs %in% rownames(aovFullFixedMod)]
      for (missTerm in missTerms) {
        fixedTermSets <- strsplit(x = fixedTerms, split = ":")
        missTermSet <- unlist(strsplit(x = missTerm, split = ":"))
        remPos <- !sapply(X = fixedTermSets, FUN = setequal, missTermSet)
        fixedTerms <- fixedTerms[remPos]
      }
    }
    ## Get rand terms that indicate zero variance components.
    ## lm reorders variables in model terms.
    ## use sets of variables in terms to compare them.
    aovTermSets <- strsplit(x = rownames(aovFullFixedMod), split = ":")
    for (randTerm in randTerms) {
      ## Convert term to set of variables in term.
      randTermSet <- unlist(strsplit(x = randTerm, split = ":"))
      ## Get position of term in ANOVA table by comparing sets.
      randTermPos <- sapply(X = aovTermSets, FUN = setequal, randTermSet)
      ## Get MSS for current term.
      MSSRandTerm <- aovFullFixedMod[randTermPos, "Mean Sq"]
      if (is.na(MSSRandTerm)) MSSRandTerm <- 0
      ## For all other terms in the ANOVA table that have the current term
      ## as a subset the MSS cannot be higher.
      ## If it is the corresponding variance component is possibly zero.
      for (i in seq_along(aovTermSets)) {
        if ((all(randTermSet %in% aovTermSets[i]) ||
             ## Always include the residual term for comparison.
             i == nrow(aovFullFixedMod)) &&
            !is.nan(aovFullFixedMod[i, "Mean Sq"]) &&
            aovFullFixedMod[i, "Mean Sq"] > MSSRandTerm) {
          warning("Mean Sum of Squares for ", randTerm, " smaller than Mean ",
                  "Sum of Squares for ", rownames(aovFullFixedMod)[i], ".\n",
                  "Possible zero variance components.\n", call. = FALSE)
        }
      }
    }
    ## Fit the fully random model and compute how much variation
    ## is explained by each of the model terms.
    ## This is stored as fullRandVC and included in the output to create
    ## a nice summary.
    engine = "lme4"
    if (engine == "lme4") {
      ## Construct input for full random model.
      fullRandTxt <- paste0("`", trait, "`~",
                            paste(paste0("(1|", c(fixedTerms, randTerms), ")"),
                                  collapse = "+"))
      fullRandMod <- lme4::lmer(formula(fullRandTxt), data = TDTot)
      fullRandVC <- as.data.frame(lme4::VarCorr(fullRandMod))
      rownames(fullRandVC) <- fullRandVC[["grp"]]
      rownames(fullRandVC)[nrow(fullRandVC)] <- "residuals"
      vcovTot <- sum(fullRandVC[["vcov"]])
      fullRandVC[["vcovPerc"]] <- fullRandVC[["vcov"]] / vcovTot
    }
    ## Reorder rows and vars within terms in rownames to match orginal
    ## function call.
    VCVars <- sapply(X = strsplit(x = rownames(fullRandVC)[-nrow(fullRandVC)],
                                  split = ":"), FUN = function(var) {
                                    paste0(sort(var), collapse = "_")
                                  })
    randVars <- sapply(X = strsplit(x = modTerms, split = ":"),
                       FUN = function(var) {
                         paste0(sort(var), collapse = "_")
                       })
    fullRandVC <- fullRandVC[c(match(randVars, VCVars), nrow(fullRandVC)),
                             colnames(fullRandVC) %in%
                               c("vcov", "stdError", "vcovPerc")]
    ## Create tables for diagnostics.
    diagTabs <- lapply(X = fixedTerms, FUN = function(fixedTerm) {
      fixedVars <- unlist(strsplit(x = fixedTerm, split = ":"))
      fixedVarLevs <- lapply(X = fixedVars, FUN = function(fixedVar) {
        unique(TDTot[[fixedVar]])
      })
      fullTab <- expand.grid(c(list(unique(TDTot[["genotype"]])), fixedVarLevs),
                             KEEP.OUT.ATTRS = FALSE)
      missTab <- fullTab[!interaction(fullTab) %in%
                           interaction(TDTot[c("genotype", fixedVars)]), ]
      colnames(missTab) <- c("genotype", fixedVars)
      return(missTab)
    })
    ## Create the full fixed part of the model as a character.
    ## This is identical for lme4 and asreml so only needs to be done once.
    fixedTxt <- paste0("`", trait, "`~", paste(fixedTerms, collapse = "+"))
    if (engine == "lme4") {
      randTxt <- paste(paste0("(1|", randTerms, ")"), collapse = "+")
      formTxt <- paste(fixedTxt, "+", randTxt)
      ## Fit the actual model.
      mr <- lme4::lmer(formula(formTxt), data = TDTot, weights = TDTot[["wt"]])
      ## Construct STA object.
    } 
    
    ## Create output.
    
    nestingFactor = NULL
    locationYear = NULL
    regionLocationYear = NULL
    
    res <- createVarComp(fitMod = mr, modDat = TDTot, trait = trait,
                         nestingFactor = nestingFactor, useLocYear = locationYear,
                         useRegionLocYear = regionLocationYear,
                         fullRandVC = fullRandVC,
                         aovFullFixedMod = aovFullFixedMod, engine = engine,
                         diagTabs = diagTabs)
    summary(res)
    dropsVarComp <- res
    df <- cbind(df, dropsVarComp$fullRandVC$vcovPerc)
    df_herit <- cbind(df_herit,herit(dropsVarComp))
    
  }
  
  colnames(df)<-traits 
  rownames(df) <- c("trial", "genotype","genotype:trial", "trial:repId","residuals")
  
  colnames(df_herit) <- traits 
  
  write.table(df_herit,"GxE_and_Herit/Herit_statgenGxE_collapsedE_FT_Meso.txt",sep='\t',quote=F,row.names=F)
  
  long <- df %>% 
    tibble::rownames_to_column(var = "Row") %>%
    pivot_longer(
      cols = `dpp`, #Changed from HSW`:`Starch.db 
      names_to = "trait",
      values_to = "value"
    ) %>% arrange(Row,value)
  
  long$trait <- factor(long$trait, levels = traits
  )
  
  long$Row <- factor(long$Row, levels = c("trial", "genotype", "genotype:trial", "trial:repId", "residuals")
  )
  
  long1 <- long %>%
    mutate(var1 = recode(Row, 'trial' = 'Environment','genotype' = 'Genotype', 'genotype:trial' = 'Genotype:Environment', 'trial:repId' = 'Environment:Rep', 'residuals'='Residuals'))
  
  write.table(long,"GxE_and_Herit/VarComps_statgenGxE_collapsedE_FT_Meso.txt",sep='\t',quote=F,row.names=F)
  
  ggplot(long, aes(fill=Row, y=value, x=trait, label=ifelse(value>0.001, round(value,digits = 3),""))) + 
    geom_bar(position="fill", stat="identity") +
    scale_x_discrete(limits=c("dpp"), #Changed from "HSW", "Protein.db", "Starch.db", "Fat.db", "Ash.db" 
                     labels=c("FT_dpp"))+ #Changed from "HSW", "Protein", "Starch","Fat","Ash"
    theme_minimal() +
    ylab("Proportion Variance Explained") +
    theme(axis.text.x = element_text(  color = "black"), #angle = 45, hjust=1
          axis.text.y = element_text(color="black"),
          panel.grid.major.x = element_blank()) + labs(fill = "")+
                          geom_text(size = 4, position = position_stack(vjust = 0.5)) + theme(axis.text = element_text(size = 15)) +theme(axis.title = element_text(size = 15))
  
    ggsave("GxE_and_Herit/GxE_FT_Meso.tiff", units="in", width=10, height=10, dpi=300, compression = 'lzw')
    ggsave("GxE_and_Herit/GxE_FT_Meso.pdf", units="in", width=10, height=10, dpi=300)
    
    # Commented 10/7/26, looks to have been for lima but vc_plot object does not exist
    #
    # ggplot(vc_plot, aes(x = trait, y = Percent, fill = component)) +
    # geom_bar(stat = "identity", width = 0.75) +
    # scale_fill_manual(
    #   values = c("E" = "#66c2a5", "G" = "#fc8d62", "G×E" = "#8da0cb", "R" = "#bdbdbd"),
    #   name = "Component") +
    # scale_x_discrete(limits=c("HSW", "Ash.db", "Protein.db", "Fat.db", "Starch.db"), 
    #                  labels=c("HSW", "Ash","Protein", "Fat","Starch"))+
    # labs(
    #   title = "Lima Variance Components Admixed",
    #   y = "% Variance Explained",
    #   x = "Trait") +
    # theme_minimal(base_size = 22) +
    # theme(
    #   axis.text.x = element_text(  color = "black"), #angle = 45, hjust=1
    #   axis.text.y = element_text(color="black"),
    #   panel.grid.major.x = element_blank()
    # )