#______________________________________________________________________________#
# MET variance component analysis — combined decomposition figure -------------
#
# APPROACH (adapted from CHD's version):
# All six analytical contexts fit the SAME model — random = ~ Env + Line.ID +
# Env:Line.ID + at(Env, reps):Rep, residual = dsum(... | Env) — and differ only
# in which rows of the data are supplied. Keeping Env as the dsum splitter in
# every context is the key: each residual section is always exactly one real
# site-year with its own intact Bed x Tier grid, so the STA-identified spatial
# structures apply everywhere without any level remapping. Selecting the right
# structures for a subset is then just an intersect() against the environments
# present.
#
# Within a location, Env is equivalent to Year; within a year, Env is
# equivalent to Location. Components are relabelled per context (E / Y / L)
# only at the plotting stage.
#
# CHD notes 9/27/26: This script was copied from ASReml_GxE_VarDecomp_
# StackBarChart_residual_error_for_across_ELY_FINAL_Lima_withHerit.R; changes
# made in this script, compared to that script: 1) switching back to Env rather
# than Field so fewer adaptations needed from cowpea manuscript; 2) added Meso,
# Andean, and admixed panels
#______________________________________________________________________________#
rm(list = ls())
library(asreml)
library(tidyverse)
library(cowplot)
asreml.options(drop.unused.levels = TRUE)

#______________________________________________________________________________#
# Load dataset ----------------------------------------------------------------
#______________________________________________________________________________#
setwd("C:/Users/chdiep/Documents/ResearchProjects/LIMA/SpatialModelFitting")
dat = read.csv("Flowering_Davis24_25_gxe.csv",head=T)
head(dat)
str(dat)
names(dat)
unique(dat$Field)
#Two colname changes for compatibility with this script base
colnames(dat)[which(colnames(dat) == "Field")] = "Env"
colnames(dat)[which(colnames(dat) == "LDP")] = "Line.ID"

dat$Pop = NA
dat[dat$Kgroup == "K1", "Pop"] = "Meso"
dat[dat$Kgroup == "K2", "Pop"] = "B"
dat[dat$Kgroup == "K3", "Pop"] = "Meso"
dat[dat$Kgroup == "K4", "Pop"] = "Meso"
dat[dat$Kgroup == "K5", "Pop"] = "Meso"
dat[dat$Kgroup == "K6", "Pop"] = "C"
dat[dat$Kgroup == "K7", "Pop"] = "C"
dat[dat$Kgroup == "K8", "Pop"] = "Meso"
dat[dat$Kgroup == "admixed", "Pop"] = "Admixed"

#______________________________________________________________________________#
# Grid-complete every environment ---------------------------------------------
# ar1(Bed):ar1(Tier) needs a complete rectangular Bed x Tier crossing within
# each environment. The STA/MET scripts expand-grid only three site-years, but
# dsum's implied-observation-count check shows others (e.g. CVARS) also have
# missing physical plot rows. Generalize: expand each Env over its own
# observed min:max range.
#______________________________________________________________________________#
dat$Bed  <- as.numeric(as.character(dat$Bed))
dat$Tier <- as.numeric(as.character(dat$Tier))

orig_col_order <- colnames(dat)

dat <- dat %>%
  group_by(Env) %>%
  group_modify(~ {
    grid <- expand.grid(
      Bed  = min(.x$Bed,  na.rm = TRUE):max(.x$Bed,  na.rm = TRUE),
      Tier = min(.x$Tier, na.rm = TRUE):max(.x$Tier, na.rm = TRUE)
    )
    merge(.x, grid, by = c("Bed", "Tier"), all.y = TRUE)
  }) %>%
  ungroup()

dat <- dat[, orig_col_order]   # merge() reorders by-columns to the front

dat[dat$Env == "UCD24", "Location"] = "Davis"
dat[dat$Env == "UCD25", "Location"] = "Davis"
dat[dat$Env == "UCD24", "Year"] = 2024
dat[dat$Env == "UCD25", "Year"] = 2025

#______________________________________________________________________________#
# Column indexing -------------------------------------------------------------
#______________________________________________________________________________#
non.trait.idx <- c(1:15,19:20) #CHD changed so that "dpp", "GDD", and "avg_plot_height" are the only trait cols; the latter two will be excluded in exclude_pattern
trait.idx     <- c(16:18) #same note as in line above
dat[, non.trait.idx] <- lapply(dat[, non.trait.idx], as.factor)
dat[, trait.idx]     <- lapply(dat[, trait.idx], as.numeric)

dat_Meso = dat
dat_Meso[which(!(dat_Meso$Pop == "Meso")),trait.idx] <- NA

exclude_pattern        <- "GDD|avg_plot_height"
trait.idx_for_analysis <- trait.idx[!grepl(exclude_pattern, colnames(dat)[trait.idx])]

replicated_envs <- c("UCD24","UCD25")

#______________________________________________________________________________#
# Load STA residual error structures ------------------------------------------
#______________________________________________________________________________#
# Important note: the file being loaded here contains both BIC and LRT best resid structs;
# the file suffix is only as it is because in that particular spatial model fitting run,
# LRT (in the overall spatial model fitting pipeline) was used to determine the best residual structure in spatial model fitting; but both the BIC- and
# LRT-determined best residual structures were stored in this file, which is from a run that had used BIC.
residual_info <- read.delim(
  "Model_Performance_Comparison_BestResErrorStrucs_compiled_2026-09-20__FT_usingBIC.txt",
  sep = "\t", header = TRUE
)

residual_info$Best.Res.Er.Str.LRT <- gsub(
  "~\\(id\\(Bed\\):ar1\\(Tier\\)\\)", "~id(Bed):ar1(Tier)",
  residual_info$Best.Res.Er.Str.LRT
)

structure_map <- c(
  "~ar1(Bed):ar1(Tier)" = 1,
  "~id(Bed):ar1(Tier)"  = 2,
  "~ar1(Bed):id(Tier)"  = 3,
  "~id(Bed):id(Tier)"   = 4,
  "~id(units)"          = 4
)

# trait -> list of 4 character vectors of Env names, one per structure index
ResStrucs_ov <- residual_info %>%
  dplyr::select(Env, Trait, Best.Res.Er.Str.LRT) %>%
  distinct() %>%
  mutate(StructureIndex = structure_map[Best.Res.Er.Str.LRT]) %>%
  group_by(Trait, StructureIndex) %>%
  summarise(Envs = list(Env), .groups = "drop") %>%
  complete(Trait, StructureIndex = 1:4, fill = list(Envs = list(character()))) %>%
  arrange(Trait, StructureIndex) %>%
  group_by(Trait) %>%
  summarise(StructList = list(map(1:4, ~ Envs[StructureIndex == .x][[1]])),
            .groups = "drop") %>%
  deframe()

#______________________________________________________________________________#
# Trait keep list and short labels ---------------------------------------------
#______________________________________________________________________________#
trait_labels <- c(
  "dpp" = "dpp"
)

#______________________________________________________________________________#
# Core model fitter -----------------------------------------------------------
# Always splits dsum by Env; caller supplies the data subset.
#______________________________________________________________________________#
fit_vc_model <- function(dat_sub, this.trait) {
  
  dat_sub      <- droplevels(dat_sub)
  dat_sub$y    <- dat_sub[[this.trait]]
  envs_present <- as.character(unique(dat_sub$Env))
  
  # dsum needs each Env contiguous and Bed/Tier ordered within it
  dat_sub <- dat_sub[order(dat_sub$Env, dat_sub$Bed, dat_sub$Tier), ]
  
  # --- Residual structures: just intersect against envs present -------------
  strucs <- if (this.trait %in% names(ResStrucs_ov)) {
    lapply(ResStrucs_ov[[this.trait]], intersect, envs_present)
  } else NULL
  
  if (!is.null(strucs)) {
    # Every Env must land in exactly one section. Environments the STA
    # skipped would otherwise be unassigned and trigger a cryptic 'names' attribute error.
    orphans <- setdiff(envs_present, unlist(strucs))
    if (length(orphans) > 0) {
      message("    no STA structure for: ", paste(orphans, collapse = ", "),
              " -> id(Bed):id(Tier)")
      strucs[[4]] <- c(strucs[[4]], orphans)
    }
    if (all(lengths(strucs) == 0)) strucs <- NULL
  }
  
  if (!is.null(strucs)) {
    # deparse() defaults to width.cutoff = 60 and silently splits long
    # objects across multiple lines, which paste0() then vectorizes into a
    # malformed formula. Force wide output and collapse to one string.
    resid_form <- as.formula(paste0(
      "~ dsum(~ ar1(Bed):ar1(Tier) + id(Bed):ar1(Tier) + ",
      "ar1(Bed):id(Tier) + id(Bed):id(Tier) | Env, levels = ",
      paste(deparse(strucs, width.cutoff = 500L), collapse = " "), ")"
    ))
  } else {
    # at() is not permitted in a residual formula; heterogeneous residual
    # variance by group must also go through dsum().
    resid_form <- ~ dsum(~ id(units) | Env)
  }
  
  # --- Random formula -------------------------------------------------------
  rep_envs <- intersect(replicated_envs, envs_present)
  random_form <- if (length(rep_envs) > 0) {
    as.formula(paste0(
      "~ Env + Line.ID + Env:Line.ID + at(Env, c(",
      paste0("'", rep_envs, "'", collapse = ", "), ")):Rep"
    ))
  } else {
    ~ Env + Line.ID + Env:Line.ID
  }
  
  mod <- tryCatch(
    asreml(fixed = y ~ 1, random = random_form, residual = resid_form,
           data = dat_sub, na.action = na.method(y = "include", x = "include"),
           workspace = "6gb", maxit = 100),
    error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL }
  )
  if (is.null(mod)) return(NULL)
  mod <- update(mod)
  
  vc <- summary(mod)$varcomp
  data.frame(component = rownames(vc),
             variance  = vc[, "component"],
             bound     = vc[, "bound"],
             stringsAsFactors = FALSE)
}

#______________________________________________________________________________#
# Tidy raw varcomp output into the five plotted components --------------------
# Rep terms are summed; residual sections are AVERAGED (not summed) so that R
# represents a typical single-environment error. Summing would count residual
# once per environment while G is counted once, inflating R and shrinking the
# apparent contribution of G, E and GxE.
#______________________________________________________________________________#
tidy_vc <- function(vc_df) {
  
  # ar1 sections emit correlation parameters (bounded -1..1) alongside
  # variances. They are not variance components and must not enter the mean.
  vc_df <- vc_df[!grepl("cor$|!cor", vc_df$component), ]
  
  main    <- c("Env", "Line.ID", "Env:Line.ID")
  is_rep  <- grepl("^at\\(Env.*\\):Rep", vc_df$component)
  is_main <- vc_df$component %in% main
  is_res  <- !is_main & !is_rep
  
  out <- vc_df[is_main, c("component", "variance", "bound")]
  
  if (any(is_rep)) {
    out <- rbind(out, data.frame(component = "Env:Rep",
                                 variance  = sum(vc_df$variance[is_rep]),
                                 bound     = NA_character_))
  }
  if (any(is_res)) {
    out <- rbind(out, data.frame(component = "units!R",
                                 variance  = mean(vc_df$variance[is_res], na.rm = TRUE),
                                 bound     = NA_character_))
  }
  out
}

#______________________________________________________________________________#
# Analysis contexts -----------------------------------------------------------
# label = what the Env-based terms are called in this context's panel.
#______________________________________________________________________________#
contexts <- list(
  #list(panel = "Across All", label = "E", min_envs = 1, #CHD changed min env to 1 (though will likely replicate Davis panel)
  #     subset = function(d) d), #CHD commented Across all as within Davis is sufficient for FT
  ## CHD commented next two lines, not conducting within-year across-loc analysis for limas
  # list(panel = "2022",       label = "L", min_envs = 2,
  #      subset = function(d) filter(d, as.character(Year) == "2022")),
  # list(panel = "2023",       label = "L", min_envs = 2,
  #      subset = function(d) filter(d, as.character(Year) == "2023")),
  ## CHD added next three lists (Meso, Andean [k=2 only, not including k6 and k7], and admixed)
  list(panel = "Meso",       label = "E", min_envs = 2,
        subset = function(d) filter(d, as.character(Pop) == "Meso")),
  #list(panel = "Andean",       label = "A", min_envs = 2, #CHD commented, sample size small
  #      subset = function(d) filter(d, as.character(Pop) == "B")),
  #list(panel = "Admixed",       label = "Adm", min_envs = 2, #CHD commented, sample size small
  #     subset = function(d) filter(d, as.character(Pop) == "Admixed")),
  list(panel = "Davis",      label = "Y", min_envs = 1, #CHD changed min_envs to 1
        subset = function(d) filter(d, Location == "Davis")) #CHD removed comma since the two lines below are now empty
  ## Commented next two lists, not applicable for FT
  # list(panel = "Central.Ferry",    label = "Y", min_envs = 2,
  #      subset = function(d) filter(d, Location == "Central.Ferry")),
  # list(panel = "Coachella",    label = "Y", min_envs = 2,
  #      subset = function(d) filter(d, Location == "Coachella"))
)

#______________________________________________________________________________#
# Run all contexts x traits ---------------------------------------------------
#______________________________________________________________________________#
vc_results   <- list()
boundary_log <- list()

for (ctx in contexts) {
  
  if(ctx$panel == "Meso"){
    dat_ctx <- dat_Meso %>% droplevels()
  }else{
  dat_ctx <- ctx$subset(dat) %>% droplevels()
  }
  
  for (col.idx in trait.idx_for_analysis) {
    this.trait <- colnames(dat_ctx)[col.idx] #CHD added _ctx suffix
    key <- paste(ctx$panel, this.trait, sep = "__")
    cat("\n[", ctx$panel, "] ", this.trait, sep = "")
    
    envs_with_data <- dat_ctx %>%
      group_by(Env) %>% 
      summarise(n = sum(!is.na(.data[[this.trait]])), .groups = "drop") %>%
      filter(n > 0) %>% pull(Env) %>% as.character()
    
    if (length(envs_with_data) < ctx$min_envs) { message("  skipping"); next }
    
    dat_trait <- dat_ctx %>% filter(as.character(Env) %in% envs_with_data)
    
    vc_raw <- fit_vc_model(dat_trait, this.trait)
    if (is.null(vc_raw)) next
    
    bad <- vc_raw %>% filter(bound %in% c("B", "F"))
    if (nrow(bad) > 0)
      boundary_log[[key]] <- bad %>% mutate(panel = ctx$panel, trait = this.trait)
    
    vc_results[[key]] <- tidy_vc(vc_raw) %>%
      mutate(panel = ctx$panel,
             label = ctx$label,
             trait = this.trait,
             Percent = variance / sum(variance) * 100)
  }
}

#______________________________________________________________________________#
# Assemble, relabel components -------------------------------------------------
#______________________________________________________________________________#
vc_all <- bind_rows(vc_results) %>%
  filter(trait %in% names(trait_labels)) %>%
  mutate(
    trait     = unname(trait_labels[trait]),
    component = case_when(
      component == "Line.ID"     ~ "G",
      component == "Env:Line.ID" ~ paste0("G\u00d7", label),
      component == "Env"         ~ label,
      component == "Env:Rep"     ~ paste0(label, ":Rep"),
      component == "units!R"     ~ "R", #CHD note: it appears in line 235 resid from best error structure is being ascribed to what is called units, even if not using a units resid struct, so just a matter of nomenclature
      TRUE                       ~ component
    )
  )

setwd("GxE_and_Herit")

if (length(boundary_log) > 0) {
  write.csv(bind_rows(boundary_log), "VC_boundary_log_combined_LRT_FT.csv", row.names = FALSE)
  cat("\nBoundary log written:", nrow(bind_rows(boundary_log)), "parameter(s).\n")
} else cat("\nNo boundary parameters detected.\n")

write.csv(vc_all, "VC_components_all_contexts_LRT_FT.csv", row.names = FALSE)

#______________________________________________________________________________#
# Figure ----------------------------------------------------------------------
#______________________________________________________________________________#
vc_theme <- theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y        = element_text(size = 7),
    axis.title.y       = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 8, hjust = 0.5,
                                      margin = margin(b = 2))
  )

component_colors <- c(
  "G" = "#fc8d62",
  "G\u00d7Y" = "#8da0cb", "G\u00d7E" = "#8da0cb", "G\u00d7L" = "#8da0cb",
  "Y" = "#66c2a5", "E" = "#66c2a5", "L" = "#66c2a5",
  "Y:Rep" = "#e78ac3", "E:Rep" = "#e78ac3", "L:Rep" = "#e78ac3",
  "R" = "#bdbdbd"
)

global_trait_order <- c("dpp")

make_panel <- function(panel_name, title_str) {
  
  lab <- vc_all$label[vc_all$panel == panel_name][1]
  lv  <- c("G", paste0("G\u00d7", lab), lab, paste0(lab, ":Rep"), "R")
  
  df <- vc_all %>% filter(panel == panel_name)
  
  # Pad traits absent from this panel so the x-axis stays aligned
  missing <- setdiff(global_trait_order, unique(df$trait))
  if (length(missing) > 0) {
    df <- bind_rows(df, expand.grid(trait = missing, component = lv,
                                    stringsAsFactors = FALSE) %>%
                      mutate(Percent = 0, variance = 0, panel = panel_name))
  }
  
  # Factor levels applied AFTER bind_rows: bind_rows coerces a factor +
  # character conflict to character, silently dropping level order and
  # leaving ggplot to stack alphabetically.
  df$component <- factor(df$component, levels = lv)
  df$trait     <- factor(df$trait,     levels = global_trait_order)
  
  ggplot(df, aes(x = trait, y = Percent, fill = component)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_text(aes(label = ifelse(Percent >= 5, round(Percent, 0), "")),
              position = position_stack(vjust = 0.5),
              size = 2.0, colour = "white", fontface = "bold") +
    scale_fill_manual(values = component_colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(y = "% Variance", x = NULL, title = title_str) +
    vc_theme
}

#p_all     <- make_panel("Across All", "Across all environments") #CHD commented
# p_2022    <- make_panel("2022",       "2022 (across locations)") #CHD commented, not conducting within-year analysis for limas
# p_2023    <- make_panel("2023",       "2023 (across locations)") #CHD commented, not conducting within-year analysis for limas
p_Meso    <- make_panel("Meso",       "Mesoamerican (n=248)")
#p_Andean    <- make_panel("Andean",       "Andean (n=47)") #CHD commented, sample size small
#p_Admixed    <- make_panel("Admixed",       "Admixed (n=56)") #CHD commented, sample size small
p_davis   <- make_panel("Davis",      "Davis, CA (across years)")
#p_cf <- make_panel("Central.Ferry",    "Central Ferry, WA (across years)") #CHD commented, not applicable for FT
#p_thermal <- make_panel("Coachella",    "CVARS, CA (across years)") #CHD commented, not applicable for FT

legend_df <- data.frame(
  component = factor(c("G", "G\u00d7E/Y/L", "E/Y/L", "E/Y/L:Rep", "R"),
                     levels = c("G", "G\u00d7E/Y/L", "E/Y/L", "E/Y/L:Rep", "R")),
  x = 1:5, y = 1
)

shared_legend <- get_legend(
  ggplot(legend_df, aes(x, y, fill = component)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("G" = "#fc8d62", "G\u00d7E/Y/L" = "#8da0cb",
                 "E/Y/L" = "#66c2a5", "E/Y/L:Rep" = "#e78ac3", "R" = "#bdbdbd"),
      name = "Component") +
    theme(legend.position = "right",
          legend.title    = element_text(face = "bold", size = 8),
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.35, "cm"))
)

left_col  <- plot_grid(p_davis, ncol = 1, align = "v", axis = "lr", #rm'ed p_all, p_Andean, p_Admixed
                       labels = c("a"), label_size = 8, #rm'ed "b","c","d"
                       label_fontface = "bold") #CHD removed p_2022, p_2023
right_col <- plot_grid(p_Meso, ncol = 1, align = "v", #rm'ed p_cf, p_thermal
                       axis = "lr", labels = c("b"), label_size = 8, #changed to b rather than e through g
                       label_fontface = "bold")

final_fig <- plot_grid(plot_grid(left_col, right_col, ncol = 2),
                       shared_legend, ncol = 2, rel_widths = c(1, 0.14))

ggsave("Fig_VC_Combined_2panel.png", final_fig, #CHD changed 6 to 2
       width = 174, height = 174, units = "mm", dpi = 300)

## Here to end was pasted and adapted from "Heritability_w_GxE_patch.R" from "ExactScriptsFromCowpea"
# =============================================================================
# 5. PER-TRAIT: OUTLIER DETECTION -> HERITABILITY
# =============================================================================
herit.across    <- NULL
OutliersCombDFs <- NULL

#' Build a dsum residual formula for the environments actually present.
#' Because Env is always the splitter, selecting structures for a subset is
#' just an intersect — no level remapping is ever needed.
build_resid_formula <- function(this.trait, envs_present) {
  
  strucs <- if (this.trait %in% names(ResStrucs_ov)) {
    lapply(ResStrucs_ov[[this.trait]], intersect, envs_present)
  } else NULL
  
  if (!is.null(strucs)) {
    # Every Env must be assigned to exactly one section. Environments the STA
    # skipped would otherwise be
    # unassigned and trigger a cryptic 'names' attribute error.
    orphans <- setdiff(envs_present, unlist(strucs))
    if (length(orphans) > 0) strucs[[4]] <- c(strucs[[4]], orphans)
    if (all(lengths(strucs) == 0)) strucs <- NULL
  }
  
  if (is.null(strucs)) return(~ dsum(~ id(units) | Env))
  
  # deparse() defaults to width.cutoff = 60 and silently splits long objects
  # over multiple lines, which paste0() then vectorises into a malformed
  # formula. Force wide output and collapse to a single string.
  as.formula(paste0(
    "~ dsum(~ ar1(Bed):ar1(Tier) + id(Bed):ar1(Tier) + ",
    "ar1(Bed):id(Tier) + id(Bed):id(Tier) | Env, levels = ",
    paste(deparse(strucs, width.cutoff = 500L), collapse = " "), ")"
  ))
}

#this 'contexts' for loop was introduced 10/20/26 (around the trait loop) to estim herit in both all genos and Mesoamerican genos
for (ctx2 in contexts[1:2]) {
  
  # added this if/else statement 10/20/26 to estim herit in both all genos and Mesoamerican genos
  if(ctx2$panel == "Meso"){
    dat_ctx2 <- dat_Meso %>% droplevels()
  }else{
    dat_ctx2 <- ctx2$subset(dat) %>% droplevels()
  }
for (col.idy in trait.idx_for_analysis) {
  
  this.trait <- colnames(dat_ctx2)[col.idy] #changed from dat to dat_ctx2 10/20/26 to estim herit in all genos and Mesoamerican genos
  cat("\n\n===== ", this.trait, " =====\n", sep = "")
  
  envs_present <- as.character(unique(dat_ctx2$Env)) #changed from dat to dat_ctx2 10/20/26 to estim herit in all genos and Mesoamerican genos
  rep_envs   <- intersect(replicated_envs, envs_present)
  resid_form <- build_resid_formula(this.trait, envs_present)
  
  random_form <- if (length(rep_envs) > 0) {
    as.formula(paste0("~ Env + Line.ID + Env:Line.ID + at(Env, c(",
                      paste0("'", rep_envs, "'", collapse = ", "), ")):Rep"))
  } else {
    ~ Env + Line.ID + Env:Line.ID
  }
  
  # ---------------------------------------------------------------------------
  # 5b. Heritability model
  # ---------------------------------------------------------------------------
  MET_herit <- tryCatch({
    m <- asreml(fixed = as.formula(paste(this.trait, "~ 1")),
                random = random_form, residual = resid_form,
                na.action = na.method(y = "include", x = "include"),
                data = dat_ctx2, workspace = "6.0gb", maxit = 100) # CHD changed from dat.clean, no MET outliers were detected; #changed from dat to dat.ctx2 10/20/26 to estim herit in all genos and Mesoamerican genos
    update(m)
  }, error = function(e) {
    message("  H2 model failed: ", conditionMessage(e)); NULL
  })
  
  if (is.null(MET_herit)) next
  
  # --- Variance components ---------------------------------------------------
  vc      <- summary(MET_herit)$varcomp
  vc$term <- rownames(vc)
  
  # ar1 sections emit correlation parameters (bounded -1..1) alongside
  # variances; they are not variance components and must be excluded.
  is_cor  <- grepl("cor$|!cor", vc$term)
  is_rep  <- grepl("^at\\(Env.*\\):Rep", vc$term)
  is_main <- vc$term %in% c("Env", "Line.ID", "Env:Line.ID")
  is_res  <- !is_cor & !is_rep & !is_main
  
  getvc <- function(nm) if (nm %in% vc$term) vc$component[vc$term == nm] else NA_real_
  
  gvar  <- getvc("Line.ID")
  gevar <- getvc("Env:Line.ID")
  # Average, not sum: any single observation is subject to one environment's
  # residual, so the mean represents a typical within-environment error.
  rvar  <- if (any(is_res)) mean(vc$component[is_res], na.rm = TRUE) else NA_real_
  
  # --- Effective replication -------------------------------------------------
  rep_tab <- dat_ctx2 %>% # CHD changed from dat.clean, no MET outliers were detected; #changed from dat to dat.ctx2 10/20/26 to estim herit in all genos and Mesoamerican genos
    filter(!is.na(.data[[this.trait]])) %>%
    group_by(Env, Line.ID) %>% summarise(n = n(), .groups = "drop") %>%
    group_by(Env) %>% summarise(mean_reps = mean(n), .groups = "drop")
  
  n_e <- nrow(rep_tab)
  n_r <- mean(rep_tab$mean_reps)   # approximate: CVARS is unreplicated
  
  # --- Cullis H2 -------------------------------------------------------------
  avsed <- tryCatch(
    predict(MET_herit, classify = "Line.ID",
            workspace = "5gb", pworkspace = "5gb")$avsed,
    error = function(e) NA_real_
  )
  
  H2_Cullis <- if (!is.na(gvar) && !is.na(avsed) && gvar > 0) {
    1 - (avsed^2 / (2 * gvar))
  } else NA_real_
  
  # --- Classical line-mean H2 ------------------------------------------------
  H2_classical <- if (!is.na(gvar) && gvar > 0) {
    gvar / (gvar +
              ifelse(is.na(gevar), 0, gevar) / n_e +
              ifelse(is.na(rvar),  0, rvar)  / (n_e * n_r))
  } else NA_real_
  
  herit.across <- rbind(herit.across, data.frame(
    panel        = ctx2$panel, #added 10/20/26 to estimate herit for both all genos and Mesoamerican genos
    Trait        = this.trait,
    H2_Cullis    = H2_Cullis,
    H2_classical = H2_classical,
    Vg           = gvar,
    Vge          = gevar,
    Vresid_mean  = rvar,
    n_envs       = n_e,
    mean_reps    = round(n_r, 2),
    avsed        = avsed,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  H2 Cullis = %.3f | H2 classical = %.3f  (%d envs, %.2f reps)\n",
              H2_Cullis, H2_classical, n_e, n_r))
}
}#end for ctx2 in contexts; added 10/20/26 to estimate herit for both all genos and Mesoamerican genos

# =============================================================================
# 6. WRITE RESULTS
# =============================================================================
write.table(
  herit.across,
  paste0("All_Envs_H2_Cullis_and_Classical_LRT_FT_AllGenosAndMeso.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n\nDone. Traits with estimates:", nrow(herit.across), "\n")
print(herit.across[, c("Trait", "H2_Cullis", "H2_classical", "n_envs")])
