# Functions to compute Generalized Additive Models

# Load libraries
library(tidyverse)
library(mgcv)
library(patchwork)
library(ggnewscale)
library(scales)

# Load palettes
source("scripts/00-palettes.R")


#### Compute Generalized Additive Models (GAMs) ####
# Define function to format response and explanatory variables for constructing GAMs
formatGAMVars <- function(lakeid, response_data, response_var) {
  # Format univariate response variable
  response <- response_data %>%
    filter(lake_id == lakeid)
  
  # Join all environmental data
  env_all <- watertemp_byyear_interpolated %>%
    full_join(sdd_byyear_interpolated, by = c("lake_id", "year")) %>%
    full_join(chemistry_byyear_interpolated, by = c("lake_id", "year")) %>%
    full_join(climate_eccc_byyear %>%
                rename(precip_eccc = precip), by = c("lake_id", "year")) %>%
    full_join(climate_ela_byyear %>%
                rename(precip_ela = precip), by = c("lake_id", "year")) %>%
    filter(lake_id == lakeid) %>%
    arrange(year)
  
  # Join all response and explanatory variables
  variables <- response %>%
    full_join(env_all, by = c("lake_id", "year")) %>%
    filter(!is.na(eval(as.name(paste(response_var))))) %>%
    arrange(year)
  
  # List all response and explanatory variables
  gam_variables <- list(env_all = env_all,
                        variables = variables)
  
  return(gam_variables)
}

# Define function to plot relationships between response and explanatory variables
plotVars <- function(lakeid, response_data, response_var) {
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  # Plot relationships between response and explanatory variables
  variables %>%
    pivot_longer(all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year"))]),
                 names_to = "var", values_to = "value", values_drop_na = TRUE) %>%
    ggplot() +
    facet_wrap(~var, scales = "free_x", strip.position = "bottom") +
    geom_point(aes(x = value, y = eval(as.name(paste(response_var))),
                   colour = lakeid),
               alpha = 0.7) +
    labs(y = as.name(paste(response_var))) +
    ggtitle(lakeid) +
    scale_colour_manual(values = palette_lake, guide = "none") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.grid = element_blank(),
          panel.background = element_blank())
}


#### Fit a GAM for each explanatory variable ####
# Define function to fit linear models and GAMs for each explanatory variable
modelIndVar <- function(lakeid, response_data, response_var) {
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  # Specify which variables to discard
  vars_discard <- c("minairtemp", "maxairtemp", "precip_ela",
                    "epilimnion_ph", "profundal_ph",
                    "epilimnion_oxygen", "profundal_oxygen")
  
  # Create vector of all environmental variables
  vars_all <- names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_discard))]
  
  models_ind <- tibble(NULL)
  for (i in 1:length(vars_all)) {
    dataset_tmp <- unique(response_data$dataset)
    
    # Extract a single explanatory variable
    var_tmp <- vars_all[i]
    
    variables_tmp <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", var_tmp))])) %>%
      drop_na()
    
    # Fit linear model using a single explanatory variable
    lm_tmp <- gam(as.formula(paste(response_var, "~", var_tmp)),
                  data = variables_tmp, method = "REML", select = TRUE)
    
    if (exists("lm_tmp")) {
      lm_table_tmp <- summary(lm_tmp)$p.table %>%
        as_tibble(rownames = "explanatory_var") %>%
        rename(lm_estimate = Estimate,
               lm_std_error = `Std. Error`,
               lm_t_value = `t value`,
               lm_p_value = `Pr(>|t|)`) %>%
        filter(explanatory_var == var_tmp) %>%
        mutate(lm_r2 = summary(lm_tmp)$r.sq,
               lm_dev_expl = summary(lm_tmp)$dev.expl,
               lm_aic = AIC(lm_tmp))
    }
    
    # Fit GAM using a single explanatory variable
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", gsub("$", ")", gsub("^", "s(", var_tmp)))),
                       data = variables_tmp, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      smooth_table_tmp <- summary(gam_tmp)$s.table %>%
        as_tibble(rownames = "explanatory_var") %>%
        mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
        rename(gam_edf = edf,
               gam_red_df = Ref.df,
               gam_f = `F`,
               gam_p_value = `p-value`) %>%
        mutate(gam_r2 = summary(gam_tmp)$r.sq,
               gam_dev_expl = summary(gam_tmp)$dev.expl,
               gam_aic = AIC(gam_tmp))
    }
    
    # Combine linear model and GAM outputs
    if (exists("lm_table_tmp") & exists("smooth_table_tmp")) {
      models_tmp <- lm_table_tmp %>%
        left_join(smooth_table_tmp, by = "explanatory_var") %>%
        mutate(dataset = dataset_tmp,
               response_var = response_var,
               lake_id = lakeid) %>%
        relocate(dataset, 1) %>%
        relocate(lake_id, .after = dataset) %>%
        relocate(response_var, .before = explanatory_var)
    } else if (exists("lm_table_tmp") & !exists("smooth_table_tmp")) {
      models_tmp <- lm_table_tmp %>%
        mutate(dataset = dataset_tmp,
               response_var = response_var,
               lake_id = lakeid) %>%
        relocate(dataset, 1) %>%
        relocate(lake_id, .after = dataset) %>%
        relocate(response_var, .before = explanatory_var)
    } else if (!exists("lm_table_tmp") & exists("smooth_table_tmp")) {
      models_tmp <- smooth_table_tmp %>%
        mutate(dataset = dataset_tmp,
               response_var = response_var,
               lake_id = lakeid) %>%
        relocate(dataset, 1) %>%
        relocate(lake_id, .after = dataset) %>%
        relocate(response_var, .before = explanatory_var)
    }
    
    models_ind <- models_ind %>%
      bind_rows(models_tmp)
    
    rm(list = ls(pattern = "_tmp"))
  }
  
  return(models_ind)
}

# Define function to fit linear models and GAMs for each explanatory variable across all lakes
modelIndVarAll <- function(response_data, response_var) {
  # Compute models for all lakes
  models_ind_all <- bind_rows(modelIndVar("L226N", response_data, response_var),
                              modelIndVar("L226S", response_data, response_var),
                              modelIndVar("L227", response_data, response_var),
                              modelIndVar("L224", response_data, response_var),
                              modelIndVar("L373", response_data, response_var))
  
  return(models_ind_all)
}


#### Fit a GAM for each category of explanatory variable ####
# Categorize variables by type (climate or limnology) and water column stratum
explanatory_vars_categories <- tibble(explanatory_category = c("climate",
                                                               "limno_epilimnion",
                                                               "limno_profundal"),
                                      vars = c(toString(c("meanairtemp",
                                                          "precip_eccc",
                                                          "epilimnion_temp",
                                                          "hypolimnion_temp")),
                                               toString(c("sdd",
                                                          "epilimnion_tn",
                                                          "epilimnion_tp"#,
                                                          #"epilimnion_ph",
                                                          #"epilimnion_oxygen"
                                                          )),
                                               toString(c("profundal_tn",
                                                          "profundal_tp"#,
                                                          #"profundal_ph",
                                                          #"profundal_oxygen"
                                                          ))))

explanatory_vars <- explanatory_vars_categories %>%
  mutate(var = strsplit(vars, ", ")) %>%
  unnest(var) %>%
  pull(var)

# Define function to fit GAMs with linear or smooth terms from the same environmental category
computeGAMEnv <- function(lakeid, response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  # Evaluate which variables should be modeled with linear or smooth terms
  model_eval <- modelIndVarAll(response_data, response_var) %>%
    mutate(model_winner = case_when(!is.na(lm_aic) & is.na(gam_aic) ~ "lm",
                                    is.na(lm_aic) & !is.na(gam_aic) ~ "gam",
                                    
                                    lm_aic < gam_aic & abs(lm_aic - gam_aic) >= 4 ~ "lm",
                                    gam_aic < lm_aic & abs(gam_aic - lm_aic) >= 4 ~ "gam",
                                    
                                    gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                    gam_edf >= 2 ~ "gam",
                                    
                                    lm_p_value < 0.05 & gam_p_value >= 0.05 ~ "lm",
                                    lm_p_value >= 0.05 & gam_p_value < 0.05 ~ "gam",
                                    
                                    TRUE ~ "gam")) %>%
    filter(lake_id == lakeid)
  
  # Retain variables with at least one statistically significant model
  model_eval_signif <- model_eval %>%
    filter(lm_p_value < 0.05 | gam_p_value < 0.05)
  
  # Define category-specific variables
  explanatory_categories <- unique(explanatory_vars_categories$explanatory_category)
  
  gam_summaries <- tibble(NULL)
  for (i in 1:length(explanatory_categories)) {
    category_tmp <- explanatory_categories[i]
    
    category_vars_tmp <- explanatory_vars_categories$vars[which(explanatory_vars_categories$explanatory_category == category_tmp)]
    category_vars_tmp <- unlist(strsplit(category_vars_tmp, ", "))
    
    category_variables <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", category_vars_tmp))])) %>%
      drop_na()
    
    model_signif_category <- model_eval_signif %>%
      filter(explanatory_var %in% category_vars_tmp)
    
    if (nrow(model_signif_category) > 0) {
      # Create vector of explanatory variables based on linear or smoothed terms
      terms_category <- model_signif_category %>%
        mutate(term = case_when(model_winner == "lm" ~ explanatory_var,
                                model_winner == "gam" ~ gsub("$", ")", gsub("^", "s(", explanatory_var)))) %>%
        pull(term)
      
      # Fit GAM using linear or smoothed terms
      gam_category <- gam(as.formula(paste(response_var, "~", paste(terms_category, collapse = " + "))),
                          data = category_variables, method = "REML", select = TRUE)
      
      # Extract variables with significant terms
      if (all(!grepl("s\\(", terms_category))) {  # Linear terms only
        gam_category_vars_signif_lm <- summary(gam_category)$p.table %>%
          as_tibble(rownames = "explanatory_var") %>%
          rename(lm_estimate = Estimate,
                 lm_std_error = `Std. Error`,
                 lm_t_value = `t value`,
                 lm_p_value = `Pr(>|t|)`) %>%
          filter(explanatory_var != "(Intercept)") %>%
          filter(lm_p_value < 0.05) %>%
          pull(explanatory_var)
        
        gam_category_vars_signif <- gam_category_vars_signif_lm
      } else if (all(grepl("s\\(", terms_category))) {  # Smooth terms only
        gam_category_vars_signif_smooth <- summary(gam_category)$s.table %>%
          as_tibble(rownames = "explanatory_var") %>%
          mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
          rename(gam_edf = edf,
                 gam_red_df = Ref.df,
                 gam_f = `F`,
                 gam_p_value = `p-value`) %>%
          filter(gam_p_value < 0.05) %>%
          pull(explanatory_var)
        
        gam_category_vars_signif <- gam_category_vars_signif_smooth
      } else {  # Combination of linear and smooth terms
        gam_category_vars_signif_lm <- summary(gam_category)$p.table %>%
          as_tibble(rownames = "explanatory_var") %>%
          rename(lm_estimate = Estimate,
                 lm_std_error = `Std. Error`,
                 lm_t_value = `t value`,
                 lm_p_value = `Pr(>|t|)`) %>%
          filter(explanatory_var != "(Intercept)") %>%
          filter(lm_p_value < 0.05) %>%
          pull(explanatory_var)
        
        gam_category_vars_signif_smooth <- summary(gam_category)$s.table %>%
          as_tibble(rownames = "explanatory_var") %>%
          mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
          rename(gam_edf = edf,
                 gam_red_df = Ref.df,
                 gam_f = `F`,
                 gam_p_value = `p-value`) %>%
          filter(gam_p_value < 0.05) %>%
          pull(explanatory_var)
        
        gam_category_vars_signif <- c(gam_category_vars_signif_lm, gam_category_vars_signif_smooth)
      }
      
      # Fit GAM using only variables that were significant in full GAM
      terms_signif_category <- model_signif_category %>%
        filter(explanatory_var %in% gam_category_vars_signif) %>%
        mutate(term = case_when(model_winner == "lm" ~ explanatory_var,
                                model_winner == "gam" ~ gsub("$", ")", gsub("^", "s(", explanatory_var)))) %>%
        pull(term)
      
      vars_signif_category <- str_remove_all(str_remove_all(terms_signif_category, "^s\\("), "\\)")
      
      category_variables_signif <- variables %>%
        dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_signif_category))])) %>%
        drop_na()
      
      if (length(terms_signif_category) > 0) {
        gam_signif_category <- gam(as.formula(paste(response_var, "~", paste(terms_signif_category, collapse = " + "))),
                                   data = category_variables_signif, method = "REML", select = TRUE)
        # Summarize all GAMs
        gam_summary_tmp <- tibble(dataset = dataset_tmp,
                                  lake_id = lakeid,
                                  response_var = response_var,
                                  explanatory_category = category_tmp,
                                  
                                  gam_full_vars = str_remove_all(str_remove_all(toString(terms_category), "s\\("), "\\)"),
                                  gam_full_terms = toString(terms_category),
                                  gam_full_r2 = summary(gam_category)$r.sq,
                                  gam_full_dev_expl = summary(gam_category)$dev.expl,
                                  gam_full_aic = AIC(gam_category),
                                  
                                  gam_signif_vars = str_remove_all(str_remove_all(toString(terms_signif_category), "s\\("), "\\)"),
                                  gam_signif_terms = toString(terms_signif_category),
                                  gam_signif_r2 = summary(gam_signif_category)$r.sq,
                                  gam_signif_dev_expl = summary(gam_signif_category)$dev.expl,
                                  gam_signif_aic = AIC(gam_signif_category))
      } else {
        # Summarize full GAM only
        gam_summary_tmp <- tibble(dataset = dataset_tmp,
                                  lake_id = lakeid,
                                  response_var = response_var,
                                  explanatory_category = category_tmp,
                                  
                                  gam_full_vars = str_remove_all(str_remove_all(toString(terms_category), "s\\("), "\\)"),
                                  gam_full_terms = toString(terms_category),
                                  gam_full_r2 = summary(gam_category)$r.sq,
                                  gam_full_dev_expl = summary(gam_category)$dev.expl,
                                  gam_full_aic = AIC(gam_category))
      }
      
      gam_summaries <- gam_summaries %>%
        bind_rows(gam_summary_tmp)
    }
  }
  
  return(gam_summaries)
}

# Define function to fit GAMs with linear or smooth terms from the same environmental category across all lakes
computeGAMEnvAll <- function(response_data, response_var) {
  gam_summaries_all <- bind_rows(computeGAMEnv("L226N", response_data, response_var),
                                 computeGAMEnv("L226S", response_data, response_var),
                                 computeGAMEnv("L227", response_data, response_var),
                                 computeGAMEnv("L224", response_data, response_var),
                                 computeGAMEnv("L373", response_data, response_var))
  
  return(gam_summaries_all)
}

# Define function to plot deviance explained of GAM based on explanatory variables from the same category
plotGAMEnvAll <- function(response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  gam_summaries_tmp <- computeGAMEnvAll(response_data, response_var)
  
  barplot_tmp <- gam_summaries_tmp %>%
    ggplot() +
    facet_wrap(~explanatory_category, nrow = 1, strip.position = "bottom") +
    geom_bar(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                 y = gam_full_dev_expl * 100,
                 fill = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373"))),
             position = position_dodge(preserve = "single"), stat = "identity") +
    scale_fill_manual(values = palette_lake, name = "Lake") +
    labs(y = "Deviance explained (%)") +
    ylim(c(0, 100)) +
    ggtitle(str_c(dataset_tmp, response_var, sep = " ")) +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom")
  
  heatmap_tmp <- gam_summaries_tmp %>%
    mutate(gam_full_var = strsplit(gam_full_vars, ", ")) %>%
    unnest(gam_full_var) %>%
    mutate(signif = case_when(str_detect(gam_signif_vars, gam_full_var) ~ "",
                              TRUE ~ "NS")) %>%
    ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               y = factor(gam_full_var, levels = rev(explanatory_vars)))) +
    geom_tile(fill = "#E0E0E0", colour = "white") +
    geom_text(aes(label = signif)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw() %+replace%
    theme(axis.title = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank())
  
  barplot_tmp + heatmap_tmp
}


#### Fit a GAM for explanatory variables across different categories ####
# Define function to fit GAMs with explanatory variable across environmental categories
computeGAM <- function(lakeid, response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  gam_summaries_tmp <- computeGAMEnv(lakeid, response_data, response_var)
  
  if ("gam_signif_terms" %in% colnames(gam_summaries_tmp)) {
    # Create vector of explanatory variables across environmental categories
    terms_env <- gam_summaries_tmp %>%
      pull(gam_signif_terms)
    terms_env <- unlist(strsplit(terms_env, ", "))
    terms_env <- terms_env[!is.na(terms_env)]
    
    vars_env <- str_remove_all(str_remove_all(terms_env, "^s\\("), "\\)")
    
    env_variables <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_env))])) %>%
      drop_na()
    
    # Fit GAM using explanatory variables across environmental categories
    gam_env <- gam(as.formula(paste(response_var, "~", paste(terms_env, collapse = " + "))),
                   data = env_variables, method = "REML", select = TRUE)
    
    # Extract significant terms
    if (all(!grepl("s\\(", terms_env))) {  # Linear terms only
      gam_env_terms_signif_lm <- summary(gam_env)$p.table %>%
        as_tibble(rownames = "term") %>%
        rename(lm_estimate = Estimate,
               lm_std_error = `Std. Error`,
               lm_t_value = `t value`,
               lm_p_value = `Pr(>|t|)`) %>%
        filter(term != "(Intercept)") %>%
        filter(lm_p_value < 0.05) %>%
        pull(term)
      
      gam_env_terms_signif <- gam_env_terms_signif_lm
    } else if (all(grepl("s\\(", terms_env))) {  # Smooth terms only
      gam_env_terms_signif_smooth <- summary(gam_env)$s.table %>%
        as_tibble(rownames = "term") %>%
        rename(gam_edf = edf,
               gam_red_df = Ref.df,
               gam_f = `F`,
               gam_p_value = `p-value`) %>%
        filter(gam_p_value < 0.05) %>%
        pull(term)
      
      gam_env_terms_signif <- gam_env_terms_signif_smooth
    } else {  # Combination of linear and smooth terms
      gam_env_terms_signif_lm <- summary(gam_env)$p.table %>%
        as_tibble(rownames = "term") %>%
        rename(lm_estimate = Estimate,
               lm_std_error = `Std. Error`,
               lm_t_value = `t value`,
               lm_p_value = `Pr(>|t|)`) %>%
        filter(term != "(Intercept)") %>%
        filter(lm_p_value < 0.05) %>%
        pull(term)
      
      gam_env_terms_signif_smooth <- summary(gam_env)$s.table %>%
        as_tibble(rownames = "term") %>%
        rename(gam_edf = edf,
               gam_red_df = Ref.df,
               gam_f = `F`,
               gam_p_value = `p-value`) %>%
        filter(gam_p_value < 0.05) %>%
        pull(term)
      
      gam_env_terms_signif <- c(gam_env_terms_signif_lm, gam_env_terms_signif_smooth)
    }
    
    # Fit GAM using only variables that were significant in full GAM
    terms_signif_env <- gam_env_terms_signif
    vars_signif_env <- str_remove_all(str_remove_all(terms_signif_env, "^s\\("), "\\)")
    
    env_variables_signif <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_signif_env))])) %>%
      drop_na()
    
    if (length(terms_signif_env) > 0) {
      gam_signif_env <- gam(as.formula(paste(response_var, "~", paste(terms_signif_env, collapse = " + "))),
                            data = env_variables_signif, method = "REML", select = TRUE)
      # Summarize all GAMs
      gam_summary_tmp <- tibble(dataset = dataset_tmp,
                                lake_id = lakeid,
                                response_var = response_var,
                                
                                gam_full_vars = str_remove_all(str_remove_all(toString(terms_env), "s\\("), "\\)"),
                                gam_full_terms = toString(terms_env),
                                gam_full_r2 = summary(gam_env)$r.sq,
                                gam_full_dev_expl = summary(gam_env)$dev.expl,
                                gam_full_aic = AIC(gam_env),
                                
                                gam_signif_vars = str_remove_all(str_remove_all(toString(terms_signif_env), "s\\("), "\\)"),
                                gam_signif_terms = toString(terms_signif_env),
                                gam_signif_r2 = summary(gam_signif_env)$r.sq,
                                gam_signif_dev_expl = summary(gam_signif_env)$dev.expl,
                                gam_signif_aic = AIC(gam_signif_env))
    } else {
      # Summarize full GAM only
      gam_summary_tmp <- tibble(dataset = dataset_tmp,
                                lake_id = lakeid,
                                response_var = response_var,
                                
                                gam_full_vars = str_remove_all(str_remove_all(toString(terms_env), "s\\("), "\\)"),
                                gam_full_terms = toString(terms_env),
                                gam_full_r2 = summary(gam_env)$r.sq,
                                gam_full_dev_expl = summary(gam_env)$dev.expl,
                                gam_full_aic = AIC(gam_env))
    }
    
    return(gam_summary_tmp)
  }
}

# Define function to fit GAMs with explanatory variable across environmental categories across all lakes
computeGAMAll <- function(response_data, response_var) {
  gam_summaries_all <- bind_rows(computeGAM("L226N", response_data, response_var),
                                 computeGAM("L226S", response_data, response_var),
                                 computeGAM("L227", response_data, response_var),
                                 computeGAM("L224", response_data, response_var),
                                 computeGAM("L373", response_data, response_var))
  
  return(gam_summaries_all)
}

# Define function to plot GAM partial dependencies
plotGAMPartial <- function(lakeid, response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  gam_summaries_tmp <- computeGAM(lakeid, response_data, response_var)
  
  if (!any(class(gam_summaries_tmp) == "NULL")) {
    # Create vector of explanatory variables across environmental categories
    terms_full <- gam_summaries_tmp$gam_full_terms
    terms_full <- unlist(strsplit(terms_full, ", "))
    terms_full <- terms_full[!is.na(terms_full)]
    
    vars_full <- str_remove_all(str_remove_all(terms_full, "^s\\("), "\\)")
    
    full_variables <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_full))])) %>%
      drop_na()
    
    # Fit GAM using explanatory variables across environmental categories
    gam_full <- gam(as.formula(paste(response_var, "~", paste(terms_full, collapse = " + "))),
                    data = full_variables, method = "REML", select = TRUE)
    
    # Identify significant terms
    if (all(!grepl("s\\(", terms_full))) {  # Linear terms only
      gam_full_terms_signif_lm <- summary(gam_full)$p.table %>%
        as_tibble(rownames = "term") %>%
        rename(lm_estimate = Estimate,
               lm_std_error = `Std. Error`,
               lm_t_value = `t value`,
               lm_p_value = `Pr(>|t|)`) %>%
        filter(term != "(Intercept)") %>%
        filter(lm_p_value < 0.05) %>%
        pull(term)
      
      gam_full_terms_signif <- gam_full_terms_signif_lm
    } else if (all(grepl("s\\(", terms_full))) {  # Smooth terms only
      gam_full_terms_signif_smooth <- summary(gam_full)$s.table %>%
        as_tibble(rownames = "term") %>%
        rename(gam_edf = edf,
               gam_red_df = Ref.df,
               gam_f = `F`,
               gam_p_value = `p-value`) %>%
        filter(gam_p_value < 0.05) %>%
        pull(term)
      
      gam_full_terms_signif <- gam_full_terms_signif_smooth
    } else {  # Combination of linear and smooth terms
      gam_full_terms_signif_lm <- summary(gam_full)$p.table %>%
        as_tibble(rownames = "term") %>%
        rename(lm_estimate = Estimate,
               lm_std_error = `Std. Error`,
               lm_t_value = `t value`,
               lm_p_value = `Pr(>|t|)`) %>%
        filter(term != "(Intercept)") %>%
        filter(lm_p_value < 0.05) %>%
        pull(term)
      
      gam_full_terms_signif_smooth <- summary(gam_full)$s.table %>%
        as_tibble(rownames = "term") %>%
        rename(gam_edf = edf,
               gam_red_df = Ref.df,
               gam_f = `F`,
               gam_p_value = `p-value`) %>%
        filter(gam_p_value < 0.05) %>%
        pull(term)
      
      gam_full_terms_signif <- c(gam_full_terms_signif_lm, gam_full_terms_signif_smooth)
    }
    
    if (any(grepl("^s\\(", gam_full_terms_signif))) {
      # Format residuals
      residuals <- full_variables %>%
        add_partial_residuals(gam_full) %>%  # Add residuals for smoothed terms only
        pivot_longer(all_of(vars_full), names_to = "var", values_to = "var_value") %>%
        pivot_longer(starts_with("s("), names_to = "smooth", values_to = "smooth_value") %>%
        filter(str_detect(smooth, var))
      
      # Plot partial dependencies for smooth terms
      smooth_estimates(gam_full) %>%
        add_confint() %>%
        pivot_longer(all_of(unique(residuals$var)), names_to = "var", values_to = "value", values_drop_na = TRUE) %>%
        mutate(signif = case_when(smooth %in% gam_full_terms_signif ~ TRUE,
                                  TRUE ~ FALSE)) %>%
        ggplot() +
        facet_wrap(~var, scales = "free", strip.position = "bottom") +
        geom_rug(aes(x = var_value),
                 data = residuals,
                 sides = "b", length = grid::unit(0.02, "npc")) +
        geom_ribbon(aes(x = value, ymin = lower_ci, ymax = upper_ci),
                    alpha = 0.2) +
        geom_line(aes(x = value, y = est, linetype = signif),
                  lwd = 1) +
        geom_point(aes(x = var_value, y = smooth_value, colour = lakeid),
                   data = residuals,
                   alpha = 0.7, stroke = 0) +
        labs(y = "Partial effect") +
        ggtitle(str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
        scale_colour_manual(values = palette_lake, guide = "none") +
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
        theme_bw() %+replace%
        theme(axis.title.x = element_blank(),
              axis.text = element_text(colour = "black"),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(colour = "black"),
              panel.grid = element_blank(),
              panel.background = element_blank())
    }
  }
}


#### Fit GAMs based on nutrients, climate, and nutrients + climate ####
# Specify nutrient and climate variables
vars_nutrients <- c("epilimnion_tp", "epilimnion_tn")
vars_climate <- c("meanairtemp")

# Define function to compute GAMs based on nutrients, climate, and nutrients + climate
computeGAMNutriClim <- function(lakeid, response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  variables_nutrients <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_nutrients))])) %>%
    drop_na()
  
  variables_climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_climate))])) %>%
    drop_na()
  
  variables_nutrients.climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("lake_id", "year", vars_nutrients, vars_climate))])) %>%
    drop_na()
  
  # Compute and tweak GAMs based on nutrients
  n_knots <- 10
  while (!exists("gam_nutrients1") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub("$", paste0(", k = ", n_knots, ")"), gsub("^", "s(", vars_nutrients)), collapse = " + "))),
                       data = variables_nutrients, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_nutrients1 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!exists("gam_nutrients1")) {
    gam_nutrients1 <- gam(as.formula(paste(response_var, "~", paste(vars_nutrients, collapse = " + "))),
                          data = variables_nutrients, method = "REML", select = TRUE)
  }
  
  if (!is.null(summary(gam_nutrients1)$s.table)) {
    gam_nutrients1_summary <- summary(gam_nutrients1)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(term_winner = case_when(model_winner == "lm" ~ explanatory_var,
                                     model_winner == "gam" ~ str_c("s(", explanatory_var, ")")))
    
    gam_nutrients1_term_winners <- unique(gam_nutrients1_summary$term_winner)
  } else if (is.null(summary(gam_nutrients1)$s.table)) {
    gam_nutrients1_term_winners <- vars_nutrients
  }
  
  n_knots <- 10
  while (!exists("gam_nutrients2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), gam_nutrients1_term_winners), collapse = " + "))),
                       data = variables_nutrients, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_nutrients2 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!is.null(summary(gam_nutrients2)$p.table) & is.null(summary(gam_nutrients2)$s.table)) {
    gam_nutrients2_summary_linear <- summary(gam_nutrients2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_nutrients2_vars_linear_signif <- gam_nutrients2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients2_vars_p_value <- unique(gam_nutrients2_summary_linear$var_signif)
    gam_nutrients2_vars_signif <- gam_nutrients2_vars_linear_signif
    
  } else if (!is.null(summary(gam_nutrients2)$s.table) & is.null(summary(gam_nutrients2)$p.table)) {
    gam_nutrients2_summary_smooth <- summary(gam_nutrients2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_nutrients2_vars_smooth_signif <- gam_nutrients2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients2_vars_p_value <- unique(gam_nutrients2_summary_smooth$var_signif)
    gam_nutrients2_vars_signif <- gam_nutrients2_vars_smooth_signif
    
  } else if (!is.null(summary(gam_nutrients2)$p.table) & !is.null(summary(gam_nutrients2)$s.table)) {
    gam_nutrients2_summary_linear <- summary(gam_nutrients2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_nutrients2_vars_linear_signif <- gam_nutrients2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients2_summary_smooth <- summary(gam_nutrients2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_nutrients2_vars_smooth_signif <- gam_nutrients2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients2_vars_p_value <- c(unique(gam_nutrients2_summary_linear$var_signif), unique(gam_nutrients2_summary_smooth$var_signif))
    gam_nutrients2_vars_signif <- c(gam_nutrients2_vars_linear_signif, gam_nutrients2_vars_smooth_signif)
  }
  
  gam_nutrients2_summary <- tibble(dataset = dataset_tmp,
                                   lake_id = lakeid,
                                   response_var = response_var,
                                   explanatory_category = "nutrients",
                                   
                                   gam_vars = toString(vars_nutrients),
                                   gam_terms = toString(gam_nutrients1_term_winners),
                                   gam_vars_p_value = toString(gam_nutrients2_vars_p_value),
                                   gam_r2 = summary(gam_nutrients2)$r.sq,
                                   gam_dev_expl = summary(gam_nutrients2)$dev.expl,
                                   gam_aic = AIC(gam_nutrients2))
  
  gam_predictions_nutrients2 <- predict.gam(gam_nutrients2, se.fit = TRUE, type = "terms")
  
  gam_predictions_nutrients2_fit <- gam_predictions_nutrients2$fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("fit.", .x))
  
  gam_predictions_nutrients2_sefit <- gam_predictions_nutrients2$se.fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("sefit.", .x))
  
  variables_nutrients_fit2 <- bind_cols(variables_nutrients,
                                        gam_predictions_nutrients2_fit,
                                        gam_predictions_nutrients2_sefit) %>%
    pivot_longer(any_of(vars_nutrients), names_to = "var", values_to = "var_value") %>%
    pivot_longer(starts_with("fit."), names_to = "fit", values_to = "fit_value") %>%
    mutate(fit = str_remove(fit, "fit.")) %>%
    pivot_longer(starts_with("sefit."), names_to = "sefit", values_to = "sefit_value") %>%
    mutate(sefit = str_remove(sefit, "sefit.")) %>%
    filter(var == fit) %>%
    filter(var == sefit) %>%
    dplyr::select(-c(fit, sefit)) %>%
    mutate(signif = case_when(var %in% gam_nutrients2_vars_signif ~ TRUE,
                              TRUE ~ FALSE)) %>%
    mutate(sefit_min = fit_value - sefit_value,
           sefit_max = fit_value + sefit_value)
  
  gam_nutrients2_contribplot <- variables_nutrients_fit2 %>%
    ggplot() +
    geom_rect(aes(xmin = start_year, xmax = end_year,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations %>%
                mutate(start_year = year(start_date),
                       end_year = year(end_date) + 1) %>%
                filter(lake_id == lakeid)) +
    scale_fill_manual(values = palette_manipulation, guide = "none") +
    new_scale_fill() +
    geom_hline(yintercept = 0) +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_nutrients_fit2 %>%
    #               filter(var == "epilimnion_tn"),
    #             fill = "orange", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_nutrients_fit2 %>%
                 filter(var == "epilimnion_tn"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Oranges",
                         direction = 1,
                         trans = "sqrt",
                         limits = c(min(chemistry_byyear_interpolated$epilimnion_tn), max(chemistry_byyear_interpolated$epilimnion_tn)),
                         #breaks = c(0, 1, 2),
                         name = "TN (mg/L)") +
    new_scale_fill() +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_nutrients_fit2 %>%
    #               filter(var == "epilimnion_tp"),
    #             fill = "green", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_nutrients_fit2 %>%
                 filter(var == "epilimnion_tp"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Greens",
                         direction = 1,
                         trans = "log10",
                         limits = c(min(chemistry_byyear_interpolated$epilimnion_tp), max(chemistry_byyear_interpolated$epilimnion_tp)),
                         #breaks = c(1, 10, 100, 1000, 5000),
                         labels = comma,
                         name = "TP (\u00b5g/L)") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_size(name = response_var,
               limits = c(-1.1, 1.1),
               range = c(0.2, 3)) +
    ylim(c(-1.5, 1.5)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = "Effect") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0"))
  
  # Compute and tweak GAMs based on climate
  n_knots <- 10
  while (!exists("gam_climate1") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub("$", paste0(", k = ", n_knots, ")"), gsub("^", "s(", vars_climate)), collapse = " + "))),
                       data = variables_climate, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_climate1 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!exists("gam_climate1")) {
    gam_climate1 <- gam(as.formula(paste(response_var, "~", paste(vars_climate, collapse = " + "))),
                        data = variables_climate, method = "REML", select = TRUE)
  }
  
  if (!is.null(summary(gam_climate1)$s.table)) {
    gam_climate1_summary <- summary(gam_climate1)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(term_winner = case_when(model_winner == "lm" ~ explanatory_var,
                                     model_winner == "gam" ~ str_c("s(", explanatory_var, ")")))
    
    gam_climate1_term_winners <- unique(gam_climate1_summary$term_winner)
  } else if (is.null(summary(gam_climate1)$s.table)) {
    gam_climate1_term_winners <- vars_climate
  }
  
  n_knots <- 10
  while (!exists("gam_climate2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), gam_climate1_term_winners), collapse = " + "))),
                       data = variables_climate, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_climate2 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!is.null(summary(gam_climate2)$p.table) & is.null(summary(gam_climate2)$s.table)) {
    gam_climate2_summary_linear <- summary(gam_climate2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_climate2_vars_linear_signif <- gam_climate2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_climate2_vars_p_value <- unique(gam_climate2_summary_linear$var_signif)
    gam_climate2_vars_signif <- gam_climate2_vars_linear_signif
    
  } else if (!is.null(summary(gam_climate2)$s.table) & is.null(summary(gam_climate2)$p.table)) {
    gam_climate2_summary_smooth <- summary(gam_climate2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_climate2_vars_smooth_signif <- gam_climate2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_climate2_vars_p_value <- unique(gam_climate2_summary_smooth$var_signif)
    gam_climate2_vars_signif <- gam_climate2_vars_smooth_signif
    
  } else if (!is.null(summary(gam_climate2)$p.table) & !is.null(summary(gam_climate2)$s.table)) {
    gam_climate2_summary_linear <- summary(gam_climate2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_climate2_vars_linear_signif <- gam_climate2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_climate2_summary_smooth <- summary(gam_climate2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_climate2_vars_smooth_signif <- gam_climate2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_climate2_vars_p_value <- c(unique(gam_climate2_summary_linear$var_signif), unique(gam_climate2_summary_smooth$var_signif))
    gam_climate2_vars_signif <- c(gam_climate2_vars_linear_signif, gam_climate2_vars_smooth_signif)
  }
  
  gam_climate2_summary <- tibble(dataset = dataset_tmp,
                                 lake_id = lakeid,
                                 response_var = response_var,
                                 explanatory_category = "climate",
                                 
                                 gam_vars = toString(vars_climate),
                                 gam_terms = toString(gam_climate1_term_winners),
                                 gam_vars_p_value = toString(gam_climate2_vars_p_value),
                                 gam_r2 = summary(gam_climate2)$r.sq,
                                 gam_dev_expl = summary(gam_climate2)$dev.expl,
                                 gam_aic = AIC(gam_climate2))
  
  gam_predictions_climate2 <- predict.gam(gam_climate2, se.fit = TRUE, type = "terms")
  
  gam_predictions_climate2_fit <- gam_predictions_climate2$fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("fit.", .x))
  
  gam_predictions_climate2_sefit <- gam_predictions_climate2$se.fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("sefit.", .x))
  
  variables_climate_fit2 <- bind_cols(variables_climate,
                                      gam_predictions_climate2_fit,
                                      gam_predictions_climate2_sefit) %>%
    pivot_longer(any_of(vars_climate), names_to = "var", values_to = "var_value") %>%
    pivot_longer(starts_with("fit."), names_to = "fit", values_to = "fit_value") %>%
    mutate(fit = str_remove(fit, "fit.")) %>%
    pivot_longer(starts_with("sefit."), names_to = "sefit", values_to = "sefit_value") %>%
    mutate(sefit = str_remove(sefit, "sefit.")) %>%
    filter(var == fit) %>%
    filter(var == sefit) %>%
    dplyr::select(-c(fit, sefit)) %>%
    mutate(signif = case_when(var %in% gam_climate2_vars_signif ~ TRUE,
                              TRUE ~ FALSE)) %>%
    mutate(sefit_min = fit_value - sefit_value,
           sefit_max = fit_value + sefit_value)
  
  gam_climate2_contribplot <- variables_climate_fit2 %>%
    ggplot() +
    geom_rect(aes(xmin = start_year, xmax = end_year,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations %>%
                mutate(start_year = year(start_date),
                       end_year = year(end_date) + 1) %>%
                filter(lake_id == lakeid)) +
    scale_fill_manual(values = palette_manipulation, guide = "none") +
    new_scale_fill() +
    geom_hline(yintercept = 0) +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_climate_fit2 %>%
    #               filter(var == "meanairtemp"),
    #             fill = "steelblue", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              variables_climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_climate_fit2 %>%
                 filter(var == "meanairtemp"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "RdBu",
                         direction = -1,
                         limits = c(min(climate_eccc_byyear$meanairtemp), max(climate_eccc_byyear$meanairtemp)),
                         #trans = "sqrt",
                         #breaks = c(0, 1, 2),
                         name = "Temp (\u00B0C)") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_size(name = response_var,
               limits = c(-1.1, 1.1),
               range = c(0.2, 3)) +
    ylim(c(-1.5, 1.5)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(y = "Effect") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0"))
  
  # Compute and tweak GAMs based on nutrients and climate
  n_knots <- 10
  while (!exists("gam_nutrients.climate1") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), c(gam_nutrients1_term_winners, gam_climate1_term_winners)), collapse = " + "))),
                       data = variables_nutrients.climate, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_nutrients.climate1 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!exists("gam_nutrients.climate1")) {
    if (exists("gam_nutrients1_summary") & exists("gam_climate1_summary")) {
      gam_nutrients.climate1 <- gam(as.formula(paste(response_var, "~", paste(str_remove_all(str_remove_all(c(unique(gam_nutrients1_summary$term_winner), unique(gam_climate1_summary$term_winner)), "^s\\("), "\\)"), collapse = " + "))),
                                    data = variables_nutrients.climate, method = "REML", select = TRUE)
    } else if (!exists("gam_nutrients1_summary") & exists("gam_climate1_summary")) {
      gam_nutrients.climate1 <- gam(as.formula(paste(response_var, "~", paste(str_remove_all(str_remove_all(unique(gam_climate1_summary$term_winner), "^s\\("), "\\)"), collapse = " + "))),
                                    data = variables_nutrients.climate, method = "REML", select = TRUE)
    }
  }
  
  if (!is.null(summary(gam_nutrients.climate1)$p.table) & is.null(summary(gam_nutrients.climate1)$s.table)) {
    gam_nutrients.climate1_summary_linear <- summary(gam_nutrients.climate1)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)")
    
    gam_nutrients.climate1_summary_vars <- unique(gam_nutrients.climate1_summary_linear$explanatory_var)
  } else if (!is.null(summary(gam_nutrients.climate1)$s.table) & is.null(summary(gam_nutrients.climate1)$p.table)) {
    gam_nutrients.climate1_summary_smooth <- summary(gam_nutrients.climate1)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(term_winner = case_when(model_winner == "lm" ~ explanatory_var,
                                     model_winner == "gam" ~ str_c("s(", explanatory_var, ")")))
    
    gam_nutrients.climate1_summary_vars <- unique(gam_nutrients.climate1_summary_smooth$term_winner)
  } else if (!is.null(summary(gam_nutrients.climate1)$p.table) & !is.null(summary(gam_nutrients.climate1)$s.table)) {
    gam_nutrients.climate1_summary_linear <- summary(gam_nutrients.climate1)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)")
    
    gam_nutrients.climate1_summary_smooth <- summary(gam_nutrients.climate1)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(term_winner = case_when(model_winner == "lm" ~ explanatory_var,
                                     model_winner == "gam" ~ str_c("s(", explanatory_var, ")")))
    
    gam_nutrients.climate1_summary_vars <- c(unique(gam_nutrients.climate1_summary_linear$explanatory_var), unique(gam_nutrients.climate1_summary_smooth$term_winner))
  }
  
  n_knots <- 10
  while (!exists("gam_nutrients.climate2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), gam_nutrients.climate1_summary_vars), collapse = " + "))),
                       data = variables_nutrients.climate, method = "REML", select = TRUE), silent = TRUE)
    
    if (exists("gam_tmp")) {
      gam_nutrients.climate2 <- gam_tmp
      rm(gam_tmp)
      
    } else {
      n_knots <- n_knots - 1
      print(paste0("Now trying ", n_knots, " knots"))
    }
  }
  
  if (!is.null(summary(gam_nutrients.climate2)$p.table) & is.null(summary(gam_nutrients.climate2)$s.table)) {
    gam_nutrients.climate2_summary_linear <- summary(gam_nutrients.climate2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_nutrients.climate2_vars_linear_signif <- gam_nutrients.climate2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients.climate2_vars_p_value <- unique(gam_nutrients.climate2_summary_linear$var_signif)
    gam_nutrients.climate2_vars_signif <- gam_nutrients.climate2_vars_linear_signif
    
  } else if (!is.null(summary(gam_nutrients.climate2)$s.table) & is.null(summary(gam_nutrients.climate2)$p.table)) {
    gam_nutrients.climate2_summary_smooth <- summary(gam_nutrients.climate2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(model_winner = case_when(gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
                                      TRUE ~ "gam")) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_nutrients.climate2_vars_smooth_signif <- gam_nutrients.climate2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients.climate2_vars_p_value <- unique(gam_nutrients.climate2_summary_smooth$var_signif)
    gam_nutrients.climate2_vars_p_value <- c(unique(gam_nutrients.climate2_summary_linear$var_signif), unique(gam_nutrients.climate2_summary_smooth$var_signif))
    
  } else if (!is.null(summary(gam_nutrients.climate2)$p.table) & !is.null(summary(gam_nutrients.climate2)$s.table)) {
    gam_nutrients.climate2_summary_linear <- summary(gam_nutrients.climate2)$p.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      rename(lm_estimate = Estimate,
             lm_std_error = `Std. Error`,
             lm_t_value = `t value`,
             lm_p_value = `Pr(>|t|)`) %>%
      filter(explanatory_var != "(Intercept)") %>%
      mutate(var_signif = str_c(explanatory_var, lm_p_value, sep = "|"))
    
    gam_nutrients.climate2_vars_linear_signif <- gam_nutrients.climate2_summary_linear %>%
      filter(lm_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients.climate2_summary_smooth <- summary(gam_nutrients.climate2)$s.table %>%
      as_tibble(rownames = "explanatory_var") %>%
      mutate(explanatory_var = str_remove_all(str_remove_all(explanatory_var, "^s\\("), "\\)")) %>%
      rename(gam_edf = edf,
             gam_red_df = Ref.df,
             gam_f = `F`,
             gam_p_value = `p-value`) %>%
      mutate(var_signif = str_c(explanatory_var, gam_p_value, sep = "|"))
    
    gam_nutrients.climate2_vars_smooth_signif <- gam_nutrients.climate2_summary_smooth %>%
      filter(gam_p_value < 0.05) %>%
      pull(explanatory_var)
    
    gam_nutrients.climate2_vars_p_value <- c(unique(gam_nutrients.climate2_summary_linear$var_signif), unique(gam_nutrients.climate2_summary_smooth$var_signif))
    gam_nutrients.climate2_vars_signif <- c(gam_nutrients.climate2_vars_linear_signif, gam_nutrients.climate2_vars_smooth_signif)
  }
  
  gam_nutrients.climate2_summary <- tibble(dataset = dataset_tmp,
                                           lake_id = lakeid,
                                           response_var = response_var,
                                           explanatory_category = "nutrients.climate",
                                           
                                           gam_vars = toString(c(vars_nutrients, vars_climate)),
                                           gam_terms = toString(gam_nutrients.climate1_summary_vars),
                                           gam_vars_p_value = toString(gam_nutrients.climate2_vars_p_value),
                                           gam_r2 = summary(gam_nutrients.climate2)$r.sq,
                                           gam_dev_expl = summary(gam_nutrients.climate2)$dev.expl,
                                           gam_aic = AIC(gam_nutrients.climate2))
  
  gam_predictions_nutrients.climate2 <- predict.gam(gam_nutrients.climate2, se.fit = TRUE, type = "terms")
  
  gam_predictions_nutrients.climate2_fit <- gam_predictions_nutrients.climate2$fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("fit.", .x))
  
  gam_predictions_nutrients.climate2_sefit <- gam_predictions_nutrients.climate2$se.fit %>%
    as_tibble() %>%
    rename_with(~gsub("\\)$", "", gsub("^s\\(", "", .x))) %>%
    rename_with(~str_c("sefit.", .x))
  
  variables_nutrients.climate_fit2 <- bind_cols(variables_nutrients.climate,
                                                gam_predictions_nutrients.climate2_fit,
                                                gam_predictions_nutrients.climate2_sefit) %>%
    pivot_longer(any_of(c(vars_nutrients, vars_climate)), names_to = "var", values_to = "var_value") %>%
    pivot_longer(starts_with("fit."), names_to = "fit", values_to = "fit_value") %>%
    mutate(fit = str_remove(fit, "fit.")) %>%
    pivot_longer(starts_with("sefit."), names_to = "sefit", values_to = "sefit_value") %>%
    mutate(sefit = str_remove(sefit, "sefit.")) %>%
    filter(var == fit) %>%
    filter(var == sefit) %>%
    dplyr::select(-c(fit, sefit)) %>%
    mutate(signif = case_when(var %in% gam_nutrients.climate2_vars_signif ~ TRUE,
                              TRUE ~ FALSE)) %>%
    mutate(sefit_min = fit_value - sefit_value,
           sefit_max = fit_value + sefit_value)
  
  gam_nutrients.climate2_contribplot <- variables_nutrients.climate_fit2 %>%
    ggplot() +
    geom_rect(aes(xmin = start_year, xmax = end_year,
                  ymin = -Inf, ymax = Inf,
                  fill = event),
              data = ela_manipulations %>%
                mutate(start_year = year(start_date),
                       end_year = year(end_date) + 1) %>%
                filter(lake_id == lakeid)) +
    scale_fill_manual(values = palette_manipulation, guide = "none") +
    new_scale_fill() +
    geom_hline(yintercept = 0) +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_nutrients.climate_fit2 %>%
    #               filter(var == "epilimnion_tn"),
    #             fill = "orange", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_nutrients.climate_fit2 %>%
                 filter(var == "epilimnion_tn"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Oranges",
                         direction = 1,
                         trans = "sqrt",
                         limits = c(min(chemistry_byyear_interpolated$epilimnion_tn), max(chemistry_byyear_interpolated$epilimnion_tn)),
                         #breaks = c(0, 1, 2),
                         name = "TN (mg/L)") +
    new_scale_fill() +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_nutrients.climate_fit2 %>%
    #               filter(var == "epilimnion_tp"),
    #             fill = "green", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_nutrients.climate_fit2 %>%
                 filter(var == "epilimnion_tp"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "Greens",
                         direction = 1,
                         trans = "log10",
                         limits = c(min(chemistry_byyear_interpolated$epilimnion_tp), max(chemistry_byyear_interpolated$epilimnion_tp)),
                         #breaks = c(1, 10, 100, 1000, 5000),
                         labels = comma,
                         name = "TP (\u00b5g/L)") +
    new_scale_fill() +
    # geom_ribbon(aes(ymin = sefit_min,
    #                 ymax = sefit_max),
    #             variables_nutrients.climate_fit2 %>%
    #               filter(var == "meanairtemp"),
    #             fill = "steelblue", alpha = 0.2) +
    geom_line(aes(x = year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = year, y = fit_value,
                   fill = var_value,
                   size = eval(as.name(paste(response_var)))),
               data = variables_nutrients.climate_fit2 %>%
                 filter(var == "meanairtemp"),
               colour = "black", pch = 21, stroke = 0.1) +
    scale_fill_distiller(palette = "RdBu",
                         direction = -1,
                         limits = c(min(climate_eccc_byyear$meanairtemp), max(climate_eccc_byyear$meanairtemp)),
                         #trans = "sqrt",
                         #breaks = c(0, 1, 2),
                         name = "Temp (\u00B0C)") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_size(name = response_var,
               limits = c(-1.1, 1.1),
               range = c(0.2, 3)) +
    scale_x_continuous(expand = c(0,0)) +
    ylim(c(-1.5, 1.5)) +
    labs(y = "Effect") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0"))
  
  # Combine contribution plots
  if (exists("gam_nutrients2_summary") & exists("gam_climate2_summary") & exists("gam_nutrients.climate2_summary")) {
    gam_contribplot_all <- gam_nutrients2_contribplot +
      gam_climate2_contribplot +
      gam_nutrients.climate2_contribplot +
      plot_annotation(title = str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
      plot_layout(guides = "collect", widths = c(1, 1, 1))
    
    # Summarize all GAMs
    gam_summaries <- bind_rows(gam_nutrients2_summary,
                               gam_climate2_summary,
                               gam_nutrients.climate2_summary)
    
  } else if (!exists("gam_nutrients2_summary") & exists("gam_climate2_summary") & exists("gam_nutrients.climate2_summary")) {
    gam_contribplot_all <- plot_spacer() +
      gam_climate2_contribplot +
      gam_nutrients.climate2_contribplot +
      plot_annotation(title = str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
      plot_layout(guides = "collect", widths = c(1, 1, 1))
    
    # Summarize all GAMs
    gam_summaries <- bind_rows(gam_climate2_summary,
                               gam_nutrients.climate2_summary)
  }
  
  gam_all <- list(summaries = gam_summaries,
                  plot = gam_contribplot_all)
  
  return(gam_all)
}

# Define function to fit GAMs with explanatory variable across environmental categories across all lakes
computeGAMNutriClimAll <- function(response_data, response_var) {
  gam_summaries_all <- bind_rows(computeGAMNutriClim("L226N", response_data, response_var)$summaries,
                                 computeGAMNutriClim("L226S", response_data, response_var)$summaries,
                                 computeGAMNutriClim("L227", response_data, response_var)$summaries,
                                 computeGAMNutriClim("L224", response_data, response_var)$summaries,
                                 computeGAMNutriClim("L373", response_data, response_var)$summaries)
  
  return(gam_summaries_all)
}

# Define function to plot contributions of nutrients and climate from GAMs
plotGAMNutriClimContrib <- function(response_data, response_var) {
  computeGAMNutriClim("L226N", response_data, response_var)$plot /
    computeGAMNutriClim("L226S", response_data, response_var)$plot /
    computeGAMNutriClim("L227", response_data, response_var)$plot /
    computeGAMNutriClim("L224", response_data, response_var)$plot /
    computeGAMNutriClim("L373", response_data, response_var)$plot +
    plot_layout(guides = "collect")
}
