# Generalized Additive Models of sediment assemblages
# https://r.qcbs.ca/workshops/r-workshop-08/

# setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
# setwd("~/Dropbox/projects/ela18s/")

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
  env_all <- syncEnv(lakeid, watertemp_byyear_interpolated, "epilimnion_temp") %>%
    full_join(syncEnv(lakeid, watertemp_byyear_interpolated, "metalimnion_temp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, watertemp_byyear_interpolated, "hypolimnion_temp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, sdd_byyear_interpolated, "sdd"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "epilimnion_tn"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "epilimnion_tp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "epilimnion_ph"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "epilimnion_oxygen"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "profundal_tn"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "profundal_tp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "profundal_ph"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, chemistry_byyear_interpolated, "profundal_oxygen"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_eccc_byyear, "meanairtemp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_eccc_byyear, "precip") %>%
                rename(precip_eccc = precip), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_ela_byyear, "maxairtemp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_ela_byyear, "minairtemp"), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_ela_byyear, "precip") %>%
                rename(precip_ela = precip), by = "layer_sync") %>%
    full_join(syncEnv(lakeid, climate_ela_byyear, "ice_duration"), by = "layer_sync") %>%
    dplyr::select(-starts_with("nyear"))
  
  # Join all response and explanatory variables
  variables <- response %>%
    full_join(env_all, by = c("layer_order" = "layer_sync")) %>%
    filter(!is.na(eval(as.name(paste(response_var))))) %>%
    arrange(midpt)
  
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
    pivot_longer(all_of(names(env_all)[which(names(env_all) != "layer_sync")]),
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
          strip.text = element_text(colour = "black"),
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
  vars_all <- names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_discard))]
  
  models_ind <- tibble(NULL)
  for (i in 1:length(vars_all)) {
    dataset_tmp <- unique(response_data$dataset)
    
    # Extract a single explanatory variable
    var_tmp <- vars_all[i]
    
    variables_tmp <- variables %>%
      dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", var_tmp))])) %>%
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

# # Combine results from linear models and GAMs for each explanatory variable
# models_ind_asvs_pco_all <- bind_rows(modelIndVarAll(pcoa %>%
#                                                       filter(dataset == "Sediment ASVs"), "PCo1"),
#                                      modelIndVarAll(pcoa %>%
#                                                       filter(dataset == "Sediment ASVs"), "PCo2"),
#                                      
#                                      modelIndVarAll(pcoa %>%
#                                                       filter(dataset == "Phototroph ASVs"), "PCo1"),
#                                      modelIndVarAll(pcoa %>%
#                                                       filter(dataset == "Phototroph ASVs"), "PCo2"))
# 
# # Visualize results from linear models and GAMs for each explanatory variable
# (models_ind_asvs_pco_all_heatmap <- models_ind_asvs_pco_all %>%
#     mutate(model_winner = case_when(!is.na(lm_aic) & is.na(gam_aic) ~ "lm",
#                                     is.na(lm_aic) & !is.na(gam_aic) ~ "gam",
#                                     
#                                     lm_aic < gam_aic & abs(lm_aic - gam_aic) >= 4 ~ "lm",
#                                     gam_aic < lm_aic & abs(gam_aic - lm_aic) >= 4 ~ "gam",
#                                     
#                                     gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
#                                     gam_edf >= 2 ~ "gam",
#                                     
#                                     lm_p_value < 0.05 & gam_p_value >= 0.05 ~ "lm",
#                                     lm_p_value >= 0.05 & gam_p_value < 0.05 ~ "gam",
#                                     
#                                     TRUE ~ "gam")) %>%
#     mutate(p_value = case_when(model_winner == "lm" ~ lm_p_value,
#                                model_winner == "gam" ~ gam_p_value),
#            dev_expl = case_when(model_winner == "lm" ~ lm_dev_expl,
#                                 model_winner == "gam" ~ gam_dev_expl)) %>%
#     mutate(signif_level = case_when(p_value >= 0.05 ~ "",
#                                     p_value < 0.05 & p_value > 0.01 ~ "*",
#                                     p_value <= 0.01 & p_value > 0.001 ~ "**",
#                                     p_value <= 0.001 ~ "***")) %>%
#     mutate(model_winner_abbr = case_when(model_winner == "lm" ~ "L",
#                                          model_winner == "gam" ~ "S"),
#            dataset_response_var = str_c(dataset, response_var, sep = " ")) %>%
#     ggplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
#                y = factor(explanatory_var, levels = rev(unique(models_ind_asvs_pco_all$explanatory_var))))) +
#     facet_wrap(~dataset_response_var, nrow = 1) +
#     geom_tile(aes(fill = dev_expl * 100), colour = "white") +
#     geom_text(aes(label = str_c(model_winner_abbr, "\n", signif_level))) +
#     scale_fill_gradient2(low = "red", mid = "white", high = "dodgerblue",
#                          midpoint = 0, na.value = "white", limits = c(0,100),
#                          name = "Deviance explained (%)") +
#     scale_x_discrete(expand = c(0,0)) +
#     scale_y_discrete(expand = c(0,0)) +
#     theme_bw() %+replace%
#     theme(axis.text = element_text(colour = "black"),
#           axis.title = element_blank(),
#           axis.ticks.x = element_blank(),
#           strip.background = element_blank(),
#           strip.text = element_text(colour = "black"),
#           panel.grid = element_blank(),
#           legend.position = "bottom"))
# #ggsave("figures/ela18s_pco_models_ind_heatmap.pdf", models_ind_asvs_pco_all_heatmap, width = 14, height = 8, units = "in", device = cairo_pdf)

# # Plot curves for all models based on individual variables
# dataset_tmp <- unique(response_data$dataset)
# 
# # Format response and explanatory variables
# env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
# variables <- formatGAMVars(lakeid, response_data, response_var)$variables
# 
# # Evaluate which variables should be modeled with linear or smooth terms
# model_eval <- modelIndVarAll(response_data, response_var) %>%
#   mutate(model_winner = case_when(!is.na(lm_aic) & is.na(gam_aic) ~ "lm",
#                                   is.na(lm_aic) & !is.na(gam_aic) ~ "gam",
#                                   
#                                   lm_aic < gam_aic & abs(lm_aic - gam_aic) >= 4 ~ "lm",
#                                   gam_aic < lm_aic & abs(gam_aic - lm_aic) >= 4 ~ "gam",
#                                   
#                                   gam_edf >= 0.5 & gam_edf <= 1.5 ~ "lm",
#                                   gam_edf >= 2 ~ "gam",
#                                   
#                                   lm_p_value < 0.05 & gam_p_value >= 0.05 ~ "lm",
#                                   lm_p_value >= 0.05 & gam_p_value < 0.05 ~ "gam",
#                                   
#                                   TRUE ~ "gam")) %>%
#   mutate(p_value = case_when(model_winner == "lm" ~ lm_p_value,
#                              model_winner == "gam" ~ gam_p_value),
#          dev_expl = case_when(model_winner == "lm" ~ lm_dev_expl,
#                               model_winner == "gam" ~ gam_dev_expl)) %>%
#   mutate(signif = case_when(p_value < 0.05 ~ TRUE,
#                             TRUE ~ FALSE),
#          signif_level = case_when(p_value >= 0.05 ~ "",
#                                   p_value < 0.05 & p_value > 0.01 ~ "*",
#                                   p_value <= 0.01 & p_value > 0.001 ~ "**",
#                                   p_value <= 0.001 ~ "***")) %>%
#   mutate(model_winner_abbr = case_when(model_winner == "lm" ~ "L",
#                                        model_winner == "gam" ~ "S"),
#          dataset_response_var = str_c(dataset, response_var, sep = " ")) %>%
#   filter(lake_id == lakeid)
# 
# smooth_plots <- list()
# for (i in 1:length(unique(model_eval$explanatory_var))) {
#   explanatory_var_tmp <- model_eval$explanatory_var[i]
#   model_winner_tmp <- model_eval$model_winner[i]
#   signif_tmp <- model_eval$signif[i]
#   variables_tmp <- variables %>%
#     dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", explanatory_var_tmp))])) %>%
#     drop_na()
#   
#   if (model_winner_tmp == "lm") {
#     plot_tmp <- variables_tmp %>%
#       mutate(signif = signif_tmp) %>%
#       ggplot(aes(x = eval(as.name(paste(explanatory_var_tmp))), y = eval(as.name(paste(response_var))))) +
#       #geom_rug(sides = "b", length = grid::unit(0.02, "npc")) +
#       geom_point(aes(colour = lakeid)) +
#       stat_smooth(aes(linetype = signif), method = "lm", lwd = 1, colour = "black") +
#       labs(x = explanatory_var_tmp, y = response_var) +
#       scale_colour_manual(values = palette_lake, guide = "none") +
#       scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
#       theme_bw() %+replace%
#       theme(axis.text = element_text(colour = "black"),
#             strip.background = element_blank(),
#             strip.placement = "outside",
#             strip.text = element_text(colour = "black"),
#             panel.grid = element_blank(),
#             panel.background = element_blank())
#   } else if (model_winner_tmp == "gam") {
#     # Fit GAM
#     gam_tmp <- gam(as.formula(paste(response_var, "~", paste0("s(", explanatory_var_tmp, ")"))),
#                    data = variables_tmp, method = "REML", select = TRUE)
#     
#     # Format residuals
#     residuals <- variables_tmp %>%
#       add_partial_residuals(gam_tmp) %>%  # Add residuals for smoothed terms only
#       pivot_longer(all_of(explanatory_var_tmp), names_to = "var", values_to = "var_value") %>%
#       pivot_longer(starts_with("s("), names_to = "smooth", values_to = "smooth_value") %>%
#       filter(str_detect(smooth, var))
#     
#     # Plot partial dependencies for smooth terms
#     plot_tmp <- smooth_estimates(gam_tmp) %>%
#       add_confint() %>%
#       pivot_longer(all_of(unique(residuals$var)), names_to = "var", values_to = "value", values_drop_na = TRUE) %>%
#       mutate(signif = case_when(signif_tmp ~ TRUE,
#                                 TRUE ~ FALSE)) %>%
#       ggplot() +
#       geom_rug(aes(x = var_value),
#                data = residuals,
#                sides = "b", length = grid::unit(0.02, "npc")) +
#       geom_ribbon(aes(x = value, ymin = lower_ci, ymax = upper_ci),
#                   alpha = 0.2) +
#       geom_line(aes(x = value, y = est, linetype = signif),
#                 lwd = 1) +
#       geom_point(aes(x = var_value, y = smooth_value, colour = lakeid),
#                  data = residuals,
#                  alpha = 0.7, stroke = 0) +
#       labs(x = explanatory_var_tmp,
#            y = response_var) +
#       #ggtitle(str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
#       scale_colour_manual(values = palette_lake, guide = "none") +
#       scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
#       theme_bw() %+replace%
#       theme(axis.text = element_text(colour = "black"),
#             strip.background = element_blank(),
#             strip.placement = "outside",
#             strip.text = element_text(colour = "black"),
#             panel.grid = element_blank(),
#             panel.background = element_blank())
#   }
#   
#   smooth_plots[[i]] <- plot_tmp
# }
# 
# gam_plot <- wrap_plots(smooth_plots) +
#   plot_layout(nrow = 1) +
#   plot_annotation(title = str_c(lakeid, dataset_tmp, response_var, sep = " "))
# 
# return(gam_plot)


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
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_nutrients))])) %>%
    drop_na()
  
  variables_climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_climate))])) %>%
    drop_na()
  
  variables_nutrients.climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_nutrients, vars_climate))])) %>%
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
  
  n_knots <- 10
  while (!exists("gam_nutrients2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), unique(gam_nutrients1_summary$term_winner)), collapse = " + "))),
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
                                   gam_terms = toString(gam_nutrients1_summary$term_winner),
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
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
  
  n_knots <- 10
  while (!exists("gam_climate2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), unique(gam_climate1_summary$term_winner)), collapse = " + "))),
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
                                 gam_terms = toString(gam_climate1_summary$term_winner),
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(1900, 2020, 10)) +
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
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), c(unique(gam_nutrients1_summary$term_winner), unique(gam_climate1_summary$term_winner))), collapse = " + "))),
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
    gam_nutrients.climate1 <- gam(as.formula(paste(response_var, "~", paste(str_remove_all(str_remove_all(c(unique(gam_nutrients1_summary$term_winner), unique(gam_climate1_summary$term_winner)), "^s\\("), "\\)"), collapse = " + "))),
                                  data = variables_nutrients.climate, method = "REML", select = TRUE)
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
    labs(y = "Effect") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0"))
  
  # Combine contribution plots
  gam_contribplot_all <- gam_nutrients2_contribplot +
    gam_climate2_contribplot +
    gam_nutrients.climate2_contribplot +
    plot_annotation(title = str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
    plot_layout(guides = "collect", widths = c(1, 2, 1))
  
  # Summarize all GAMs
  gam_summaries <- bind_rows(gam_nutrients2_summary,
                             gam_climate2_summary,
                             gam_nutrients.climate2_summary)
  
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

# Define function plot contributions of nutrients and climate from GAMs
plotGAMNutriClimContrib <- function(response_data, response_var) {
  computeGAMNutriClim("L226N", response_data, response_var)$plot /
    computeGAMNutriClim("L226S", response_data, response_var)$plot /
    computeGAMNutriClim("L227", response_data, response_var)$plot /
    computeGAMNutriClim("L224", response_data, response_var)$plot /
    computeGAMNutriClim("L373", response_data, response_var)$plot +
    plot_layout(guides = "collect")
}


#### Fit GAMs based on response/predictor variables within monitoring periods ####
# Define function to compute GAMs based on nutrients, climate, and nutrients + climate
# within monitoring period
computeGAMNutriClimMonitor <- function(lakeid, response_data, response_var) {
  dataset_tmp <- unique(response_data$dataset)
  
  # Format response and explanatory variables
  env_all <- formatGAMVars(lakeid, response_data, response_var)$env_all
  variables <- formatGAMVars(lakeid, response_data, response_var)$variables
  
  # Filter response and explanatory variables within TP/TN monitoring period
  env_all <- env_all %>%
    filter(!is.na(epilimnion_tp) & !is.na(epilimnion_tn))
  
  variables <- variables %>%
    filter(!is.na(epilimnion_tp) & !is.na(epilimnion_tn))
  
  variables_nutrients <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_nutrients))])) %>%
    drop_na()
  
  variables_climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_climate))])) %>%
    drop_na()
  
  variables_nutrients.climate <- variables %>%
    dplyr::select(-all_of(names(env_all)[which(!names(env_all) %in% c("layer_sync", vars_nutrients, vars_climate))])) %>%
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
  
  n_knots <- 10
  while (!exists("gam_nutrients2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), unique(gam_nutrients1_summary$term_winner)), collapse = " + "))),
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
                                   gam_terms = toString(gam_nutrients1_summary$term_winner),
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
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
  
  n_knots <- 10
  while (!exists("gam_climate2") & n_knots >= 1) {
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), unique(gam_climate1_summary$term_winner)), collapse = " + "))),
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
                                 gam_terms = toString(gam_climate1_summary$term_winner),
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(1900, 2020, 10)) +
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
    try(gam_tmp <- gam(as.formula(paste(response_var, "~", paste(gsub(")$", paste0(", k = ", n_knots, ")"), c(unique(gam_nutrients1_summary$term_winner), unique(gam_climate1_summary$term_winner))), collapse = " + "))),
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
    gam_nutrients.climate1 <- gam(as.formula(paste(response_var, "~", paste(str_remove_all(str_remove_all(c(unique(gam_nutrients1_summary$term_winner), unique(gam_climate1_summary$term_winner)), "^s\\("), "\\)"), collapse = " + "))),
                                  data = variables_nutrients.climate, method = "REML", select = TRUE)
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tn"),
              colour = "orange") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "epilimnion_tp"),
              colour = "green") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    geom_line(aes(x = midpt_year, y = fit_value, linetype = signif),
              data = variables_nutrients.climate_fit2 %>%
                filter(var == "meanairtemp"),
              colour = "steelblue") +
    geom_point(aes(x = midpt_year, y = fit_value,
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
    ylim(c(-4, 4)) +
    labs(y = "Effect") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#00B0F0"))
  
  # Combine contribution plots
  gam_contribplot_all <- gam_nutrients2_contribplot +
    gam_climate2_contribplot +
    gam_nutrients.climate2_contribplot +
    plot_annotation(title = str_c(lakeid, dataset_tmp, response_var, sep = " ")) +
    plot_layout(guides = "collect", widths = c(1, 2, 1))
  
  # Summarize all GAMs
  gam_summaries <- bind_rows(gam_nutrients2_summary,
                             gam_climate2_summary,
                             gam_nutrients.climate2_summary)
  
  gam_all <- list(summaries = gam_summaries,
                  plot = gam_contribplot_all)
  
  return(gam_all)
}

# Define function to fit GAMs with explanatory variable across environmental categories across all lakes
# within monitoring period
computeGAMNutriClimAllMonitor <- function(response_data, response_var) {
  gam_summaries_all <- bind_rows(computeGAMNutriClimMonitor("L226N", response_data, response_var)$summaries,
                                 computeGAMNutriClimMonitor("L226S", response_data, response_var)$summaries,
                                 computeGAMNutriClimMonitor("L227", response_data, response_var)$summaries,
                                 computeGAMNutriClimMonitor("L224", response_data, response_var)$summaries,
                                 computeGAMNutriClimMonitor("L373", response_data, response_var)$summaries)
  
  return(gam_summaries_all)
}
