# Taxonomic and sequence composition

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")
setwd("~/Library/CloudStorage/Dropbox-Personal/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Visualize taxonomic composition ####
# Sediment taxonomic composition
(tax_barplots_sediments <- sediment_melt %>%
    mutate(subdivision = case_when(is.na(subdivision) & division == "Opisthokonta" ~ "Unclassified Opisthokonta",
                                   is.na(subdivision) & division == "Alveolata" ~ "Unclassified Alveolata",
                                   is.na(subdivision) & is.na(division) ~ "Unclassified subdivision",
                                   !is.na(subdivision) & !subdivision %in% tax_order$subdivision[which(tax_order$relseqs >= 0.1)] ~ "Other subdivision",
                                   TRUE ~ subdivision)) %>%
    group_by(lake_id, upperpt, lowerpt, interval, midpt_year, subdivision) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               nrow = 1) +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
                 fill = factor(subdivision, levels = c(tax_order$subdivision[which(tax_order$subdivision %in% names(palette_subdivision))],
                                                       "Other subdivision", "Unclassified Alveolata", "Unclassified Opisthokonta", "Unclassified subdivision"))),
             stat = "identity") +
    scale_fill_manual(values = palette_subdivision, na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.1,
                  label = str_c(round(midpt_year, 0))),
              data = sediment_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "black") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0)) +
    coord_flip(ylim = c(0, 1), clip = "off") +
    labs(x = "Depth interval (cm)",
         y = "Relative sequence abundance",
         fill = "Sequence abundance (%)") +
    theme_bw() %+replace%
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_sediment_barplots.pdf", tax_barplots_sediments, width = 16, height = 10, units = "in", device = cairo_pdf)


#### Visualize trophic function composition ####
(trophic_barplots_sediments <- sediment_melt %>%
    mutate(trophic_mode = str_replace(trophic_mode, "saprotrophs", "heterotrophs")) %>%
    group_by(lake_id, upperpt, lowerpt, midpt_year, trophic_mode) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    mutate(interval = paste0(format(upperpt, nsmall = 1), " \U2012 ", format(lowerpt, nsmall = 1))) %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               nrow = 1) +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
                 fill = factor(trophic_mode, levels = rev(names(palette_function)))),
             stat = "identity") +
    scale_fill_manual(values = rev(palette_function), na.value = "#EEEEEE") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.1,
                  label = str_c(round(midpt_year, 0))),
              data = sediment_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "black") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0)) +
    coord_flip(ylim = c(0, 1), clip = "off") +
    labs(x = "Depth interval (cm)",
         y = "Relative sequence abundance",
         fill = "Sequence abundance (%)") +
    theme_bw() %+replace%
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_trophic_functions_barplots_sediments.pdf", trophic_barplots_sediments, width = 16, height = 10, units = "in", device = cairo_pdf)


#### Visualize sequence abundance by sample ####
# Total number of sequences per samples
sediment_nseqs <- sediment_melt %>%
  group_by(lake_id, upperpt, lowerpt, midpt_year) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup()

sediment_nseqs_means <- sediment_nseqs %>%
  group_by(lake_id) %>%
  summarize(mean_nseqs = mean(nseqs),
            sd_nseqs = sd(nseqs)) %>%
  ungroup()

(nseqs_barplots_sediments <- sediment_nseqs %>%
    mutate(interval = paste0(format(upperpt, nsmall = 1), " \U2012 ", format(lowerpt, nsmall = 1))) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = nseqs),
             fill = "grey85",
             stat = "identity") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -5500,
                  label = str_c(round(midpt_year, 0))),
              data = sediment_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "black") +
    geom_hline(aes(yintercept = mean_nseqs),
               data = sediment_nseqs_means,
               colour = "red", linetype = "dashed") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0),
                       label = scales::comma) +
    coord_flip(clip = "off") +
    labs(x = "Depth interval (cm)",
         y = "Number of sequences") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_sediment_nseqs_barplots.pdf", nseqs_barplots_sediments, width = 16, height = 10, units = "in", device = cairo_pdf)
