# Algal taxa composition

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Assess trophic functional group ASV/sequence compositions ####
# Calculate total n sequences per sample
sediment_nseqs <- sediment_melt %>%
  group_by(sample_id) %>%
  summarize(total_nseqs = sum(nseqs)) %>%
  ungroup()

# Calculate n ASVs assigned to trophic functional groups
nasvs_trophic <- sediment_melt %>%
  select(-sample_id, -nseqs) %>%
  distinct(asv_code, .keep_all = TRUE) %>%
  group_by(trophic_mode) %>%
  tally(name = "n_asvs") %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/length(unique(sediment_melt$asv_code)) * 100)
sum(nasvs_trophic$n_asvs) == length(unique(sediment_melt$asv_code))  # Should evaluate to TRUE
sum(nasvs_trophic$pct_asvs) == 100  # Should evaluate to TRUE

# Calculate n sequences assigned to trophic functional groups
nseqs_trophic <- sediment_melt %>%
  select(-sample_id) %>%
  group_by(trophic_mode) %>%
  summarize(n_seqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = n_seqs/sum(sediment_melt$nseqs) * 100)
sum(nseqs_trophic$n_seqs) == sum(sediment_melt$nseqs)  # Should evaluate to TRUE
sum(nseqs_trophic$pct_seqs) == 100  # Should evaluate to TRUE

# Combine n ASVs and n sequences data for trophic functional groups
nasvs_nseqs_trophic <- nasvs_trophic %>%
  left_join(nseqs_trophic, by = "trophic_mode") %>%
  mutate(pct_seqs = -pct_seqs) %>%
  pivot_longer(!trophic_mode, names_to = c(".value", "var"), names_pattern = "(.*)_(.*)")

# Create back-to-back bar plot of n ASVs and n sequences assigned to taxonomic divisions
(nasvs_nseqs_trophic_barplot <- nasvs_nseqs_trophic %>%
    ggplot() +
    geom_bar(aes(x = factor(trophic_mode, levels = rev(names(palette_function))),
                 y = pct,
                 fill = factor(trophic_mode, levels = rev(names(palette_function)))),
             stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = palette_function, na.value = "black") +
    scale_y_continuous(labels = abs) +
    labs(y = "Sequences (%) | ASVs (%)",
         fill = "Trophic group") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right"))
#ggsave("figures/ela18s_nasvs_nseqs_trophic_functions_barplot.pdf", nasvs_nseqs_trophic_barplot, width = 8, height = 6, units = "in", device = cairo_pdf)


#### Examine phototroph diversity ####
# Summarize phototroph diversity by lake and ASV
phototrophs_melt %>%
  group_by(lake_id, supergroup, division, subdivision, class, order, family, genus) %>%
  summarize(nseqs = sum(nseqs))

phototrophs_melt %>%
  group_by(lake_id, supergroup, division, subdivision, class) %>%
  summarize(nseqs = sum(nseqs))

phototrophs_nseqs <- phototrophs_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

# Order phototrophic taxa by decreasing total relative sequence abundance
phototroph_order <- phototrophs_melt %>%
  group_by(supergroup, division, subdivision, class) %>%
  summarize(relseqs = sum(relseqs)) %>%
  ungroup() %>%
  arrange(-relseqs) %>%
  filter(!is.na(class))

phototroph_supergroup_order <- phototroph_order %>%
  distinct(supergroup) %>%
  mutate(supergroup_order = row_number())

phototroph_division_order <- phototroph_order %>%
  distinct(division) %>%
  mutate(division_order = row_number())

phototroph_subdivision_order <- phototroph_order %>%
  distinct(subdivision) %>%
  mutate(subdivision_order = row_number())

phototroph_class_order <- phototroph_order %>%
  distinct(class) %>%
  mutate(class_order = row_number())

phototroph_order <- phototroph_order %>%
  left_join(phototroph_supergroup_order, by = "supergroup") %>%
  left_join(phototroph_division_order, by = "division") %>%
  left_join(phototroph_subdivision_order, by = "subdivision") %>%
  left_join(phototroph_class_order, by = "class") %>%
  arrange(supergroup_order, division_order, subdivision_order, class_order)

phototrophs_top <- phototroph_order %>%
  filter(relseqs >= 0.01)

# Plot phototroph stratigraphies
(phototrophs_barplots <- samples_unique %>%
    full_join(phototrophs_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    left_join(phototrophs_nseqs, by = "sample_id") %>%
    mutate(relseqs = nseqs/nseqs_total) %>%
    mutate(class = case_when(class %in% phototrophs_top$class ~ class,
                             TRUE ~ "Other")) %>%
    group_by(lake_id, upperpt, lowerpt, interval, midpt_year, class) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               nrow = 1) +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
                 fill = factor(class, levels = phototroph_order$class)),
             stat = "identity") +
    scale_fill_manual(values = palette_phototroph, na.value = "black") +
    #scale_fill_manual(values = rainbow(length(phototrophs_top$class)), na.value = "black") +
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
         fill = "Class") +
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
#ggsave("figures/ela18s_phototroph_class_barplots.pdf", phototrophs_barplots, width = 16, height = 10, units = "in", device = cairo_pdf)


#### Examine diatom diversity ####
# Summarize diatom diversity by lake and ASV
diatoms_melt %>%
  group_by(lake_id, class, order, family, genus, species, asv_code) %>%
  summarize(nseqs = sum(nseqs))

diatoms_nseqs <- diatoms_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

# Visualize diatom genus diversity in sediment stratigraphies
(diatoms_genus_nseqs_barplot <- samples_unique %>%
    full_join(diatoms_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = nseqs,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(diatoms_melt$genus)[which(!is.na(unique(diatoms_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -170,
                  label = str_c(round(midpt_year, 0))),
              data = diatoms_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, max(diatoms_nseqs$nseqs_total), by = 500),
                       label = scales::comma) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Number of sequences",
         fill = "Diatom genus") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_genus_nseqs_barplots.pdf", diatoms_genus_nseqs_barplot, width = 16, height = 10, units = "in", device = cairo_pdf)

(diatoms_genus_relseqs_barplot <- samples_unique %>%
    full_join(diatoms_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs * 100,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(diatoms_melt$genus)[which(!is.na(unique(diatoms_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.4,
                  label = str_c(round(midpt_year, 0))),
              data = diatoms_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, 100, by = 1)) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Sequences (%)",
         fill = "Diatom genus") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_genus_relseqs_barplots.pdf", diatoms_genus_relseqs_barplot, width = 16, height = 10, units = "in", device = cairo_pdf)

# Visualize diatom family diversity in sediment stratigraphies
(diatoms_family_relseqs_barplot <- samples_unique %>%
    full_join(diatoms_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs * 100,
                 fill = family), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(diatoms_melt$family)[which(!is.na(unique(diatoms_melt$family)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.4,
                  label = str_c(round(midpt_year, 0))),
              data = diatoms_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, 100, by = 1)) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Sequences (%)",
         fill = "diatom family") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_family_relseqs_barplots.pdf", diatoms_family_relseqs_barplot, width = 16, height = 10, units = "in", device = cairo_pdf)


#### Examine chrysophyte diversity ####
# Summarize Chrysophyceae diversity by lake and ASV
chrysophyceae_melt %>%
  group_by(lake_id, class, order, family, genus, species, asv_code) %>%
  summarize(nseqs = sum(nseqs))

chrysophyceae_nseqs <- chrysophyceae_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

# Visualize Chrysophyceae genus diversity in sediment stratigraphies
(chrysophyceae_genus_nseqs_barplot <- samples_unique %>%
    full_join(chrysophyceae_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = nseqs,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(chrysophyceae_melt$genus)[which(!is.na(unique(chrysophyceae_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -6000,
                  label = str_c(round(midpt_year, 0))),
              data = chrysophyceae_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, max(chrysophyceae_nseqs$nseqs_total), by = 10000),
                       label = scales::comma) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Number of sequences",
         fill = "Chrysophyceae genus") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_chrysophyceae_genus_nseqs_barplots.pdf", chrysophyceae_genus_nseqs_barplot, width = 16, height = 10, units = "in", device = cairo_pdf)

(chrysophyceae_genus_relseqs_barplot <- samples_unique %>%
    full_join(chrysophyceae_melt, by = c("sample_id", "lake_id", "interval", "upperpt")) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs * 100,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(chrysophyceae_melt$genus)[which(!is.na(unique(chrysophyceae_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -10,
                  label = str_c(round(midpt_year, 0))),
              data = chrysophyceae_melt %>%
                distinct(lake_id, interval, upperpt, midpt_year),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Sequences (%)",
         fill = "Chrysophyceae genus") +
    theme_bw() %+replace%
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(2.3, "lines"),
          axis.line = element_line(colour = "black"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_chrysophyceae_genus_relseqs_barplots.pdf", chrysophyceae_genus_relseqs_barplot, width = 16, height = 10, units = "in", device = cairo_pdf)


#### Examine mixotroph composition ####
# Recalculate relative sequences
phototrophs_nseqs <- phototrophs_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

mixotrophs_melt <- phototrophs_melt %>%
  dplyr::select(-relseqs) %>%
  left_join(phototrophs_nseqs, by = "sample_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  select(-nseqs_total) %>%
  group_by(sample_id, lake_id,
           upperpt, lowerpt, midpt,
           midpt_year, upperpt_year, lowerpt_year,
           layer_order, interval,
           trophic_mode) %>%
  summarize(relseqs = sum(relseqs)) %>%
  ungroup() %>%
  filter(trophic_mode == "mixotrophs" | trophic_mode == "phototrophs/mixotrophs")

mixotrophs_melt %>%
  ggplot() +
  facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
             scales = "free_x", space = "free_x") +
  geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs * 100,
               fill = trophic_mode), stat = "identity") +
  geom_text(aes(x = reorder(interval, -upperpt),
                y = -1,
                label = str_c(round(midpt_year, 0))),
            data = chrysophyceae_melt %>%
              distinct(lake_id, interval, upperpt, midpt_year),
            colour = "grey30") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  coord_flip(clip = "off") +
  labs(x = "Sediment depth (cm)",
       y = "Sequences (%)",
       fill = "Chrysophyceae genus") +
  theme_bw() %+replace%
  theme(axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(2.3, "lines"),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"))


#### Count algal ASVs and sequences ####
length(unique(phototrophs_melt$asv_code))  # n ASVs
sum(phototrophs_melt$nseqs)  # n sequences
