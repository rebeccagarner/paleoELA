# Taxonomic composition

# Load libraries
library(tidyverse)
library(patchwork)

# Load palettes
source("scripts/00-palettes.R")

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Assess taxonomic group ASV/sequence compositions ####
# Calculate total n sequences per sample
sediment_nseqs <- sediment_melt %>%
  group_by(sample_id) %>%
  summarize(total_nseqs = sum(nseqs)) %>%
  ungroup()

# Calculate n ASVs assigned to taxonomic subdivisions
nasvs_subdivision <- sediment_melt %>%
  select(-sample_id, -nseqs) %>%
  distinct(asv_code, .keep_all = TRUE) %>%
  group_by(supergroup, division, subdivision) %>%
  tally(name = "n_asvs") %>%
  ungroup() %>%
  mutate(pct_asvs = n_asvs/length(unique(sediment_melt$asv_code)) * 100)
sum(nasvs_subdivision$n_asvs) == length(unique(sediment_melt$asv_code))  # Should evaluate to TRUE
sum(nasvs_subdivision$pct_asvs) == 100  # Should evaluate to TRUE

# Calculate n sequences assigned to taxonomic subdivisions
nseqs_subdivision <- sediment_melt %>%
  select(-sample_id) %>%
  group_by(supergroup, division, subdivision) %>%
  summarize(n_seqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = n_seqs/sum(sediment_melt$nseqs) * 100)
sum(nseqs_subdivision$n_seqs) == sum(sediment_melt$nseqs)  # Should evaluate to TRUE
sum(nseqs_subdivision$pct_seqs) == 100  # Should evaluate to TRUE

# Combine n ASVs and n sequences data for taxonomic divisions
nasvs_nseqs_subdivision <- nasvs_subdivision %>%
  left_join(nseqs_subdivision, by = c("supergroup", "division", "subdivision")) %>%
  mutate(pct_seqs = -pct_seqs) %>%
  pivot_longer(!c(supergroup, division, subdivision), names_to = c(".value", "var"), names_pattern = "(.*)_(.*)")

# Create back-to-back bar plot of n ASVs and n sequences assigned to taxonomic divisions
(nasvs_nseqs_subdivision_barplot <- nasvs_nseqs_subdivision %>%
    # mutate(supergroup = case_when(is.na(supergroup) ~ "Unclassified supergroup",
    #                               TRUE ~ supergroup)) %>%
    ggplot() +
    geom_bar(aes(x = division,
                 y = pct,
                 fill = subdivision),
             stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = palette_subdivision, na.value = "black") +
    scale_y_continuous(labels = abs) +
    labs(y = "Composition (%)",
         fill = "Taxonomy") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          legend.position = "right"))
#ggsave("figures/ela18s_nasvs_nseqs_subdivision_barplot.pdf", nasvs_nseqs_subdivision_barplot, "pdf", width = 8, height = 6, units = "in")


#### Target classical paleolimnological proxies ####
# Summarize diatoms diversity by lake and ASV
diatoms_melt %>%
  group_by(lake_id, class, order, family, genus, species, asv_code) %>%
  summarize(nseqs = sum(nseqs))

diatoms_nseqs <- diatoms_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup()

# Visualize diatoms genus diversity in sediment stratigraphies
(diatoms_genus_barplot <- diatoms_melt %>%
    dplyr::select(-relseqs) %>%
    left_join(diatoms_nseqs, by = "sample_id") %>%
    mutate(relseqs = nseqs/nseqs_total) %>%
    group_by(lake_id, upperpt, lowerpt, midpt_year, interval, genus) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               nrow = 1) +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(diatoms_melt$genus)[which(!is.na(unique(diatoms_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.17,
                  label = str_c(round(midpt_year, 0))),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    #scale_y_continuous(expand = c(0,0)) +
    coord_flip(ylim = c(0, 1), clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Relative sequence abundance",
         fill = "Taxonomy") +
    theme_bw() %+replace%
    theme(panel.spacing = unit(2.3, "lines"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_genus_barplots.pdf", diatoms_genus_barplot, "pdf", width = 16, height = 10, units = "in")

(diatoms_genus_nseqs_barplot <- diatoms_melt %>%
    mutate(interval = paste0(format(upperpt, nsmall = 1), " - ", format(lowerpt, nsmall = 1))) %>%
    ggplot() +
    facet_grid(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               scales = "free_x", space = "free_x") +
    geom_bar(aes(x = reorder(interval, -upperpt), y = nseqs,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(diatoms_melt$genus)[which(!is.na(unique(diatoms_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.17,
                  label = str_c(round(year, 0))),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(breaks = seq(0, max(diatoms_nseqs$nseqs_total), by = 500)) +
    coord_flip(clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Number of sequences",
         fill = "Taxonomy") +
    theme_bw() %+replace%
    theme(panel.spacing = unit(2.3, "lines"),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_genus_nseqs_barplots.pdf", diatoms_genus_nseqs_barplot, "pdf", width = 16, height = 10, units = "in")

# Visualize diatoms family diversity in sediment stratigraphies
(diatoms_family_barplot <- diatoms_melt %>%
  dplyr::select(-relseqs) %>%
  left_join(diatoms_nseqs, by = "sample_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  group_by(lake_id, upperpt, lowerpt, year, family) %>%
  summarize(relseqs = sum(relseqs)) %>%
  ungroup() %>%
  mutate(interval = paste0(format(upperpt, nsmall = 1), " - ", format(lowerpt, nsmall = 1))) %>%
  ggplot() +
  facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
             nrow = 1) +
  geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
               fill = family), stat = "identity") +
  scale_fill_manual(values = rainbow(length(unique(diatoms_melt$family)[which(!is.na(unique(diatoms_melt$family)))])), na.value = "black") +
  geom_text(aes(x = reorder(interval, -upperpt),
                y = -0.17,
                label = str_c(round(year, 0))),
            colour = "grey30") +
  scale_x_discrete(position = "top") +
  #scale_y_continuous(expand = c(0,0)) +
  coord_flip(ylim = c(0, 1), clip = "off") +
  labs(x = "Sediment depth (cm)",
       y = "Relative sequence abundance",
       fill = "Taxonomy") +
  theme_bw() %+replace%
  theme(panel.spacing = unit(2.3, "lines"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_diatoms_family_barplots.pdf", diatoms_family_barplot, "pdf", width = 16, height = 10, units = "in")

# Extract Chrysophyceae ASVs
chrysophyceae_melt <- sediment_melt %>%
  filter(class == "Chrysophyceae")

# Summarize Chrysophyceae diversity by lake and ASV
chrysophyceae_melt %>%
  group_by(lake_id, class, order, family, genus, species, asv_code) %>%
  summarize(nseqs = sum(nseqs))

chrysophyceae_nseqs <- chrysophyceae_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs)) %>%
  ungroup() %>%
  left_join(sediment_nseqs, by = "sample_id") %>%
  mutate(pct_seqs = nseqs_total/total_nseqs * 100) %>%
  dplyr::select(-total_nseqs)

# Visualize chrysophyceae genus diversity in sediment stratigraphies
(chrysophyceae_genus_barplot <- chrysophyceae_melt %>%
    dplyr::select(-relseqs) %>%
    left_join(chrysophyceae_nseqs, by = "sample_id") %>%
    mutate(relseqs = nseqs/nseqs_total) %>%
    group_by(lake_id, upperpt, lowerpt, year, genus) %>%
    summarize(relseqs = sum(relseqs)) %>%
    ungroup() %>%
    mutate(interval = paste0(format(upperpt, nsmall = 1), " - ", format(lowerpt, nsmall = 1))) %>%
    ggplot() +
    facet_wrap(~factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
               nrow = 1) +
    geom_bar(aes(x = reorder(interval, -upperpt), y = relseqs,
                 fill = genus), stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(chrysophyceae_melt$genus)[which(!is.na(unique(chrysophyceae_melt$genus)))])), na.value = "black") +
    geom_text(aes(x = reorder(interval, -upperpt),
                  y = -0.17,
                  label = str_c(round(year, 0))),
              colour = "grey30") +
    scale_x_discrete(position = "top") +
    #scale_y_continuous(expand = c(0,0)) +
    coord_flip(ylim = c(0, 1), clip = "off") +
    labs(x = "Sediment depth (cm)",
         y = "Relative sequence abundance",
         fill = "Taxonomy") +
    theme_bw() %+replace%
    theme(panel.spacing = unit(2.3, "lines"),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")))
#ggsave("figures/ela18s_chrysophyceae_genus_barplots.pdf", chrysophyceae_genus_barplot, "pdf", width = 16, height = 10, units = "in")


#### Assess unclassified fraction of assemblages ####
# Calculate unclassified percentage of sequences per sample
pct_seqs_unclassified <- sediment_melt %>%
  filter(is.na(subdivision)) %>%
  group_by(sample_id) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  left_join(sediment_nseqs, by = "sample_id") %>%
  mutate(pct_unclassified = nseqs/total_nseqs * 100) %>%
  mutate(lake_id = str_remove_all(sample_id, "_.*"))

# Plot comparison of unclassified percentage of sequences per sample
(pct_seqs_unclassified_boxplots <- pct_seqs_unclassified %>%
    ggplot() +
    geom_boxplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     y = pct_unclassified, fill = lake_id)) +
    scale_fill_manual(values = palette_lake) +
    labs(y = "Unclassified sequences (%)") +
    theme_bw() %+replace%
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank()))

# Plot comparison of unclassified percentage of sequences per supergroup per sample
(pct_seqs_unclassified_supergroups_boxplots <- sediment_melt %>%
    filter(is.na(subdivision)) %>%
    group_by(sample_id, supergroup, division) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ungroup() %>%
    left_join(sediment_nseqs, by = "sample_id") %>%
    mutate(pct_unclassified = nseqs/total_nseqs * 100) %>%
    mutate(lake_id = str_remove_all(sample_id, "_.*")) %>%
    ggplot() +
    geom_boxplot(aes(x = factor(lake_id, levels = c("L226N", "L226S", "L227", "L224", "L373")),
                     y = pct_unclassified, fill = supergroup)) +
    #scale_fill_manual(values = palette_supergroup) +
    labs(y = "Unclassified sequences (%)") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          panel.grid = element_blank()))

# Combine plots of unclassified percentage of sequences
(pct_seqs_unclassified_boxplots_all <- pct_seqs_unclassified_boxplots /
    pct_seqs_unclassified_supergroups_boxplots)
#ggsave("figures/ela18s_pct_sequences_unclassified.pdf", pct_seqs_unclassified_boxplots_all, "pdf", width = 6, height = 6, units = "in")

# Identify abundant unclassified ASVs
sediment_melt %>%
  filter(is.na(subdivision)) %>%
  group_by(asv_code, sequence, domain, supergroup, division) %>%
  summarize(nseqs = sum(nseqs)) %>%
  arrange(-nseqs)


#### Assess groups considered missing from sediments ####
sum(sediment_melt$relseqs)  # Should evaluate to n of samples
sediment_melt %>%
  group_by(subdivision) %>%
  summarize(nseqs = sum(nseqs)) %>%
  ungroup() %>%
  mutate(pct_seqs = nseqs/sum(sediment_melt$nseqs) * 100)
