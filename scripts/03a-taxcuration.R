# Taxonomy curation

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)
library(seqinr)
library(ggpubr)
library(ape)
library(ggtree)
library(lubridate)

# Load palettes
source("scripts/00-palettes.R")


#### Load and format ASV data ####
# Import melted ASV sequence count and taxonomy data
asv_melt <- read_tsv("output/melted/ela18s_all_melt.tsv", col_names = TRUE,
                     col_types = cols(date = col_date()))
length(unique(asv_melt$asv_code))

# Format taxonomy information
taxonomy <- asv_melt %>%
  distinct(asv_code, sequence, domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
  arrange(asv_code)


#### Visualize tree taxonomy ####
# # Import tree
# tree_all <- read.tree("output/trees/ela18s_rmseqs_trimprimers_gtr.tree")
# 
# # Join tree and taxonomy
# ggtree_all <- ggtree(tree_all, layout = "rectangular", size = 0.05)
# 
# ggtree_all$data <- ggtree_all$data %>%
#   left_join(taxonomy, by = c("label" = "asv_code"))
# 
# # Visualize tree with taxonomy labels
# (tree_taxonomy <- ggtree_all +
#     geom_tippoint(data = ggtree_all$data,
#                   aes(colour = division), size = 0.5) +
#     geom_tiplab(data = ggtree_all$data,
#                 aes(label = paste(supergroup, division, class, order, family, genus, sep = "|")), size = 0.2) +
#     scale_colour_manual(values = palette_division, na.value = "black") +
#     theme(legend.position = "none"))
# #ggsave("figures/ela18s_all_gtr.pdf", tree_taxonomy, width = 12, height = 100, units = "in", limitsize = FALSE)
# 
# (tree_taxonomy_macrobes <- ggtree_all +
#     geom_tippoint(data = ggtree_all$data %>%
#                     filter(division %in% c("Metazoa", "Fungi") | class == "Embryophyceae"),
#                   aes(colour = division), size = 0.5) +
#     geom_tiplab(data = ggtree_all$data,
#                 aes(label = paste(supergroup, division, class, order, family, genus, sep = "|")), size = 0.2) +
#     scale_colour_manual(values = palette_division, na.value = "black") +
#     theme(legend.position = "none"))
# #ggsave("figures/ela18s_all_gtr_macrobes.pdf", tree_taxonomy_macrobes, width = 12, height = 100, units = "in", limitsize = FALSE)


#### Extract microeukaryotes and fungi ####
# Separate Metazoa and land plants from other eukaryotes (protists and fungi)
# Identify ASVs assigned to Metazoa
metazoa_asvs <- taxonomy %>%
  filter(subdivision == "Metazoa") %>%
  pull(asv_code)

# Identify ASVs assigned to land plants
plant_asvs <- taxonomy %>%
  filter(class == "Embryophyceae") %>%
  pull(asv_code)

# Identify ASVs not assigned at the domain level
unassigned_domain_asvs <- taxonomy %>%
  filter(is.na(domain)) %>%
  pull(asv_code)

# Remove ASVs assigned to non-target taxa or not assigned at the domain level
microeuks_melt <- asv_melt %>%
  filter(!asv_code %in% metazoa_asvs & !asv_code %in% plant_asvs & !asv_code %in% unassigned_domain_asvs)
length(unique(microeuks_melt$asv_code))

# Calculate relative abundance of sequences
samples_nseqs <- microeuks_melt %>%
  group_by(sample_id) %>%
  summarize(nseqs_total = sum(nseqs))

microeuks_melt <- microeuks_melt %>%
  left_join(samples_nseqs, by = "sample_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  select(-nseqs_total)
sum(microeuks_melt$relseqs) == length(unique(microeuks_melt$sample_id))  # Should evaluate to TRUE

# Write melted microeukaryotes dataset to file
# microeuks_melt %>%
#   write_tsv("output/melted/ela18s_microeuks_melt.tsv", col_names = TRUE)

# Write microeukaryote sequences to a fasta file
# Remove primer sequences
microeuks_seqs <- microeuks_melt %>%
  distinct(asv_code, sequence)
# write.fasta(sequences = as.list(microeuks_seqs$sequence),
#             names = microeuks_seqs$asv_code,
#             file.out = "output/fasta/ela18s_microeuks.fasta")

# # Import tree
# tree_microeuks <- read.tree("output/trees/ela18s_microeuks_gtr.tree")
# 
# # Join tree and taxonomy
# ggtree_microeuks <- ggtree(tree_microeuks, layout = "rectangular", size = 0.05)
# 
# ggtree_microeuks$data <- ggtree_microeuks$data %>%
#   left_join(taxonomy, by = c("label" = "asv_code"))
# 
# # Visualize tree with taxonomy labels
# (tree_taxonomy_microeuks <- ggtree_microeuks +
#     geom_tippoint(data = ggtree_microeuks$data,
#                   aes(colour = division), size = 0.5) +
#     geom_tiplab(data = ggtree_microeuks$data,
#                 aes(label = paste(supergroup, division, class, order, family, genus, sep = "|")), size = 0.2) +
#     scale_colour_manual(values = palette_division, na.value = "black") +
#     theme(legend.position = "none"))
# #ggsave("figures/ela18s_microeuks_gtr.pdf", tree_taxonomy_microeuks, width = 12, height = 100, units = "in", limitsize = FALSE)


#### Tax group ASV/sequence compositions (pie charts) ####
# Pie chart n ASVs assigned to metazoa/fungi/plants/other
(taxgroup_asvs_piechart <- asv_melt %>%
   select(-sample_id, -nseqs) %>%
   distinct(asv_code, .keep_all = TRUE) %>%
   mutate(tax_group = case_when(subdivision == "Metazoa" ~ "Metazoa",
                                subdivision == "Fungi" ~ "Fungi",
                                class == "Embryophyceae" ~ "Plants",
                                TRUE ~ "Other (protists)")) %>%
   group_by(tax_group) %>%
   tally(name = "nasvs") %>%
   ggplot(aes(x = "", y = nasvs, fill = tax_group)) +
   geom_bar(stat = "identity") +
   coord_polar("y") +
   geom_text(aes(label = str_c(nasvs, " (", round(nasvs/(length(unique(asv_melt$asv_code)))*100, 1), "%)"), x = 1),
             position = position_stack(vjust = 0.5),
             size = 4) +
   labs(fill = "Taxonomy") +
   scale_fill_manual(values = palette_taxgroup) +
   theme_void())

# Pie chart n sequences assigned to metazoa/fungi/plants/other
(taxgroup_seqs_piechart <- asv_melt %>%
    select(-sample_id) %>%
    mutate(tax_group = case_when(subdivision == "Metazoa" ~ "Metazoa",
                                 subdivision == "Fungi" ~ "Fungi",
                                 class == "Embryophyceae" ~ "Plants",
                                 TRUE ~ "Other (protists)")) %>%
    group_by(tax_group) %>%
    summarize(sum_nseqs = sum(nseqs)) %>%
    ggplot(aes(x = "", y = sum_nseqs, fill = tax_group)) +
    geom_bar(stat = "identity") +
    coord_polar("y") +
    geom_text(aes(label = str_c(sum_nseqs, " (", round(sum_nseqs/sum(asv_melt$nseqs)*100, 1), "%)"), x = 1),
              position = position_stack(vjust = 0.5),
              size = 4) +
    scale_fill_manual(values = palette_taxgroup) +
    theme_void())

# Combine pie charts n ASVs and n sequences assigned to metazoa/fungi/plants/other
(taxgroup_piecharts <- ggarrange(taxgroup_asvs_piechart, taxgroup_seqs_piechart,
                                 labels = c("ASVs", "Sequences"),
                                 common.legend = TRUE, legend = "right"))
#ggsave("figures/ela18s_taxgroup_piecharts.pdf", taxgroup_piecharts, width = 10, height = 4, units = "in")
