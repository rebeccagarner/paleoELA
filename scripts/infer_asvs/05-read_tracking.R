# Track reads through ASV table creation pipeline

# Load libraries
library(tidyverse)


#### Import and format data ####
# Import read tracking table
read_tracking <- read_tsv("output/dada2/ela18s_readtracking.tsv", col_names = TRUE) %>%
  rename(sample_id = X1) %>%
  mutate(sample_id = str_remove(sample_id, "ELA-")) %>%
  mutate(sample_type = case_when(grepl("\\d\\d\\d\\d\\d\\d\\d\\d", sample_id) ~ "filter",
                                 TRUE ~ "sediment"))


### Visualize read tracking ####
# Plot number of sequences at each step of the pipeline
(read_tracking_lineplot <- read_tracking %>%
    mutate(dada_mean = (dada_f + dada_r)/2) %>%
    select(-dada_f, -dada_r, -final_perc_reads_retained) %>%
    pivot_longer(!c(sample_id, sample_type), names_to = "step", values_to = "nseqs") %>%
    ggplot() +
    geom_line(aes(x = factor(step, levels = c("dada2_input", "filtered", "dada_mean", "merged", "nonchim")),
                  y = nseqs,
                  colour = sample_type,
                  group = sample_id)) +
    labs(y = "Number of sequences",
         colour = "Sample type") +
    theme_bw() %+replace%
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/ela18s_read_tracking_lineplot.pdf", read_tracking_lineplot, width = 8, height = 6, units = "in")

# Plot histogram of number of samples vs. number of sequences
(read_tracking_histogram <- read_tracking %>%
    mutate(dada_mean = (dada_f + dada_r)/2) %>%
    select(-dada_f, -dada_r, -final_perc_reads_retained) %>%
    pivot_longer(!c(sample_id, sample_type), names_to = "step", values_to = "nseqs") %>%
    ggplot() +
    facet_wrap(~factor(step, levels = c("dada2_input", "filtered", "dada_mean", "merged", "nonchim")), ncol = 1) +
    geom_histogram(aes(x = nseqs,
                       y = ..count..,
                       fill = sample_type)) +
    labs(x = "Number of sequences",
         y = "Number of samples",
         fill = "Sample type") +
    theme_bw() %+replace%
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/ela18s_read_tracking_histogram.pdf", read_tracking_histogram, width = 4, height = 10, units = "in")

(read_tracking_pct_seqs_histogram <- read_tracking %>%
    ggplot() +
    geom_histogram(aes(x = final_perc_reads_retained,
                       y = ..count..,
                       fill = sample_type)) +
    labs(x = "Sequences retained (%)",
         y = "Number of samples",
         fill = "Sample type") +
    theme_bw() %+replace%
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
#ggsave("figures/ela18s_read_tracking_pct_seqs_histogram.pdf", read_tracking_pct_seqs_histogram, width = 5, height = 4, units = "in")
