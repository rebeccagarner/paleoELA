# ENA read manifest files (filters)

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load sediment data
source("scripts/03e-sediment_data.R")


#### Import and format data ####
# Import and format sample accessions
samples <- read_tsv("output/ena/Webin-accessions-2023-11-24T19_23_34.681Z.txt", col_names = TRUE) %>%
  rename(sample_accession = ACCESSION,
         sample_id = ALIAS) %>%
  select(sample_id, sample_accession) %>%
  filter(grepl("^ERS", sample_accession))

# Import and format read files
reads <- read_tsv("data/sequencing/ela18s_read_files_filters.txt", col_names = "file_name") %>%
  mutate(read_sample = str_remove(file_name, "_R[:digit:].fastq.gz.*")) %>%
  mutate(file_type = case_when(grepl("_R1.fastq.gz$", file_name) ~ "fwd_reads",
                               grepl("_R2.fastq.gz$", file_name) ~ "rev_reads",
                               grepl("_R1.fastq.gz.md5$", file_name) ~ "fwd_reads_md5",
                               grepl("_R2.fastq.gz.md5$", file_name) ~ "rev_reads_md5")) %>%
  pivot_wider(names_from = file_type, values_from = file_name)

metadata_format <- metadata_format <- asv_melt %>%
  filter(sample_type == "filter") %>%
  distinct(sample_id, lake_id, date, depth_unit) %>%
  mutate(sample_id = str_c("ELA-", str_replace_all(sample_id, "_", "-"))) %>%
  mutate(depth = as.numeric(str_remove(depth_unit, "m$")))

reads_samples <- reads %>%
  left_join(metadata_format, by = c("read_sample" = "sample_id")) %>%
  left_join(samples, by = c("read_sample" = "sample_id"))


#### Write manifest files ####
# Define path for saving manifest files
path <- "output/ena/read_manifests/filters/"

for (i in 1:nrow(reads_samples)) {
  reads_samples_tmp <- reads_samples[i,]
  
  manifest <- tibble(STUDY = "ERP155436",  #Study accession or unique name (alias)
                     SAMPLE = reads_samples_tmp$sample_accession,  #Sample accession or unique name (alias)
                     NAME = str_c("IISD-ELA ", reads_samples_tmp$lake_id, " ", reads_samples_tmp$depth_unit, " sampled on ", reads_samples_tmp$date),  #Unique experiment name
                     PLATFORM = "ILLUMINA",  #See permitted values. Not needed if INSTRUMENT is provided.
                     INSTRUMENT = "Illumina MiSeq",  #See permitted values
                     INSERT_SIZE = "225",  #Insert size for paired reads
                     LIBRARY_NAME = reads_samples_tmp$read_sample,  #Library name (optional)
                     LIBRARY_SOURCE = "METAGENOMIC",  #See permitted values
                     LIBRARY_SELECTION = "PCR",  #See permitted values
                     LIBRARY_STRATEGY = "AMPLICON",  #See permitted values
                     DESCRIPTION = "18S rRNA gene V7 region amplicons",  #free text library description (optional)
                     FASTQ_fwd = reads_samples_tmp$fwd_reads,
                     FASTQ_rev = reads_samples_tmp$rev_reads) %>%
    pivot_longer(everything()) %>%
    mutate(name = case_when(grepl("^FASTQ", name) ~ "FASTQ",
                            TRUE ~ name))
  
  manifest %>%
    write_tsv(paste0(path, "manifest-", reads_samples_tmp$read_sample, ".txt"), col_names = FALSE)
}
