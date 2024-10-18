# ENA read manifest files

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

# Format sample IDs in read files
reads <- reads %>%
  mutate(sample_id = str_remove(read_sample, "ELA-")) %>%
  mutate(lake_id = str_remove(sample_id, "_.*")) %>%
  mutate(depth_unit = str_remove(sample_id, paste0(lake_id, "_"))) %>%
  separate(depth_unit, c("d1", "d2", "d3"), "-", remove = FALSE) %>%
  mutate(d3 = str_remove(d3, "cm")) %>%
  mutate(upperpt = case_when(d1 == "5" & d2 == "5" & d3 == "5" ~ 5,
                             d1 == "5" & d2 == "5" & d3 == "6" ~ 5.5,
                             d1 == d2 ~ as.numeric(d1),
                             TRUE ~ as.numeric(paste0(d1, ".", d2))),
         lowerpt = case_when(d1 == "5" & d2 == "5" & d3 == "5" ~ 5.5,
                             d1 == "5" & d2 == "5" & d3 == "6" ~ 6,
                             d1 == d2 ~ as.numeric(paste0(d2, ".", d3)),
                             TRUE ~ as.numeric(d3))) %>%
  select(-c(depth_unit, d1, d2, d3)) %>%
  mutate(depth_unit = str_c(upperpt, "-", lowerpt, "cm")) %>%
  mutate(sample_id = str_c(lake_id, "_", depth_unit)) %>%
  select(sample_id, read_sample, starts_with(c("fwd", "rev")))

metadata_format <- metadata %>%
  mutate(upperpt_str = case_when(grepl(".5", upperpt) ~ as.character(upperpt),
                                 TRUE ~ str_c(upperpt, ".0")),
         lowerpt_str = case_when(grepl(".5", lowerpt) ~ as.character(lowerpt),
                                 TRUE ~ str_c(lowerpt, ".0"))) %>%
  mutate(sample_id2 = str_c(lake_id, "_", upperpt_str, "-", lowerpt_str, "cm"))

reads_samples <- reads %>%
  left_join(metadata_format, by = "sample_id") %>%
  left_join(samples, by = c("sample_id2" = "sample_id"))


#### Write manifest files ####
# Define path for saving manifest files
path <- "output/ena/read_manifests/"

for (i in 1:nrow(reads_samples)) {
  reads_samples_tmp <- reads_samples[i,]
  
  manifest <- tibble(STUDY = "ERP131591",  #Study accession or unique name (alias)
                     SAMPLE = reads_samples_tmp$sample_accession,  #Sample accession or unique name (alias)
                     NAME = str_c("IISD-ELA ", reads_samples_tmp$lake_id, " ", reads_samples_tmp$upperpt_str, "-", reads_samples_tmp$lowerpt_str, " cm"),  #Unique experiment name
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
