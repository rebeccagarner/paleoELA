# Functional diversity

setwd("C:/Users/Gandalf/Dropbox/projects/ela18s/")
setwd("~/Dropbox/projects/ela18s/")

# Load libraries
library(tidyverse)

# Load palettes
source("scripts/00-palettes.R")


#### Import and format data ####
# Import melted ASV data
asv_melt <- read_tsv("output/melted/ela18s_microeuks_melt.tsv", col_names = TRUE,
                     col_types = cols(date = col_date()))

# Import pre-formatted FUNGuild database (downloaded 2023/04/21)
funguild <- read_tsv("data/function/funguild_db.tsv", col_names = TRUE) %>%
  dplyr::select(taxon, trophicMode, guild) %>%
  filter(str_count(taxon, " ") <= 1) %>%
  distinct() %>%
  arrange(taxon)


#### Link taxonomy and function ####
# Create taxonomy table
taxonomy <- asv_melt %>%
  distinct(asv_code, domain, supergroup, division, subdivision, class, order, family, genus, species)

# Create non-redundant taxonomy table
taxonomy_nr <- taxonomy %>%
  distinct(domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
  arrange(domain, supergroup, division, subdivision, class, order, family, genus, species)

# Link taxonomy and trophic function summarized from Adl et al. 2019 "Revisions to
# the Classification, Nomenclature, and Diversity of Eukaryotes"
taxonomy_function <- taxonomy_nr %>%
  mutate(trophic_mode = case_when(
    # Amoebozoa
    species == "Difflugia_bacillariarum" ~ "mixotrophs",  # To verify
    supergroup == "Amoebozoa" ~ "heterotrophs",
    
    # Apusomonada
    division == "Apusomonada" ~ "heterotrophs",
    
    # Breviatea
    division == "Breviatea" ~ "heterotrophs",
    
    # Opisthokonta
    genus == "Archaeorhizomyces" ~ "heterotrophs",
    genus == "Kionochaeta" ~ "saprotrophs",
    genus == "Teratosphaeria" ~ "parasites",
    order == "Gromochytriales" ~ "parasites",
    class == "Aphelidiomycota" ~ "parasites",
    class == "Blastocladiomycota" ~ "heterotrophs",
    class == "Chytridiomycota" ~ "saprotrophs/parasites",
    class == "Opisthosporidia" ~ "parasites",
    class == "Rozellomycota" ~ "parasites",
    subdivision == "Choanoflagellata" ~ "heterotrophs",
    division == "Opisthokonta" ~ "heterotrophs",
    
    supergroup == "Obazoa" ~ "heterotrophs",
    
    # Archaeplastida
    genus == "Helicosporidium" ~ "parasites",
    subdivision == "Chlorophyta" ~ "phototrophs",
    subdivision == "Streptophyta" ~ "phototrophs",
    supergroup == "Archaeplastida" ~ "phototrophs",
    
    # Cryptista
    genus == "Cryptomonas" ~ "mixotrophs",
    genus == "Geminigera" ~ "mixotrophs",
    genus == "Goniomonas" ~ "heterotrophs",
    genus == "Komma" ~ "mixotrophs",
    genus == "Plagioselmis" ~ "mixotrophs",
    division == "Kathablepharidacea" ~ "heterotrophs",
    
    # Haptista
    genus == "Diacronema"~ "mixotrophs",
    genus == "Chrysochromulina" ~ "mixotrophs",
    genus == "Prymnesium" ~ "mixotrophs",
    division == "Centroplasthelida" ~ "heterotrophs",
    
    # Alveolata
    genus == "Alphamonas" ~ "heterotrophs",
    genus == "Anurofeca" ~ "parasites",
    genus == "Apocalathium" ~ "phototrophs",
    genus == "Asulcocephalium" ~ "phototrophs",
    genus == "Biecheleria" ~ "heterotrophs/mixotrophs",
    genus == "Blastodinium" ~ "mixotrophs",
    genus == "Borghiella" ~ "phototrophs",
    genus == "Ceratium" ~ "phototrophs/mixotrophs/heterotrophs",
    genus == "Chromera" ~ "phototrophs",
    genus == "Chytriodinium" ~ "parasites",
    genus == "Cladocopium" ~ "phototrophs",
    genus == "Gymnodinium" ~ "phototrophs/mixotrophs",
    genus == "Gyrodinium" ~ "heterotrophs",
    genus == "Ichthyophthirius" ~ "parasites",
    genus == "Karenia" ~ "heterotrophs/mixotrophs",
    genus == "Nusuttodinium" ~ "heterotrophs/mixotrophs",
    genus == "Paramoebidium" ~ "parasites",
    genus == "Peridinium" ~ "phototrophs",
    genus == "Tetrahymena" ~ "parasites",
    family == "Bursariomorphida" ~ "heterotrophs",
    family == "Symbiodiniaceae" ~ "phototrophs",
    class == "Phyllopharyngea" ~ "parasites",
    subdivision == "Apicomplexa" ~ "parasites",
    subdivision == "Ciliophora" ~ "heterotrophs",
    subdivision == "Perkinsea" ~ "parasites",
    
    # Stramenopiles
    genus == "Acrispumella" ~ "heterotrophs",
    genus == "Dinobryon" ~ "mixotrophs",
    genus == "Epipyxis" ~ "phototrophs",
    genus == "Hydrurus" ~ "phototrophs",
    order == "Chromulinales" ~ "phototrophs",
    order == "Hyphochytriales" ~ "saprotrophs",
    order == "Paraphysomonadales" ~ "heterotrophs",
    order == "Synurales" ~ "phototrophs",
    class == "Bacillariophyceae" ~ "phototrophs",
    class == "Bolidophyceae" ~ "phototrophs",
    class == "Chrysophyceae" ~ "phototrophs",
    class == "Coscinodiscophyceae" ~ "phototrophs",
    class == "Dictyochophyceae" ~ "phototrophs",
    class == "Eustigmatophyceae" ~ "phototrophs",
    class == "Mediophyceae" ~ "phototrophs",
    class == "Peronosporomycetes" ~ "parasites",
    class == "Pirsoniales" ~ "parasites",
    class == "Phaeothamniophyceae" ~ "phototrophs",
    class == "Raphidophyceae" ~ "phototrophs",
    class == "Xanthophyceae" ~ "phototrophs",
    subdivision == "Bigyra" ~ "heterotrophs",
    
    # Rhizaria
    genus == "Aquavolon" ~ "heterotrophs",
    order == "Cercomonadida" ~ "heterotrophs",
    order == "Euglyphida" ~ "heterotrophs",
    order == "Glissomonadida" ~ "heterotrophs",
    order == "Novel-clade-10" ~ "heterotrophs",
    order == "Pansomonadida" ~ "heterotrophs",
    order == "Paracercomonadida" ~ "heterotrophs",
    order == "Plasmodiophorida" ~ "parasites",
    order == "Thaumatomonadida" ~ "heterotrophs",
    order == "Vampyrellida" ~ "heterotrophs",
    
    # Telonemia
    division == "Telonemia" ~ "heterotrophs",
    
    # Excavata
    genus == "Hexamita" ~ "parasites",
    order == "Aphagea" ~ "heterotrophs",
    order == "Euglenophyceae" ~ "phototrophs",
    order == "Hibberdiales" ~ "phototrophs",
    order == "Trypanosomatida" ~ "parasites",
    class == "Heterolobosea" ~ "heterotrophs",
    class == "Kinetoplastea" ~ "heterotrophs",
    division == "Metamonada" ~ "heterotrophs",
    
    # CRuMs
    division == "Rigifilida" ~ "heterotrophs",
    
    # Other
    division == "Ancyromonadida" ~ "heterotrophs"
  ))

# Annotate fungi at species or genus level based on FUNGuild database
taxonomy_function$guild <- NA_character_
for (i in 1:nrow(taxonomy_function)) {
  
  if (taxonomy_function$subdivision[i] == "Fungi" & !is.na(taxonomy_function$genus[i])) {
    taxonomy_tmp <- taxonomy_function[i,] %>%
      dplyr::select(genus, species) %>%
      mutate(species = str_replace_all(species, "_", " ")) %>%
      mutate(species = case_when(grepl(" sp\\.$", species) ~ NA_character_,
                                 TRUE ~ species))
    
    species_tmp <- taxonomy_function$species[i]
    genus_tmp <- taxonomy_function$genus[i]
    
    print(paste0("Now on row: ", i, ", genus: ", genus_tmp, ", species: ", species_tmp))
    
    if (species_tmp %in% funguild$taxon) {
      taxonomy_function$trophic_mode[i] <- funguild$trophicMode[which(funguild$taxon == species_tmp)]
      taxonomy_function$guild[i] <- funguild$guild[which(funguild$taxon == species_tmp)]
    } else if (genus_tmp %in% funguild$taxon) {
      taxonomy_function$trophic_mode[i] <- funguild$trophicMode[which(funguild$taxon == genus_tmp)]
      taxonomy_function$guild[i] <- funguild$guild[which(funguild$taxon == genus_tmp)]
    }
    
  }
}

# Harmonize FUNGuild annotations with Adl et al. 
taxonomy_function <- taxonomy_function %>%
  mutate(trophic_mode = str_trim(trophic_mode, "both")) %>%
  mutate(trophic_mode = str_to_lower(trophic_mode)) %>%
  mutate(trophic_mode = str_replace_all(trophic_mode, "-", "s/")) %>%
  mutate(trophic_mode = case_when(!grepl("s$", trophic_mode) ~ str_replace(trophic_mode, "$", "s"),
                                  TRUE ~ trophic_mode))


#### Assess functional diversity ####
# Join trophic function information with melted ASV data
asv_melt_function <- asv_melt %>%
  left_join(taxonomy_function %>%
              dplyr::select(-guild), by = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"))

# Write melted ASV and function data to file
# asv_melt_function %>%
#   write_tsv("output/melted/ela18s_microeuks_melt_function.tsv", col_names = TRUE)

# Evaluate ASV data without functional annotations (in priority order of relative abundance)
(not_annotated <- asv_melt_function %>%
    filter(is.na(trophic_mode)) %>%
    group_by(domain, supergroup, division, subdivision, class, order, family, genus, species) %>%
    summarize(nseqs = sum(nseqs)) %>%
    ungroup() %>%
    arrange(-nseqs))

# Assess functional diversity assignments in sediment samples
sediment_melt_function <- asv_melt_function %>%
  filter(sample_type == "sediment")

length(unique(sediment_melt_function$sample_id))
length(unique(sediment_melt_function$asv_code))
sum(sediment_melt_function$nseqs)

sediment_melt_function %>%
  mutate(trophic_assigned = case_when(!is.na(trophic_mode) ~ TRUE,
                                      TRUE ~ FALSE)) %>%
  group_by(trophic_assigned) %>%
  summarize(nasvs = length(unique(asv_code)),
            nseqs = sum(nseqs))

