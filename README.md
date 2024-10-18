# paleoELA
Scripts associated with the manuscript

---

Garner RE, Taranu ZE, Higgins SN, Paterson MJ, Gregory-Eaves I, Walsh DA. Eutrophication and warming drive algal community shifts in synchronized time series of experimental lakes.

 ---

## scripts/

- *infer_asvs/*
    - **00-workflow_18s.txt** Overview of workflow for 18S rRNA amplicon sequence variant (ASV) inference.
    - **01-rename_fastq.R** Modify raw fastq file names.
    - **02-create_sample_list.R** Create list of samples.
    - **03-cutadapt.sh** Trim primers and sequencing adapters in Cutadapt.
    - **04-dada2.R** Run DADA2 pipeline on Cutadapt-trimmed reads.
    - **05-read_tracking.R** Track read loss through the ASV table creation pipeline.
- **00-clean_pr2.R** Clean and dereplicate Protist Ribosomal Reference database (PR2) v. 5.0.0 eukaryotic taxonomy.
- **00-ena_read_manifest.R** Generate European Nucleotide Archive (ENA) read manifest files.
- **00-ena_sample_template.R** Populate template for ENA sample submission.
- **00-00-ena_webin_reads.sh** Submit read manifest files via ENA Webin Client.
- **00-palettes.R** Define palettes for figures.
- **00a-dating.R** Consolidate sediment dating information.
- **00b_map.R** Map region and lakes.
- **02-curate_asvs.R** Curate ASV dataset.
- **03a-taxcuration.R** Retain ASVs assigned to protists and fungi.
- **03b-alignment_tree.txt** Align ASVs to Silva reference database in SINA.
- **03d-trophic_function.R** Linked functional assignments to ASVs by taxon.
- **03e-sediment_data.R** Parse sediment ASV metadata including sediment intervals and age-dates.
- **04-stratigraphy.R** Evaluate sediment ASV assemblage taxonomic, functional, and sequence composition.
- **05-pca.R** Perform principal component analyses on sediment ASV assemblages.
- **07-algae.R** Assess algal taxonomic composition of sediment samples.
- **07-nasvs_nseqs_taxonomy.R** Assess sediment taxonomic group ASV/sequence compositions.
- **08-ela_env_data.R** Define functions to summarize Experimental Lakes Area (ELA) environmental data.
- **08a-ela_climate_data.R** Summarize ELA climate data.
- **08d-ela_chemistry_data.R** Summarize ELA nutrient data.
- **09-eccc_climate.R** Summarize Environment Canada historical climate data.
- **09-ela_phytoplankton_data.R** Format ELA algal monitoring data.
- **10-ela_phytoplankton_biomass.R** Summarize ELA algal biomass data.
- **10-ela_phytoplankton_species.R** Summarize ELA algal species cell count, density, and biomass data.
- **11-ela_monitoring.R** Consolidate ELA monitoring data.
- **12-pca_phytoplankton.R** Perform principal component analyses on monitoring assemblages.
- **18-rv.R** Calculat RV multivariate correlations between monitoring and sediment assemblages.
- **21-gam_devexpl_plots.R** Plot generalized additive model (GAM) results.
- **21-gam_functions_phytoplankton.R** Define functions to compute GAMs on monitoring assemblages.
- **21-gam_functions_seddna.R** Define functions to compute GAMs on sediment assemblages.
- **21a-gam_phytoplankton.R** Compute GAMs on monitoring assemblages.
- **21b-gam_seddna.R** Compute GAMs on sediment assemblages.
- **21c-gam_seddna_monitor.R** Compute GAMs on sediment assemblages truncated to monitoring time frames.
- **23-plot_nutriclim.R** Plot nutrient and air temperature trends.
- **24-pc_loadings.R** Extract top principal component taxon loadings.
- **25-ela_taxonomy_plots.R** Plot algal community time series.
- **26-pc1.R** Plot first principal component axis.
