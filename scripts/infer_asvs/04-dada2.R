# Run DADA2 pipeline on Cutadapt-trimmed reads
# For ELA full-core sediments and 2018 GF/C filters 18S rRNA gene amplicon sequencing data

library(dada2)
packageVersion("dada2")

setwd("~/projects/ela18s/")


# Identify Cutadapt-trimmed samples and read files
(samples <- scan("samples", what = "character"))

(forward_reads <- paste0(samples, "_R1_trimmed.fastq.gz"))
(reverse_reads <- paste0(samples, "_R2_trimmed.fastq.gz"))

# Prepare to filter reads
filtered_forward_reads <- paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fastq.gz")

# Define forward and reverse read lengths (based on quality plots)
forward_read_length <- 200
reverse_read_length <- 210

# Filter reads
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE = c(2, 2), rm.phix = TRUE, minLen = 160, truncLen = c(forward_read_length, reverse_read_length), compress = TRUE, multithread = TRUE)

# Learn error rates
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)
#plotErrors(err_forward_reads, nominalQ = TRUE)
#plotErrors(err_reverse_reads, nominalQ = TRUE)

# Dereplicate reads
derep_forward <- derepFastq(filtered_forward_reads, verbose = TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = TRUE)
names(derep_reverse) <- samples

# Infer ASVs
dada_forward <- dada(derep_forward, err = err_forward_reads, pool = TRUE, multithread = TRUE)
dada_reverse <- dada(derep_reverse, err = err_reverse_reads, pool = TRUE, multithread = TRUE)

# Merge forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang = TRUE, minOverlap = 30)

# Create sequence table
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)

# Remove chimerae
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
sum(seqtab.nochim)/sum(seqtab)
save(seqtab.nochim, file = "ela18s_seqtab_nochim.rda")

# Track read loss through the pipeline
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names = samples, dada2_input = filtered_out[,1], filtered = filtered_out[,2], dada_f = sapply(dada_forward, getN), dada_r = sapply(dada_reverse, getN), merged = sapply(merged_amplicons, getN), nonchim = rowSums(seqtab.nochim), final_perc_reads_retained = round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "ela18s_readtracking.tsv", quote = FALSE, sep = "\t", col.names = NA)

# Assign taxonomy against PR2 v. 4.12.0
# taxonomy <- assignTaxonomy(seqtab.nochim, "~/projects/lp2017_amplicons/ref_dbs/pr2_version_4.12.0_18S_dada2.fasta.gz", taxLevels = c("kingdom","supergroup","division","class","order","family","genus","species"), multithread = TRUE)
# save(taxonomy, file = "ela18s_taxonomy_pr2v412.rda")

# Assign taxonomy against PR2 v. 5.0.0
taxonomy <- assignTaxonomy(seqtab.nochim, "~/Desktop/ref_dbs/dada2_taxonomy/pr2_version_5.0.0_SSU_dada2.fasta.gz", taxLevels = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"), tryRC = TRUE, multithread = TRUE)
save(taxonomy, file = "ela18s_taxonomy_pr2v5_tryrc.rda")

# Assign ASV codes to sequences
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# Create ASV fasta file
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ela18s_asvs.fasta")
