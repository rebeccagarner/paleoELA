# R script for changing raw fastq file names
rm(list = ls())
(myfiles <- list.files(pattern = "fastq.gz"))
gsub(".*\\.L", "L", myfiles)
file.rename(from = myfiles, to = gsub(".*\\.L", "L", myfiles))
list.files()
