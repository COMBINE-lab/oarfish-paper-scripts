#! /nfshomes/zzare/miniconda3/envs/coverage-pattern/bin/Rscript

## Libraries
library(R.utils)
library(PolyAtailor)
library(Biostrings)
library(conflicted)

#https://github.com/BMILAB/PolyAtailor/blob/f380be097592608857d9b8fb7f70dbf2bfddf799/inst/doc/PolyAtailor_batchCompare_example.Rmd
conflict_prefer("Position", "ggplot2")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("rbind", "base")
conflict_prefer("cbind", "base")
conflict_prefer("strsplit", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("reduce", "IRanges")
conflict_prefer("geom_bar", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("compose", "purrr")
conflict_prefer("last", "dplyr")
conflict_prefer("simplify", "purrr")
conflict_prefer("%>%", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("extract", "tidyr")

## the args part is taken from this link: https://github.com/IARCbioinfo/R-tricks
## Collect arguments
args <- commandArgs(TRUE)

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)


## unzip the .fastq.gz file
#gz_file_path <- args$input  # Path to the .fastq.gz file
#unzipped_file_path <- args$output  # Path where the .fastq file will be saved

# Unzip the .fastq.gz file to a specified location
#gunzip(gz_file_path, destname = unzipped_file_path, remove = FALSE)

# Use the unzipped .fastq file
#https://rdrr.io/github/XHWUlab/PolyAtailor/man/tailScan.html
#fastqfile <- unzipped_file_path
fastqfile <- args$input
tailDF <- tailScan(fastqfile, findUmi = FALSE, resultpath = args$output, samplename = args$samplename, mapping = FALSE)
head(tailDF)
output_path <- file.path(args$output, paste0(args$samplename, "_tailDF_len.tsv"))
write.table(tailDF, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

q(save = "no")
