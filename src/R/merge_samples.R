#!/usr/bin/env Rscript

# Import libraries
library(tidyverse)
library(dada2)
library(xlsx)
library(phyloseq)

# User supplied parameters
args <- commandArgs(trailingOnly = TRUE)

# Data directory
data <- args[0]
out_file <- args[1]

# List consensus files
file_names <- list.files(
  path = data,
  pattern = "_cons_otu.csv",
  full.names = TRUE
)

sample_names <- list.files(
  path = data,
  pattern = "_cons_otu.csv",
  full.names = FALSE
)

# loads all sequence tables in a list with data.frames
cons_tables <- lapply(fileNames, read.csv)

# name the dataframes within the list
names(cons_tables) <- gsub(
  pattern = "\\_cons_otu.csv$",
  replacement = "",
  x = sampleNames
)

# set colnames suitable for dada2
colnames <- c("sequence", "abundance")
cons_tables <- lapply(cons_tables, setNames, colnames)

# Generate final sequence table
seqtab_final <- dada2::makeSequenceTable(
  samples = cons_tables
)

# Save sequence table to fasta file
dada2::uniquesToFasta(
  unqs = seqtab_final,
  fout = file.path(data, out_file),
  ids = colnames(seqtab_final)
)
