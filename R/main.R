print("main.R: Starting...")
# Package names
packages <- c("parallelly", "parallel", "minpack.lm", "BiocManager", "curl", "gtools")
# Install packages that aren't already there
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
  BiocManager::install("variancePartition")
}
print("main.R: Loading packages...")
# Load packages
invisible(lapply(packages, library, character.only = TRUE))
library(variancePartition)
print("main.R: Correcting Datasets...")
cohort_data<-read.csv('./Data/cohort_data.csv', row.names = 1)
# Select for subjects in the Schizophrenia condition, "SZ", and control condition, "C"
header  <- scan("./Data/Project_1005-rawCounts.csv", nlines = 2, what = character())
expression_data<-read.csv("./Data/Project_1005-rawCounts.csv", sep='\t', skip=2, header=FALSE,row.names=1)
library(curl)
# annotations <- curl("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2024-04-01.txt")
# annotations <- read.delim(annotations, header=T, skip=0, sep = "\t")
# print(head(annotations))

# colnames(expression_data) <- c(rownames(cohort_data))
# rows_of_expression_data <- annotations$ensembl_gene_id[match(rownames(expression_data), annotations$entrez_id)]
# expression_data$Gene <- rows_of_expression_data
# write.csv(expression_data, file = "./Data/Project_1005-rawCounts-annotated.csv", row.names = FALSE)
# Generate a new empty expression data dataframe to contain filtered expression data by Schizophrenia condition group.
# Get columns of expression data
expression_data_subjects<-rownames(cohort_data)
# Create an empty vector for edited expression_data_subjects vector
print("main.R: Initializing External Functions...")
# source("./R/variance_partitioning_analysis.R")
# system.time(source("./R/RhythmicityCode.R"))
source('./R/Curve_Drawing.R')
source('./R/fitSinCurve.R')
source('./R/generate_jakob_plots.R')
# source("./R/generate_p_values.R")
# print("main.R: Running ./R/DE.R...")
# source("./R/DE.R")
print("main.R: Routine Complete.\nSee Results in the `./Results/` folder.")

