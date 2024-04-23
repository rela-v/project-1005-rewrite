print("main.R: Starting...")
# Package names
packages <- c("parallelly", "parallel", "minpack.lm", "BiocManager")
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
relevant_conditions<-c("SZ","C")
cohort_data.n<-cohort_data[cohort_data$DDx %in% relevant_conditions,]
write.csv(cohort_data.n, "./Data/revised_cohort_data.csv")
# Number of rows in the filtered cohort data (number of subjects)
revised_cohort_data<-read.csv('./Data/revised_cohort_data.csv', row.names = 1)
n_c<-nrow(revised_cohort_data)
# Complete expression data containing differential gene expression data per subject.
expression_data<-read.csv("./Data/Project_1005-rawCounts.csv", sep='\t')
colnames(expression_data) <- rownames(cohort_data)
# Generate a new empty expression data dataframe to contain filtered expression data by Schizophrenia condition group.
# Get columns of expression data
expression_data_subjects<-rownames(cohort_data)
# Create an empty vector for edited expression_data_subjects vector
print("main.R: Initializing External Functions...")
source("./R/variance_partitioning_analysis.R")
system.time(source("./R/RhythmicityCode.R"))
# print("main.R: Running ./R/DE.R...")
# source("./R/DE.R")
print("main.R: Routine Complete.\nSee Results in the `./Results/` folder.")

