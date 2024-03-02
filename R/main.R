print("main.R: Starting...")
# Package names
packages <- c("parallelly", "parallel", "minpack.lm", "lme4", "BiocManager", "qvalue")
# Install packages that aren't already there
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("qvalue")
}
print("main.R: Loading packages...")
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

print("main.R: Correcting Datasets...")
cohort_data<-read.csv('./Data/cohort_data.csv', row.names = 1)
# Select for subjects in the Schizophrenia condition, "SZ", and control condition, "C"
relevant_conditions<-c("SZ","C")
cohort_data.n<-cohort_data[cohort_data$DDx %in% relevant_conditions,]
write.csv(cohort_data.n, "./Data/revised_cohort_data.csv")
# Number of rows in the filtered cohort data (number of subjects)
cohort_data<-read.csv('./Data/revised_cohort_data.csv', row.names = 1)
n_c<-nrow(cohort_data)
# Complete expression data containing differential gene expression data per subject.
expression_data<-read.csv("./Data/expression_data.csv", row.names = 1)
# Generate a new empty expression data dataframe to contain filtered expression data by Schizophrenia condition group.
# Get columns of expression data
expression_data_subjects<-colnames(expression_data)
# Create an empty vector for edited expression_data_subjects vector
new.dat.ls <- c()
for (i in 1:ncol(expression_data)) {
  new.dat<-strsplit(x=strsplit(expression_data_subjects[i],split='_')[[1]][2],split='S')[[1]][2]
  new.dat.ls<- c(new.dat.ls,new.dat)
}
# new.dat.ls: "162" "158" "136", etc.
# Create an 'subject_mappings' vector containing the indices of each subject
subject_mappings<-match(row.names(cohort_data), new.dat.ls)
# subject_mappings: `52 51 45`, etc.
# fix expression data

expression_data_new<-expression_data[,c(subject_mappings)]
write.csv(expression_data_new, "./Data/revised_expression_data.csv")
print("main.R: Installing/verifying packages...")

print("main.R: Initializing External Functions...")
source("./R/fitSinCurve.R")
source("./R/Curve_Drawing.R")
source("./R/bestModelSelection.R")
print("main.R: Running RhythmicityCode.R...")
source("./R/RhythmicityCode.R")
# print("main.R: Running ./R/DE.R...")
# source("./R/DE.R")
print("main.R: Routine Complete.\nSee Results in the `./Results/` folder.")
