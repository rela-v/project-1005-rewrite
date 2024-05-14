# variance_partitioning_analysis.R

# All variables:
# age
# Sex
# DDx
# psychosis
# suicide
# PMI
# pH
# Alcohol
# substance.abuse.severity
# smoking.history
# Race
# ZT
# onset
# Hemisphere
# Brain.weight
# duration
# lifetime.antipschotics..fluphenazine.eq

non_na_antipsychotics_idx <- which(is.na(cohort_data$lifetime.antipsychotics..fluphenazine.eq.)==FALSE)
cohort_data$lifetime.antipsychotics..fluphenazine.eq.[non_na_antipsychotics_idx] <- log(cohort_data$lifetime.antipsychotics..fluphenazine.eq.)
infinite_antipsychotics_idx <- which(!is.finite(cohort_data$lifetime.antipsychotics..fluphenazine.eq.))
cohort_data$lifetime.antipsychotics..fluphenazine.eq.[infinite_antipsychotics_idx] <- 0
cohort_data$substance.abuse.severity <- as.character(cohort_data$substance.abuse.severity)
formula_varpart_SZ <- ~ age + (1 | Sex) + (1 | suicide) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. +
  (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere) +
  brain.weight
formula_varpart_BD <- ~ age + (1 | Sex) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere)
formula_varpart_MDD <- ~ age + (1 | Sex) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + ZT + (1 | hemisphere)

formula_control <- ~ (1 | Race) + ZT + (1 | Sex) + age + (1 | smoking.history) + pH + brain.weight + PMI

# cordata <- canCorPairs(form, cohort_data)
# plotCorrMatrix(cordata)

# Get variance partitioning plot

# Create functtion to grab metadata and expression data for a particular condition along with control condition

grab_metadata_expression_data <- function(relevant_condition) {
  filtered_cohort_data <- cohort_data[(cohort_data$DDx %in% relevant_condition),]
  filtered_expression_data <- expression_data[colnames(expression_data) %in% rownames(filtered_cohort_data)]
  return(list(filtered_cohort_data, filtered_expression_data))
}

print_column_levels <- function(dataframe) {
  for (col in colnames(dataframe)) {
    levels_count <- length(unique(dataframe[[col]]))
    cat("Column:", col, "\tLevels:", levels_count, "\n")
  }
}


library(parallelly)
library(parallel)
multicoreParam <- MulticoreParam(workers = availableCores())

plotVarPartition_for_disease <- function(condition) {
  formula <- switch(condition,
                    "SZ" = formula_varpart_SZ,
                    "BD" = formula_varpart_BD,
                    "MDD" = formula_varpart_MDD,
                    "C" = formula_control)
  data <- grab_metadata_expression_data(condition)
  print_column_levels(data[[1]])
  varPart <- fitExtractVarPartModel(data[[2]], formula, data[[1]], BPPARAM= multicoreParam)
  vp <- sortCols(varPart)
  percentbarplot <- plotPercentBars(vp[1:10, ], ggtitle="Variance Explained by Top 10 Variables")
  varpartplot <- plotVarPart(vp, ggtitle=paste("Variance Explained by Top 10 Variables in ", condition))
  dir.create(paste('./Results/variance_partitioning_plots/', condition, sep=''), recursive=TRUE)
  ggsave(paste("./Results/variance_partitioning_plots/", condition, "/variance_partitioningbarplot_", condition, ".png", sep=''), plot=percentbarplot)
  ggsave(paste("./Results/variance_partitioning_plots/", condition, "/variance_partitioningviolinplot_", condition, ".png", sep=''), plot=varpartplot)
  save(varPart, file=paste("./Results/variance_partitioning_plots/", condition, "/variance_partitioning_", condition, ".RData", sep=''))
}

for (condition in c("SZ", "BD", "MDD", "C")) {
  print(paste0("Running variance partitioning for ", condition))
  plotVarPartition_for_disease(condition)
}


