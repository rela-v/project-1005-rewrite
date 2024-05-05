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
formula_varpart_SZ <- ~ age + (1 | Sex) + (1 | DDx) + (1 | psychosis) + (1 | suicide) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. +
  (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere) +
  brain.weight
formula_varpart_BD_MDD <- ~ age + (1 | Sex) + (1 | DDx) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere)

cordata <- canCorPairs(form, cohort_data)
plotCorrMatrix(cordata)

# Get variance partitioning plot

# Create functtion to grab metadata and expression data for a particular condition along with control condition

grab_metadata_expression_data <- function(relevant_condition) {
  filtered_cohort_data <- cohort_data[(cohort_data$DDx %in% relevant_condition | cohort_data$DDx %in% 'C'),]
  filtered_expression_data <- expression_data[colnames(expression_data) %in% rownames(filtered_cohort_data)]
  return(list(filtered_cohort_data, filtered_expression_data))
}

get_levels_per_variable <- function(cohort_data_filtered) {
  for(col in colnames(cohort_data_filtered)) {
    string <- paste("Levels for the particular variable ", col, " are: ", length(unique(cohort_data_filtered[,col])))
    print(string)
  }
  return('done')
}


library(parallelly)
library(parallel)
multicoreParam <- MulticoreParam(workers = availableCores())
data <- grab_metadata_expression_data("MDD")
plotVarPartition_for_disease <- function(condition) {
  varPart <- fitExtractVarPartModel(data[[2]], formula_varpart_BD_MDD, data[[1]], BPPARAM= multicoreParam)
  vp <- sortCols(varPart)
  percentbarplot <- plotPercentBars(vp[1:10, ], ggtitle="Variance Explained by Top 10 Variables")
  varpartplot <- plotVarPart(vp, ggtitle=paste("Variance Explained by Top 10 Variables in ", condition))
  dir.create(paste('./Results/', condition, sep=''), recursive=TRUE)
  ggsave(paste("./Results/", condition, "/variance_partitioningbarplot_", condition, ".png", sep=''), plot=percentbarplot)
  ggsave(paste("./Results/", condition, "/variance_partitioningviolinplot_", condition, ".png", sep=''), plot=varpartplot)
  save(varPart, file=paste("./Results/", condition, "/variance_partitioning_", condition, ".RData", sep=''))
}

for (condition in c("MDD")) {
  plotVarPartition_for_disease(condition)
}


