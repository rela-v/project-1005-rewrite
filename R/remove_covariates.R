# Removing covariates from data using data from the variance partitioning analysis

# Create formulae of covariates of interest

formula_varpart_SZ <- ~ age + (1 | Sex) + (1 | suicide) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. +
  (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere) +
  brain.weight
formula_varpart_BD <- ~ age + (1 | Sex) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere)
formula_varpart_MDD <- ~ age + (1 | Sex) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + ZT + (1 | hemisphere)
formula_control <- ~ (1 | Race) + ZT + (1 | Sex) + age + (1 | smoking.history) + pH + brain.weight + PMI

formula_wo_sub_SZ  <- ~ (1 | Sex) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. + (1 | ethanol.severity) + (1 | Race) + ZT + (1 | hemisphere)
formula_wo_sub_MDD <- ~ (1 | Sex) + (1 | smoking.history) + ZT + (1 | hemisphere)
formula_wo_sub_BPD <- ~ (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere)
formula_wo_sub_control <- ~ ZT + (1 | Sex) + (1 | smoking.history)


formula_BPD_sub <- ~ (1 | ethanol.severity) + PMI + pH + age + (1 | Sex) 
formula_MDD_sub <- ~ (1 | substance.abuse.severity) + age + PMI + pH + (1 | ethanol.severity)
formula_SZ_sub <- ~ (1 | substance.abuse.severity) + (1 | smoking.history) + age + brain.weight + (1 | suicide) 
formula_control_sub <- ~ PMI + brain.weight + age + pH + (1 | Race)

library(lme4)

print_column_levels <- function(dataframe) {
  for (col in colnames(dataframe)) {
    levels_count <- length(unique(dataframe[[col]]))
    cat("Column:", col, "\tLevels:", levels_count, "\n")
  }
}

grab_metadata_expression_data <- function(relevant_condition) {
  filtered_cohort_data <- cohort_data[(cohort_data$DDx %in% relevant_condition),]
  filtered_expression_data <- expression_data[colnames(expression_data) %in% rownames(filtered_cohort_data)]
  return(list(filtered_cohort_data, filtered_expression_data))
}

multicoreParam <- MulticoreParam(workers = availableCores())
# Remove covariates from data
remove_covariates <- function(condition) {
  prelim_data <- grab_metadata_expression_data(condition)
  metadata <- prelim_data[[1]]
  
  metadata <- sapply(metadata, function(x) {
                       if (any(is.na(x))) {
                         x[is.na(x)] <- 0
                       } else if (any(is.nan(x))) {
                         x[is.nan(x)] <- 0
                       } else if (any(is.infinite(x))) {
                         x[is.infinite(x)] <- 0
                       } else {
                         x
                       } 
                })
  metadata <- as.data.frame(metadata)
  metadata[c("ZT", "age", "onset", "duration", "lifetime.antipsychotics..fluphenazine.eq.", "pH", "brain.weight", "PMI")] <- lapply(metadata[c("ZT", "age", "onset", "duration", "lifetime.antipsychotics..fluphenazine.eq.", "pH", "brain.weight", "PMI")], as.numeric)
  metadata$lifetime.antipsychotics..fluphenazine.eq. <- log10(metadata$lifetime.antipsychotics..fluphenazine.eq. + 1)
  print(metadata$lifetime.antipsychotics..fluphenazine.eq.)
  expr <- prelim_data[[2]]
  print_column_levels(metadata)
  # Load the current varPart model
  formula <- switch(condition,
                    "SZ" = formula_varpart_SZ,
                    "BD" = formula_varpart_BD,
                    "MDD" = formula_varpart_MDD,
                    "C" = formula_control)
  sub_form <- switch(condition,
                     "SZ" = formula_SZ_sub,
                     "BD" = formula_BPD_sub,
                     "MDD" = formula_MDD_sub,
                     "C" = formula_control_sub)
  formula_wo_sub <- switch(condition,
                           "SZ" = formula_wo_sub_SZ,
                           "BD" = formula_wo_sub_BPD,
                           "MDD" = formula_wo_sub_MDD,
                           "C" = formula_wo_sub_control)

  
  # Fit the variance partitioning model
  modelfit <- fitVarPartModel(expr, formula=sub_form, data=metadata, BPPARAM= multicoreParam) 
  residual  <- residuals(modelfit, expr)
  save(residual, file=paste("./Results/variance_partitioning_plots/", condition, "/residualmodelfit_", condition, ".RData", sep=''))
  return(residual)
} 
# BD, MDD, SZ, C
for (condition in c("C")) {
  print(paste0("Beginning analysis on ", condition))
  remove_covariates(condition)
}
