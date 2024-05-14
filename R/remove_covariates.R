# Removing covariates from data using data from the variance partitioning analysis

# Create formulae of covariates of interest
formula_SZ <- ~ age + (1 | Sex) + (1 | DDx) + (1 | psychosis) + (1 | suicide) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. +
  (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere) +
  brain.weight
formula_BPD_MDD <- ~ age + (1 | Sex) + (1 | DDx) + PMI + pH + (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + (1 | hemisphere)
formula_control <- ~ (1 | Race) + ZT + (1 | Sex) + age + (1 | smoking.history) + pH + brain.weight + PMI

formula_BPD_sub <- ~ age + PMI + pH + (1 | ethanol.severity) + (1 | hemisphere) 
formula_MDD_sub <- ~ PMI + age + pH + (1 | ethanol.severity) + (1 | Race)
formula_SZ_sub <- ~ lifetime.antipsychotics..fluphenazine.eq. + PMI + brain.weight + pH + age 
formula_control_sub <- ~ PMI + brain.weight + age + pH + (1 | Race)

library(lme4)

print_column_levels <- function(dataframe) {
  for (col in colnames(dataframe)) {
    levels_count <- length(unique(dataframe[[col]]))
    cat("Column:", col, "\tLevels:", levels_count, "\n")
  }
}


# Remove covariates from data
remove_covariates <- function(condition) {
  prelim_data <- grab_metadata_expression_data(condition)
  metadata <- prelim_data[[1]]
  expr <- prelim_data[[2]]
  print_column_levels(metadata)
  # Load the current varPart model
  formula <- switch(condition,
                    "SZ" = formula_SZ,
                    "BP" = formula_BPD_MDD,
                    "MDD" = formula_BPD_MDD,
                    "C" = formula_control)
  sub_form <- switch(condition,
                     "SZ" = formula_SZ_sub,
                     "BP" = formula_BPD_sub,
                     "MDD" = formula_MDD_sub,
                     "C" = formula_control_sub)

  load(paste0('./Results/variance_partitioning_plots/', condition, '/variance_partitioning_', condition, '.RData'))
  
  # Fit the variance partitioning model
  print("db1")
  init <- fitVarPartModel(exprObj=expr, formula= formula, data=metadata, BPPARAM= multicoreParam)
  print("db2")
  res <- fitVarPartModel(exprObj=expr, formula=sub_form, data=metadata, BPPARAM= multicoreParam)
  print("db3")
  final_models <- init-res
  # Get the residuals
  print("db4")
  final_model <- varPart-res 
  print("db5")
  variancepart <- calcVarPart(final_model)
  print("db6")
  vp <- sortCols(final_model)
  percentbarplot <- plotPercentBars(vp[1:10, ], ggtitle="Variance Explained by Top 10 Variables")
  varpartplot <- plotVarPart(vp, ggtitle=paste("Variance Explained by Top 10 Variables in ", condition))
  ggsave(paste0("./Results/variance_partitioning_plots/", condition, "/corrrected_variance_partitioningbarplot_", condition, ".png", sep=''), plot=percentbarplot)
  ggsave(paste0("./Results/variance_partitioning_plots/", condition, "/corrected_variance_partitioningviolinplot_", condition, ".png", sep=''), plot=varpartplot)
  save(final_model, file=paste("./Results/variance_partitioning_plots/", condition, "/corrected_variance_partitioning_", condition, ".RData", sep=''))
  return(final_model)
} 

for (condition in c("BP", "MDD", "SZ", "C")) {
  print(paste0("Beginning analysis on ", condition))
  remove_covariates(condition)
}

