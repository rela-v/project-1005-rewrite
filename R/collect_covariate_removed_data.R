# collect_covariate_removed_data.R

generate_file_name  <- function(condition) {
  # 'Results/variance_partitioning_plots/BD/residualmodelfit_BD.RData'
  filename <- paste0("Results/variance_partitioning_plots/", condition, "/residualmodelfit_", condition, ".RData")
}

grab_residuals <- function(condition) {
  filename <- generate_file_name(condition)
  load(filename)
  write.csv(residual, paste0("Results/residuals/", condition, "/residualdata_", condition, ".csv"))
  rm(residual)
}
