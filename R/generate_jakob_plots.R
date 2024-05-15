# Generate Plots for Jakob Hartmann

# Names of rows in gene expression matrix
expression_data_of_interest <- read.csv("./Data/Project_1005-rawCounts-annotated.csv")
cohort_data <- read.csv("./Data/cohort_data.csv", row.names = 1)
colnames(expression_data_of_interest)  <- c(rownames(cohort_data), 'Gene')
Symbols <- expression_data_of_interest$Gene
# Observed dataframe of genes and fitted curve parameters
# observed_para_c <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
# Time of death (clinical)
grab_metadata_expression_data <- function(relevant_condition) {
  filtered_cohort_data <- cohort_data[(cohort_data$DDx %in% relevant_condition),]
  filtered_expression_data <- expression_data_of_interest[colnames(expression_data_of_interest) %in% rownames(filtered_cohort_data)]
  return(list(filtered_cohort_data, filtered_expression_data))
}

multicoreParam <- MulticoreParam(workers = availableCores())
data_C <- grab_metadata_expression_data("C")
data_SZ <- grab_metadata_expression_data("SZ")
data_BD <- grab_metadata_expression_data("BD")
data_MDD <- grab_metadata_expression_data("MDD")

cohort_data_C <- data_C[[1]]
cohort_data_SZ <- data_SZ[[1]]
cohort_data_BD <- data_BD[[1]]
cohort_data_MDD <- data_MDD[[1]]

expression_data_C <- data_C[[2]][15150,]
expression_data_SZ <- data_SZ[[2]][15150,]
expression_data_BD <- data_BD[[2]][15150,]
expression_data_MDD <- data_MDD[[2]][15150,]


set.seed(123)
# 1->number of rows: fit the sine curves, assign the output to each row.
`nls.lm` <- minpack.lm::nls.lm
# SKA2 Index  <- 15150


#' generate_observed_parameters()
#' Generate the parameters for the fitted sinusoid observed in the expressed genes over time 
#' @returns Dataframe containing "A" (amplitude), "phase" (phase), "offset" (x offset), "peak" (peak value), and "R2" (correlation).
#' @param index represents the index of the gene to be analyzed.
#' @param num_expressed_genes represents the number of genes there are to be analyzed in the dataframe.
#' @param gene_expression_dataframe represents the dataframe used as input for the analysis.
#' @param zeitgeber_times represents the Zeitgeber Time, or the relative timescale used for the circadian analysis.
#' @examples
#' out_pars <- generate_observed_parameters(index=1, num_expressed_genes=10345, gene_expression_dataframe=df_with_gene_data, zeitgeber_time=df_with_ZT_data
print('RhythmicityCode.R: Defining generate_observed_parameters function...')
generate_observed_parameters <- function(index, condition, num_expressed_genes, gene_expression_dataframe, zeitgeber_times) {
    covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", condition, "/residualmodelfit_", condition, ".RData", sep='')
  load(covariate_adjusted_observations_filename)
  out <- fitSinCurve(xx=as.numeric(zeitgeber_times), observed=as.numeric(residual))
  out_row <- data.frame(A=out$A, phase=out$phase, offset=out$offset, peak=out$peak, R2=out$R2) 
  return(out_row)
}

# Generate observed parameters for the SKA2 gene for each condition
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for control...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "C", "/residualmodelfit_", "C", ".RData", sep='')
load(covariate_adjusted_observations_filename)
observed_para_c_C <- generate_observed_parameters(index=1, condition="C", num_expressed_genes=nrow(residual), gene_expression_dataframe=residual, zeitgeber_times=cohort_data_C$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for SZ...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "SZ", "/residualmodelfit_", "SZ", ".RData", sep='')
load(covariate_adjusted_observations_filename)
observed_para_c_SZ <- generate_observed_parameters(index=1, condition="SZ", num_expressed_genes=nrow(residual), gene_expression_dataframe=residual, zeitgeber_times=cohort_data_SZ$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for BD...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "BD", "/residualmodelfit_", "BD", ".RData", sep='')
load(covariate_adjusted_observations_filename)
observed_para_c_BD <- generate_observed_parameters(index=1, condition="BD", num_expressed_genes=nrow(residual), gene_expression_dataframe=residual, zeitgeber_times=cohort_data_BD$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for MDD...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "MDD", "/residualmodelfit_", "MDD", ".RData", sep='')
load(covariate_adjusted_observations_filename)
observed_para_c_MDD <- generate_observed_parameters(index=1, condition="MDD", num_expressed_genes=nrow(residual), gene_expression_dataframe=residual, zeitgeber_times=cohort_data_MDD$ZT)


print('RhythmicityCode.R: Defining generate_circadian_drawings function...')
# Generate scatter plots
gene_names <- c("SKA2")
number_of_top_genes <- nrow(top.control)
system("mkdir -p ./Results/Rhythmicity/C/PDF")

print("Done")
dir.create('./Results/Rhythmicity/PDF/Control')
fileName_C <- paste('./Results/Rhythmicity/PDF/Control', '_SKA2_', '.png', sep='')
dir.create('./Results/Rhythmicity/PDF/SZ')
fileName_SZ <- paste('./Results/Rhythmicity/PDF/SZ', '_SKA2_', '.png', sep='')
dir.create('./Results/Rhythmicity/PDF/BD')
fileName_BD <- paste('./Results/Rhythmicity/PDF/BD', '_SKA2_', '.png', sep='')
dir.create('./Results/Rhythmicity/PDF/MDD')
fileName_MDD <- paste('./Results/Rhythmicity/PDF/MDD', '_SKA2_', '.png', sep='')

print('RhythmicityCode.R: Generating circadian drawings for SKA2 gene for control...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "C", "/residualmodelfit_", "C", ".RData", sep='')
load(covariate_adjusted_observations_filename)
print(residual)
png(fileName_C)
circadianDrawing(tod=cohort_data_C$ZT, expr=residual, apar=observed_para_c_C, labels=cohort_data_C$suicide, specInfo=special_information)
dev.off()
print('RhythmicityCode.R: Generating circadian drawings for SKA2 gene for SZ...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "SZ", "/residualmodelfit_", "SZ", ".RData", sep='')
load(covariate_adjusted_observations_filename)
png(fileName_SZ)
circadianDrawing(tod=cohort_data_SZ$ZT, expr=residual, apar=observed_para_c_SZ, labels=cohort_data_SZ$suicide, specInfo=special_information)
dev.off()
print('RhythmicityCode.R: Generating circadian drawings for SKA2 gene for BD...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "BD", "/residualmodelfit_", "BD", ".RData", sep='')
load(covariate_adjusted_observations_filename)
png(fileName_BD)
circadianDrawing(tod=cohort_data_BD$ZT, expr=residual, apar=observed_para_c_BD, labels=cohort_data_BD$suicide, specInfo=special_information)
dev.off()
print('RhythmicityCode.R: Generating circadian drawings for SKA2 gene for MDD...')
covariate_adjusted_observations_filename <- paste("./Results/variance_partitioning_plots/", "MDD", "/residualmodelfit_", "MDD", ".RData", sep='')
load(covariate_adjusted_observations_filename)
png(fileName_MDD)
circadianDrawing(tod=cohort_data_MDD$ZT, expr=residual, apar=observed_para_c_MDD, labels=cohort_data_MDD$suicide, specInfo=special_information)
dev.off()


print('generate_jakob_plots.R: Routine complete.')



