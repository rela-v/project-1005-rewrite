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
generate_observed_parameters <- function(index, num_expressed_genes, gene_expression_dataframe, zeitgeber_times) {
  out <- fitSinCurve(xx=as.numeric(zeitgeber_times), observed=as.numeric(gene_expression_dataframe[index,]))
  model_sum <- out[[2]]
  print(model_sum)
  out <- out[[1]]
  out_row <- data.frame(Symbols=Symbols[index], A=out$A, phase=out$phase, offset=out$offset, peak=out$peak, R2=out$R2)
  return(out_row)
}

# Generate observed parameters for the SKA2 gene for each condition
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for control...')
observed_para_c_C <- generate_observed_parameters(index=1, num_expressed_genes=nrow(expression_data_C), gene_expression_dataframe=expression_data_C, zeitgeber_times=cohort_data_C$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for SZ...')
observed_para_c_SZ <- generate_observed_parameters(index=1, num_expressed_genes=nrow(expression_data_SZ), gene_expression_dataframe=expression_data_SZ, zeitgeber_times=cohort_data_SZ$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for BD...')
observed_para_c_BD <- generate_observed_parameters(index=1, num_expressed_genes=nrow(expression_data_BD), gene_expression_dataframe=expression_data_BD, zeitgeber_times=cohort_data_BD$ZT)
print('RhythmicityCode.R: Generating observed parameters for SKA2 gene for MDD...')
observed_para_c_MDD <- generate_observed_parameters(index=1, num_expressed_genes=nrow(expression_data_MDD), gene_expression_dataframe=expression_data_MDD, zeitgeber_times=cohort_data_MDD$ZT)

library(gtools)
# Generate permutations of the observed parameters
print('RhythmicityCode.R: Defining generate_permutations function...')
generate_permutations <- function(observed_parameters, num_permutations) {
  combinations <- data.frame(matrix(ncol=5, nrow=num_permutations))
  colnames(combinations) <- c("A", "phase", "offset", "peak", "R2")
  for(i in 1:num_combinations) {
    combinations[i,] <- gtools::combinations(, m, tVec2, repeats.allowed = TRUE)
  }
  return(permutations)
}

# Generate permutations for the observed parameters
print('RhythmicityCode.R: Generating permutations for the observed parameters...')
num_permutations <- 1000
permutations_C <- generate_permutations(expression_data_C, num_permutations)
permutations_SZ <- generate_permutations(expression_data_SZ, num_permutations)
permutations_BD <- generate_permutations(expression_data_BD, num_permutations)
permutations_MDD <- generate_permutations(expression_data_MDD, num_permutations)

# Generate null parameters for the observed parameters
print('RhythmicityCode.R: Defining generate_null_parameters function...')
generate_null_parameters <- function(permutations, model_fn) {
  null_parameters <- data.frame(matrix(ncol=1, nrow=nrow(permutations)))
  for(i in 1:nrow(permutations)) {
    null_parameters[i,] <- model_fn(permutations[i,])$R2
  }
  return(null_parameters)
}

# Generate null parameters for the observed parameters
print('RhythmicityCode.R: Generating null parameters for the observed parameters...')
null_parameters_C <- generate_null_parameters(permutations_C, fitSinCurve)
null_parameters_SZ <- generate_null_parameters(permutations_SZ, fitSinCurve)
null_parameters_BD <- generate_null_parameters(permutations_BD, fitSinCurve)
null_parameters_MDD <- generate_null_parameters(permutations_MDD, fitSinCurve)

# Generate p-values for the observed parameters
print('RhythmicityCode.R: Defining generate_p_values function...')
generate_p_values <- function(observed_parameters, permuted_parameters) {

  p_value <- (permuted_parameters[which(permuted_parameters >= observed_parameters)])/length(permuted_parameters)
  return(p_value)
}

# Generate p-values for the observed parameters
print('RhythmicityCode.R: Generating p-values for the observed parameters...')
R2_permuted_values_C <- generate_p_values(observed_para_c_C$R2, null_parameters_C)
R2_permuted_values_SZ <- generate_p_values(observed_para_c_SZ$R2, null_parameters_SZ)
R2_permuted_values_BD <- generate_p_values(observed_para_c_BD$R2, null_parameters_BD)
R2_permuted_values_MDD <- generate_p_values(observed_para_c_MDD$R2, null_parameters_MDD)

# Generate p-values for the observed parameters
print('RhythmicityCode.R: Defining generate_circadian_drawings function...')

p_value_C <- length((which(R2_permuted_values_C >= observed_para_c_C$R2)==TRUE))/length(R2_permuted_values_C)
p_value_SZ <- length((which(R2_permuted_values_SZ >= observed_para_c_SZ$R2)==TRUE))/length(R2_permuted_values_SZ)
p_value_BD <- length((which(R2_permuted_values_BD >= observed_para_c_BD$R2)==TRUE))/length(R2_permuted_values_BD)
p_value_MDD <- length((which(R2_permuted_values_MDD >= observed_para_c_MDD$R2)==TRUE))/length(R2_permuted_values_MDD)

print('RhythmicityCode.R: Defining generate_circadian_drawings function...')
# Generate scatter plots
gene_names <- c("SKA2")
number_of_top_genes <- nrow(top.control)
data_labels <- cohort_data$suicide
special_information <- "Control"
system("mkdir -p ./Results/Rhythmicity/C/PDF")

#' generate_circadian_drawings() 
#' Generate the drawings of the observed circadian rhythmicity of the genes in question
#' @returns NULL, produces PDFs in the './Results/Rhythmicity/PDF/' folder.
#' @param index represents the number of the gene that will be analyzed.
#' @param top_control_dataframe represents the dataframe of the most rhythmic genes analyzed (based on R2).
#' @param circadian_drawing_function represents the function used to produce the drawings. Customizable in './R/Curve_Drawing.R.'
#' @param cohort_data_dataframe represents the dataframe of cohort metadata used in the input.
#' @param gene_expression_dataframe represents the dataframe of gene expression values.
#' @param gene_names represents the names of the genes contained in the gene_expression_dataframe.
#' @param data_labels represents the special data labels included in the circadian_drawing_function.
#' @param special_information represents, in this case, the group of the genes being analyzed.
#' @examples
#' generate_circadian_drawings(index=1, top_control_dataframe=top_control_genes_df, circadian_drawing_function=circadianDrawing_fn, cohort_data_dataframe=cohort_df, gene_expression_dataframe=gene_expression_df, gene_names=genes, data_labels=labels, special_information="Control") 

generate_circadian_drawings <- function(index, top_control_dataframe, circadian_drawing_function, zeitgeber_times, gene_expression_dataframe, gene_names, data_labels, special_information) {
  a_gene <- "SKA2" 
  fileName <- paste('./Results/Rhythmicity/PDF/Top_Control', a_gene, '.pdf', sep='')
  pdf(fileName)
  circadian_drawing_function(tod=zeitgeber_times, expr=gene_expression_dataframe[gene_names[index],], apar=top_control_dataframe[index,], labels=data_labels, specInfo=special_information)
  dev.off()
}
print("Done")
dir.create('./Results/Rhythmicity/PDF/Control')
fileName_C <- paste('./Results/Rhythmicity/PDF/Control', '_SKA2_', '.pdf', sep='')
dir.create('./Results/Rhythmicity/PDF/SZ')
fileName_SZ <- paste('./Results/Rhythmicity/PDF/SZ', '_SKA2_', '.pdf', sep='')
dir.create('./Results/Rhythmicity/PDF/BD')
fileName_BD <- paste('./Results/Rhythmicity/PDF/BD', '_SKA2_', '.pdf', sep='')
dir.create('./Results/Rhythmicity/PDF/MDD')
fileName_MDD <- paste('./Results/Rhythmicity/PDF/MDD', '_SKA2_', '.pdf', sep='')
pdf(fileName_C)
C_ZT_drawing <- circadianDrawing(tod=cohort_data_C$ZT, expr=expression_data_C[1,], apar=observed_para_c_C, labels=cohort_data_C$suicide, specInfo="Control")
dev.off()
pdf(fileName_SZ)
SZ_ZT_drawing <- circadianDrawing(tod=cohort_data_SZ$ZT, expr=expression_data_SZ[1,], apar=observed_para_c_SZ, labels=cohort_data_SZ$suicide, specInfo="SZ")
dev.off()
pdf(fileName_BD)
BD_ZT_drawing <- circadianDrawing(tod=cohort_data_BD$ZT, expr=expression_data_BD[1,], apar=observed_para_c_BD, labels=cohort_data_BD$suicide, specInfo="BD")
dev.off()
pdf(fileName_MDD)
MDD_ZT_drawing <- circadianDrawing(tod=cohort_data_MDD$ZT, expr=expression_data_MDD[1,], apar=observed_para_c_MDD, labels=cohort_data_MDD$suicide, specInfo="MDD")
dev.off()

print('generate_jakob_plots.R: Routine complete.')


