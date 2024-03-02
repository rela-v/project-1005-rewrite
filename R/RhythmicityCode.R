#' This is the RhythmicityCode.R file, which calculates rhythmicity in genes.
print('RhythmicityCode.R: Starting...')
library('minpack.lm')
source('./R/fitSinCurve.R')
cohort_data<-read.csv('./Data/revised_cohort_data.csv', row.names=1)
gene_expression_dataframe<-read.csv('./Data/revised_expression_data.csv', row.names=1)
n.cores <- availableCores()
num_expressed_genes <- nrow(gene_expression_dataframe)

# Names of rows in gene expression matrix
Symbols<-row.names(gene_expression_dataframe)
# Observed dataframe of genes and fitted curve parameters
# observed_para_c <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
# Time of death (clinical)
zeitgeber_times<-cohort_data$ZT
# 1->number of rows: fit the sine curves, assign the output to each row.
`nls.lm` <- minpack.lm::nls.lm

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
  out_row <- data.frame(Symbols=Symbols[index], A=out$A, phase=out$phase, offset=out$offset, peak=out$peak, R2=out$R2)
  return(out_row)
}

print('RhythmicityCode.R: Running generate_observed_parameters function in parallel...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing generate_observed_parameters
clusterExport(cl, c("generate_observed_parameters", "fitSinCurve", "nls.lm", "parStartVal", "Symbols", "num_expressed_genes", "gene_expression_dataframe", "zeitgeber_times"))
observed_para_c <- do.call(rbind,parLapply(cl, X=1:num_expressed_genes, fun=generate_observed_parameters, num_expressed_genes, gene_expression_dataframe, zeitgeber_times))
stopCluster(cl)

observed_para_c_sorted <- observed_para_c[order(observed_para_c$R2, decreasing=TRUE),]

print('RhythmicityCode.R: Defining generate_circadian_drawings function...')
# Generate scatter plots
top.control <- observed_para_c_sorted[1:4,]
gene_names <- as.numeric(row.names(top.control))
number_of_top_genes <- nrow(top.control)
data_labels <- cohort_data$suicide
special_information <- "Control"
system("mkdir -p ./Results/Rhythmicity/PDF")

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
  a_gene <- as.character(gene_names[index])
  fileName <- paste('./Results/Rhythmicity/PDF/Top_Control', a_gene, '.pdf', sep='')
  pdf(fileName)
  circadian_drawing_function(tod=zeitgeber_times, expr=gene_expression_dataframe[gene_names[index],], apar=top_control_dataframe[index,], labels=data_labels, specInfo=special_information)
  dev.off()
}

print('RhythmicityCode.R: Running generate_circadian_drawings function in parallel...')
cl <- makeClusterPSOCK(n.cores)
source('./R/Curve_Drawing.R')
# Parallelizing generate_circadian_drawings from range 1:nrow(top.control)
clusterExport(cl, c("generate_circadian_drawings", "number_of_top_genes", "top.control", "circadianDrawing", "zeitgeber_times", "gene_expression_dataframe", "gene_names", "data_labels", "special_information"))
parLapply(cl, X=1:number_of_top_genes, fun=generate_circadian_drawings, top_control_dataframe=top.control, circadian_drawing_function=circadianDrawing, zeitgeber_times=zeitgeber_times, gene_expression_dataframe=gene_expression_dataframe, gene_names=gene_names, data_labels=data_labels, special_information=special_information)
stopCluster(cl)

