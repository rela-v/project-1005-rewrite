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

print('RhythmicityCode.R: Initializing create_shuffled_TOD function...')
system('mkdir -p ./nullFolder')
setwd('./nullFolder/')
groupName <- 'control'
permutations <- 10 #Alter to 1000 permutations upon completion.

shuffleTOD_df <- data.frame(matrix(ncol=0, nrow=length(zeitgeber_times)))

#' create_shuffled_TOD()
#' A function to produce shuffled times of death (TOD) for the null distributions.
#' @returns dataframe, containing shuffled zeitgeber times.
#' @param index represents in this case which permutation of the shuffle is being generated.
#' @param zeitgeber_times represents the dataframe containing the zeitgeber times of death of the subjects.
#' @examples
#' create_shuffled_TOD(index=1, out_df=some_output_dataframe, zeitgeber_times=some_dataframe_with_TOD
create_shuffled_TOD <- function(index, zeitgeber_times) {
  set.seed(index)
  return(data.frame(sample(zeitgeber_times)))
}
print('RhythmicityCode.R: Running create_shuffled_TOD function in parallel...')

cl <- makeClusterPSOCK(n.cores)
# Parallelized create_shuffled_TOD function in the range from 1:permutations
clusterExport(cl, c("create_shuffled_TOD", "permutations", "zeitgeber_times"))
shuffleTOD_df <- do.call(cbind,parLapply(cl, X=1:permutations, fun=create_shuffled_TOD, zeitgeber_times))
stopCluster(cl)

print('RhythmicityCode.R: Initializing generate_null_data_filenames function...')
#' generate_null_data_filenames()
#' Generate the filenames for the null_data files. 
#' @returns string, containing the filename of the null_data at position 'index'.
#' @param index represents which permutation the null_data is being produced for.
#' @param groupName represents the group (e.g. Control, Old, Young, etc.) that these files are being produced for.
#' @examples
#' generate_null_data_filenames(index=1, groupName="Control")
generate_null_data_filenames <- function(index, groupName) {
  return(paste('null_', groupName, '_', index, '.rdata', sep=''))
}

print('RhythmicityCode.R: Running generate_null_data_filenames in parallel...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing execution of the generate_null_data_filenames function over the range 1:permutations
clusterExport(cl, c("generate_null_data_filenames", "groupName"))
null_para_files <- do.call(rbind, parLapply(cl, X=1:permutations, fun=generate_null_data_filenames, groupName))
stopCluster(cl)

print('RhythmicityCode.R: Initializing generate_null_data_parameters function...')
#' generate_null_data_parameters()
#' Generate the null parameters for each permutation's shuffled time of deaths. 
#' @returns dataframe, containing a row of sinusoid parameters (amplitude, phase, offset, peak, and R2).
#' @param index represents which gene the sinusoid parameters are being generated for.
#' @param permutation represents which permutation of the time of death shuffling is being used to generate parameters.
#' @param gene_expression_dataframe represents the dataframe with gene expression data inside
#' @param curve_fitting_function represents the function that is being used to generate the null_parameters.
#' @examples
#' output_row <- generate_null_data_parameters(index=1, permutation=2, shuffled_TOD_dataframe=shuffleTOD_df, gene_expression_dataframe=gene_expression_df, curve_fitting_function=fitSinCurve)
generate_null_data_parameters <- function(index, permutation, shuffled_TOD_dataframe, gene_expression_dataframe, curve_fitting_function) {
  out <- curve_fitting_function(xx=shuffled_TOD_dataframe[,permutation], observed=unlist(gene_expression_dataframe[index,]))
  out_row <- data.frame(out$A, out$phase, out$offset, out$peak, out$R2)
  return(out_row)
}

print('RhythmicityCode.R: Running generate_null_data_parameters function in parallel...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing the generate_null_data_parameters function across the range 1:num_expressed_genes to append to dataframe null_pare
clusterExport(cl, c("generate_null_data_parameters", "parStartVal", "shuffleTOD_df", "nls.lm", "gene_expression_dataframe", "fitSinCurve"))
for(permutation in 1:permutations) {
  clusterExport(cl, "permutation")
  null_pare <- do.call(rbind, parLapply(cl, X=1:num_expressed_genes, fun=generate_null_data_parameters, permutations, shuffleTOD_df, gene_expression_dataframe, fitSinCurve))
  colnames(null_pare) <- c("A", "phase", "offset", "peak", "R2")
  save(null_pare, file=null_para_files[permutation,])
}
stopCluster(cl)

print('RhythmicityCode.R: Initializing collect_variable_from_permutations function...')
null_pare_A <- data.frame(matrix(0, num_expressed_genes, permutations))
null_pare_phase <- data.frame(matrix(0, num_expressed_genes, permutations))
null_pare_offset <- data.frame(matrix(0, num_expressed_genes, permutations))
null_pare_peak <- data.frame(matrix(0, num_expressed_genes, permutations))
null_pare_R2 <- data.frame(matrix(0, num_expressed_genes, permutations))

# Method: 
# 1) parallelize getting all of the null_pare files
# 2) get all of the null_pare_A, phase, offset, peak, R2...
# 3) cbind all A together, all phase together, etc.

#' collect_variable_from_permutations() 
#' Function for collecting a particular variable from all permutation dataframes 
#' @returns dataframe, containing the specific column of the dataframe.
#' @param permutation represents which gene is being analyzed.
#' @param variable represents which variable is being grabbed.
#' @examples
#' collect_variable_from_permutations(1, 'A')
collect_variable_from_permutations <- function(permutation, variable) {
  null_pare <- get(load(null_para_files[permutation]))
  return(null_pare[variable])
}
print('RhythmicityCode.R: Running collect_variable_from_permutations function...')
cl <- makeClusterPSOCK(n.cores)
# Parallelize the collect_variable_from_permutations function across the range 1:permutations.

for(variable in c('A', 'phase', 'offset', 'peak', 'R2')) {
  print(variable)
  clusterExport(cl, c("variable", "collect_variable_from_permutations", "null_para_files", "permutations"))
  assign(paste('null_para', variable, sep='_'), do.call(cbind, parLapply(cl, X=1:permutations, fun=collect_variable_from_permutations, variable)),envir=.GlobalEnv)
}
stopCluster(cl)
null_para <- list(null_para_A=null_para_A, null_para_phase=null_para_phase, null_para_offset=null_para_offset, null_para_peak=null_para_peak, null_para_R2=null_para_R2)
null_para_file <- paste('null_', groupName, '.rdata', sep='')
save(null_para, file=null_para_file)
setwd('..')
para_R2_pool <- cbind(observed_para_c$R2, null_para$null_para_R2)
R2Rank_para <- wilcox.test(para_R2_pool[,1], para_R2_pool[,2])$p.value
observed_para_c$pvalue <- R2Rank_para
observed_para_c$qvalue <- p.adjust(observed_para_c$pvalue, 'BH')
observed_para_c_sorted <- observed_para_c[order(observed_para_c$pvalue),]
write.csv(observed_para_c_sorted, "./Results/Rhythmicity/Example_result.csv")

# Rhythmicity gain/loss

#Split data into two groups: For example we will split the controls into younger and older 
old_index<-which(cohort_data$age>=50)
young_index<-which(cohort_data$age<50)

expr.old<-gene_expression_dataframe[,old_index]
expr.young<-gene_expression_dataframe[,young_index]

tod_o<-cohort_data$ZT[old_index]
tod_y<-cohort_data$ZT[young_index]

observed_para_y <- data.frame(genes=Symbols,A=numeric(num_expressed_genes), phase=numeric(num_expressed_genes), offset=numeric(num_expressed_genes), peak=numeric(num_expressed_genes), R2=numeric(num_expressed_genes))
observed_para_o <- data.frame(genes=Symbols,A=numeric(num_expressed_genes), phase=numeric(num_expressed_genes), offset=numeric(num_expressed_genes), peak=numeric(num_expressed_genes), R2=numeric(num_expressed_genes))

print('RhythmicityCode.R: Running generate_observed_parameters function for old cohort in parallel...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing generate_observed_parameters
clusterExport(cl, c("generate_observed_parameters", "fitSinCurve", "nls.lm", "parStartVal", "Symbols", "num_expressed_genes", "expr.old", "tod_o"))
observed_para_c <- do.call(rbind,parLapply(cl, X=1:num_expressed_genes, fun=generate_observed_parameters, num_expressed_genes, expr.old, tod_o))
stopCluster(cl)

print('RhythmicityCode.R: Running generate_observed_parameters function for young cohort in parallel...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing generate_observed_parameters
clusterExport(cl, c("generate_observed_parameters", "fitSinCurve", "nls.lm", "parStartVal", "Symbols", "num_expressed_genes", "expr.young", "tod_y"))
observed_para_c <- do.call(rbind,parLapply(cl, X=1:num_expressed_genes, fun=generate_observed_parameters, num_expressed_genes, expr.young, tod_y))
stopCluster(cl)
