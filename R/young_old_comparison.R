# young_old_comparison.R

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
# Generate Null Distributions: 
setwd('./nullFolder')
library(doParallel)
groupName <- 'old'
thisData <- expr.old

###REPEAT NULL GENERATION FOR Old###

print('RhythmicityCode.R: Initializing create_shuffled_TOD function for old cohort...')
groupName <- 'old'
shuffleTOD_df <- data.frame(matrix(ncol=0, nrow=length(tod_o)))

print('RhythmicityCode.R: Running create_shuffled_TOD function in parallel...')

cl <- makeClusterPSOCK(n.cores)
# Parallelized create_shuffled_TOD function in the range from 1:permutations for old cohort
clusterExport(cl, c("create_shuffled_TOD", "permutations", "tod_o"))
shuffleTOD_df <- do.call(cbind,parLapply(cl, X=1:permutations, fun=create_shuffled_TOD, tod_o))
stopCluster(cl)

print('RhythmicityCode.R: Running generate_null_data_filenames in parallel for old cohort...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing execution of the generate_null_data_filenames function over the range 1:permutations
clusterExport(cl, c("generate_null_data_filenames", "groupName"))
null_para_files <- do.call(rbind, parLapply(cl, X=1:permutations, fun=generate_null_data_filenames, groupName))
stopCluster(cl)

print('RhythmicityCode.R: Running generate_null_data_parameters function in parallel for old cohort...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing the generate_null_data_parameters function across the range 1:num_expressed_genes to append to dataframe null_pare
clusterExport(cl, c("generate_null_data_parameters", "parStartVal", "shuffleTOD_df", "nls.lm", "expr.old", "fitSinCurve"))
for(permutation in 1:permutations) {
  clusterExport(cl, "permutation")
  null_pare <- do.call(rbind, parLapply(cl, X=1:num_expressed_genes, fun=generate_null_data_parameters, permutations, shuffleTOD_df, expr.old, fitSinCurve))
  colnames(null_pare) <- c("A", "phase", "offset", "peak", "R2")
  save(null_pare, file=null_para_files[permutation,])
}
stopCluster(cl)

print('RhythmicityCode.R: Running collect_variable_from_permutations function for old cohort...')
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

###REPEAT NULL GENERATION FOR Young###
groupName <- 'young'
print('RhythmicityCode.R: Running create_shuffled_TOD function in parallel for young cohort...')

cl <- makeClusterPSOCK(n.cores)
# Parallelized create_shuffled_TOD function in the range from 1:permutations
clusterExport(cl, c("create_shuffled_TOD", "permutations", "tod_y"))
shuffleTOD_df <- do.call(cbind,parLapply(cl, X=1:permutations, fun=create_shuffled_TOD, tod_y))
stopCluster(cl)

print('RhythmicityCode.R: Running generate_null_data_filenames in parallel for young cohort...')
cl <- makeClusterPSOCK(n.cores)
# Parallelizing execution of the generate_null_data_filenames function over the range 1:permutations
clusterExport(cl, c("generate_null_data_filenames", "groupName"))
null_para_files <- do.call(rbind, parLapply(cl, X=1:permutations, fun=generate_null_data_filenames, groupName))
stopCluster(cl)

print('RhythmicityCode.R: Running generate_null_data_parameters function in parallel for young cohort...')
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

print('RhythmicityCode.R: Running collect_variable_from_permutations function for young cohort...')
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

R2change <- observed_para_o$R2 - observed_para_y$R2 
Ashift <- abs(abs(observed_para_o$A) - abs(observed_para_y$A))
phaseshift0 <- observed_para_o$phase - observed_para_y$phase
phaseshift <- pmin(abs(phaseshift0), 24 - abs(phaseshift0))
intshift <- abs(observed_para_o$offset - observed_para_y$offset)

R2changeNULL <- R2change
AshiftNULL <- Ashift
phaseshiftNULL <- phaseshift
intshiftNULL <- intshift

print('RhythmicityCode.R: Running collect_variable_from_permutations function...')
cl <- makeClusterPSOCK(n.cores)
# Parallelize the collect_variable_from_permutations function across the range 1:permutations.

for(variable in c('A', 'phase', 'offset', 'peak', 'R2')) {
  print(variable)
  clusterExport(cl, c("variable", "collect_variable_from_permutations", "null_para_files", "permutations"))
  assign(paste('null_para', variable, sep='_'), do.call(cbind, parLapply(cl, X=1:permutations, fun=collect_variable_from_permutations, variable)),envir=.GlobalEnv)
}




stopCluster(cl)

setwd('..')
