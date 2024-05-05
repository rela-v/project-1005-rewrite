setwd('./nullFolder')
# Get null_para_files for young and old,
#' compareGroups()
#' Comparison of variable between two groups 
#' @returns list of values of groups 1 and 2 fo\ the variable in question.
#' @param permutation represents the number of the permutation in question.
#' @param group1 represents the first group in the comparison.
#' @param group2 represents the second group in the comparison.
#' @param variable represents the variable to be grabbed.
#' @examples
#' compareGroups(1, 'old', 'young', variable)
compareGroups <- function(permutation, group1, group2, variable) {
  # Load groups 1 and 2
  group1_ret <- get(load(paste('./null_', group1, '_', permutation, '.rdata', sep='')))[variable]
  group2_ret <- get(load(paste('./null_', group2, '_', permutation, '.rdata', sep='')))[variable] 
  # Logic for p-value generation by Wiloxon rank sum test
  group_pvalues <- wilcox.test(as.numeric(unlist(group1_ret)), as.numeric(unlist(group2_ret)), alternative = "two.sided")$p.value
  return(group_pvalues)
}

cl <- makeClusterPSOCK(n.cores)
for(variable in c('R2', 'A', 'phase', 'offset')) {
  print(variable)
  clusterExport(cl, c("variable", "compareGroups", "permutations"))
  assign(paste(variable, 'shiftNULL', sep=''), do.call(cbind, parLapply(cl, X=1:permutations, fun=compareGroups, group1='old', group2='young', variable=variable)))
}
stopCluster(cl)

# Calculate p-values with two-sample independent t-test

p <- nrow(observed_para_o)
R2gainPvalue <- 1-rank(R2shiftNULL)[1:p]/length(R2shiftNULL) 
R2losePvalue <- rank(R2shiftNULL)[1:p]/length(R2shiftNULL)
AshiftPvalue <- 1 - rank(AshiftNULL)[1:p]/length(AshiftNULL)
phasePvalue <- 1 - rank(phaseshiftNULL)[1:p]/length(phaseshiftNULL)
offsetshiftPvalue <- 1 - rank(offsetshiftNULL)[1:p]/length(offsetshiftNULL)
result2<-data.frame(cbind(R2gainPvalue, R2losePvalue, AshiftPvalue, phasePvalue, offsetshiftPvalue))
row.names(result2)<-observed_para_o$genes
result2_sorted<-result2[order(result2$R2gainPvalue, decreasing = FALSE), ]
setwd("../")
write.csv(result2_sorted, "./Results/Rhythmicity/Example_Result2.csv")

print('RhythmicityCode.R: Routine complete.')
