print("")
library('minpack.lm')
source('./R/fitSinCurve.R')
cohort_data<-read.csv('./Data/revised_cohort_data.csv', row.names=1)
expression_data<-read.csv('./Data/revised_expression_data.csv', row.names=1)
n.cores <- availableCores()
# Number of rows in gene expression matrix
n_e<-nrow(expression_data)
print(n_e)
# Names of rows in gene expression matrix
Symbols<-row.names(expression_data)
# Observed dataframe of genes and fitted curve parameters
# observed_para_c <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
# Time of death (clinical)
ZT<-cohort_data$ZT
# 1->number of rows: fit the sine curves, assign the output to each row.
`nls.lm` <- minpack.lm::nls.lm
function_1 <- function(index, data_0, data_1, data_2) {
  out <- fitSinCurve(xx=as.numeric(ZT),observed=as.numeric(expression_data[index,]))
  out_row <- data.frame(Symbols=rownames(expression_data[index,]), A=out$A, phase=out$phase, out=out$offset, peak=out$peak, R2=out$R2)
  return(out_row)
}
cl <- makeClusterPSOCK(n.cores)
clusterExport(cl, c("function_1", "fitSinCurve", "n_e", "expression_data", "ZT", "nls.lm", "parStartVal"))
observed_para_c <- do.call(rbind,parLapply(cl, X=1:n_e, fun=function_1, data_0=n_e, data_1=expression_data, data_2=ZT))
stopCluster(cl)

# 1->number of rows: fit the sine curves, assign the output to each row.
observed_para_c_sorted<-observed_para_c[order(observed_para_c$R2, decreasing = TRUE),]

# Generate Scatter plots 
top.control<-observed_para_c_sorted[1:4,]
row_index<-as.numeric(row.names(top.control))
n_t <- nrow(top.control)
labels1<-cohort_data$suicide
specInfo<-"Control"
system("mkdir -p ./Results/Rhythmicity/PDF")

cl <- makeClusterPSOCK(n.cores)

# foreach(i=1:nrow(top.control)) %dopar% {
function_2 <- function(index, data_0, function_0, data_1, data_2, data_3, data_4, data_5) {
  agene <- as.character(data_0$genes[index])
  fileName <- paste('./Results/Rhythmicity/PDF/Top_Control',agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=data_1$ZT, expr=data_2[data_3[index],], apar=data_0[index,],labels=data_4, specInfo=data_5)
  dev.off()
}
clusterExport(cl, c("function_2", "top.control", "circadianDrawing", "cohort_data", "expression_data", "row_index", "labels1", "specInfo", "n_t"))
parLapply(cl, X=1:n_t, fun=function_2, data_0=top.control, function_0=circadianDrawing, data_1=cohort_data, data_2=expression_data, data_3=row_index, data_4=labels1, data_5=specInfo)
stopCluster(cl)
##Generate p/q values 
system("mkdir -p ./nullFolder")
setwd("./nullFolder/")
groupName <- 'control'
B<-10 #we use 1000 permutations in the paper, but to save time and computation power,
# the results in this folder reflect 10 permutations
# shuffleTOD_df <- data.frame(matrix(ncol=0, nrow=length(ZT)))
# index = the index used to identify which "run" this is
# out_df = the output dataframe
# data_0 = ZT
#cl <- makeClusterPSOCK(n.cores)
#function_3 <- function(index, out_df, data_0) {
#  set.seed(index)
#  return(data.frame(sample(data_0)))
#}
#clusterExport(cl, c("function_3", "shuffleTOD_df", "B", "ZT"))
#shuffleTOD_df <- do.call(cbind,parLapply(cl, X=1:B, fun=function_3, out_df=shuffleTOD_df, data_0=ZT))
#stopCluster(cl)

# index = index of run
# data_0 = shuffleTOD_df
# fn_0 = fitSinCurve
# data_1 = expression_data

#function_4 <- function(index, super_index, data_0, data_1, fn_0) {
#  out <- fn_0(xx=data_0[,super_index], observed=unlist(data_1[index,]))
#  out_row <- data.frame(out$A, out$phase, out$offset, out$peak, out$R2)
#  return(out_row)
#}
#function_4b <- function(index, data_0) {
#  return(paste('null_', data_0, '_', index, '.rdata', sep=''))
#}
#
#
#cl <- makeClusterPSOCK(n.cores)
#clusterExport(cl, c("B", "function_4b", "groupName"))
#null_para_files <- do.call(rbind, parLapply(cl, X=1:B, fun=function_4b, data_0=groupName))
#colnames(null_pare) <- c("A", "phase", "offset", "peak", "R2")
#stopCluster(cl)
#print(head(null_para_files,10))

#for(b in 1:B) {
#  cl <- makeClusterPSOCK(n.cores)
#  clusterExport(cl, c("b", "nls.lm", "n_e", "shuffleTOD_df", "fitSinCurve", "function_4", "expression_data", "groupName", "parStartVal"))
#  null_pare <- do.call(rbind, parLapply(cl, X=1:n_e, fun=function_4, super_index=b, data_0=shuffleTOD_df, data_1=expression_data, fn_0=fitSinCurve))
#  colnames(null_pare) <- c("A", "phase", "offset", "peak", "R2")
#  save(null_pare, file=null_para_files[b])
#  stopCluster(cl)
#  print(b)
#}

 
#null_pare_A <- data.frame(matrix(0, n_e, B))
#null_pare_phase <- data.frame(matrix(0, n_e, B))
#null_pare_offset <- data.frame(matrix(0, n_e, B))
#null_pare_peak <- data.frame(matrix(0, n_e, B))
#null_pare_R2 <- data.frame(matrix(0, n_e, B))


# foreach(i=1:B) %dopar% {
# index = index of run
# data_0 = groupName
#function_5 <- function(index, data_0) {
#  print(index)
#  null_pare <- get(load(null_para_files[index])) # Load
#  return(null_pare[data_0])
#}
#cl <- makeClusterPSOCK(n.cores)
#clusterExport(cl, c("B", "function_5", "null_para_files"))
#for(i in c('A', 'phase', 'offset', 'peak', 'R2')) {
#  assign(paste('null_pare_', i, sep=''), do.call(cbind,parLapply(cl, X=1:B, fun=function_5, data_0=i)))
#}
#stopCluster(cl)
#null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
#                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
#null_para_file <- paste('null_',groupName,'.rdata',sep='')
#save(null_para,file=null_para_file)
#setwd('..')
#null_para <- get(load("./null_control.rdata"))
#para_R2_pool <- cbind(observed_para_c$R2,null_para$null_para_R2)
#R2Rank_para <- 1 - (rank(para_R2_pool)[1:length(observed_para_c$R2)] - 0.5)/length(para_R2_pool)
#observed_para_c$pvalue <- R2Rank_para
#observed_para_c$qvalue <- p.adjust(observed_para_c$pvalue, 'BH')
#observe_para_c_sorted<-observed_para_c[order(observed_para_c$pvalue),]
#setwd("../")
#write.csv(observe_para_c_sorted, "./Results/Rhythmicity/Example_result.csv")
#
##Rhythmicity gain/loss example 
#Split data into two groups: For example we will split the controls into younger and older 
#old_index<-which(cohort_data$age>=50)
#young_index<-which(cohort_data$age<50)
#
#expr.old<-expression_data[,old_index]
#expr.young<-expression_data[,young_index]
#
#tod_o<-cohort_data$ZT[old_index]
#tod_y<-cohort_data$ZT[young_index]

# observed_para_y <- data.frame(genes=Symbols,A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
# 
# foreach(i=1:n_e) %dopar% {
#   out <- fitSinCurve(xx=tod_y,observed=as.numeric(expr.young[i,]))
#   observed_para_y[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
#   if(i%%1000==0) print(i)
# }

# observed_para_o <- data.frame(genes=Symbols,A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))

# foreach(i=1:n_e) %dopar% {
#   out <- fitSinCurve(xx=tod_o,observed=as.numeric(expr.old[i,]))
#   observed_para_o[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
#   if(i%%1000==0) print(i)
# }
# 
# #Generate Null Distributions: 
# setwd('./nullFolder')
# library(doParallel)
# groupName <- 'old'
# thisData <- expr.old
# B<-10
# result <- foreach(b = 1:B) %dopar% {
#  set.seed(b)
#   print(b)	
#   library(minpack.lm)
#   source('../R/fitSinCurve.R')
#   
#   null_pare <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
#   null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
# 
#  shuffleTOD <- sample(tod_o)
#   
#  for (i in 1:n_e) {
#     out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
#     null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
#   }
#   save(null_pare,file=null_para_file)	
# }
# 
# null_pare_A <- matrix(0,n_e,B)
# null_pare_phase <- matrix(0,n_e,B)
# null_pare_offset <- matrix(0,n_e,B)
# null_pare_peak <- matrix(0,n_e,B)
# null_pare_R2 <- matrix(0,n_e,B)
# 
# foreach(b=1:B) %dopar% {
#   set.seed(b)
#   print(b)
#   file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
#   
#   load(file11)
#   null_pare_A[,b] <- null_pare$A
#   null_pare_phase[,b] <- null_pare$phase
#   null_pare_offset[,b] <- null_pare$offset
#   null_pare_peak[,b] <- null_pare$peak
#   null_pare_R2[,b] <- null_pare$R2		
# }
# null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
#                   null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
# 
# null_para_file <- paste('null_',groupName,'.rdata',sep='')
# save(null_para,file=null_para_file)
# 
# ###REPEAT NULL GENERATION FOR Young###
# groupName <- 'young'
# thisData <- expr.young
# B<-10
# result <- foreach(b = 1:B) %dopar% {
#   set.seed(b)
#   print(b)	
#   library(minpack.lm)
#   source('../R/fitSinCurve.R')
#   
#   null_pare <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
#   null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
#   shuffleTOD <- sample(tod_y)
#   
#   for (i in 1:n_e) {
#     out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
#     null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
#   }		
#   save(null_pare,file=null_para_file)	
# }
# 
# null_pare_A <- matrix(0,n_e,B)
# null_pare_phase <- matrix(0,n_e,B)
# null_pare_offset <- matrix(0,n_e,B)
# null_pare_peak <- matrix(0,n_e,B)
# null_pare_R2 <- matrix(0,n_e,B)
# 
# foreach(b=1:B) %dopar% {
#   set.seed(b)
#   print(b)
#   file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
 #  
#   load(file11)
 #  null_pare_A[,b] <- null_pare$A
#   null_pare_phase[,b] <- null_pare$phase
#  null_pare_offset[,b] <- null_pare$offset
#   null_pare_peak[,b] <- null_pare$peak
#   null_pare_R2[,b] <- null_pare$R2		
# }
# 
# null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
 #                 null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
# 
# null_para_file <- paste('null_',groupName,'.rdata',sep='')
# save(null_para,file=null_para_file)
# 
# R2change <- observed_para_o$R2 - observed_para_y$R2 
# Ashift <- abs(abs(observed_para_o$A) - abs(observed_para_y$A))
# phaseshift0 <- observed_para_o$phase - observed_para_y$phase
# phaseshift <- pmin(abs(phaseshift0), 24 - abs(phaseshift0))
# intshift <- abs(observed_para_o$offset - observed_para_y$offset)
# 
# R2changeNULL <- R2change
# AshiftNULL <- Ashift
# phaseshiftNULL <- phaseshift
# intshiftNULL <- intshift
# 
# B <- 10
# foreach(b=1:B) %dopar% {
#   set.seed(b)
#   if(b%%50 == 0) print(b)
#   file_old <- paste0("./",'null_old_',b,'.rdata')
#   file_young <- paste0("./",'null_young_',b,'.rdata')
#   null_para_old <- get(load(file_old))
#   null_para_young <- get(load(file_young))
#   
#   aR2change <- null_para_old$R2 - null_para_young$R2
#   aAshift <- abs(abs(null_para_old$A) - abs(null_para_young$A))
#   aphaseshift <- pmin(abs(null_para_old$phase - null_para_young$phase),24 - abs(null_para_old$phase - null_para_young$phase))
#   aintshift <- abs(null_para_old$offset - null_para_young$offset)
#   
#   R2changeNULL <- c(R2changeNULL,aR2change)
#   AshiftNULL <- c(AshiftNULL,aAshift)
#   phaseshiftNULL <- c(phaseshiftNULL,aphaseshift)
#   intshiftNULL <- c(intshiftNULL,aintshift)
# }

# p <- nrow(observed_para_o)

# R2gainPvalue <- 1-rank(R2changeNULL)[1:p]/length(R2changeNULL) 
# R2losePvalue <- rank(R2changeNULL)[1:p]/length(R2changeNULL) 
# AshiftPvalue <- 1 - rank(AshiftNULL)[1:p]/length(AshiftNULL)
# phasePvalue <- 1 - rank(phaseshiftNULL)[1:p]/length(phaseshiftNULL)
# intshiftPvalue <- 1 - rank(intshiftNULL)[1:p]/length(intshiftNULL)

# result2<-data.frame(cbind(R2gainPvalue, R2losePvalue, AshiftPvalue, phasePvalue, intshiftPvalue))
# row.names(result2)<-observed_para_o$genes
# result2_sorted<-result2[order(result2$R2gainPvalue, decreasing = FALSE), ]
# setwd("../")
# write.csv(result2_sorted, "./Results/Rhythmicity/Example_Result2.csv")
#print("RhythmicityCode.R: Protocol complete.")
