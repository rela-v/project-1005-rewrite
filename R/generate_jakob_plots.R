# Generate Plots for Jakob Hartmann

install.packages("lmtest")
library(lmtest)
# Names of rows in gene expression matrix
expression_data_of_interest <- read.csv("./Data/Project_1005-rawCounts-annotated.csv")
cohort_data <- read.csv("./Data/cohort_data.csv", row.names = 1)
colnames(expression_data_of_interest)  <- c(rownames(cohort_data), 'Gene')
Symbols <- expression_data_of_interest$Gene

### Utility functions

# Observed dataframe of genes and fitted curve parameters
# observed_para_c <- data.frame(A=numeric(n_e), phase=numeric(n_e), offset=numeric(n_e), peak=numeric(n_e), R2=numeric(n_e))
# Time of death (clinical)
grab_metadata_expression_data <- function(relevant_condition) {
  filtered_cohort_data <- cohort_data[(cohort_data$DDx %in% relevant_condition),]
  filtered_expression_data <- expression_data_of_interest[colnames(expression_data_of_interest) %in% rownames(filtered_cohort_data)]
  return(list(filtered_cohort_data, filtered_expression_data))
}

generate_observed_parameters <- function(index, condition, num_expressed_genes, gene_expression_dataframe, zeitgeber_times) {
  observed_para <- data.frame(A=numeric(num_expressed_genes), phase=numeric(num_expressed_genes), offset=numeric(num_expressed_genes), peak=numeric(num_expressed_genes), R2=numeric(num_expressed_genes))
  observed_para[index,] <- fitSinCurve(xx=as.numeric(zeitgeber_times), observed=as.numeric(gene_expression_dataframe))
  out_row <- data.frame(A=observed_para$A[index], phase=observed_para$phase[index], offset=observed_para$offset[index], peak=observed_para$peak[index], R2=observed_para$R2[index])
  return(out_row)
}

#' generate_null_parameters()
#' Generate the null parameters for the fitted sinusoid in the expressed genes over shuffled TOD
#' @returns Dataframe containing "A" (amplitude), "phase" (phase), "offset" (x offset), "peak" (peak value), and "R2" (correlation).
#' @param index represents the index of the gene to be analyzed.
#' @param num_expressed_genes represents the number of genes there are to be analyzed in the dataframe.
#' @param gene_expression_dataframe represents the dataframe used as input for the analysis.
#' @param zeitgeber_times represents the Zeitgeber Time, or the relative timescale used for the circadian analysis.
#' @examples
#' out_pars <- generate_null_parameters(index=1, num_expressed_genes=10345, gene_expression_dataframe=df_with_gene_data, zeitgeber_time=df_with_ZT_data)
print('RhythmicityCode.R: Defining generate_null_parameters function...')
generate_null_parameters <- function(index, condition, num_expressed_genes, gene_expression_dataframe, zeitgeber_times) {
  shuffled_zeitgeber_times <- sample(zeitgeber_times, length(cohort_data_C$ZT), replace=TRUE)
  out <- lm(gene_expression_dataframe[index,] ~ shuffled_zeitgeber_times)
  out_row <- summary(out)$r.squared
  return(out_row)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

get_lrtest <- function(residuals, zeitgeber_times) {
  lrtest_results <- list()
  lrtest_results <- lrtest(lm(range01(residuals) ~ zeitgeber_times), lm(asin(range01(residuals)) ~ zeitgeber_times))
  lrtest_results_pvalue <- pchisq(lrtest_results$Chisq[2], df=lrtest_results$`#Df`[2], lower.tail=FALSE)
  return(lrtest_results_pvalue)
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

expression_data_C <- unlist(data_C[[2]][15150,])
expression_data_SZ <- unlist(data_SZ[[2]][15150,])
expression_data_BD <- unlist(data_BD[[2]][15150,])
expression_data_MDD <- unlist(data_MDD[[2]][15150,])


# List of tasks:
# 1. have the SKA2 expression values for each subject in excel
# 2. export them in a higher resolution
# 3. methods section
# 4. send the figure without covariate effects

# SKA2 expression values for each subject before removing covariates
ska2_control_expression <- expression_data_C
dir.create("./Results/C")
write.csv(ska2_control_expression, file="./Results/C/ska2_control_expression.csv")
ska2_sz_expression <- expression_data_SZ
dir.create("./Results/SZ")
write.csv(ska2_sz_expression, file="./Results/SZ/ska2_sz_expression.csv")
ska2_bd_expression <- expression_data_BD
dir.create("./Results/BD")
write.csv(ska2_bd_expression, file="./Results/BD/ska2_bd_expression.csv")
ska2_mdd_expression <- expression_data_MDD
dir.create("./Results/MDD")
write.csv(ska2_mdd_expression, file="./Results/MDD/ska2_mdd_expression.csv")

# SKA2 expression values for each subject after removing covariates
# File path template:   save(residual, file=paste("./results/variance_partitioning_plots/", condition, "/residualmodelfit_", condition, ".rdata", sep=''))
load("./results/variance_partitioning_plots/C/residualmodelfit_C.rdata")
ska2_control_expression_no_covariates <- residual[1,]
write.csv(ska2_control_expression_no_covariates, file="./Results/C/ska2_control_expression_no_covariates.csv")
load("./results/variance_partitioning_plots/SZ/residualmodelfit_SZ.rdata")
ska2_sz_expression_no_covariates <- residual[1,]
write.csv(ska2_sz_expression_no_covariates, file="./Results/SZ/ska2_sz_expression_no_covariates.csv")
load("./results/variance_partitioning_plots/BD/residualmodelfit_BD.rdata")
ska2_bd_expression_no_covariates <- residual[1,]
write.csv(ska2_bd_expression_no_covariates, file="./Results/BD/ska2_bd_expression_no_covariates.csv")
load("./results/variance_partitioning_plots/MDD/residualmodelfit_MDD.rdata")
ska2_mdd_expression_no_covariates <- residual[1,]
write.csv(ska2_mdd_expression_no_covariates, file="./Results/MDD/ska2_mdd_expression_no_covariates.csv")

# Generate plots for Jakob in high resolution for SKA2 expression values for control group
png("./Results/C/ska2_control_diurnal_expression.png", width=800, height=800)
# circadiandrawing <- function(tod, expr, apar, labels, specinfo=null){	
apar_C_w <- generate_observed_parameters(index=1, condition="C", num_expressed_genes=1, gene_expression_dataframe=expression_data_C, zeitgeber_times=cohort_data_C$ZT)
apar_SZ_w <- generate_observed_parameters(index=1, condition="SZ", num_expressed_genes=1, gene_expression_dataframe=expression_data_SZ, zeitgeber_times=cohort_data_SZ$ZT)
apar_BD_w <- generate_observed_parameters(index=1, condition="BD", num_expressed_genes=1, gene_expression_dataframe=expression_data_BD, zeitgeber_times=cohort_data_BD$ZT)
apar_MDD_w <- generate_observed_parameters(index=1, condition="MDD", num_expressed_genes=1, gene_expression_dataframe=expression_data_MDD, zeitgeber_times=cohort_data_MDD$ZT)
xlim <- c(min(c(cohort_data_C$ZT, cohort_data_SZ$ZT, cohort_data_BD$ZT, cohort_data_MDD$ZT)), max(c(cohort_data_C$ZT, cohort_data_SZ$ZT, cohort_data_BD$ZT, cohort_data_MDD$ZT)))
ylim <- c(min(c(expression_data_C, expression_data_SZ, expression_data_BD, expression_data_MDD))-50, max(c(expression_data_C, expression_data_SZ, expression_data_BD, expression_data_MDD))+50)
circadianDrawing(tod=cohort_data_C$ZT, expr=expression_data_C, specInfo=paste("Diurnal Expression of SKA2 in Control Group"),
                 apar=apar_C_w, labels=cohort_data_C$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for SZ group
png("./Results/SZ/ska2_sz_diurnal_expression.png", width=800, height=800)

circadianDrawing(tod=cohort_data_SZ$ZT, expr=expression_data_SZ, specInfo=paste("Diurnal Expression of SKA2 in Schizophrenia Group"),
                 apar=apar_SZ_w, labels=cohort_data_SZ$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for BD group
png("./Results/BD/ska2_bd_diurnal_expression.png", width=800, height=800)

circadianDrawing(tod=cohort_data_BD$ZT, expr=expression_data_BD, specInfo=paste("Diurnal Expression of SKA2 in Bipolar Disorder Group"),
                 apar=apar_BD_w, labels=cohort_data_BD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for MDD group
png("./Results/MDD/ska2_mdd_diurnal_expression.png", width=800, height=800)

circadianDrawing(tod=cohort_data_MDD$ZT, expr=expression_data_MDD, specInfo=paste("Diurnal Expression of SKA2 in Major Depressive Disorder Group"),
                 apar=apar_MDD_w, labels=cohort_data_MDD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for control group after removing covariates
png("./Results/C/ska2_control_diurnal_expression_no_covariates.png", width=800, height=800)
apar_C_n <- generate_observed_parameters(index=1, condition="C", num_expressed_genes=1, gene_expression_dataframe=ska2_control_expression_no_covariates, zeitgeber_times=cohort_data_C$ZT)
apar_SZ_n <- generate_observed_parameters(index=1, condition="SZ", num_expressed_genes=1, gene_expression_dataframe=ska2_sz_expression_no_covariates, zeitgeber_times=cohort_data_SZ$ZT)
apar_BD_n <- generate_observed_parameters(index=1, condition="BD", num_expressed_genes=1, gene_expression_dataframe=ska2_bd_expression_no_covariates, zeitgeber_times=cohort_data_BD$ZT)
apar_MDD_n <- generate_observed_parameters(index=1, condition="MDD", num_expressed_genes=1, gene_expression_dataframe=ska2_mdd_expression_no_covariates, zeitgeber_times=cohort_data_MDD$ZT)
xlim <- c(-6, 18)
ylim <- c(min(c(ska2_control_expression_no_covariates, ska2_sz_expression_no_covariates, ska2_bd_expression_no_covariates, ska2_mdd_expression_no_covariates))-50, max(c(ska2_control_expression_no_covariates, ska2_sz_expression_no_covariates, ska2_bd_expression_no_covariates, ska2_mdd_expression_no_covariates))+50)
circadianDrawing(tod=cohort_data_C$ZT, expr=ska2_control_expression_no_covariates, paste("Diurnal Expression of SKA2 in Control Group After Removing Covariates"),
                 apar=apar_C_n, labels=cohort_data_C$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for SZ group after removing covariates
png("./Results/SZ/ska2_sz_diurnal_expression_no_covariates.png", width=800, height=800)

circadianDrawing(tod=cohort_data_SZ$ZT, expr=ska2_sz_expression_no_covariates, paste("Diurnal Expression of SKA2 in Schizophrenia Group After Removing Covariates"),
                 apar=apar_SZ_n, labels=cohort_data_SZ$DDx, xlim=xlim, ylim=ylim)

dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for BD group after removing covariates
png("./Results/BD/ska2_bd_diurnal_expression_no_covariates.png", width=800, height=800)

circadianDrawing(tod=cohort_data_BD$ZT, expr=ska2_bd_expression_no_covariates, paste("Diurnal Expression of SKA2 in Bipolar Disorder Group After Removing Covariates"),
                 apar=apar_BD_n, labels=cohort_data_BD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for MDD group after removing covariates
png("./Results/MDD/ska2_mdd_diurnal_expression_no_covariates.png", width=800, height=800)

circadianDrawing(tod=cohort_data_MDD$ZT, expr=ska2_mdd_expression_no_covariates, paste("Diurnal Expression of SKA2 in Major Depressive Disorder Group After Removing Covariates"),
                 apar=apar_MDD_n, labels=cohort_data_MDD$DDx, xlim=xlim, ylim=ylim)
dev.off()

print(paste("R2 for SKA2 in Control Group: ", apar_C_w$R2, "R2 for SKA2 in Schizophrenia Group: ", apar_SZ_w$R2, "R2 for SKA2 in Bipolar Disorder Group: ", apar_BD_w$R2, "R2 for SKA2 in Major Depressive Disorder Group: ", apar_MDD_w$R2))
print(paste("p-value for SKA2 in Control Group: ", get_lrtest(unlist(ska2_control_expression), cohort_data_C$ZT), "p-value for SKA2 in Schizophrenia Group: ", get_lrtest(unlist(ska2_sz_expression), cohort_data_SZ$ZT), "p-value for SKA2 in Bipolar Disorder Group: ", get_lrtest(unlist(ska2_bd_expression), cohort_data_BD$ZT), "p-value for SKA2 in Major Depressive Disorder Group: ", get_lrtest(unlist(ska2_mdd_expression), cohort_data_MDD$ZT)))

print(paste("R2 for SKA2 in Control Group After Removing Covariates: ", apar_C_n$R2, "R2 for SKA2 in Schizophrenia Group After Removing Covariates: ", apar_SZ_n$R2, "R2 for SKA2 in Bipolar Disorder Group After Removing Covariates: ", apar_BD_n$R2, "R2 for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", apar_MDD_n$R2))
print(paste("p-value for SKA2 in Control Group After Removing Covariates: ", get_lrtest(unlist(ska2_control_expression_no_covariates), cohort_data_C$ZT), "p-value for SKA2 in Schizophrenia Group After Removing Covariates: ", get_lrtest(unlist(ska2_sz_expression_no_covariates), cohort_data_SZ$ZT), "p-value for SKA2 in Bipolar Disorder Group After Removing Covariates: ", get_lrtest(unlist(ska2_bd_expression_no_covariates), cohort_data_BD$ZT), "p-value for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", get_lrtest(unlist(ska2_mdd_expression_no_covariates), cohort_data_MDD$ZT)))
