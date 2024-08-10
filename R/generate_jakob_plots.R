# Generate Plots for Jakob Hartmann

install.packages("lmtest")
library(lmtest)
library(readxl)
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

generate_observed_parameters <- function(gene_expression_dataframe, zeitgeber_times) {
  observed_para <- fitSinCurve(xx=as.numeric(zeitgeber_times), observed=as.numeric(gene_expression_dataframe))
  out_row <- data.frame(A=observed_para$A, phase=observed_para$phase, offset=observed_para$offset, peak=observed_para$peak, R2=observed_para$R2)
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
# 1. Create function to shuffle zeitgeber times
# 2. Get null parameters overall for all 4 conditions
# 3. Get p-value for each condition from null distribution of parameters vs. observed parameters

get_shuffled_data <- function(expression_data, zeitgeber_times) {
  shuffled_zeitgeber_times <- sample(zeitgeber_times, length(expression_data), replace=TRUE)
  return(shuffled_zeitgeber_times)
}

t_expression_data <- c(expression_data_C, expression_data_SZ, expression_data_BD, expression_data_MDD)
t_ZT <- c(cohort_data_C$ZT, cohort_data_SZ$ZT, cohort_data_BD$ZT, cohort_data_MDD$ZT)

get_pvalue <- function(expression_data_o, ZT_o, expression_data_n=t_expression_data, ZT_n=t_ZT) {
  n_e <- 10000
  null_parameters <- do.call(rbind, lapply(1:n_e, function(x) generate_observed_parameters(gene_expression_dataframe=expression_data_n, zeitgeber_times=get_shuffled_data(expression_data_n,ZT_n))))
  observed_parameters <- generate_observed_parameters(gene_expression_dataframe=expression_data_o, zeitgeber_times=ZT_o)
  p_value <- length(which(null_parameters$R2 > observed_parameters$R2))/n_e 
  return(p_value)
}

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
pdf("./Results/C/ska2_control_diurnal_expression.pdf")
# circadiandrawing <- function(tod, expr, apar, labels, specinfo=null){	
apar_C_w <- generate_observed_parameters(gene_expression_dataframe=expression_data_C, zeitgeber_times=cohort_data_C$ZT)
apar_SZ_w <- generate_observed_parameters(gene_expression_dataframe=expression_data_SZ, zeitgeber_times=cohort_data_SZ$ZT)
apar_BD_w <- generate_observed_parameters(gene_expression_dataframe=expression_data_BD, zeitgeber_times=cohort_data_BD$ZT)
apar_MDD_w <- generate_observed_parameters(gene_expression_dataframe=expression_data_MDD, zeitgeber_times=cohort_data_MDD$ZT)
xlim <- c(min(c(cohort_data_C$ZT, cohort_data_SZ$ZT, cohort_data_BD$ZT, cohort_data_MDD$ZT)), max(c(cohort_data_C$ZT, cohort_data_SZ$ZT, cohort_data_BD$ZT, cohort_data_MDD$ZT)))
ylim <- c(min(c(expression_data_C, expression_data_SZ, expression_data_BD, expression_data_MDD))-50, max(c(expression_data_C, expression_data_SZ, expression_data_BD, expression_data_MDD))+50)
circadianDrawing(tod=cohort_data_C$ZT, expr=expression_data_C, specInfo=paste("Diurnal Expression of SKA2 in Control Group"),
                 apar=apar_C_w, labels=cohort_data_C$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for SZ group
pdf("./Results/SZ/ska2_sz_diurnal_expression.pdf")

circadianDrawing(tod=cohort_data_SZ$ZT, expr=expression_data_SZ, specInfo=paste("Diurnal Expression of SKA2 in Schizophrenia Group"),
                 apar=apar_SZ_w, labels=cohort_data_SZ$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for BD group
pdf("./Results/BD/ska2_bd_diurnal_expression.pdf")

circadianDrawing(tod=cohort_data_BD$ZT, expr=expression_data_BD, specInfo=paste("Diurnal Expression of SKA2 in Bipolar Disorder Group"),
                 apar=apar_BD_w, labels=cohort_data_BD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for MDD group
pdf("./Results/MDD/ska2_mdd_diurnal_expression.pdf")

circadianDrawing(tod=cohort_data_MDD$ZT, expr=expression_data_MDD, specInfo=paste("Diurnal Expression of SKA2 in Major Depressive Disorder Group"),
                 apar=apar_MDD_w, labels=cohort_data_MDD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for control group after removing covariates
pdf("./Results/C/ska2_control_diurnal_expression_no_covariates.pdf")
apar_C_n <- generate_observed_parameters(gene_expression_dataframe=ska2_control_expression_no_covariates, zeitgeber_times=cohort_data_C$ZT)
apar_SZ_n <- generate_observed_parameters(gene_expression_dataframe=ska2_sz_expression_no_covariates, zeitgeber_times=cohort_data_SZ$ZT)
apar_BD_n <- generate_observed_parameters(gene_expression_dataframe=ska2_bd_expression_no_covariates, zeitgeber_times=cohort_data_BD$ZT)
apar_MDD_n <- generate_observed_parameters(gene_expression_dataframe=ska2_mdd_expression_no_covariates, zeitgeber_times=cohort_data_MDD$ZT)
xlim <- c(-6, 18)
ylim <- c(min(c(ska2_control_expression_no_covariates, ska2_sz_expression_no_covariates, ska2_bd_expression_no_covariates, ska2_mdd_expression_no_covariates))-50, max(c(ska2_control_expression_no_covariates, ska2_sz_expression_no_covariates, ska2_bd_expression_no_covariates, ska2_mdd_expression_no_covariates))+50)
circadianDrawing(tod=cohort_data_C$ZT, expr=ska2_control_expression_no_covariates, paste("Diurnal Expression of SKA2 in Control Group After Removing Covariates"),
                 apar=apar_C_n, labels=cohort_data_C$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for SZ group after removing covariates
pdf("./Results/SZ/ska2_sz_diurnal_expression_no_covariates.pdf")

circadianDrawing(tod=cohort_data_SZ$ZT, expr=ska2_sz_expression_no_covariates, paste("Diurnal Expression of SKA2 in Schizophrenia Group After Removing Covariates"),
                 apar=apar_SZ_n, labels=cohort_data_SZ$DDx, xlim=xlim, ylim=ylim)

dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for BD group after removing covariates
pdf("./Results/BD/ska2_bd_diurnal_expression_no_covariates.pdf")

circadianDrawing(tod=cohort_data_BD$ZT, expr=ska2_bd_expression_no_covariates, paste("Diurnal Expression of SKA2 in Bipolar Disorder Group After Removing Covariates"),
                 apar=apar_BD_n, labels=cohort_data_BD$DDx, xlim=xlim, ylim=ylim)
dev.off()

# Generate plots for Jakob in high resolution for SKA2 expression values for MDD group after removing covariates
pdf("./Results/MDD/ska2_mdd_diurnal_expression_no_covariates.pdf")

circadianDrawing(tod=cohort_data_MDD$ZT, expr=ska2_mdd_expression_no_covariates, paste("Diurnal Expression of SKA2 in Major Depressive Disorder Group After Removing Covariates"),
                 apar=apar_MDD_n, labels=cohort_data_MDD$DDx, xlim=xlim, ylim=ylim)
dev.off()



print(paste("R2 for SKA2 in Control Group: ", apar_C_w$R2, "R2 for SKA2 in Schizophrenia Group: ", apar_SZ_w$R2, "R2 for SKA2 in Bipolar Disorder Group: ", apar_BD_w$R2, "R2 for SKA2 in Major Depressive Disorder Group: ", apar_MDD_w$R2))
print(paste("Amplitude for SKA2 in Control Group: ", apar_C_w$A, "Amplitude for SKA2 in Schizophrenia Group: ", apar_SZ_w$A, "Amplitude for SKA2 in Bipolar Disorder Group: ", apar_BD_w$A, "Amplitude for SKA2 in Major Depressive Disorder Group: ", apar_MDD_w$A))
print(paste("Phase for SKA2 in Control Group: ", apar_C_w$phase, "Phase for SKA2 in Schizophrenia Group: ", apar_SZ_w$phase, "Phase for SKA2 in Bipolar Disorder Group: ", apar_BD_w$phase, "Phase for SKA2 in Major Depressive Disorder Group: ", apar_MDD_w$phase))
print(paste("p-value for SKA2 in Control Group: ", get_pvalue(unlist(ska2_control_expression), cohort_data_C$ZT), "p-value for SKA2 in Schizophrenia Group: ", get_pvalue(unlist(ska2_sz_expression), cohort_data_SZ$ZT), "p-value for SKA2 in Bipolar Disorder Group: ", get_pvalue(unlist(ska2_bd_expression), cohort_data_BD$ZT), "p-value for SKA2 in Major Depressive Disorder Group: ", get_pvalue(unlist(ska2_mdd_expression), cohort_data_MDD$ZT)))

print(paste("R2 for SKA2 in Control Group After Removing Covariates: ", apar_C_n$R2, "R2 for SKA2 in Schizophrenia Group After Removing Covariates: ", apar_SZ_n$R2, "R2 for SKA2 in Bipolar Disorder Group After Removing Covariates: ", apar_BD_n$R2, "R2 for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", apar_MDD_n$R2))
print(paste("Amplitude for SKA2 in Control Group After Removing Covariates: ", apar_C_n$A, "Amplitude for SKA2 in Schizophrenia Group After Removing Covariates: ", apar_SZ_n$A, "Amplitude for SKA2 in Bipolar Disorder Group After Removing Covariates: ", apar_BD_n$A, "Amplitude for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", apar_MDD_n$A))
print(paste("Phase for SKA2 in Control Group After Removing Covariates: ", apar_C_n$phase, "Phase for SKA2 in Schizophrenia Group After Removing Covariates: ", apar_SZ_n$phase, "Phase for SKA2 in Bipolar Disorder Group After Removing Covariates: ", apar_BD_n$phase, "Phase for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", apar_MDD_n$phase))
print(paste("p-value for SKA2 in Control Group After Removing Covariates: ", get_pvalue(unlist(ska2_control_expression_no_covariates), cohort_data_C$ZT), "p-value for SKA2 in Schizophrenia Group After Removing Covariates: ", get_pvalue(unlist(ska2_sz_expression_no_covariates), cohort_data_SZ$ZT), "p-value for SKA2 in Bipolar Disorder Group After Removing Covariates: ", get_pvalue(unlist(ska2_bd_expression_no_covariates), cohort_data_BD$ZT), "p-value for SKA2 in Major Depressive Disorder Group After Removing Covariates: ", get_pvalue(unlist(ska2_mdd_expression_no_covariates), cohort_data_MDD$ZT)))

### Mouse analysis
# Load the data
data <- read_excel("./Data/mouseSKA2_HIP_updated_JH.xlsx", sheet = 1)
vehicle_data <- data[which(data[[2]] == "vehicle"),]
dex_data <- data[which(data[[2]] == "Dex"),]
lfchange <- data[[4]]
zeitgeber_times <- data[[3]]
lfchange_vehicle <- vehicle_data[[4]]
zeitgeber_times_vehicle <- vehicle_data[[3]]
lfchange_dex <- dex_data[[4]]
zeitgeber_times_dex <- dex_data[[3]]
# Generate plots for Jakob in high resolution for SKA2 expression values for mouse data
pdf("./Results/mouse_ska2_diurnal_expression_vehicle.pdf")
xlim <- c(min(zeitgeber_times), max(zeitgeber_times))
ylim <- c(min(lfchange)-sd(lfchange), max(lfchange)+sd(lfchange))
apar_mouse_vehicle <- generate_observed_parameters(gene_expression_dataframe=lfchange_vehicle, zeitgeber_times=zeitgeber_times_vehicle)
circadianDrawing(tod=zeitgeber_times_vehicle, expr=lfchange_vehicle, specInfo=paste("Diurnal Expression of SKA2 in Mouse Hippocampus: Vehicle"),
                 apar=apar_mouse_vehicle, labels=rep("Mouse", length(lfchange_vehicle)), xlim=xlim, ylim=ylim)
dev.off()

print(paste("R2 for SKA2 in Mouse Hippocampus for Vehicle: ", apar_mouse_vehicle$R2))
print(paste("Amplitude for SKA2 in Mouse Hippocampus for Vehicle: ", apar_mouse_vehicle$A))
print(paste("Phase for SKA2 in Mouse Hippocampus for Vehicle: ", apar_mouse_vehicle$phase))
print(paste("p-value for SKA2 in Mouse Hippocampus for Vehicle: ", get_pvalue(lfchange_vehicle, zeitgeber_times_vehicle)))

pdf("./Results/mouse_ska2_diurnal_expression_dex.pdf")
xlim <- c(min(zeitgeber_times), max(zeitgeber_times))
ylim <- c(min(lfchange)-sd(lfchange), max(lfchange)+sd(lfchange))
apar_mouse_dex <- generate_observed_parameters(gene_expression_dataframe=lfchange_dex, zeitgeber_times=zeitgeber_times_dex)
circadianDrawing(tod=zeitgeber_times_dex, expr=lfchange_dex, specInfo=paste("Diurnal Expression of SKA2 in Mouse Hippocampus: Dexamethasone"), apar=apar_mouse_dex, labels=rep("Mouse", length(lfchange_dex)), xlim=xlim, ylim=ylim)

dev.off()
print(paste("R2 for SKA2 in Mouse Hippocampus for Dex: ", apar_mouse_dex$R2))
print(paste("Amplitude for SKA2 in Mouse Hippocampus for Dex: ", apar_mouse_dex$A))
print(paste("Phase for SKA2 in Mouse Hippocampus for Dex: ", apar_mouse_dex$phase))
print(paste("p-value for SKA2 in Mouse Hippocampus for Dex: ", get_pvalue(lfchange_dex, zeitgeber_times_dex)))

# Generate null distribution of R2 values
n_e <- 1000
null_parameters_vehicle <- do.call(rbind, lapply(1:n_e, function(x) generate_observed_parameters(gene_expression_dataframe=lfchange_vehicle, zeitgeber_times=sample(zeitgeber_times, length(zeitgeber_times), replace=TRUE))))
null_parameters_dex <- do.call(rbind, lapply(1:n_e, function(x) generate_observed_parameters(gene_expression_dataframe=lfchange_dex, zeitgeber_times=sample(zeitgeber_times, length(zeitgeber_times), replace=TRUE))))

# Generate a binned frequency table of R2 values
vehicle_R2 <- null_parameters_vehicle$R2
dex_R2 <- null_parameters_dex$R2
vehicle_R2_freq_table <- table(cut(vehicle_R2, breaks=seq(-1, 1, by=0.1)))
dex_R2_freq_table <- table(cut(dex_R2, breaks=seq(-1,1, by=0.1)))

fisher_test <- fisher.test(vehicle_R2_freq_table, dex_R2_freq_table)
print(fisher_test)

# Same process for phase
vehicle_phase <- null_parameters_vehicle$phase
dex_phase <- null_parameters_dex$phase
vehicle_phase_freq_table <- table(cut(vehicle_phase, breaks=seq(min(c(vehicle_phase, dex_phase)), max(c(vehicle_phase, dex_phase)), by=0.1)))
dex_phase_freq_table <- table(cut(dex_phase, breaks=seq(min(c(vehicle_phase, dex_phase)), max(c(vehicle_phase, dex_phase)), by=0.1)))

fisher_test_phase <- fisher.test(vehicle_phase_freq_table, dex_phase_freq_table, simulate.p.value=TRUE)
print(fisher_test_phase)

# Same process for amplitude
vehicle_A <- null_parameters_vehicle$A
dex_A <- null_parameters_dex$A
vehicle_A_freq_table <- table(cut(vehicle_A, breaks=seq(min(c(vehicle_A, dex_A)), max(c(vehicle_A, dex_A)), by=0.001)))
dex_A_freq_table <- table(cut(dex_A, breaks=seq(min(c(vehicle_A, dex_A)), max(c(vehicle_A, dex_A)), by=0.001)))

fisher_test_A <- fisher.test(vehicle_A_freq_table, dex_A_freq_table, simulate.p.value=TRUE)
print(fisher_test_A)

### Generate parameters for amygdala data: p-value for phase, amplitude, and R2, as well as comparison p-value

# ska2_control_expression_no_covariates
# ska2_bd_expression_no_covariates
apar_C <- generate_observed_parameters(gene_expression_dataframe=ska2_control_expression_no_covariates, zeitgeber_times=cohort_data_C$ZT)
apar_BD <- generate_observed_parameters(gene_expression_dataframe=ska2_bd_expression_no_covariates, zeitgeber_times=cohort_data_BD$ZT)
print(paste("phase for control: ", apar_C$phase, ", phase for BD: ", apar_BD$phase))
print(paste("amplitude for control: ", apar_C$A, ", amplitude for BD: ", apar_BD$A))
print(paste("R2 for control: ", apar_C$R2, ", R2 for BD: ", apar_BD$R2))
print(paste("p-value for control R2: ", get_pvalue(unlist(ska2_control_expression_no_covariates), cohort_data_C$ZT), ", p-value for BD R2: ", get_pvalue(unlist(ska2_bd_expression_no_covariates), cohort_data_BD$ZT)))

# Generate comparison p-values for BD vs control
# Generate binned frequency table of R2 values
# For example:
# vehicle_R2 <- null_parameters_vehicle$R2
# dex_R2 <- null_parameters_dex$R2
# vehicle_R2_freq_table <- table(cut(vehicle_R2, breaks=seq(-1, 1, by=0.1)))
# dex_R2_freq_table <- table(cut(dex_R2, breaks=seq(-1,1, by=0.1)))

# fisher_test <- fisher.test(vehicle_R2_freq_table, dex_R2_freq_table)
# print(fisher_test)

# So for control vs. BD in phase:
n_e <- 1000
null_parameters_control <- do.call(rbind, lapply(1:n_e, function(x) generate_observed_parameters(gene_expression_dataframe=ska2_control_expression_no_covariates, zeitgeber_times=sample(cohort_data_C$ZT, length(cohort_data_C$ZT), replace=TRUE))))
null_parameters_bd <- do.call(rbind, lapply(1:n_e, function(x) generate_observed_parameters(gene_expression_dataframe=ska2_bd_expression_no_covariates, zeitgeber_times=sample(cohort_data_BD$ZT, length(cohort_data_BD$ZT), replace=TRUE))))

control_phase <- null_parameters_control$phase
bd_phase <- null_parameters_bd$phase
control_phase_freq_table <- table(cut(control_phase, breaks=seq(min(c(control_phase, bd_phase)), max(c(control_phase, bd_phase)), by=0.1)))
bd_phase_freq_table <- table(cut(bd_phase, breaks=seq(min(c(control_phase, bd_phase)), max(c(control_phase, bd_phase)), by=0.1)))

fisher_test_phase <- fisher.test(control_phase_freq_table, bd_phase_freq_table, simulate.p.value=TRUE)
print("Phase comparison p-value")
print(fisher_test_phase)

# For amplitude
control_A <- null_parameters_control$A
bd_A <- null_parameters_bd$A

control_A_freq_table <- table(cut(control_A, breaks=seq(min(c(control_A, bd_A)), max(c(control_A, bd_A)), by=0.001)))
bd_A_freq_table <- table(cut(bd_A, breaks=seq(min(c(control_A, bd_A)), max(c(control_A, bd_A)), by=0.001)))

fisher_test_A <- fisher.test(control_A_freq_table, bd_A_freq_table, simulate.p.value=TRUE)
print("Amplitude comparison p-value")
print(fisher_test_A)

# For R2
control_R2 <- null_parameters_control$R2
bd_R2 <- null_parameters_bd$R2

control_R2_freq_table <- table(cut(control_R2, breaks=seq(-1, 1, by=0.1)))
bd_R2_freq_table <- table(cut(bd_R2, breaks=seq(-1, 1, by=0.1)))

fisher_test_R2 <- fisher.test(control_R2_freq_table, bd_R2_freq_table, simulate.p.value=TRUE)
print("R2 comparison p-value")
print(fisher_test_R2)


