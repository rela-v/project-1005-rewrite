##DE Analysis 
expression_data<-read.csv("./Data/revised_expression_data.csv", row.names = 1)
cohort_data<-read.csv('./Data/revised_cohort_data.csv', row.names = 1)
n_c<-nrow(cohort_data)
n_e<-nrow(expression_data)
Symbols<-row.names(expression_data)

ZT<-cohort_data$ZT
night.index.neg<-which(ZT < 0)
night.index.pos<-which(ZT<=18 & ZT>12)
night.index<-c(night.index.neg,night.index.pos) #32
morn.index<-which(ZT>=0 & ZT<=12) #60
tod_index<-rep(1,n_c)
tod_index<-replace(tod_index, morn.index, 0)
cohort_data$tod_index<-tod_index

covariates0<-c("Sex","PMI", "age", "psychosis", "suicide")
covariables3 <- cohort_data[,covariates0]
UniqueID <- cohort_data$hemisphere
VariableListOne = combn(covariates0,1)
VariableListTwo = combn(covariates0,2)
Diagnosis<-as.factor(tod_index)

set.seed(15213)
L3startTime = Sys.time()
sfInit(parallel=TRUE,type="SOCK", cpus=24) 
GeneIndeces = 1:n_e
sfExport("expression_data") 
sfExport("covariables3") 
sfExport("Diagnosis") 
sfExport("UniqueID") 
sfExport("VariableListOne") 
sfExport("VariableListTwo") 
sfExport("GeneIndeces") 

indRIM<-c()
sfExport("bestModelSelection") ##Very long, just running first 100 genes 
indRIM<-sfLapply(1:10,function(x) try(bestModelSelection(x,expression_data,covariables3,Diagnosis,UniqueID,VariableListOne,VariableListTwo),silent=TRUE ) ) 
sfStop()

#Permutation to correct for p-value bias 
rawdata8<-expression_data[1:10,]
B<-10
#Corrected 
system("mkdir -p ./Results/DE/NullData")
foreach(b=1:B) %dopar% {
  set.seed(b)
  print(b)
  rawdata8_b <- rawdata8[,sample(ncol(rawdata8))]
  result_PVL3_SZ <- replicate(nrow(rawdata8_b),list())
  names(result_PVL3_SZ) <- rownames(rawdata8_b)
  sfExport("rawdata8_b")
  result_PVL3_SZ<-sfLapply(1:nrow(rawdata8_b),function(x) try(bestModelSelection(x,rawdata8_b,covariables3,Diagnosis,UniqueID,VariableListOne,VariableListTwo),silent=TRUE ) ) 
  afile <- paste('./Results/DE/NullData/DE_null_',b,'.rdata',sep='')
  save(result_PVL3_SZ,file=afile)
}

RIM.NULL<-c()
foreach(b=1:10) %dopar% {
  set.seed(b)
  print(b)
  aNull<-get(load(paste("./Results/DE/NullData/DE_null_", b,".rdata", sep = "")))
  RIM.NULL[[b]]<-aNull
}

unlistRIM<-unlist(RIM.NULL)
names.pval<-grep("lrt.pvalue", names(unlistRIM))
pval.null<-unname(as.numeric(unlistRIM[names.pval]))
unlistRIM2<-unlist(indRIM)
names.pval<-grep("lrt.pvalue", names(unlistRIM2))
pval.biased<-unname(as.numeric(unlistRIM2[names.pval]))

p.corr<-c()
numerator<-c()
denominator<-length(pval.null)
foreach(i=1:length(pval.biased)) %dopar% {
  set.seed(i)
  numerator[i]<-length(which(pval.null<pval.biased[i])) #empirical p valuse calculation
  p.corr[i]<-numerator[i]/denominator
  print(i) 
}

# bh.q<-p.adjust(p.corr, "BH") 
# library(qvalue)
# q<-qvalue(p.corr)$qvalues
result3<-data.frame(cbind(pval.biased)) #, p.corr, bh.q, q))
row.names(result3)<-Symbols[1:10]
write.csv(result3, "./Results/DE/Example_result3.csv")
