# Load the data and make sure the rownames and colnames are matching
expression_data <- read.csv("./Data/Project_1005-rawCounts.csv", sep='\t', skip=2, header=FALSE,row.names=1)
cohort_data <- read.csv("./Data/cohort_data.csv", row.names = 1)

# Make sure the rownames and colnames are matching
colnames(expression_data) <- c(rownames(cohort_data))
