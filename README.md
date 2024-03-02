## Circadian Analysis [REPRODUCTION]: Rhythmicity detection and DE analysis 
The purpose of this repository is to reproduce results from Dr. Cahill and colleagues.

# Required R packages: 
- minpack.lm
- Snowfall
- doParallel
- Lme4
- qvalue (bioconductor) <br/>
*This code has been tested on R.4.3.2

# Installation Instructions: 
1. Clone the repository
```
git clone https://github.com/rela-v/Circadian-Analysis-Re
```
2. Open RStudio and run "Circadian-Analysis-Re/R/main.R"
3. See Console for messages including errors.
4. See "Circadian-Analysis-Re/Results" for compiled results.
# Clinical Data:
A subset of expression, and clinical data for 104 control subjects is available in csv form in the rhythmicity folder. 
Clinical variables of interest include:
- corrected time of death (TOD) in ZT. 
- Hospital site (Mount Sinai vs Pitt)
- Diagnosis (control vs schizophrenia)
- Gender 
- Ethnicity

# Expression data:  
- 13914 genes (all genes that were used in the full analysis in the paper)
- Data has already been filtered according to method discussed in the paper supplement 
- Expression units are in log2 CPM 

# Demo:
- fitSinCurve.R and Curve_Drawing.R are external functions required to run the main analysis in RhythmicityCode.R 
- BestModelSelection.R is an external function required to run the main DE.R analysis
- Run all external functions, so that they are in the R environment 
- Run code as is in RhythmicityCode.R 
- Run code as is in DE.R <br/>
*NOTE: Each method requires permutation. Permutations were set to 10 to save time and computation space. Permuted files are included in the DE (DE NullData) and Rhythmicity (nullFolder) folders to save reviewers time. If permutation files are recreated by the reviewer, output will be slightly different due to the randomness of permutation. 

# Expected output:
- Results/Rhythmicity/Example_result.csv is the example output of curve fitting parameters and p value. 
- Results/Rhythmicity/PDF plots of top circadian genes 
- Results/Rhythmicity/Example_result2.csv is the example output for the shift in rhythmicity analysis <br/>
*NOTE: Because only control subjects are used in the example data, two groups were created using age to demonstrate gain and loss of rhythmicity 

- Results/DE/Example_result3.csv example of DE output including original p-value from the likelihood ratio test, corrected p-value, Benjamini Hochberg corrected p-value, and Storey's q value. <br/>
*NOTE: Because only control subjects are used in the example data, effect of a binary indicator for time of day (morning vs night) is analyzed. 

Estimated run time: 
30min-1hr

# Output:

- Control, ENSG00000104907.12: ; peak = −2 -> TRMT1 tRNA methyltransferase 1
- Control, ENSG00000184271.17: ; peak = −5 -> POU Domain, Class 6, Transcription Factor 1
- Control, ENSG00000231680.1: ; peak = 6 -> LINC02723 Gene - Long Intergenic Non-Protein Coding RNA
- Control, ENSG00000232006.9: ; peak = −4 -> [PROMOTER] HECW1 Gene - HECT, C2 And WW Domain Containing E3 Ubiquitin Protein Ligase 1
