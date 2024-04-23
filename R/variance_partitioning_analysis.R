# variance_partitioning_analysis.R

# All variables:
# age
# Sex
# DDx
# psychosis
# suicide
# PMI
# pH
# Alcohol
# substance.abuse.severity
# smoking.history
# Race
# ZT
# onset
# Hemisphere
# Brain.weight
# duration
# lifetime.antipschotics..fluphenazine.eq

form <- ~ age + (1 | Sex) + (1 | DDx) + (1 | psychosis) + (1 | suicide) + PMI + pH + lifetime.antipsychotics..fluphenazine.eq. +
  (1 | ethanol.severity) + (1 | substance.abuse.severity) + (1 | smoking.history) + (1 | Race) + ZT + onset + (1 | hemisphere) +
  brain.weight

cordata <- canCorPairs(form, cohort_data)
plotCorrMatrix(cordata)
