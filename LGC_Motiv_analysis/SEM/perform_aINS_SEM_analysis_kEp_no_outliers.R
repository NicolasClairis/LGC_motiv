## notes regarding fit of the data: 
# CFI and TLI should be close to 1
# RMSEA should be <= 0.08 (reasonable approximate fit) or better <= 0.05 (close-fit) while if >= 0.1
# check https://stats.oarc.ucla.edu/r/seminars/rsem/ for more details also


# download library allowing to read matlab files
#install.packages("R.matlab")
library(R.matlab)

#install.packages("lavaan")
library(lavaan)
#install.packages('semPlot')
library(semPlot)

# read the data prepared for SEM
setwd("P:/boulot/postdoc_CarmenSandi/results/SEM")
data_for_SEM <- readMat("bloodLac_aINSLac_MRS_aINS_GLM235_EpEm_Ech_40subs_kEp_no_outliers.mat")

# rename variables
blood_Lac = data_for_SEM$mtrx.kEp.no.outliers[,1]
aINS_Lac = data_for_SEM$mtrx.kEp.no.outliers[,2]
aINS_fMRI = data_for_SEM$mtrx.kEp.no.outliers[,3]
kEp = data_for_SEM$mtrx.kEp.no.outliers[,4]
#create data frame
df = data.frame(blood_Lac, aINS_Lac, aINS_fMRI, kEp)

##############################################################################
## test 1: kEp full path

# define the SEM model to try (check https://bookdown.org/jdholster1/idsr/structural-equation-modeling.html for more details as well)
# Specify model, can also write it in a single line
sem_model <- '
aINS_Lac ~ blood_Lac
aINS_fMRI ~ aINS_Lac + blood_Lac
kEp ~ aINS_fMRI + aINS_Lac + blood_Lac
'

# Estimate our model
fit_model1 <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model1, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
# details <- summary(fit_model1, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
# make a nice display
semPaths(fit_model1, "std", edge.label.cex = 1.0, curvePivot = TRUE)


##############################################################################
## test 2: blood=>aINS=>fMRI=>bhv + blood => fMRI

## kEp
sem_model <- '
aINS_Lac ~ blood_Lac
aINS_fMRI ~ aINS_Lac + blood_Lac
kEp ~ aINS_fMRI
'
fit_model2 <- sem(sem_model,data=df)
summary(fit_model2, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model2, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model2, "std", edge.label.cex = 1.0, curvePivot = TRUE)


##############################################################################
## test 3: blood=>aINS=>fMRI=>bhv

## kEp
sem_model <- '
aINS_Lac ~ blood_Lac
aINS_fMRI ~ aINS_Lac
kEp ~ aINS_fMRI + aINS_Lac
'
fit_model3 <- sem(sem_model,data=df)
summary(fit_model3, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model3, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model3, "std", edge.label.cex = 1.0, curvePivot = TRUE)


##############################################################################
## test 4: blood=>aINS=>fMRI=>bhv (no direct path)

## kEp
sem_model <- '
aINS_Lac ~ blood_Lac
aINS_fMRI ~ aINS_Lac
kEp ~ aINS_fMRI
'
fit_model4 <- sem(sem_model,data=df)
summary(fit_model4, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model4, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model4, "std", edge.label.cex = 1.0, curvePivot = TRUE)

##############################################################################
## test 5: aINS=>fMRI=>bhv (basic mediation without blood)

## THE
sem_model <- '
aINS_fMRI ~ aINS_Lac
kEp ~ aINS_fMRI + aINS_Lac
'
fit_model5 <- sem(sem_model,data=df)
summary(fit_model5, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model5, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model5, "std", edge.label.cex = 1.0, curvePivot = TRUE)