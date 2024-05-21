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
# old model
#data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM200_EpEm_Ech_58subs_PHE_no_outliers.mat")
#data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM232_EpEm_Ech_59subs_PHE_no_outliers.mat")
#data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM235_EpEm_Ech_59subs_PHE_no_outliers.mat")
#data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM242_EpEm_Ech_59subs_PHE_no_outliers.mat")
#data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM244_EpEm_Ech_58subs_PHE_no_outliers.mat")
# current model
data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM263_EpEm_Ech_59subs_PHE_no_outliers.mat")

# rename variables
blood_Lac = data_for_SEM$mtrx.PHE.no.outliers[,1]
dmPFC_Lac = data_for_SEM$mtrx.PHE.no.outliers[,2]
dmPFC_fMRI = data_for_SEM$mtrx.PHE.no.outliers[,3]
PHE = data_for_SEM$mtrx.PHE.no.outliers[,4]
#create data frame
df = data.frame(blood_Lac, dmPFC_Lac, dmPFC_fMRI, PHE)

##############################################################################
## test 1: PHE full path

# define the SEM model to try (check https://bookdown.org/jdholster1/idsr/structural-equation-modeling.html for more details as well)
# Specify model, can also write it in a single line
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI ~ dmPFC_Lac + blood_Lac
PHE ~ dmPFC_fMRI + dmPFC_Lac + blood_Lac
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
## test 2: blood=>dmPFC=>fMRI=>bhv + blood => fMRI

## PHE
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI ~ dmPFC_Lac + blood_Lac
PHE ~ dmPFC_fMRI
'
fit_model2 <- sem(sem_model,data=df)
summary(fit_model2, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model2, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model2, "std", edge.label.cex = 1.0, curvePivot = TRUE)


##############################################################################
## test 3: blood=>dmPFC=>fMRI=>bhv

## PHE
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI ~ dmPFC_Lac
PHE ~ dmPFC_fMRI + dmPFC_Lac
'
fit_model3 <- sem(sem_model,data=df)
summary(fit_model3, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model3, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model3, "std", edge.label.cex = 1.0, curvePivot = TRUE)


##############################################################################
## test 4: blood=>dmPFC=>fMRI=>bhv (no direct path)

## PHE
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI ~ dmPFC_Lac
PHE ~ dmPFC_fMRI
'
fit_model4 <- sem(sem_model,data=df)
summary(fit_model4, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
#details <- summary(fit_model4, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
semPaths(fit_model4, "std", edge.label.cex = 1.0, curvePivot = TRUE)