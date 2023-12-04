
# download library allowing to read matlab files
#install.packages("R.matlab")
library(R.matlab)

#install.packages("lavaan")
library(lavaan)
#install.packages('semPlot')
library(semPlot)

# read the data prepared for SEM
setwd("P:/boulot/postdoc_CarmenSandi/results/SEM")
data_for_SEM <- readMat("bloodLac_dmPFCLac_MRS_dmPFC_GLM200_EpEm_Ech_63subs.mat")
mtrx_SEM = data_for_SEM$mtrx.THE

# rename variables
blood_Lac = data_for_SEM$mtrx.THE[,1]
dmPFC_Lac = data_for_SEM$mtrx.THE[,2]
dmPFC_fMRI_GLM200 = data_for_SEM$mtrx.THE[,3]
THE = data_for_SEM$mtrx.THE[,4]
PHE = data_for_SEM$mtrx.THE[,5]
MHE = data_for_SEM$mtrx.THE[,6]
kEp = data_for_SEM$mtrx.THE[,7]
kEm = data_for_SEM$mtrx.THE[,8]

## PHE
#create data frame
df = data.frame(blood_Lac, dmPFC_Lac, dmPFC_fMRI_GLM200, PHE)

# define the SEM model to try (check https://bookdown.org/jdholster1/idsr/structural-equation-modeling.html for more details as well)
# Specify model, can also write it in a single line
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI_GLM200 ~ dmPFC_Lac + blood_Lac
PHE ~ dmPFC_fMRI_GLM200 + dmPFC_Lac + blood_Lac
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# make a nice display
semPaths(fit, "std", edge.label.cex = 1.0, curvePivot = TRUE)


## kEp
#create data frame
df = data.frame(blood_Lac, dmPFC_Lac, dmPFC_fMRI_GLM200, kEp)

# define the SEM model to try (check https://bookdown.org/jdholster1/idsr/structural-equation-modeling.html for more details as well)
# Specify model, can also write it in a single line
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI_GLM200 ~ dmPFC_Lac + blood_Lac
kEp ~ dmPFC_fMRI_GLM200 + dmPFC_Lac + blood_Lac
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# make a nice display
semPaths(fit, "std", edge.label.cex = 1.0, curvePivot = TRUE)


## MHE
#create data frame
df = data.frame(blood_Lac, dmPFC_Lac, dmPFC_fMRI_GLM200, MHE)

# define the SEM model to try (check https://bookdown.org/jdholster1/idsr/structural-equation-modeling.html for more details as well)
# Specify model, can also write it in a single line
sem_model <- '
dmPFC_Lac ~ blood_Lac
dmPFC_fMRI_GLM200 ~ dmPFC_Lac + blood_Lac
MHE ~ dmPFC_fMRI_GLM200 + dmPFC_Lac + blood_Lac
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# make a nice display
semPaths(fit, "std", edge.label.cex = 1.0, curvePivot = TRUE)