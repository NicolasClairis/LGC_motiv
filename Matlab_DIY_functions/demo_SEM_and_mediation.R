# Demo on SEM to compare with Matlab home-made mediation script

# workind directory
setwd("P:/boulot/postdoc_CarmenSandi/results/SEM")

## load the data
demo1 <- readMat("demo_SEM_1.mat")
mtrx_SEM1 = demo1$mtrx.THE
demo2 <- readMat("demo_SEM_2.mat")
mtrx_SEM2 = demo2$mtrx.THE
demo3 <- readMat("demo_SEM_3.mat")
mtrx_SEM3 = demo3$mtrx.THE

################################################################################
## DEMO 0: compare mediation and SEM
# load variables
X1 = demo1$mtrx.THE[,2]
M = demo1$mtrx.THE[,3]
Y = demo1$mtrx.THE[,4]
#create data frame
df = data.frame(X1, M, Y)

# Specify model
sem_model <- '
M ~ X1
Y ~ M + X1
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
# details <- summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
# make a nice display
semPaths(fit_model, "std", edge.label.cex = 1.0, curvePivot = TRUE)


################################################################################
## DEMO 1: only indirect path
# load variables
X0 = demo1$mtrx.THE[,1]
X1 = demo1$mtrx.THE[,2]
M = demo1$mtrx.THE[,3]
Y = demo1$mtrx.THE[,4]

#create data frame
df = data.frame(X0, X1, M, Y)

# Specify model
sem_model <- '
X1 ~ X0
M ~ X1
Y ~ M
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
# details <- summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
# make a nice display
semPaths(fit_model, "std", edge.label.cex = 1.0, curvePivot = TRUE)

############################################################################################
## DEMO 2: indirect + direct path (X1=>M)
# load variables
X0 = demo2$mtrx.THE[,1]
X1 = demo2$mtrx.THE[,2]
M = demo2$mtrx.THE[,3]
Y = demo2$mtrx.THE[,4]

#create data frame
df = data.frame(X0, X1, M, Y)

# Specify model
sem_model <- '
X1 ~ X0
M ~ X1
Y ~ M + X1
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
# details <- summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
# make a nice display
semPaths(fit_model, "std", edge.label.cex = 1.0, curvePivot = TRUE)

#######################
## DEMO 3: indirect + direct path (X1=>M and X0=>M)
# load variables
X0 = demo3$mtrx.THE[,1]
X1 = demo3$mtrx.THE[,2]
M = demo3$mtrx.THE[,3]
Y = demo3$mtrx.THE[,4]

#create data frame
df = data.frame(X0, X1, M, Y)

# Specify model
sem_model <- '
X1 ~ X0
M ~ X1
Y ~ M + X0 + X1
'

# Estimate our model
fit_model <- sem(sem_model,data=df)
# Summarize the results since it is degree of fred = 0, it is just identified, so fit.measures won't give more
summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
# in case you need to extract p.value very low use the following lines:
# details <- summary(fit_model, fit.measures=TRUE, standardized = TRUE, rsquare=TRUE)
#details$pe
# make a nice display
semPaths(fit_model, "std", edge.label.cex = 1.0, curvePivot = TRUE)