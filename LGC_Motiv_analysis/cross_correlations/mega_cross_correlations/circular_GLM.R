# adaptation of circular_plot.R and circular_rplot_bis.R but instead of looking at correlations like these scripts,
# circular_GLM.R will look at the betas after correcting for age and sex and, instead of filtering on r, 
#circular_GLM.R will filter on p.values.

# 0. install required libraries
# install.packages("circlize")
# install.packages("tidyverse") # install excel file


# 1. Loading libraries
library(circlize) # necessary to create the plots with circlize
library(readxl) # necessary to read excel files
library(stringr) # ?

# 2. Reading input file
# setwd('C:/Users/rdrotari/OneDrive - NESTLE/Projects/CircularPlots/ClairisEPFL/'); # Ricardo path
setwd("P:/boulot/postdoc_CarmenSandi/results/mega_correlation_matrix") # Nicolas path
betasList = read_xlsx('GLM_table_signed_betas_75subs.xlsx')

# 3. Selecting relevant variables and organizing them in groups
varList = unique(betasList$var1)

# shortening the variable names by removing the general category
# general
galList0 = varList[startsWith(varList,'gal') & varList != 'gal_age' & varList != 'gal_sex'] # selection of variables associated to general + ignore age and sex which are just NaN and removed from the analysis (no meaning since they are already in the explaining variables of the GLM)
galgroup0 = which(startsWith(galList0,'gal_')) # definition of general group
galList1 <- gsub("gal_","",galList0) # shortening the variable names
galList1 <- gsub("hexaco_","",galList1) # shortening the variable names
galList1 <- gsub("prevDay_min_avg_sleep","delta sleep",galList1) # shortening the variable names
galList2 = galList1
for (var in galList0){ # replacing short names in betasList [dataframe with value of correlation coefficient]
  ivar = which(galList0 == var)
  w1 = which(betasList$var1 == var)
  betasList$var1[w1] = galList2[ivar]
  w2 = which(betasList$var2 == var)
  betasList$var2[w2] = galList2[ivar]
}

# behavior: questionnaires + behavioral task
behaviorList0 = varList[startsWith(varList,'behavior_questionnaire') | startsWith(varList,'behavior_task')] # selection of variables associated to behavior
# questionnaires
behaviorList1 = sapply(behaviorList0,function(x) gsub('behavior_questionnaires_','',x))
# behavioral groups
behaviorgroup1 = which(startsWith(behaviorList1,'stress_anxiety_')) # definition of group 1 in behavior: stress/anxiety
behaviorgroup2 = which(startsWith(behaviorList1,'dominance_')) # definition of group 2 in behavior: dominance
behaviorgroup3 = which(startsWith(behaviorList1,'motivation_')) # definition of group 3 in behavior: motivation
behaviorgroup4 = which(startsWith(behaviorList1,'behavior_task_')) # definition of group 4: behavioral task
# task variables
behaviorList1 = sapply(behaviorList1,function(x) gsub('behavior_task_choices_','',x))
behaviorList1 = sapply(behaviorList1,function(x) gsub('behavior_task_prm_','',x))
# shortening questionnaire variable names
behaviorList2 = behaviorList1 # shortening the variable names
# stress/anxiety questionnaires
behaviorList2[behaviorgroup1] = sapply(behaviorList2[behaviorgroup1],function(x) gsub('stress_anxiety_','',x))
behaviorList2 <- gsub("CTQ_sexA","CTQsa",behaviorList2)
behaviorList2 <- gsub("CTQ_physicalA","CTQpa",behaviorList2)
behaviorList2 <- gsub("CTQ_physicalN","CTQpn",behaviorList2)
behaviorList2 <- gsub("CTQ_minDenial","CTQmd",behaviorList2)
behaviorList2 <- gsub("CTQ_emotionalA","CTQea",behaviorList2)
behaviorList2 <- gsub("CTQ_emotionalN","CTQen",behaviorList2)
# dominance questionnaires
behaviorList2[behaviorgroup2] = sapply(behaviorList2[behaviorgroup2],function(x) gsub('dominance_','',x))
behaviorList2 <- gsub("CI_contentiousness","CIc",behaviorList2)
behaviorList2 <- gsub("CI_enjCompet","CIec",behaviorList2)
# motivation questionnaires
behaviorList2[behaviorgroup3] = sapply(behaviorList2[behaviorgroup3],function(x) gsub('motivation_','',x))
behaviorList2 <- gsub("Lars_e_ActionInit","Lars_e_AI",behaviorList2)
behaviorList2 <- gsub("Lars_e_AI_Init","Lars_e_AIi",behaviorList2)
behaviorList2 <- gsub("Lars_e_AI_EverydayProd","Lars_e_AIep",behaviorList2)
behaviorList2 <- gsub("Lars_e_IntellectCuriosity","Lars_e_IC",behaviorList2)
behaviorList2 <- gsub("Lars_e_IC_Novelty","Lars_e_ICn",behaviorList2)
behaviorList2 <- gsub("Lars_e_IC_Social","Lars_e_ICs",behaviorList2)
behaviorList2 <- gsub("Lars_e_IC_Interest","Lars_e_ICi",behaviorList2)
behaviorList2 <- gsub("Lars_e_IC_Motiv","Lars_e_ICm",behaviorList2)
behaviorList2 <- gsub("Lars_e_SelfAwareness","Lars_e_SA",behaviorList2)
behaviorList2 <- gsub("Lars_e_EmotResp","Lars_e_ER",behaviorList2)
behaviorList2 <- gsub("MPSTEFS_mental","MPSTEFSm",behaviorList2)
behaviorList2 <- gsub("MPSTEFS_physical","MPSTEFSp",behaviorList2)
for (var in behaviorList0){ # replacing short names in betasList [dataframe with value of correlation coefficient]
  ivar = which(behaviorList0 == var)
  w1 = which(betasList$var1 == var)
  betasList$var1[w1] = behaviorList2[ivar]
  w2 = which(betasList$var2 == var)
  betasList$var2[w2] = behaviorList2[ivar]
}

# brain metabolites
brainList0 = varList[startsWith(varList,'brain')] # selection of variables associated to brain
brainList1 = sapply(brainList0,function(x) gsub('brainM_','',x)) # shortening the variable names
brainList1 <- gsub("Gln_div_Glu","Gln/Glu",brainList1)
braingroup1 = which(startsWith(brainList1,'dmPFC_')) # definition of group 1 in brain: dmPFC
braingroup2 = which(startsWith(brainList1,'aIns_')) # definition of group 2 in brain: aIns
brainList1 <- gsub("dmPFC_","d_",brainList1) # shortening the variable names removing dmPFC now that we know
brainList1 <- gsub("aIns_","ai_",brainList1) # shortening the variable names removing aIns now that we know
brainList2 = brainList1
for (var in brainList0){ # replacing short names in betasList [dataframe with value of correlation coefficient]
  ivar = which(brainList0 == var)
  w1 = which(betasList$var1 == var)
  betasList$var1[w1] = brainList2[ivar]
  w2 = which(betasList$var2 == var)
  betasList$var2[w2] = brainList2[ivar]
}
# plasma, saliva and whole-blood metabolites
# plasma
plasmaList0 = varList[startsWith(varList,'plasma')] # selection of variables associated to plasma
plasmaList1 = sapply(plasmaList0,function(x) gsub('plasma_','',x)) # shortening the variable names
# whole-blood
bloodList0 = varList[startsWith(varList,'wholeB')] # selection of variables associated to whole blood
bloodList0 = bloodList0[str_count(bloodList0,'_') == 1] # removing not interesting variables (based on Nestlé)
bloodList0 = bloodList0[bloodList0 != 'wholeB_tNAD']
bloodList1 = sapply(bloodList0,function(x) gsub('wholeB_','',x)) # shortening the variable names
bloodList1 <- gsub("NAD_div_NADH","NAD/NADH",bloodList1)
bloodList1 <- gsub("NADP_div_NADPH","NADP/NADPH",bloodList1)
# salivary metabolites
salivaList0 = varList[startsWith(varList,'salivary')] # selection of variables associated to saliva
salivaList1 = sapply(salivaList0,function(x) gsub('salivary_','',x)) # shortening the variable names
salivaList1 <- gsub('TESTO_div_CORT','TESTO/CORT',salivaList1)
# store info in global circulation list variables
circulationList0 = c(plasmaList0,bloodList0,salivaList0) # creating the list of variable associated to circulation (= plasma + whole blood + saliva)
circulationList1 = c(plasmaList1,bloodList1,salivaList1) # creating the list of variable associated to circulation (= plasma + whole blood + saliva)
# identify groups for colors
circulationgroup1 = which(startsWith(circulationList1,'aa_')) # definition of group 1 in circulation: amino-acids
circulationgroup2 = which(startsWith(circulationList1,'Lac')) # definition of group 2 in circulation: lactic acid
circulationgroup3 = which(startsWith(circulationList1,'fa_')) # definition of group 3 in circulation: fatty acids
circulationgroup4 = which(startsWith(circulationList0,'wholeB_')) # definition of group 4 in circulation: NAD derivatives
circulationgroup5 = which(startsWith(circulationList0,'salivary_')) # salivary measures (testosterone, cortisol, interleukins)
# shorten variable names
circulationList2 = circulationList1 # shortening the variable names
circulationList2[circulationgroup1] = sapply(circulationList1[circulationgroup1],function(x) gsub('aa_','',x))
circulationList2[circulationgroup3] = sapply(circulationList1[circulationgroup3],function(x) gsub('fa_','',x))
for (var in circulationList0){ # replacing short names in betasList [dataframe with value of correlation coefficient]
  ivar = which(circulationList0 == var)
  w1 = which(betasList$var1 == var)
  betasList$var1[w1] = circulationList2[ivar]
  w2 = which(betasList$var2 == var)
  betasList$var2[w2] = circulationList2[ivar]
}
# define list of variables to include
goodList = c(galList2,behaviorList2,brainList2,circulationList2) # creation of the "good" list of variables (= behavior + brain + circulation + saliva)

# 3. Creation of the beta matrix and the corresponding p.Value matrix
betasMat = matrix(nrow = length(goodList), ncol = length(goodList)) # initialization
pvalMat = matrix(nrow = length(goodList), ncol = length(goodList)) # initialization
rownames(betasMat) = goodList #definition of the row names (= variable names)
colnames(betasMat) = goodList #definition of the columns names (= variable names)
rownames(pvalMat) = goodList #definition of the row names (= variable names)
colnames(pvalMat) = goodList #definition of the columns names (= variable names)
for (ir in goodList){ # retrieving of correlation coefficients and put in the matrix
  print(ir)
  for (ic in goodList){
    betasMat[ir,ic] = betasList$b_var[betasList$var1 == ir & betasList$var2 == ic]
    pvalMat[ir,ic] = betasList$pvalue_var[betasList$var1 == ir & betasList$var2 == ic]
  }
}


# 4. Function to determine the color of the link according to the valence of the beta (distinguish between positive (in purple) and negative (in grey) correlations)
colorLink <- function(beta,bMin,bMax){   # Function to determine the color of the link according to the value of the correlation coefficient
  if (beta > 0){col0 = '#c51b8a'} else {col0 = '#000000'}
  if (abs(beta) > bMax) {
    tr = floor((abs(beta)-bMin)/(bMax-bMin)*256)
  }  else if (abs(beta) < bMin) {# use min color if abs(beta) < bMin
    tr = 0
  }  else {# use max color if abs(beta) > bMax
    tr = 255
  }
  if (tr < 16){
    col0 = paste(col0,'0',as.character(as.hexmode(tr)),sep='')
  } else {
    col0 = paste(col0,as.character(as.hexmode(tr)),sep='')
  }
  return(col0)
}

# 5. BASIC DEFINITION OF THE CIRCLE PLOT
# CHANGE HERE THE VALUES FOR THE THRESHOLD ON THE CORRELATION COEFFICIENT AND THE LINE WIDTH OF THE LINKS
wdt = 1.5 # line width of the link in the plot
pvalThresh = 0.01 # Threshold for the p.value of the correlation to filter

# 6. Creation of the circular plot
circos.clear()
circos.par$track.height = 0.05 # height of the circular sector
circos.par$circle.margin = 0.35 # margin of the figure
circos.par$start.degree = 45 # start point (measure in degrees): the full list starts from behavior variables
circos.par(cell.padding = c(0, 0))
circos.initialize(colnames(betasMat), xlim = c(-1,1),)

#color for general variables
bgcol0 = rep(NA,length(galList2)) # definition of the colors for circular sectors
# general variables (sex, age, sleep, etc.)
bgcol0[galgroup0] = '#ffffcc'

# color for behavior
bgcol1 = rep(NA,length(behaviorList2)) # definition of the colors for circular sectors
# green colors for behavior variables: different shades stands for different subgroups
bgcol1[behaviorgroup1] = '#238b45' # stress/anxiety questionnaires
bgcol1[behaviorgroup2] = '#66c2a4' # dominance questionnaires
bgcol1[behaviorgroup3] = '#b2e2e2' # motivation questionnaires
bgcol1[behaviorgroup4] = '#edf8fb' # motivational task 
# color for brain metabolites
bgcol2 = rep(NA,length(brainList2))
#blue colors for brain variables: different shades stands for dACC vs aIns
bgcol2[braingroup1] = '#2b8cbe' # dmPFC/dACC
bgcol2[braingroup2] = '#a6bddb' # anterior insula
# color for circulatory variables
bgcol3 = rep(NA,length(circulationList2))
#red colors for circulation variables: different shades stands for different subgroups
bgcol3[circulationgroup1] = '#a50f15' # plasma aa
bgcol3[circulationgroup2] = '#de2d26' # plasma fa
bgcol3[circulationgroup3] = '#fb6a4a' # plasma Lac
bgcol3[circulationgroup4] = '#fcae91' # whole-blood NADomics
bgcol3[circulationgroup5] = '#fee5d9' # saliva: TESTO/CORT/IL
bgcol = c(bgcol0,bgcol1,bgcol2,bgcol3)

# 7. Drawing circular sectors and writing variable names
circos.track(ylim = c(0, 1), track.height = mm_h(4),
             panel.fun = function(x, y) {}, bg.border = NA)
circos.track(ylim = c(-1,1),bg.col=bgcol,
             panel.fun = function(x, y) {circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2] + mm_y(2),CELL_META$sector.index,
                                                     facing = 'clockwise', cex = 1, adj = c(0,0.5))})

# 8. Creation of the matrix "links" with pairs of correlated variable (pval > rThres) + (at least 2 links between different categories)
nvar = nrow(betasMat) #number of variables
groupVar = rep(NA,nvar)
# define big groups
groupVar[which(goodList %in% galList2)] = 'general'
groupVar[which(goodList %in% behaviorList2)] = 'behavior'
groupVar[which(goodList %in% brainList2)] = 'brain'
groupVar[which(goodList %in% circulationList2)] = 'circulation'
# extract number of links that are above the pvalThresh threshold defined for betas
nlinks = (sum(pvalMat<=pvalThresh)-nvar)/2 #divide by 2 because all links appear twice in the matrix (but no idea what the nvar is doing)
links = matrix(NA,nrow = nlinks,ncol=6) # 6 columns: 1) var1 name, (2) var2 name, (3) beta, (4) pvalue, (5) category var1, (6) category var2
w = which(pvalMat<=pvalThresh)
il = 0
for (iw in w){
  i1 = floor((iw-1)/nvar)+1
  i2 = (iw-1)%%nvar+1
  if (i1>i2){
    il = il+1
    beta = betasList$b_var[betasList$var1 == rownames(betasMat)[i1] & betasList$var2 == colnames(betasMat)[i2]]
    pval = betasList$pvalue_var[betasList$var1 == rownames(betasMat)[i1] & betasList$var2 == colnames(betasMat)[i2]]
    w1 = which(goodList == rownames(betasMat)[i1])
    w2 = which(goodList == colnames(betasMat)[i2])
    links[il,] = c(rownames(betasMat)[i1],colnames(betasMat)[i2],beta,pval,groupVar[w1],groupVar[w2])
  }
}

# filter links to only keep those that bridge one item with two different other groups
#create variable grower to identify the number of links between categories for each item which has a significant correlation with another
allLinkVars = unique(c(links[,1],links[,2]))
grower <- matrix(nrow = length(allLinkVars), ncol = 3)
grower <- as.data.frame(grower)
grower$V1 <- allLinkVars

for (i in 1:nrow(grower)) { # loop through variables that display a significant correlation with another
  # extract category of the current variable
  groupCurrentVar1 = unique(links[links[,1] == grower[i,1],5])
  groupCurrentVar2 = unique(links[links[,2] == grower[i,1],6])
  groupCurrentVar <- unique(c(groupCurrentVar1, groupCurrentVar2))
  # extract category of all the variables correlated to it
  groupCorrelatedVar1 = unique(links[links[,1] == grower[i,1],6])
  groupCorrelatedVar2 = unique(links[links[,2] == grower[i,1],5])
  groupCorrelatedVars <- unique(c(groupCorrelatedVar1, groupCorrelatedVar2))
  grower[i,2] <- length(groupCorrelatedVars) # extract number of categories correlated to current item (potentially including same category as well)
  # remove the count for the category belonging to the same dimension
  if (groupCurrentVar %in% groupCorrelatedVars){
    grower[i,3] <- length(groupCorrelatedVars) - 1
  }
  else {
    grower[i,3] <- length(groupCorrelatedVars)
  }
}
# filter variables that correlate with at least 2 different categories
grower <- grower[grower$V3 > 1, ]

if (nrow(grower) > 0) { # filter case where there are no links to plot
  links2 = links[links[,1] %in% grower[,1], ] # remove links from any item that links only within the same category or doesn't have at least 2 different links to different categories
  
  # 9. Drawing links2 in the circular plot
  # define range for color range to be used for links
  bMinScale = 0.2
  bMaxScale = 2
  # draw the links between different categories
  for (iLink in 1:nrow(links2)){ #loop through correlations
    w1 = which(goodList == links2[iLink,1])
    w2 = which(goodList == links2[iLink,2])
    if (groupVar[w1] != groupVar[w2]){  #plot only link between variables of different groups
      circos.link(sector.index1 = links2[iLink,1], c(-wdt,wdt), sector.index2 = links2[iLink,2], c(-wdt,wdt),col = colorLink(as.numeric(links2[iLink,3]),bMinScale,bMaxScale))
    }
  }
}

# 10. plotting legend
legend(1.2,1,legend=c('General'),fill = c('#ffffcc'),cex=1)
legend(1.2,-0.1,legend = c('Q stress/anxiety','Q dominance','Q motivation','B motivation'),title='Behavior', fill = c('#238b45','#66c2a4','#b2e2e2','#edf8fb'), cex = 1)
legend(-2,-0.8,legend = c('dmPFC/dACC','aIns'),title='Brain', fill = c('#2b8cbe','#a6bddb'), cex = 1)
legend(-2.3,1.3,legend = c('plasma amino-acids','plasma lactate','plasma fatty acids','whole-blood NAD','saliva'),title='Circulation', fill = c('#a50f15','#de2d26','#fb6a4a','#fcae91','#fee5d9'), cex = 1)
# NOTE: first two arguments of the function "legend" refers to the coordinates at which the legend appear.
# Change them to move the legend in the plot

## to check specific conditions:
print(links2[links2[,6]=='behavior',]) # will display all the links involving a variable belonging to the "behavior" group