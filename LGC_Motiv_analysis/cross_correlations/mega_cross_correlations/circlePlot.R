# 1. Loading libraries
library(circlize)
library(readxl)
library(stringr)

# 2. Reading input file
corrList = read_xlsx('C:/Users/rdrotari/OneDrive - NESTLE/Projects/CircularPlots/ClairisEPFL/crosscorrel_signed_r_table_75subs.xlsx')

# 3. Selecting relevant variables and organizing them in groups
varList = unique(corrList$var1)
behaviorList0 = varList[startsWith(varList,'behavior_questionnaire')] # selection of variables associated to behavior
behaviorList1 = sapply(behaviorList0,function(x) gsub('behavior_questionnaires_','',x)) # shortening the variable names
behaviorgroup1 = which(startsWith(behaviorList1,'stress_anxiety_')) # definition of group 1 in behavior: stress/anxiety
behaviorgroup2 = which(startsWith(behaviorList1,'motivation_')) # definition of group 2 in behavior: motivation
behaviorgroup3 = which(startsWith(behaviorList1,'dominance_')) # definition of group 3 in behavior: dominance
behaviorList2 = behaviorList1 # shortening the variable names
behaviorList2[behaviorgroup1] = sapply(behaviorList1[behaviorgroup1],function(x) gsub('stress_anxiety_','',x))
behaviorList2[behaviorgroup2] = sapply(behaviorList1[behaviorgroup2],function(x) gsub('motivation_','',x))
behaviorList2[behaviorgroup3] = sapply(behaviorList1[behaviorgroup3],function(x) gsub('dominance_','',x))
for (var in behaviorList0){ # replacing short names in corrList [dataframe with value of correlation coefficient]
  ivar = which(behaviorList0 == var)
  w1 = which(corrList$var1 == var)
  corrList$var1[w1] = behaviorList2[ivar]
  w2 = which(corrList$var2 == var)
  corrList$var2[w2] = behaviorList2[ivar]
}
brainList0 = varList[startsWith(varList,'brain')] # selection of variables associated to brain
brainList1 = sapply(brainList0,function(x) gsub('brainM_','',x)) # shortening the variable names
braingroup1 = which(startsWith(brainList1,'dmPFC_')) # definition of group 1 in brain: dmPFC
braingroup2 = which(startsWith(brainList1,'aIns_')) # definition of group 2 in brain: aIns
brainList2 = brainList1
for (var in brainList0){ # replacing short names in corrList [dataframe with value of correlation coefficient]
  ivar = which(brainList0 == var)
  w1 = which(corrList$var1 == var)
  corrList$var1[w1] = brainList2[ivar]
  w2 = which(corrList$var2 == var)
  corrList$var2[w2] = brainList2[ivar]
}
plasmaList0 = varList[startsWith(varList,'plasma')] # selection of variables associated to plasma
plasmaList1 = sapply(plasmaList0,function(x) gsub('plasma_','',x)) # shortening the variable names
bloodList0 = varList[startsWith(varList,'wholeB')] # selection of variables associated to whole blood
bloodList0 = bloodList0[str_count(bloodList0,'_') == 1] # removing not interesting variables
bloodList0 = bloodList0[bloodList0 != 'wholeB_tNAD']
bloodList1 = sapply(bloodList0,function(x) gsub('wholeB_','',x)) # shortening the variable names
circulationList0 = c(plasmaList0,bloodList0) # creating the list of variable associated to circulation (= plasma + whole blood)
circulationList1 = c(plasmaList1,bloodList1) # creating the list of variable associated to circulation (= plasma + whole blood)
circulationgroup1 = which(startsWith(circulationList1,'aa_')) # definition of group 1 in circulation: aminoacids
circulationgroup2 = which(startsWith(circulationList1,'Lac')) # definition of group 2 in circulation: lactic acid
circulationgroup3 = which(startsWith(circulationList1,'fa_')) # definition of group 3 in circulation: fatty acids
circulationgroup4 = which(circulationList1 %in% bloodList1) # definition of group 4 in circulation: NAD derivatives
circulationList2 = circulationList1 # shortening the variable names
circulationList2[circulationgroup1] = sapply(circulationList1[circulationgroup1],function(x) gsub('aa_','',x))
circulationList2[circulationgroup3] = sapply(circulationList1[circulationgroup3],function(x) gsub('fa_','',x))
for (var in circulationList0){ # replacing short names in corrList [dataframe with value of correlation coefficient]
  ivar = which(circulationList0 == var)
  w1 = which(corrList$var1 == var)
  corrList$var1[w1] = circulationList2[ivar]
  w2 = which(corrList$var2 == var)
  corrList$var2[w2] = circulationList2[ivar]
}
goodList = c(behaviorList2,brainList2,circulationList2) # creation of the "good" list of variables (= behavior + brain + circulation)

# 3. Creation of the correlation matrix
corrMat = matrix(nrow = length(goodList), ncol = length(goodList)) # initialization
rownames(corrMat) = goodList #definition of the row names (= variable names)
colnames(corrMat) = goodList #definition of the columns names (= variable names)
for (ir in goodList){ # retrieving of correlation coefficients and put in the matrix
  print(ir)
  for (ic in goodList){
    corrMat[ir,ic] = corrList$r_corr[corrList$var1 == ir & corrList$var2 == ic]
  }
}


# 4. Function to determine the color of the link according to the value of the correlation coefficient
colorLink <- function(r,Thres,rMax){   # Function to determine the color of the link according to the value of the correlation coefficient
  if (r > 0){col0 = '#FF0000'} else {col0 = '#000000'}
  tr = floor((abs(r)-Thres)/(rMax-Thres)*256)
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
rThres = 0.45 # Threshold for the absolute value of the correlation coefficient

# 6. Creation of the circular plot
circos.clear()
circos.par$track.height = 0.05 # height of the circular sector
circos.par$circle.margin = 0.35 # margin of the figure
circos.par$start.degree = 45 # start point (measure in degrees): the full list starts from behavior variables
circos.initialize(colnames(corrMat), xlim = c(-1,1),) 
bgcol1 = rep(NA,length(behaviorList2)) # definition of the colors for circular sectors
bgcol1[behaviorgroup1] = '#138808FF' #green colors for behavior variables: different shades stands for different subgroups
bgcol1[behaviorgroup2] = '#13880888'
bgcol1[behaviorgroup3] = '#13880844'
bgcol2 = rep(NA,length(brainList2))
bgcol2[braingroup1] = '#FFDD00FF' #yellow colors for brain variables: different shades stands for different subgroups
bgcol2[braingroup2] = '#FFDD0066'
bgcol3 = rep(NA,length(circulationList2))
bgcol3[circulationgroup1] = '#0000FFFF' #blue colors for circulation variables: different shades stands for different subgroups
bgcol3[circulationgroup2] = '#0000FFAA'
bgcol3[circulationgroup3] = '#0000FF77'
bgcol3[circulationgroup4] = '#0000FF33'
bgcol = c(bgcol1,bgcol2,bgcol3)

# 7. Drawing circular sectors and writing variable names
circos.track(ylim = c(0, 1), track.height = mm_h(4),
             panel.fun = function(x, y) {}, bg.border = NA)
circos.track(ylim = c(-1,1),bg.col=bgcol,
             panel.fun = function(x, y) {circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2] + mm_y(2),CELL_META$sector.index,
                                                     facing = 'clockwise', cex = 1, adj = c(0,0.5))})

# 8. Creation of the matrix "links" with pairs of correlated variable (r > rThres)
nvar = nrow(corrMat)
groupVar = rep(NA,nvar)
groupVar[which(goodList %in% behaviorList2)] = 'behavior'
groupVar[which(goodList %in% brainList2)] = 'brain'
groupVar[which(goodList %in% circulationList2)] = 'circulation'
nlinks = (sum(abs(corrMat)>rThres)-nvar)/2
links = matrix(NA,nrow = nlinks,ncol=3)
w = which(abs(corrMat)>rThres)
il = 0
for (iw in w){
  i1 = floor((iw-1)/nvar)+1
  i2 = (iw-1)%%nvar+1
  if (i1>i2){
    il = il+1
    r = corrList$r_corr[corrList$var1 == rownames(corrMat)[i1] & corrList$var2 == colnames(corrMat)[i2]]
    links[il,] = c(rownames(corrMat)[i1],colnames(corrMat)[i2],r)
  }
}

# 9. Drawing links in the circular plot
for (il in 1:nrow(links)){
  w1 = which(goodList == links[il,1])
  w2 = which(goodList == links[il,2])
  if (groupVar[w1] != groupVar[w2]){  #plot only link between variables of different groups
    circos.link(sector.index1 = links[il,1], c(-wdt,wdt), sector.index2 = links[il,2], c(-wdt,wdt),col = colorLink(as.numeric(links[il,3]),rThres,0.7))
  }
}

# 10. ploitting legend
legend(1.3,-0.1,legend = c('stress/anxiety','motivation','dominance'),title='Behavior', fill = c('#138808FF','#13880888','#13880844'), cex = 1)
legend(-1.3,-0.8,legend = c('dmPCF','aIns'),title='Brain', fill = c('#FFDD00FF','#FFDD0066'), cex = 1)
legend(-1.5,0.9,legend = c('Aminoacids','Lac','Fatty Acid','NAD der.'),title='Circulation', fill = c('#0000FFFF','#0000FFAA','#0000FF77','#0000FF33'), cex = 1)
# NOTE: first two arguments of the function "legend" refers to the coordinates at which the legend appear.
# Change them to move the legend in the plot