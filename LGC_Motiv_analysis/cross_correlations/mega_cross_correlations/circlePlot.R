# 1. Loading libraries
library(circlize)
library(readxl)
library(stringr)

# Reading files, and creating correlation matrix of relevant variables (organized in groups)
corrList = read_xlsx('C:/Users/rdrotari/OneDrive - NESTLE/Projects/CircularPlots/ClairisEPFL/crosscorrel_signed_r_table_75subs.xlsx')
varList = unique(corrList$var1)
behaviorList = varList[startsWith(varList,'behavior_questionnaire')]
brainList = varList[startsWith(varList,'brain')]
plasmaList = varList[startsWith(varList,'plasma')] 
bloodList = varList[startsWith(varList,'wholeB')]
bloodList = bloodList[str_count(bloodList,'_') == 1]
bloodList = bloodList[bloodList != 'wholeB_tNAD']
circulationList = c(plasmaList,bloodList)
goodList = c(behaviorList,brainList,circulationList)
corrMat = matrix(nrow = length(goodList), ncol = length(goodList))
rownames(corrMat) = goodList
colnames(corrMat) = goodList
for (ir in goodList){
  print(ir)
  for (ic in goodList){
    corrMat[ir,ic] = corrList$r_corr[corrList$var1 == ir & corrList$var2 == ic]
  }
}


# 3. Draw circle plot
colorLink <- function(r,Thres){   # Function to determine the color of the link according to the value of the correlation coefficient
  if (r > 0){col0 = '#FF0000'} else {col0 = '#000000'}
  tr = floor((abs(r)-Thres)/(1-Thres)*256)
  if (tr < 16){
    col0 = paste(col0,'0',as.character(as.hexmode(tr)),sep='')
  } else {
    col0 = paste(col0,as.character(as.hexmode(tr)),sep='')
  }
  return(col0)
}

getBGcol <- function(colHEX,nsec){   # Function to get the alternate color of the circular sectors
  bgcol = rep(colHEX,nsec)
  trpvec = c('FF','88')
  for (i in 1:nsec){
    itr = (i-1)%%2+1
    bgcol[i] = paste(bgcol[i],trpvec[itr],sep='')
  }
  return(bgcol)
}


wdt = 1.5 # line width of the link in the plot
rThres = 0.45 # Threshold for the absolute value of the correlation coefficient

# Initialisation of the circular plot
circos.clear()
circos.par$track.height = 0.05
circos.par$circle.margin = 0.01
circos.par$start.degree = 90
circos.initialize(colnames(corrMat), xlim = c(-1,1),)
bgcol1 = getBGcol('#138808',length(behaviorList))
bgcol2 = getBGcol('#FFDD00',length(brainList))
bgcol3 = getBGcol('#0000FF',length(circulationList))
bgcol = c(bgcol1,bgcol2,bgcol3)
txtCol = c(bgcol[2:length(bgcol)],bgcol[1])
groups = rep('',nrow(corrMat))
groups[0.5*length(behaviorList)+1] = 'Behavior'
groups[0.5*length(brainList)+length(behaviorList)+1] = 'Brain'
groups[0.5*length(circulationList)+length(brainList)+length(behaviorList)] = 'Circulation'
circos.track(ylim = c(0, 1), track.height = mm_h(5),
            panel.fun = function(x, y) {circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2] + mm_y(0),groups[CELL_META$sector.numeric.index], 
                                                    facing = 'inside', font = 2,cex = 2.25, col = txtCol[CELL_META$sector.numeric.index])},bg.border = NA)
circos.track(ylim = c(0, 1), track.height = mm_h(4),
             panel.fun = function(x, y) {}, bg.border = NA)
circos.track(ylim = c(-1,1),bg.col=bgcol,
             panel.fun = function(x, y) {circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2] + mm_y(2),CELL_META$sector.numeric.index,
                                                     facing = 'clockwise', cex = 1)})

# creation of the matrix "links", which include all the information about the links in the circular plot
nvar = nrow(corrMat)
groupVar = rep(NA,nvar)
groupVar[which(goodList %in% behaviorList)] = 'behavior'
groupVar[which(goodList %in% brainList)] = 'brain'
groupVar[which(goodList %in% circulationList)] = 'circulation'
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

# loop for drawing the plot
for (il in 1:nrow(links)){
  w1 = which(goodList == links[il,1])
  w2 = which(goodList == links[il,2])
  if (groupVar[w1] != groupVar[w2]){  #plot only link between variables of different groups
    circos.link(sector.index1 = links[il,1], c(-wdt,wdt), sector.index2 = links[il,2], c(-wdt,wdt),col = colorLink(as.numeric(links[il,3]),0.1*rThres))
  }
}
