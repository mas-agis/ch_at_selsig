##plotting admixture results

library(tidyr)
setwd("D:/maulana/third_project/pca/admixture")
list.files()

#K=2
tbl2 = read.table("for_admixture.2.Q", sep="\t")
tbl2 = subset(tbl2, V1!="Gir")
tbl2 = subset(tbl2, V1!="Hereford")
tbl2 = subset(tbl2, V1!="Bohai")
tbl2 = subset(tbl2, V1!="Boran")
tbl2 = subset(tbl2, V1!="Ogaden")
tbl2 = tbl2[-c(105),] #remove BRahman outlier
tbl2 = separate(tbl2, V2, into = c("K1","K2"), sep = " ")
tbl2 = tbl2[order(tbl2$V1),]

name_bars <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- NA
  }
  return(retVec)
}
lab = name_bars(tbl2$V1)
gapping <- function(vec) {
  retGaps = rep(2.5, length(vec))
  tmp = name_bars(tbl2$V1)
  for (k in 1:length(retGaps)) {
    if ( is.na(tmp[k]))
      retGaps[k] = 0
  }
  return(retGaps)
}
gaps <- gapping(tbl2$V1)
barplot(t(as.matrix(tbl2[, 2:3])), ylab="K = 2",
        space = gaps, col = rainbow(2),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75,)

#K=3
tbl3=read.table("for_admixture.3.Q", sep="\t")
tbl3 = subset(tbl3, V1!="Gir")
tbl3 = subset(tbl3, V1!="Hereford")
tbl3 = subset(tbl3, V1!="Bohai")
tbl3 = subset(tbl3, V1!="Boran")
tbl3 = subset(tbl3, V1!="Ogaden")
tbl3 = tbl3[-c(105),] #remove BRahman outlier
tbl3 = separate(tbl3, V2, into = c("K1","K2","K3"), sep = " ")
tbl3 = tbl3[order(tbl3$V1),]
  
name_bars <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- NA
  }
  return(retVec)
}
lab = name_bars(tbl3$V1)
gapping <- function(vec) {
  retGaps = rep(2.5, length(vec))
  tmp = name_bars(tbl3$V1)
  for (k in 1:length(retGaps)) {
    if ( is.na(tmp[k]))
      retGaps[k] = 0
  }
  return(retGaps)
}
gaps <- gapping(tbl3$V1)
barplot(t(as.matrix(tbl3[, 2:4])), ylab="K = 3",
        space = gaps, col = rainbow(3),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75,)

#K=4
tbl4=read.table("for_admixture.4.Q", sep="\t")
tbl4 = subset(tbl4, V1!="Gir")
tbl4 = subset(tbl4, V1!="Hereford")
tbl4 = subset(tbl4, V1!="Bohai")
tbl4 = subset(tbl4, V1!="Boran")
tbl4 = subset(tbl4, V1!="Ogaden")
tbl4 = tbl4[-c(105),] #remove BRahman outlier
tbl4 = separate(tbl4, V2, into = c("K1","K2","K3", "K4"), sep = " ")
tbl4 = tbl4[order(tbl4$V1),]

name_bars <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- NA
  }
  return(retVec)
}
lab = name_bars(tbl4$V1)
gapping <- function(vec) {
  retGaps = rep(2.5, length(vec))
  tmp = name_bars(tbl4$V1)
  for (k in 1:length(retGaps)) {
    if ( is.na(tmp[k]))
      retGaps[k] = 0
  }
  return(retGaps)
}
gaps <- gapping(tbl4$V1)
barplot(t(as.matrix(tbl4[, 2:5])), ylab="K = 4",
        space = gaps, col = rainbow(5),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75)


#K=5
tbl5=read.table("for_admixture.5.Q", sep="\t")
tbl5 = subset(tbl5, V1!="Gir")
tbl5 = subset(tbl5, V1!="Hereford")
tbl5 = subset(tbl5, V1!="Bohai")
tbl5 = subset(tbl5, V1!="Boran")
tbl5 = subset(tbl5, V1!="Ogaden")
tbl5 = tbl5[-c(105),] #remove BRahman outlier
tbl5 = separate(tbl5, V2, into = c("K1","K2","K3", "K4", "K5"), sep = " ")
tbl5 = tbl5[order(tbl5$V1),]

name_bars <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- NA
  }
  return(retVec)
}
lab = name_bars(tbl5$V1)
gapping <- function(vec) {
  retGaps = rep(2.5, length(vec))
  tmp = name_bars(tbl5$V1)
  for (k in 1:length(retGaps)) {
    if ( is.na(tmp[k]))
      retGaps[k] = 0
  }
  return(retGaps)
}
gaps <- gapping(tbl5$V1)
barplot(t(as.matrix(tbl5[, 2:6])), ylab="K = 5",
        space = gaps, col = rainbow(5),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75)


#Multi subplots in a single r-graph
par(mfrow=c(2,1)) #defining rows and column in the figure
barplot(t(as.matrix(tbl3[, 2:4])), ylab="K = 3",
        space = gaps, col = rainbow(3),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75,)
barplot(t(as.matrix(tbl5[, 2:6])), ylab="K = 5",
        space = gaps, col = rainbow(5),
        names.arg = lab,
        las = 2,
        cex.axis = 0.1,
        cex.names = 0.75)



