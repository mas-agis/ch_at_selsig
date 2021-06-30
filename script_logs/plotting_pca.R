##plotting pca for Chinese, Austrian, and other breeds
setwd("D:/maulana/third_project/pca/")
library(ggplot2)
#load dataset
data = read.table("filtered_multi_merged_pca.eigenvec",header=T, sep="\t")
data$Breed = data$FID
data = subset(data, Breed!="Gir")
data = subset(data, Breed!="Hereford")
data = subset(data, Breed!="Bohai")
data = subset(data, Breed!="Boran")
data = subset(data, Breed!="Ogaden")
data = data[-c(105),] #remove BRahman outlier
gp <- ggplot(data,aes(x=data$PC1, y=data$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of chinese and austrian")+
  geom_point(size=3)
gp
data$Breed

##plotting admixture results
library(tidyr)
#K=3
setwd("D:/maulana/third_project/pca/admixture")
tbl=read.table("for_admixture.3.Q", sep="\t")
tbl = subset(tbl, V1!="Gir")
tbl = subset(tbl, V1!="Hereford")
tbl = subset(tbl, V1!="Bohai")
tbl = subset(tbl, V1!="Boran")
tbl = subset(tbl, V1!="Ogaden")
tbl = tbl[-c(105),] #remove BRahman outlier
tbl = separate(tbl, V2, into = c("K1","K2","K3"), sep = " ")

barplot(t(as.matrix(tbl$K1,tbl$K2,tbl$K3)), col=rainbow(3), xlab="Breed", ylab="K = 3",border=NA) #names.arg = data$Breed, 
#introgressed_uoa 