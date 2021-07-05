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

