##plotting pca for Chinese, Austrian, and other breeds
setwd("D:/maulana/third_project/pca/")
library(ggplot2)
library(dplyr)
#load dataset
data = read.table("filtered_multi_merged_pca.eigenvec",header=T, sep="\t")
data$Breed = data$FID
#adding group_id
data$group_id = data %>% group_indices(Breed)
data$group_id = as.factor(data$group_id)
#remove several breeds due to max. 30 combintaion of shape n color
data = subset(data, Breed!="Hereford")
data = subset(data, Breed!="Bohai")
data = subset(data, Breed!="Boran")
data = subset(data, Breed!="Shorthorn")
#reset rownames
rownames(data) = 1:nrow(data)
#remove brahman outlier
data = data[-c(103),] #remove BRahman outlier

##Plot PC1 and PC2 of all breeds
gp <- ggplot(data,aes(x=data$PC1, y=data$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of chinese and austrian")+
  geom_point(size=3)
gp

##Plot PC1 and PC2 of Taurine breeds 
#make copy of dataset
taurus = data
#rows with indicine breeds
indicine_index = c(6:11, 22:25, 42:45, 58:61, 91:93, 95:127, 130:165, 172:178)
#set NA for PCA1 n PC2 of indicine breeds 
taurus[indicine_index, 3] = NA
taurus[indicine_index, 4] = NA 
#plot
gp <- ggplot(taurus,aes(x=taurus$PC1, y=taurus$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(taurus$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of taurine ")+
  geom_point(size=3)
gp

##Plot PC1 and PC2 of indicine breeds 
#make copy of dataset
indicine = data
taurine_index = c(1:5, 12:21, 26:41, 46:57, 62:90, 94, 127:129, 165:171, 178:187) 
#set NA for PCA1 n PC2 of indicine breeds 
indicine[taurine_index, 3] = NA
indicine[taurine_index, 4] = NA 
#plot
gp <- ggplot(indicine,aes(x=indicine$PC1, y=indicine$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(indicine$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of taurine ")+
  geom_point(size=3)
gp

