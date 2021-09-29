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
data = data[-c(108),] #remove BRahman outlier
data = subset(data, Breed!="Ogaden")
unique(data$Breed)

##Plot PC1 and PC2 of all breeds
gp <- ggplot(data,aes(x=data$PC1, y=data$PC2, group=Breed, color=Breed, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of chinese and austrian")+
  geom_point(size=3)
gp

#plot color based on taurine/indicine/admixed data
data1 = read.table("breeds_group", col.names = c("Breed", "Subspecies"), sep = "\t")
#data1$Breed = as.factor(data1$Breed)
test = left_join(data, data1, by="Breed")
test = subset(test, Breed!="Dabieshan")
test$subspecies = as.factor(as.character(test$subspecies))
##Plot PC1 and PC2 of all breeds
gp <- ggplot(test,aes(x=test$PC1, y=test$PC2, group=Subspecies, color=Subspecies, shape=Breed)) +
  scale_shape_manual(values=1:nlevels(data$Breed)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(x="Component 1", y=" Component 2", title = "pca of chinese and austrian")+
  geom_point(size=3)
gp
