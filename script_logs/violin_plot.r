#Violin plot of FPKM
library(reshape)
setwd("D:/maulana/third_project")
data = read.table("summary_unique_genes_lenovo_&_FPKM.txt", skip = 2, fill = TRUE, na.strings = "NA")
#colnames with FPKM values
nama_col = c("V2", "V5", "V8", "V11", "V14", "V17", "V20", "V23", "V26", "V29", "V32", "V35", "V38", "V41", "V44", "V47")
#new dataset with only FPKM values
mydata = data[nama_col]
#transpose data
tes_data =  data.frame(t(mydata))
#set id column
tes_data$id = row.names(tes_data)
#melting data
mdata <- melt(tes_data, id="id") 
#remove rows with NA values
mdata = na.omit(mdata)
#set values as numeric and id as factor 
mdata$value = as.numeric(as.character(mdata$value))
mdata$id = as.factor(as.character(mdata$id))
#label names
lab =c("ihs_FL", "ihs_ME", "ihs_YB", "ihs_CHBI", "ihs_CHBI_M", "ihs_CHBI_L", 
       "iS_FL", "iS_HA", "iS_ME", "iS_YB", "iS_CHBI_M", "iS_CHBI_L", 
       "fs_FL", "fs_HA", "fs_ME", "fs_YB")
#violin plot the data
library(ggplot2)
# Basic violin plot
p <- ggplot(mdata, aes(x=id, y=value)) + 
  geom_violin()
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  scale_x_discrete(breaks= nama_col, labels= lab)

##frequency of tissue of associated genes from each test
head(data)
hist(data$V48)
