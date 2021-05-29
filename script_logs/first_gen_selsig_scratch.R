##Reading output of SelScan v.1.3.0 (iHS, nSl, iHH12)
library(dplyr)
library(ggplot2)

#isafe
getwd()
setwd("D:/maulana/third_project/isafe/atfl")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$chr <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
combined$pval = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
combined$min_log_pval = -1*log10(combined$pval)
bonf_tres = -1*log10(0.05/nrow(combined))
combined$idu = as.numeric(rownames(combined))
#plot based on min_log_pval
ggplot(combined, aes(idu, min_log_pval, colour = chr)) + 
  geom_point()
#plot based on isafe score
ggplot(combined, aes(idu, V2, colour = chr)) + 
  geom_point() 

#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$V1-1 
#keep only chr, start, and end columns
bed = select(signi_regions, "chr", "start", "V1")
write.table(bed, "atfl_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "chr", "start", "V1", "V2", "min_log_pval")
write.table(bed, "extend_atfl_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#for daf score #!!(not really good for display- scores are discrete)
mu1 = mean(combined$V3, na.rm = T)
sig1 = sd(combined$V3, na.rm = T)
combined$pval1 = pnorm(combined$V3, mu1, sig1, lower.tail = FALSE)
combined$min_log_pval1 = -1*log10(combined$pval1)
bonf_tres = -1*log10(0.05/nrow(combined))
plot(as.numeric(rownames(combined)), combined$min_log_pval1)
plot(as.numeric(rownames(combined)), combined$V3)


#iHS
setwd("D:/maulana/third_project_simul/test_ihs")
hols_29 <- read.table("chr_29.ihs.out")
tail(head(hols_29, n = 120))
hist(hols_29$V6) #unstandardized iHS
#qq-plot
plt_qq <- function(x){
  qqnorm(x, pch = 1, frame = FALSE)
  qqline(x, col = "steelblue", lwd = 2)  
}
plt_qq(hols_29$V6)
#plot unstardized value of iHS
plot(hols_29$V2,hols_29$V6,col="blue")
#empirical cumulative distribution function
x=ecdf(hols_29$V6) 
#p-value based on Gautier&Naves(2011)
hols_29$v7 <- -log10(1-2*abs(x(hols_29$V6)-0.5)) 
#Plot -log(p-values) of iHS
plot(hols_29$V2,hols_29$v7,col="blue")

#putting binning every 10000
hols_29$v8 <- as.integer(hols_29$V2/10000) 
#removing snps which occurs less than 10
hols_29 <- hols_29[hols_29$v8 %in% names(which(table(hols_29$v8)  > 9)), ]
#Summarises by bin numbers
ihs_mean <- hols_29 %>%
  group_by(v8) %>%
  summarise(ihs = mean(V6))
#empirical cumulative distribution function
x=ecdf(ihs_mean$ihs) 
#p-value based on Gautier&Naves(2011)
ihs_mean$p <- -log10(1-2*abs(x(ihs_mean$ihs)-0.5)) 
#Plot -log(p-values) of ihs
plot(ihs_mean$v8,ihs_mean$p,col="blue")

#nsl
setwd("D:/maulana/third_project_simul/nsl")
hols_29 <- read.table("chr_29.nsl.out")
head(hols_29)
summary(hols_29)
hist(hols_29$V6) #unstandardized nsl
#qq-plot
plt_qq <- function(x){
  qqnorm(x, pch = 1, frame = FALSE)
  qqline(x, col = "steelblue", lwd = 2)  
}
plt_qq(hols_29$V6)
#plot unstandardized nsl
plot(hols_29$V2,hols_29$v6,col="blue")
#histogram of unstandardized nsl
hist(hols_29$V6)
#empirical cumulative distribution function
x=ecdf(hols_29$V6) 
#p-value based on Gautier&Naves(2011)
hols_29$v7 <- -log10(1-2*abs(x(hols_29$V6)-0.5)) 
#Plot -log(p-values) of nsl
plot(hols_29$V2,hols_29$v7,col="blue")

#ihh12
setwd("D:/maulana/third_project_simul/ihh12")
hols_29 <- read.table("chr_29.ihh12.out",header = T)
head(hols_29)
hist(hols_29$ihh12) #unstandardized ihh12

#putting binning every 10000
hols_29$v8 <- as.integer(hols_29$pos/10000) 
#removing snps which occurs less than 10
hols_29 <- hols_29[hols_29$v8 %in% names(which(table(hols_29$v8)  > 9)), ]
#Summarises by bin numbers
ihh_mean <- hols_29 %>%
            group_by(v8) %>%
            summarise(ihh_12 = mean(ihh12))
#empirical cumulative distribution function
x=ecdf(ihh_mean$ihh_12) 
#p-value based on Gautier&Naves(2011)
ihh_mean$p <- -log10(1-2*abs(x(ihh_mean$ihh_12)-0.5)) 
#Plot -log(p-values) of nsl
plot(ihh_mean$v8,ihh_mean$p,col="blue")

#putting binning every 10000
hols_29$v8 <- as.integer(hols_29$pos/10000) 
unique(hols_29$v8)
akhir1 <- b %>%  mutate(Cut = cut(Pos, breaks = seq(0, max(Pos, na.rm = T), by = scan))) %>% 
  group_by(Cut,Poin,.drop = FALSE) %>%
  summarise(Ancestral_count = n()) %>%
  filter(Poin == 1)
Data[Data$c %in% names(which(table(Data$c) > 2)), ]

##trial of t-test for each snp - read more here 
#https://stackoverflow.com/questions/40708410/r-how-to-generate-a-vector-of-probabilities-normally-distributed-to-be-used-at-c
x=c(1,2,3,2,4,2,3,1,4,5,1,5,3,5,4,2,3,4,2)
mu <- mean(x)
sig <- sd(x)
null=rnorm(length(x),mean(x),sd(x))
pnorm(x, mu, sig)
y=c(4,4,4,4,4,4,4,4)
pnorm(y, mu, sig)

#trying calculating meta-ss - simulation 
#first test - ihs
z=c(6,6,6,5.2, 5.1, 5.3, 10, 6,6,6,5.2, 5.1, 5.3, 10)
data = (z - mean(z)) / sd(z)
hist(data)
mu = mean(data, na.rm = T)
sig = sd(data, na.rm = T)
pval = pnorm(data, mu, sig, lower.tail = FALSE)
pval
phi = dnorm(data, mu, sig)
phi
Zscore1 = -phi**-1 - (-phi**-1 * pval)
Zscore1
#second test
z=c(2, 1, 0.5, 3, 5, 2, 1, 2, 1, 4, 2, 1 , 2, 3)
data = (z - mean(z)) / sd(z)
hist(data)
mu = mean(data, na.rm = T)
sig = sd(data, na.rm = T)
pval = pnorm(data, mu, sig, lower.tail = FALSE)
pval
phi = dnorm(data, mu, sig)
phi
Zscore2 = -phi**-1 - (-phi**-1 * pval)
Zscore2
#third test
z=c(3, -2, -3, -1, -0.5, 2, 3, 1, -4, -2, 0.8, 1.2, 3.4, 2.7)
data = (z - mean(z)) / sd(z)
hist(data)
mu = mean(data, na.rm = T)
sig = sd(data, na.rm = T)
pval = pnorm(data, mu, sig, lower.tail = FALSE)
pval
phi = dnorm(data, mu, sig)
phi
Zscore3 = -phi**-1 - (-phi**-1 * pval)
Zscore3
#meta-ss
meta_ss = (Zscore1 + Zscore2 + Zscore3) / sqrt((1+1+1)**2)
#referring back meta-ss to standard normal distribution
meta_ss_norm = (meta_ss - mean(meta_ss)) / sd(meta_ss)
meta_ss_norm
#getting p-value for each variables in meta_ss_norm both upper-lower tails
meta_ss_pval = pnorm(meta_ss_norm, lower.tail = TRUE)
meta_ss_pval

#trying calculating meta-ss - simulation - windows bin 
#first test - ihs
z = read.csv("D:/maulana/third_project/norm_ihs/atfl/atfl_29.ihs.out.100bins.norm.10kb.windows", sep="\t", header = FALSE)[6]
hist(z$V6)
mu = mean(z$V6, na.rm = T)
sig = sd(z$V6, na.rm = T)
z$pval = pnorm(z$V6, mu, sig, lower.tail = FALSE)
head(z)
phi = dnorm(data, mu, sig)
phi
Zscore1 = -phi**-1 - (-phi**-1 * pval)
Zscore1
#second test
z=c(2, 1, 0.5, 3, 5, 2, 1, 2, 1, 4, 2, 1 , 2, 3)
data = (z - mean(z)) / sd(z)
z = read.csv("D:/maulana/third_project/nSL/atfl/atfl_29.nsl.out", sep="\t", header = FALSE)[6]
data = (z$V6 - mean(z$V6)) / sd(z$V6)
hist(data)
pval = pnorm(data, mu, sig, lower.tail = FALSE)
pval
phi = dnorm(data, mu, sig)
phi
Zscore2 = -phi**-1 - (-phi**-1 * pval)
Zscore2
#third test
z=c(3, -2, -3, -1, -0.5, 2, 3, 1, -4, -2, 0.8, 1.2, 3.4, 2.7)
data = (z - mean(z)) / sd(z)
hist(data)
pval = pnorm(data, mu, sig, lower.tail = FALSE)
pval
phi = dnorm(data, mu, sig)
phi
Zscore3 = -phi**-1 - (-phi**-1 * pval)
Zscore3
#meta-ss
meta_ss = (Zscore1 + Zscore2 + Zscore3) / sqrt((1+1+1)**2)
#referring back meta-ss to standard normal distribution
meta_ss_norm = (meta_ss - mean(meta_ss)) / sd(meta_ss)
meta_ss_norm
#getting p-value for each variables in meta_ss_norm both upper-lower tails
meta_ss_pval = pnorm(meta_ss_norm, lower.tail = TRUE)
meta_ss_pval