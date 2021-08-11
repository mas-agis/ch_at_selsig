##Reading output of SelScan v.1.3.0 (iHS, nSl, iHH12)
library(dplyr)
library(ggplot2)
library(qqman)
vignette('qqman')

#meta_ss for atfl
setwd("D:/maulana/third_project/meta_ss")
list.files()
meta_ss_atfl = read.csv("atfl_meta_ss.csv", sep='\t')
#calculate p-value
meta_ss_atfl$P = pnorm(meta_ss_atfl$meta_ss)
#rename column 'pos' as BP
names(meta_ss_atfl)[2] = "BP"
#rename column 'chr' as CHR
names(meta_ss_atfl)[1] = "CHR"
#defined order of rows as SNP names
meta_ss_atfl$SNP = as.numeric(rownames(meta_ss_atfl))
#plot using qqman
manhattan(meta_ss_atfl, suggestiveline = TRUE)

#isafe for atfl
getwd()
setwd("D:/maulana/third_project/isafe/atfl")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp/isafe")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_atfl_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chha
getwd()
setwd("D:/maulana/third_project/isafe/chha")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp/isafe")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_chha_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chme
getwd()
setwd("D:/maulana/third_project/isafe/chme")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp/isafe")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_chme_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chya
getwd()
setwd("D:/maulana/third_project/isafe/chya")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp/isafe")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_chya_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chbt
getwd()
setwd("D:/maulana/third_project/isafe/chbt")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp/isafe")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_chbt_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chbi
getwd()
setwd("D:/maulana/third_project/isafe/chbi")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
#calculating min_log_pval_and bonf_tres
bonf_tres = -log10(5e-08) #according to manhattan
combined$min_log_pval = -1*log10(combined$P)
#bonf_tres = -1*log10(0.05/nrow(combined))
#Extracting snps higher than suggestive line instead of bonf treshold!!
signi_regions = filter(combined, min_log_pval > -log10(1e-05))
signi_regions$start = signi_regions$BP-1 
setwd("D:/maulana/third_project/by_snp")
#keep chr, start, end, isafe score, and min_log_pval columns 
bed = select(signi_regions, "CHR", "start", "BP", "V2", "min_log_pval")
write.table(bed, "extend_chbi_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#isafe for chbi_low
getwd()
setwd("D:/maulana/third_project/isafe/chbi_low")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))


#isafe for chbi_med
getwd()
setwd("D:/maulana/third_project/isafe/chbi_med")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("final_", i, ".txt")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for isafe score
mu = mean(combined$V2, na.rm = T)
sig = sd(combined$V2, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V2, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))

#############################################################
##############for iHS tests##################################
#iHS for atfl
getwd()
setwd("D:/maulana/third_project/iHS/atfl")
list.files()
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("atfl_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
#manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_atfl = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_atfl, "ihs_atfl_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
#manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_atfl = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_atfl, "ihs_atfl_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chbt
getwd()
setwd("D:/maulana/third_project/iHS/chbt")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chbt_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbt = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbt, "ihs_chbt_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbt = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbt, "ihs_chbt_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chbi
getwd()
setwd("D:/maulana/third_project/iHS/chbi")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chbi_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbi = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbi, "ihs_chbi_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbi = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbi, "ihs_chbi_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chha
getwd()
setwd("D:/maulana/third_project/iHS/chha")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chha_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chha = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chha, "ihs_chha_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chha = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chha, "ihs_chha_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chme
getwd()
setwd("D:/maulana/third_project/iHS/chme")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chme_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chme = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chme, "ihs_chme_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chme = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chme, "ihs_chme_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chya
getwd()
setwd("D:/maulana/third_project/iHS/chya")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chya_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")

#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chya = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chya, "ihs_chya_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chya = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chya, "ihs_chya_both_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chbi_low
getwd()
setwd("D:/maulana/third_project/iHS/chbi_low")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chbi_low_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
  rm (temp)
}
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")
#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbi_low = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbi_low, "ihs_chbi_low_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#iHS for chbi_med
getwd()
setwd("D:/maulana/third_project/iHS/chbi_med")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chbi_med_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
  rm (temp)
}
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#setwd for saving the outputs
setwd("D:/maulana/third_project/by_snp/ihs")
#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(5e-08))
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbi_med = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbi_med, "ihs_chbi_med_right_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

############################################################################
#####Plotting Fst against CHBI##############################################
#ATFL
setwd("D:/maulana/third_project/10_kb/fst/atfl")
raw_files = list.files(pattern = "*_chbi_*")
files = grep(".fst", raw_files, value = TRUE)
rm(raw_files)
files
combined = data.frame()
for (i in files) {
  print(paste0("this is file ", i))
  temp = read.table(i, header = TRUE)
  combined = rbind(combined, temp)
  rm(temp)
}
head(combined)
hist(as.numeric(combined$MEAN_FST))
#for fst score 
mu = mean(combined$MEAN_FST, na.rm = T)
sig = sd(combined$MEAN_FST, na.rm = T)
#calculating p-value
combined$P = pnorm(combined$MEAN_FST, mu, sig, lower.tail = FALSE)
combined$snp = row.names(combined)
#plot using qqman
manhattan(combined, chr= "CHROM", bp="BIN_START", p="P", snp= "snp", suggestiveline = FALSE)
#CHHA
setwd("D:/maulana/third_project/10_kb/fst/chha")
raw_files = list.files(pattern = "*_chbi_*")
files = grep(".fst", raw_files, value = TRUE)
rm(raw_files)
files
combined = data.frame()
for (i in files) {
  print(paste0("this is file ", i))
  temp = read.table(i, header = TRUE)
  combined = rbind(combined, temp)
  rm(temp)
}
head(combined)
hist(as.numeric(combined$MEAN_FST))
#for fst score 
mu = mean(combined$MEAN_FST, na.rm = T)
sig = sd(combined$MEAN_FST, na.rm = T)
#calculating p-value
combined$P = pnorm(combined$MEAN_FST, mu, sig, lower.tail = FALSE)
combined$snp = row.names(combined)
#plot using qqman
manhattan(combined, chr= "CHROM", bp="BIN_START", p="P", snp= "snp", suggestiveline = FALSE)
#CHME
setwd("D:/maulana/third_project/10_kb/fst/chme")
raw_files = list.files(pattern = "*_chbi_*")
files = grep(".fst", raw_files, value = TRUE)
rm(raw_files)
files
combined = data.frame()
for (i in files) {
  print(paste0("this is file ", i))
  temp = read.table(i, header = TRUE)
  combined = rbind(combined, temp)
  rm(temp)
}
head(combined)
hist(as.numeric(combined$MEAN_FST))
#for fst score 
mu = mean(combined$MEAN_FST, na.rm = T)
sig = sd(combined$MEAN_FST, na.rm = T)
#calculating p-value
combined$P = pnorm(combined$MEAN_FST, mu, sig, lower.tail = FALSE)
combined$snp = row.names(combined)
#plot using qqman
manhattan(combined, chr= "CHROM", bp="BIN_START", p="P", snp= "snp", suggestiveline = FALSE)
#CHYA
setwd("D:/maulana/third_project/10_kb/fst/chya")
raw_files = list.files(pattern = "*_chbi_*")
files = grep(".fst", raw_files, value = TRUE)
rm(raw_files)
files
combined = data.frame()
for (i in files) {
  print(paste0("this is file ", i))
  temp = read.table(i, header = TRUE)
  combined = rbind(combined, temp)
  rm(temp)
}
head(combined)
hist(as.numeric(combined$MEAN_FST))
#for fst score 
mu = mean(combined$MEAN_FST, na.rm = T)
sig = sd(combined$MEAN_FST, na.rm = T)
#calculating p-value
combined$P = pnorm(combined$MEAN_FST, mu, sig, lower.tail = FALSE)
combined$snp = row.names(combined)
#plot using qqman
manhattan(combined, chr= "CHROM", bp="BIN_START", p="P", snp= "snp", suggestiveline = FALSE)

#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
#Extracting snps higher than genome-wide significant line
signi_regions = filter(combined, min_log_pval > -log10(0.05/length(row.names(combined))))
#keep chr, start, end, iHS, and min_log_pval columns 
bed_fst_chya_chbi = select(signi_regions, "CHROM", "BIN_START", "BIN_END", "MEAN_FST", "min_log_pval")
write.table(bed_fst_chya_chbi, "bed_fst_chya_chbi", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)









#empiricial cumulative distribution function
X = rnorm(100)
plot(X)
P = ecdf(X)
P(0)
P(-6)
plot(P)

#iHS for chbi
setwd("D:/maulana/third_project/iHS/chbi")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("chbi_", i, ".ihs.out")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#for iHS score 
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#calculate p-value - right tail only
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[2] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculate p-value - both tails 
combined$P = pnorm(combined$V6, mu, sig, lower.tail = TRUE)
combined$phi = pnorm(combined$P, mean(combined$P, na.rm = T), sd(combined$P, na.rm = T), lower.tail = TRUE)
head(combined)
combined$temp = (-1*(1/combined$phi)*1) - (-1*(1/combined$phi)*combined$P) 
combined$min_log_temp = -1*log10(combined$temp)

#histogram of p-value
hist()
#plot using qqman
manhattan(combined, suggestiveline = FALSE)
#calculating min_log_pval_and bonf_tres
combined$min_log_pval = -1*log10(combined$P)
bonf_tres = -1*log10(0.05/nrow(combined))
bonf_tres = -log10(5e-08)
suggestiveline = -log10(1e-05)
#Extracting snps higher than bonferroni treshold
signi_regions = filter(combined, min_log_pval > bonf_tres)
signi_regions = filter(combined, min_log_pval > suggestiveline)
signi_regions$start = signi_regions$BP-1 
#keep chr, start, end, iHS, and min_log_pval columns 
bed_chbi = select(signi_regions, "CHR", "start", "BP", "V6", "min_log_pval")
write.table(bed_chbi, "ihs_chbi_bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



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

######################################################################################################
###Scratch
#plot based on min_log_pval
ggplot(combined, aes(idu, min_log_pval, colour = chr)) + 
  geom_point() + geom_hline(yintercept=bonf_tres, linetype="dashed",color = "red")
#plot based on isafe score
ggplot(combined, aes(idu, V2, colour = chr)) + 
  geom_point() 
#manhattan plot using ggplot2
#https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)
gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

data_cum <- combined %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(V1)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)
data_cum$chr <- as.factor(data_cum$chr)
combined$chr <- as.factor(combined$chr)
data_cum
as.data.frame(data_cum)
merge(combined, as.data.frame(data_cum), by="chr")
combined <- combined %>% 
  inner_join(combined, data_cum, by = chr) %>% 
  mutate(bp_cum = V1 + bp_add)
axis_set <- combined %>% 
  group_by(chr) %>% 
  summarize(center = mean(V1))
ggplot(combined, aes(x = idu, y = min_log_pval, 
                     color = chr)) + #, size = -log10(p))) +
  geom_point() +
  geom_hline(yintercept = bonf_tres, color = "red", linetype = "dashed") + 
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) #+
scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

###Another trials for calculating the meta-SS
##atfl
#reading iHs 
setwd("D:/maulana/third_project/iHS/atfl")
data_ihs = data.frame()
for (chr in 1:29) {
  file = paste0("atfl_",chr,".ihs.out")
  temp = read.table(file, colClasses = c("NULL", "integer", "NULL", "NULL", "NULL", "character"))
  temp$chr = chr
  data_ihs = rbind(data_ihs, temp)
}
#reading iSAFE 
setwd("D:/maulana/third_project/isafe/atfl")
data_comb = data.frame()
for (chr in 1:29) {
  file = paste0("final_",chr,".txt")
  temp = read.table(file, colClasses = c("integer", "NULL", "numeric"))
  temp$chr = chr
  data_comb = rbind(data_comb, temp)
}
unique(data_comb$chr)
hist(as.numeric(data_comb$V6))

##reading normalized iHS score with non-overlapped 10Kb window
setwd("D:/maulana/third_project/10_kb/norm_ihs/atfl")
files = list.files(pattern = "*.windows")
combined = data.frame()
for (i in 1:29) {
  print(paste0("this is file ", i))
  filename = paste0("atfl_", i, ".ihs.out.100bins.norm.10kb.windows")
  temp = read.table(filename)
  temp$CHR <- i
  combined = rbind(combined, temp)
}
rm (temp)
#plot histogram
hist(combined$V6)
#filter regions having iHS score more than 3
signi_regions = filter(combined, V6 > 3)
#writing signi_regions for bed annotation
bed = select(signi_regions, "CHR", "V1", "V2", "V6")
write.table(bed, "iHS_normalized_atfl.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#run the command below in linux for annotation
#java -Xmx8g -jar snpEff.jar -i bed ARS-UCD1.2.99 ~/data/third_project/10_kb/norm_ihs/atfl/iHS_normalized_atfl.txt > ~/data/third_project/10_kb/norm_ihs/atfl/iHS_normalized_atfl_annotated
#for iHS score
mu = mean(combined$V6, na.rm = T)
sig = sd(combined$V6, na.rm = T)
#calculate p-value
combined$P = pnorm(combined$V6, mu, sig, lower.tail = FALSE)
#rename column 'V1' as BP
names(combined)[1] = "BP"
#defined order of rows as SNP names
combined$SNP = as.numeric(rownames(combined))
#plot using qqman
manhattan(combined, suggestiveline = FALSE, ylim = c(0, 30))
