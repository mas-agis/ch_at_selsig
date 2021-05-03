##Reading output of SelScan v.1.3.0 (iHS, nSl, iHH12)
library(dplyr)

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
z=c(6,6,6,5.2, 5.1, 5.3, 10)
pnorm(z, mu, sig, lower.tail = FALSE)
pnorm(z, mu, sig, lower.tail = TRUE)
diff(pnorm(z, mu, sig, lower.tail = TRUE))
a=c(1,1,1,1,1,1,1,1,1,1,1,1)
b=c(3,3,3,3,3,3,3,3)
c=c(4,4,1,1,2,2,3,3,3,3,5,5,4,6)
pnorm(z, 0, 1, lower.tail = FALSE)
c=c(4,4,-1,-1,2,-2,3,-3,-3,3,5,5,4,6)
pnorm(z, mean(c), sd(c), lower.tail = FALSE)
d=c(4,4,4,4,4,4,4,4)

setwd("D:/maulana/third_project/norm_ihs/atfl")
list.files()
ihs <- read.csv("atfl_29.ihs.out.100bins.norm.10kb.windows", sep="\t", header = FALSE)
ihs$pval <- pnorm(ihs$V6, lower.tail = FALSE, log.p = TRUE)
plot(ihs$V1, -(ihs$pval))
bonferonni <- -log10(0.05/length(ihs))
bonferonni <- -log10(0.05/1000000)

mu <- mean(ihs$V6, na.rm = TRUE)
sig <- sd(ihs$V6, na.rm = TRUE)
ihs$norm <- (ihs$V6 - mu)/sig

t.test(y, x)
t.test(y,null)
t.test(z, x)
t.test(z, null)
t.test(a, x)
t.test(a, null)
t.test(b, x)
t.test(b ,null)
t.test(c, x)
t.test(c, null)
t.test(d, null)
t.test(x, null)

val=c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)
calc <- function(x){
  1-2*abs(x-0.5)
}
d = calc(val)
