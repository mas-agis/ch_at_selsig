library(dplyr)
library(ggplot2)

#combined chr 1 from several breeds
setwd("D:/maulana/third_project/smc")
smc = read.table("combined.csv", sep=",", header = TRUE, 
                 colClasses = c("character", "character", "character", "NULL", "NULL"))
smc$label = as.factor(smc$label)
smc$x = log10(as.integer(smc$x))
smc = smc[-1,] 
smc$y = log10(as.integer(smc$y))
#plot
ggplot(data = smc, aes(x=x, y=y)) + geom_line(aes(colour=label)) + xlim(4.1, 5.5)


#smc for atfl combining all chromosomes calculation 
setwd("D:/maulana/third_project/smc/atfl")
files = list.files(pattern = "*.csv")
combined = data.frame()
for (i in files) {
  temp = read.table("plot_phased-28-atfl.csv", sep=",", header = TRUE, 
                    colClasses = c("NULL", "character", "character", "NULL", "NULL"))
  combined = rbind(combined, temp)
  rm(temp)
}
#sort rows by x column(years)
combined$x = as.integer(combined$x)
combined = combined[order(combined$x),] #2929
#discard duplicate columns
combined = distinct(combined)
combined$axis_x = log10(as.numeric(combined$x))
combined$axis_y = log10(as.numeric(combined$y))
combined$variable = as.factor("atfl")
#plot with log values of years and Ne
ggplot(data = combined, aes(x=axis_x, y=axis_y)) + geom_line(aes(colour=variable)) + xlim(3.4, 5.5)
#plot with real years and Ne doesn't work

#smoothing the Ne values before log
library(data.table)
combined$sm2 = log10(frollmean(as.numeric(combined$y), 2))
combined$variable = as.factor("atfl")
combined = combined[-1,]
#plot
ggplot(data = combined, aes(x=axis_x, y=sm2)) + geom_line(aes(colour=variable)) + xlim(3.4, 5.5)
ggplot(data = combined, aes(x=axis_x, y=sm2) + geom_line(aes(colour=variable)) + xlim(3.4, 5.5))
