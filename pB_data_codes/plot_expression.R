setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(nettools)
library(igraph)

####Load in expression data used to build each network

#data with names parsed
Final_Data <- read.delim("Final_Data_avg.csv",header=TRUE,sep=",",as.is=TRUE) 

#pull out rows for NF54 and order by hpi
NF54_data <- Final_Data[which(Final_Data$pBLine == "NF54"),]
NF54_data <- NF54_data[order(NF54_data$hpi),]
NF54_data2 <- NF54_data[,-1]
rownames(NF54_data2) <- NF54_data[,1]
NF54_data <- NF54_data2
#new matrx without metadata
n <- length(NF54_data)
NF54_data_exp <- NF54_data[, 5:n]
NF54_wgcna <- t(NF54_data_exp)


#pull out rows for PB-58 and order by hpi
PB58_data <- Final_Data[which(Final_Data$pBLine == "PB-58"),]
PB58_data <- PB58_data[order(PB58_data$hpi),]
PB58_data2 <- PB58_data[,-1]
rownames(PB58_data2) <- PB58_data[,1]
PB58_data <- PB58_data2
#new matrx without metadata
n <- length(PB58_data)
PB58_data_exp <- PB58_data[, 5:n] 
PB58_wgcna <- t(PB58_data_exp)

x <- c(2,6,10,14,18,22,26,30,34,38,42,46)

nfk13 <- NF54_data_exp$PF3D7_0305700
pbk13 <- PB58_data_exp$PF3D7_0305700

plot(x, pbk13, type = 'l', col="blue", main = "Gene 4 Expression", ylab="Expression", xlab="Time Point", lwd=3)
lines(x, nfk13, type="l", col = "darkgoldenrod3", lwd=3)

nfudp <- NF54_data_exp$PF3D7_1343600
pbudp <- PB58_data_exp$PF3D7_1340900

plot(x, pbudp, type = 'l', col="blue", main = "Gene 1 Expression", ylab="Expression", xlab="Time Point", lwd=3, ylim = c(0,5000))
lines(x, nfudp, type="l", col = "darkgoldenrod3", lwd=3)

nfunk <- NF54_data_exp$PF3D7_1343800
pbunk <- PB58_data_exp$PF3D7_1343800

plot(x, pbunk, type = 'l', col="blue", main = "Gene 5 Expression", ylab="Expression", xlab="Time Point", lwd=3, ylim = c(0,5000))
lines(x, nfunk, type="l", col = "darkgoldenrod3", lwd=3)
