setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(nettools)

####Data preperation

#data with names parsed
Final_Data <- read.delim("Final_Data_ARACNE_avg.csv",header=TRUE,sep=",",as.is=TRUE) 

#pull out rows for NF54 and order by hpi
NF54_data <- Final_Data[which(Final_Data$pBLine == "NF54"),]
NF54_data <- NF54_data[order(NF54_data$hpi),]
NF54_data2 <- NF54_data[,-1]
rownames(NF54_data2) <- NF54_data[,1]
NF54_data <- NF54_data2
#new matrx without metadata
n <- length(NF54_data)
NF54_data_exp <- NF54_data[, 5:n] 


#pull out rows for PB-58 and order by hpi
PB58_data <- Final_Data[which(Final_Data$pBLine == "PB-58"),]
PB58_data <- PB58_data[order(PB58_data$hpi),]
P58_data2 <- PB58_data[,-1]
rownames(PB58_data2) <- PB58_data[,1]
PB58_data <- PB58_data2
#new matrx without metadata
n <- length(PB58_data)
PB58_data_exp <- PB58_data[, 5:n] 



####running nettools DTWMIC
setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC")
PB35_nettools_run <- mat2adj(PB35_data_exp, method='DTWMIC')
NF54_nettools_run <- mat2adj(NF54_data_exp, method='DTWMIC')

#write.csv(PB35_nettools_run, "DTWMIC_PB35_results.csv")
#write.csv(NF54_nettools_run, "DTWMIC_NF54_results.csv")

rewire_p <- vector()
for(i in 1:5540){
  for(j in 1:5540){
    Fnf <- (1/2)*log((1+NF54_nettools_run[i,j])/(1-NF54_nettools_run[i,j]))
    Fpb <- (1/2)*log((1+PB35_nettools_run[i,j])/(1-PB35_nettools_run[i,j]))
    score <- abs(Fnf-Fpb)/sqrt((1/(12-3))+(1/(12-3)))
    p <- pnorm(score,lower.tail=FALSE)
    rewire_p <- c(rewire_p, p)
  }
}



