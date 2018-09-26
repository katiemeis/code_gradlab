setwd("~/ferdig_rotation/pB_data/Data/")

lA_Data <- read.delim("lA_Data_t_katie.csv",sep=",",header=TRUE,as.is=TRUE)
lA_Data2 <- lA_Data[, -1]
rownames(lA_Data2) <- lA_Data[,1]
lA_Data <- lA_Data2
lA_Data <- t(lA_Data)

