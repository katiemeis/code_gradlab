setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(nettools)

####Data preperation

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


#pull out rows for PB-58 and order by hpi
PB35_data <- Final_Data[which(Final_Data$pBLine == "PB-35"),]
PB35_data <- PB35_data[order(PB35_data$hpi),]
PB35_data2 <- PB35_data[,-1]
rownames(PB35_data2) <- PB35_data[,1]
PB35_data <- PB35_data2
#new matrx without metadata
n <- length(PB35_data)
PB35_data_exp <- PB35_data[, 5:n] 






####running nettools DTWMIC
#matrix or datafrme wth samples by genes
setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC")
PB35_nettools_run <- mat2adj(PB35_data_exp, method='DTWMIC')
NF54_nettools_run <- mat2adj(NF54_data_exp, method='DTWMIC')

write.csv(PB35_nettools_run, "DTWMIC_PB35_results.csv")
write.csv(NF54_nettools_run, "DTWMIC_NF54_results.csv")


#make adj matrix into edge list with value
library(reshape2)
PB35_adjmat <- read.csv("DTWMIC_PB58_results.csv")
PB35_new = melt(PB35_adjmat)
write.csv(PB35_new, "PB35_DTWMIC_edgelist.csv")

NF54_adjmat <- read.csv("DTWMIC_NF54_results.csv")
NF54_new = melt(NF54_adjmat)
write.csv(NF54_new, "NF54_DTWMIC_edgelist.csv")


#read in edgelist, make histogram, order values 
PB35_new = read.csv("PB35_DTWMIC_edgelist.csv")
PB35_new <- PB35_new[,-1]
hist(PB35_new[,3], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), xlab='MI', main = "PB58 MI Values")
ordered_PB35 <- PB35_new[with(PB35_new, order(-value)),]
PB35_top1 = ordered_PB35[1:306916,] ##top 1% of values, contains only values in bin farthest right on histogram
PB35_top5 = ordered_PB35[1:1534580,] ##top 5% of values
PB35_top10 = ordered_PB35[1:3069160,] ##top 10% of values
PB35_topbin = ordered_PB35[1:824858,] ##all values in top bin, 824858/30691600 ~ 2.69%


NF54_new = read.csv("NF54_DTWMIC_edgelist.csv")
NF54_new <- NF54_new[,-1]
hist(NF54_new[,3], breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), xlab='MI', main = "NF54 MI Values")
ordered_NF54 <- NF54_new[with(NF54_new, order(-value)),]
NF54_top1 = ordered_NF54[1:306916,] ##top 1% of values, contains only values in bin farthest right on histogram
NF54_top5 = ordered_NF54[1:1534580,]  ##top 5% of values
NF54_top10 = ordered_NF54[1:3069160,] ##top 10% of values
NF54_topbin = ordered_NF54[1:821840,] ##all values in top bin, 821840/30691600 ~ 2.68%




####aracne with nettools
setwd("~/ferdig_rotation/pB_data/my_results/aracne")
NF54_aracne <- mat2adj(NF54_data_exp, method = "ARACNE")
NF54_edges <- melt(NF54_aracne)
ordered_NF54_edges <- NF54_edges[with(NF54_edges, order(-value)),]
final_edges_NF54 <- ordered_NF54_edges[1:34874,]      ###only the nonzero connections
write.csv(final_edges_NF54, "nonzero_NF54_edges.csv")

#histogran from all MIs (including zeros), cant see the distribution 
NF54_edges <- read.csv("NF54_edges_all.csv")
NF54_edges <- NF54_edges <- NF54_edges[,-1]
hist(NF54_edges[,3], xlab='MI', main = "NF54 MI Values")

#histogram from only nonzero MIs
NF54_nonzero_edges <- read.csv("nonzero_NF54_edges.csv")
NF54_nonzero_edges <- NF54_nonzero_edges[,-1]
hist(NF54_nonzero_edges[,3], xlab='MI', main = "NF54 MI Values")

#top % of edges from nonzero
NF54_top1 <- NF54_nonzero_edges[1:349,]
NF54_top5 <- NF54_nonzero_edges[1:1744,]
NF54_top10 <-NF54_nonzero_edges[1:3487,] 


PB35_aracne <- mat2adj(PB35_data_exp, method = "ARACNE")
PB35_edges <- melt(PB35_aracne)
ordered_PB35_edges <- PB35_edges[with(PB35_edges, order(-value)),]
final_edges_PB35 <- ordered_PB35_edges[1:34864,]      ###only the nonseor connections
write.csv(final_edges_PB35, "nonzero_PB35_edges.csv")

#histogran from all MIs (including zeros), cant see the distribution 
PB35_edges <- read.csv("PB35_edges_all.csv")
PB35_edges <- PB35_edges <- PB35_edges[,-1]
hist(PB35_edges[,3], xlab='MI', main = "NF35 MI Values")

#histogram from only nonzero MIs
PB35_nonzero_edges <- read.csv("nonzero_PB35_edges.csv")
PB35_nonzero_edges <- PB35_nonzero_edges[,-1]
hist(PB35_nonzero_edges[,3], xlab='MI', main = "PB58 MI Values")

#top % of edges from nonzero
PB535_top1 <- PB35_nonzero_edges[1:349,]
PB535_top5 <- PB35_nonzero_edges[1:1743,]
PB535_top10 <-PB35_nonzero_edges[1:3486,] 


####running regular aracne
library(minet)
MIM = build.mim(NF54_data_exp, estimator = "spearman", disc = "none")
NF54_aracne_results = aracne(MIM, eps= 0)
library(Rgraphviz)
plot( as(aracne_results ,"graphNEL") )
