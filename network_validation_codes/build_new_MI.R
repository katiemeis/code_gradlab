# making new MI network from ARACNE

setwd("~/ferdig_rotation/regulon_validation/new_curation/")
new_map_data <- read.csv("GSE19468_Final.csv")
row.names(new_map_data) <- new_map_data[,1]
dat <- new_map_data[,-1]

#library(nettools)
#need matrix with samples on rows and variables on columns
#dat2 <- t(dat)
#dat2[1:5,1:5] # check, samples by genes
#MI_res <- mat2adj(dat2, method ="ARACNE")


#building with empirical probability distributions 
library(minet)
#need dataframe with samples on rows and variables on columns
dat2 <- t(dat)
dat2[1:5,1:5] # check, samples by genes
data_df <- data.frame(dat2)

MI_net <- build.mim(data_df, estimator = "mi.empirical")

graph_MI <- graph.adjacency(MI_net,weighted=T, mode="upper", diag = F)
edgesMI <- get.data.frame(graph_MI)

edgesMI_ordered <- edgesMI[order(-edgesMI$weight),]
write.csv(edgesMI_ordered, "full_networks/MIempirical_allonall_EL.csv", quote = F, row.names = F, col.names = T)



#building with Miller-Madow aymptotic bias corrected emperical estimator
library(minet)
#need dataframe with samples on rows and variables on columns
dat2 <- t(dat)
dat2[1:5,1:5] # check, samples by genes
data_df <- data.frame(dat2)

MI_net_mm <- build.mim(data_df, estimator = "mi.mm")

graph_MI_mm <- graph.adjacency(MI_net_mm,weighted=T, mode="upper", diag = F)
edgesMI_mm <- get.data.frame(graph_MI_mm)

edgesMI_ordered_mm <- edgesMI_mm[order(-edgesMI_mm$weight),]
write.csv(edgesMI_ordered_mm, "full_networks/MIempirical_allonall_EL.csv", quote = F, row.names = F, col.names = T)

