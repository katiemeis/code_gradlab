setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(parcor)

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

NF54_pcc_mat <- adalasso.net(NF54_data_exp)

adalasso_NF54 <- NF54_pcc_mat$pcor.adalasso
for(ii in 1:3705){
  adalasso_NF54[ii,ii] = 0
}
colnames(adalasso_NF54) = colnames(NF54_data_exp)
rownames(adalasso_NF54) = colnames(NF54_data_exp)
write.csv(adalasso_NF54, 'adalasso_pcc_NF54_mat.csv')

edge_list_adalasso <- graph.adjacency(adalasso_NF54, weighted=TRUE)
adalasso_df_edges <- get.data.frame(edge_list_adalasso)
write.csv(adalasso_df_edges, 'adalasso_NF54_edgelist.csv')


save(list=ls(), file='pcc_hu_workspace.RData')


#pull out rows for PB-58 and order by hpi
PB58_data <- Final_Data[which(Final_Data$pBLine == "PB-58"),]
PB58_data <- PB58_data[order(PB58_data$hpi),]
PB58_data2 <- PB58_data[,-1]
rownames(PB58_data2) <- PB58_data[,1]
PB58_data <- PB58_data2
#new matrx without metadata
n <- length(PB58_data)
PB58_data_exp <- PB58_data[, 5:n] 

PB58_pcc_mat <- adalasso.net(PB58_data_exp)

adalasso_PB58 <- PB58_pcc_mat$pcor.adalasso
for(ii in 1:3705){
  adalasso_PB58[ii,ii] = 0
}
colnames(adalasso_PB58) = colnames(PB58_data_exp)
rownames(adalasso_PB58) = colnames(PB58_data_exp)
write.csv(adalasso_PB58, 'adalasso_pcc_PB58_mat.csv')

edge_list_adalasso <- graph.adjacency(adalasso_PB58, weighted=TRUE)
adalasso_df_edges <- get.data.frame(edge_list_adalasso)
write.csv(adalasso_df_edges, 'adalasso_PB58_edgelist.csv')

