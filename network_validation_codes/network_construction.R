
#### Code for constructing the pcc and pearson networks with the Hu data #### random forest done in python


############################### lasso/pcc network construction ###########################
setwd('~/network_science/project/')
#install.packages('parcor')
library(parcor)
library(igraph)

hu_data <- read.csv('Hu_data.csv', header=TRUE)
hu_data_t <- t(hu_data)
colnames(hu_data_t) = hu_data_t[1,]
hu_data_t = hu_data_t[-1,]
class(hu_data_t) <- "numeric"

hu_pcc_mat <- adalasso.net(hu_data_t)

#regular lasso
lasso.hu <- hu_pcc_mat$pcor.lasso
for(ii in 1:3705){
  lasso.hu[ii,ii] = 0
}
colnames(lasso.hu) = colnames(hu_data_t)
rownames(lasso.hu) = colnames(hu_data_t)
write.csv(lasso.hu, 'lasso_pcc_hu_noloops.csv')
edge_list_lasso <- graph.adjacency(lasso.hu, weighted=TRUE)
lasso_df_edges <- get.data.frame(edge_list_lasso)
write.csv(lasso_df_edges, 'lasso_hu_edgelist_noloops.csv')
write.table(lasso_df_edges, 'lasso_hu_edgelist_noloops.txt')

#adaptive lasso
adalasso.hu <- hu_pcc_mat$pcor.adalasso
for(ii in 1:3705){
  adalasso.hu[ii,ii] = 0
}
colnames(adalasso.hu) = colnames(hu_data_t)
rownames(adalasso.hu) = colnames(hu_data_t)
write.csv(adalasso.hu, 'adalasso_pcc_hu_noloops.csv')
edge_list_adalasso <- graph.adjacency(adalasso.hu, weighted=TRUE)
adalasso_df_edges <- get.data.frame(edge_list_adalasso)
write.csv(adalasso_df_edges, 'adalasso_hu_edgelist_noloops.csv')
write.table(adalasso_df_edges, 'adalasso_hu_edgelist_noloops.txt')

save(list=ls(), file='pcc_hu_workspace.RData')




##################### Pearson correlation network construction ############################
hu_data <- read.csv('Hu_data.csv', header=TRUE)
hu_data_t <- t(hu_data)
colnames(hu_data_t) = hu_data_t[1,]
hu_data_t = hu_data_t[-1,]
class(hu_data_t) <- "numeric"

library(nettools)
hu_pearson_result <- mat2adj(hu_data_t, method="cor")
write.csv(hu_pearson_result, file="~/ferdig_rotation/regulon_validation/hu_pearson_network_mat.csv", quote=FALSE)

pear_net <- read.csv("../hu_pearson_network_mat.csv", as.is = T)
row.names(pear_net) <- pear_net$X
pear_net2 <- pear_net[,-1]
pear_graph <- graph.adjacency(as.matrix(pear_net2), weighted=TRUE)
pear_df_edges <- get.data.frame(pear_graph)
pear_edge_ord <- pear_df_edges[order(-pear_df_edges$weight),]


