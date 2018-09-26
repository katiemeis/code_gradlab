
setwd("~/ferdig_rotation/regulon_validation/original_nets/final_nets/")

#install.packages("MCL")
library(MCL)

#read in edgelist and turn into adjacency matrix
#I think we need symmetric matrix for it to be treated as undirected
edgelist <- read.csv("full_adalasso_renamed.csv", header=T, as.is = T)
#add in the duplicate edges
edges2 <- data.frame(edgelist$to, edgelist$from, edgelist$Weight)
colnames(edges2) <- c("from", "to", "Weight")
full_edges <- rbind(edgelist, edges2)


library(igraph)
lasso_graph <- graph.data.frame(full_edges)
summary(lasso_graph)
lasso_adj_mat <- as.matrix(as_adjacency_matrix(lasso_graph))



#mcl site says normal usage is to vary inflation between 1.2 and 5
#1.2 gives course clusters, 5 gives fine clusters
#https://micans.org/mcl/

#Vary inflation from 1.2 to 3 by 0.2, then 3 to 5 by 1

set.seed(0)
cluster_res <- mcl(x=lasso_adj_mat, addLoops = F, inflation = 1.4, max.iter =300)
node_cluster_IDs <- cluster_res$Cluster
length(unique(node_cluster_IDs)) #245 clusters
#get the cluster indexes
clust_index <- unique(node_cluster_IDs)
#get format of genes in 1st col, cluster index in 2nd col
names_vec <- row.names(con_adj_mat)
genes_cluster <- data.frame(names_vec, node_cluster_IDs)
genes_cluster$names_vec <- as.character(genes_cluster$names_vec)

#order the dataframe so cluster 0 is at top, decreasing in order
gc_ord <- genes_cluster[order(genes_cluster$node_cluster_IDs),] 

#reformat to have clusterID followed by list of genes in that cluster
gc_reform <- aggregate( .~ node_cluster_IDs, gc_ord, function(x) toString(unique(x)))



#write out only the second column - should match Qi's format of 1 cluster per line
write.table(gc_reform$names_vec, "../adalasso_clustering/lasso_MCL_i1.2_QL.csv", quote = F, row.names = F, col.names = F, sep = "")

#will need to use sed "s/, /\t/g" file.csv > file.txt to make tab delimited

#my format, used in the GO code (MCL_clustering.R which I might rename to GO_analysis or something)
write.csv(genes_cluster, "../adalasso_clustering/lasso_MCL_i1.2_KM.csv", row.names = F, col.names = T, quote=F)
