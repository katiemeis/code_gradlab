
####################################  MCL clustering on consensus network  ############################################

setwd("~/ferdig_rotation/regulon_validation/original_nets/consensus_network/")

#install.packages("MCL")
library(MCL)

#read in edgelist and turn into adjacency matrix
#I think we need symmetric matrix for it to be treated as undirected
edgelist <- read.csv("consensus_edges_duplicates.csv", header=T, as.is = T)

library(igraph)
con_graph <- graph.data.frame(edgelist)
summary(con_graph)
con_adj_mat <- as.matrix(as_adjacency_matrix(con_graph))

#mcl site says normal usage is to vary inflation between 1.2 and 5
#1.2 gives course clusters, 5 gives fine clusters
#https://micans.org/mcl/

#with default inflation 2
cluster_res <- mcl(x=con_adj_mat, addLoops = F)
node_cluster_IDs <- cluster_res$Cluster
length(unique(node_cluster_IDs)) #245 clusters
#get the cluster indexes
clust_index <- unique(node_cluster_IDs)
#get format of genes in 1st col, cluster index in 2nd col
names_vec <- row.names(con_adj_mat)
genes_cluster <- data.frame(names_vec, node_cluster_IDs)
genes_cluster$names_vec <- as.character(genes_cluster$names_vec)

#write.csv(genes_cluster, "MCL_inflation2_default/genes_to_cluster_map.csv", row.names = F, col.names = T, quote=F)

#sort by index smallest to largest
ordered_clust <- genes_cluster[order(node_cluster_IDs),]


#with inflation 4
#cluster_res_i4 <- mcl(x=con_adj_mat, inflation = 4, addLoops = F)
#node_cluster_map_i4 <- cluster_res_i4$Cluster
#length(unique(node_cluster_map_i4)) #356 clusters





###############################################  GO enrichment  ####################################################


#now lets do enrichment on each cluster
#bc3net package has a function called enrichment that will take the terms in each cluster
#compare to background using hypergeometric test
# https://www.rdocumentation.org/packages/bc3net/versions/1.0.4/topics/enrichment

#take genes in the cluster and get a set of GO terms
#what should the background be? Just the set or
library(bc3net)

#need to load in GO file and convert to named list
go_mat <- read.csv("../GO_file/PID_6_10_NEW.csv", as.is=T)
go_list <- vector("list", 699)
#for each column starting at 2, go through each row and add the gene if it maps to the GO term
for(i in 2:ncol(go_mat)){
  col_vec <- vector()
  for(j in 1:nrow(go_mat)){
   if(go_mat[j,i] == 1){
     col_vec <- c(col_vec, go_mat[j,1])
   } 
  }
  go_list[[i-1]] <- col_vec
}

names(go_list) <- colnames(go_mat)[2:ncol(go_mat)]

#we can now use go_list as our genesets variable in enrichment function
# reference should be the full set of nodes in the network, as a vector

#since we loaded the consensus network and made it a graph earlier we can use that
full_nodes <- V(con_graph)$name
#want an empty list we can put results in
enriched_list <- vector("list")

for(i in 1:length(clust_index)){
  current_cluster <- clust_index[i]
  #get list of genes in current cluster
  genes_in_c <- genes_cluster[which(genes_cluster[,2] == current_cluster),1]
  #run enrichment function
  res <- enrichment(genes = genes_in_c, reference = full_nodes, genesets = go_list, adj = "fdr")
  #get the GO terms that are enriched and put them in enriched list where the ith entry is for cluster i
  enriched_list[[i]] <- as.vector(res[which(res$padj < 0.05),1])
}

#now lets name the list so we know which cluster is which

#find the list entries that aren't empty
list_index <- vector()

for(j in 1:length(enriched_list)){
  if(length(enriched_list[[j]]) > 0){
    list_index <- c(list_index, j)
  }
}

#create a new list with only these entries
sublist <- enriched_list[list_index]

#now lets get a vector of which clusters have enriched terms
clusters_with_terms <- clust_index[list_index]
#make these the sublist list indexes
names(sublist) <- clusters_with_terms

GOs_vec <- vector()
for(k in 1:length(sublist)){
  GOs_vec <- c(GOs_vec, paste0(sublist[[k]], collapse = ","))
}

final_df <- data.frame(clusters_with_terms, GOs_vec)
#write.table(final_df, file = "MCL_inflation2_default/enriched_list.txt", sep = "\t", row.names = F, col.names = T, quote = F)







################################## EXTRAS ##############################################

## do i need to change genes_in_c and full_nodes to only include nodes with at least 1 GO term
#no, checked with the exapmle provided
data(exanet)
data(exgensets)
candidate=V(getgcc(exanet))$name
reference=V(exanet)$name
# if I get the full set of unique genes in the GO map object its 558
tun <- unique(unlist(exgensets, use.names = FALSE))
#candidate has 382 genes
#reference has 814 genes
#SO reference has genes that arent in GO map
#if we look at overlap between candidate and GO map, its 224 so not all candidate in GO map
length(intersect(candidate, tun))
#this should mean it automatically controls for these in the function
