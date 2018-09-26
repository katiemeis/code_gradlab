library(igraph)

#new network
new <- read.csv(file= "~/network_science/project/adalasso_pcc_hu_noloops_newnames.csv", header = TRUE)

new <- t(new)
colnames(new)<- new[1,]
new <- new[-1,]
write.table(new, "/User/emilyherring/Desktop/new.txt")
new.num <- mapply(new, FUN=as.numeric)
new.num <- matrix(data=new.num, nrow=3705, ncol = 3705)
row.names(new.num)<- row.names(new)
colnames(new.num)<-colnames(new)
for(i in 1:3705){
  for(j in 1:3705){
    if(new.num[i,j] !=0) {new.num[i,j] <- 1}
  }
}

new.edge <- graph.adjacency(new.num)
new.cluster <- cluster_fast_greedy(as.undirected(new.edge))
#new.cluster.betweenness <- cluster_edge_betweenness(as.undirected(new.edge))
new.communities <- communities(new.cluster)


#reduced regulon
old <- read.table(file= "~/network_science/project/regulon_reduced_newnames.csv", header = TRUE,sep = ",")
old_mat <- graph.data.frame(old)
old_mat_1 <- get.adjacency(old_mat, sparse=F)

old.edge <- graph.adjacency(old_mat_1)
old.cluster <- cluster_fast_greedy(as.undirected(old.edge))
old.communities <- communities(old.cluster)
#old.cluster.betweenness <- cluster_edge_betweenness(old.edge)
#old.communities.betweenness <- communities(old.cluster.betweenness)
