
################ plotting distributions to compare GO and cluster #################
setwd("~/ferdig_rotation/regulon_validation/original_nets/consensus_network/clustering_output/")
library(igraph)


###################### GO distribution ##########################

###lets first load in the GO file and get the dist for number of genes with GO temrs with 1 gene, 2 genes, ... n genes
GO_data <- read.csv("../../GO_file/PID_6_10_NEW.csv", as.is=T)
GO_data[1:5,1:5]
row.names(GO_data) <- GO_data[,1]
GO_data[1:5,1:5]
GO_data1 <- GO_data[,-1]
GO_data1[1:5,1:5]

###now lets take only the rows (genes) that are also in the consensus network
#read in edgelist, convert to graph, get node names
cons_el <- read.csv("../consensus_edges.csv", as.is = T)
cons_graph <- graph_from_edgelist(as.matrix(cons_el[,c(1,2)]), directed=F)
cons_nodes <- V(cons_graph)$name
#find the rows of the GO matrix that intersect the consensus nodes
GO_cons <- GO_data1[row.names(GO_data1) %in% cons_nodes,]


###now we want to get the column sums - how many genes are in GO1, GO2, GO3 ...
go_sizes <- colSums(GO_cons)
#remove any zeros because we don't care about these terms - they aren't in the network
final_go_sizes <- go_sizes[go_sizes > 0]
#check what the mx is
max(final_go_sizes) #109
hist(final_go_sizes, breaks = seq(0,109,by=1), col="cadetblue2", xlab = "Number of genes (k)", ylab = "Number of GO terms with k genes", main="GO Term Size Distribution")
length(which(final_go_sizes == 1))
length(which(final_go_sizes == 2))


############# now lets look at some of our clustering results ##############
#lets do 1.9 (most GO terms) and 3.1 (highest LOO precision)
cluster_data <- read.csv("my_format/consnet_MCL_i1.9_KM.csv")
my_table <- table(cluster_data[,2])
tail(sort(my_table), 5)
max(my_table) #235
hist(my_table, breaks = seq(0,235,by=1), col="cadetblue2", xlab = "Number of genes (k)", ylab = "Number of clusters with k genes", main="Cluster Size Distribution (MCL i=1.9)")

#now for i=3.1
cluster_data <- read.csv("my_format/consnet_MCL_i3.1_KM.csv")
my_table <- table(cluster_data[,2])
max(my_table) #814
my_table2 <- sort(my_table)
tail(sort(my_table), 5)
hist(my_table2[1:373], breaks = seq(0,39,by=1), col="cadetblue2", xlab = "Number of genes (k)", ylab = "Number of clusters with k genes", main="Cluster Size Distribution (MCL i=3.1)")


