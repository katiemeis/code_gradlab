###calculate modularity for networks
#takes in edgelist (with weight)
#also need a GO term mapping (matrix gene by GO term)

setwd("~/ferdig_rotation/regulon_validation/series_thresholds/")

#need to only look at genes associated with GO terms that have 10 or more genes
gocounts <- read.csv("PID_6_10.csv")
row.names(gocounts) <- gocounts[,1]
gocounts <- gocounts[,-1]

#count how many genes are associated with each term
gosums <- colSums(gocounts)
gosumsvec <- as.integer(gosums)

#only keep go terms with at least ten genes
go_ten_index <- which(gosumsvec > 9)
go_ten_table <- gocounts[,go_ten_index]

#delete genes that are all zeros so we only keep nodes associated with an ID n > 9
go_final_ten <- go_ten_table[rowSums(go_ten_table) > 0,]

#get the list of nodes included in GOs w/ at least 10 genes in them
#these are the only genes we want to consider
gonodes <- row.names(go_final_ten)

#read in the graph
graph1 <- read.csv("rf_150000.csv", header=T, as.is=T)

#find rows of graph1 are edges between GO terms with at least 10 genes
#only keep edges between gonodes
ind_vec <- which(graph1$V1 %in% gonodes & graph1$V2 %in% gonodes)
final_edges <- graph1[ind_vec,]
node_sums <- 


#modulatiy
# Q = (1/2m) sumc sumij [Aij - (kjki)/2m](1/OiOj)




################# extra tests 
ind_vec <- vector()
for(i in 1:dim(graph1)[1]){
  r <- graph1$V1[i]
  s <- graph1$V2[i]
  if(graph1$V1[i] %in% gonodes & graph1$V2[i] %in% gonodes){
    ind_vec <- c(ind_vec, i)
  }
  
}
library(igraph)
smat <- as.matrix(get.adjacency((graph.data.frame(graph1))))
snode <- row.names(smat)

intersect(gonodes, snode)
s <- "PF3D7_0101500"
r <- "PF3D7_0101500"
if(s %in% gonodes & r %in% gonodes){print(1)}
