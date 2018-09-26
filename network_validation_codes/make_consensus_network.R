library(igraph)
setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")

############################### edge overlap for networks ##########################
#enter edgelist file names here
A <- read.csv("pearson_75000.csv", as.is = T)
colnames(A) <- c("V1", "V2", "V3")
sym_edgesA <- data.frame(A$V2, A$V1, A$V3)
A2 <- as.data.frame(rbind(as.matrix(A), as.matrix(sym_edgesA)))

B <- read.csv("regulon_75000.csv", as.is = T)
colnames(B) <- c("V1", "V2", "V3")
sym_edgesB <- data.frame(B$V2, B$V1, B$V3)
B2 <- as.data.frame(rbind(as.matrix(B), as.matrix(sym_edgesB)))

C <- read.csv("rf_75000.csv", as.is = T)
colnames(C) <- c("V1", "V2", "V3")
sym_edgesC <- data.frame(C$V2, C$V1, C$V3)
C2 <- as.data.frame(rbind(as.matrix(C), as.matrix(sym_edgesC)))


Amatrix <- graph.data.frame(A2)
Bmatrix <- graph.data.frame(B2)
Cmatrix <- graph.data.frame(C2)

ABintersect <- intersection(Amatrix, Bmatrix)
BCintersect <- intersection(Bmatrix, Cmatrix)
ACintersect <- intersection(Amatrix, Cmatrix)
ABCintersect <- intersection(Amatrix, Bmatrix, Cmatrix)

summary(Amatrix)
summary(Bmatrix)
summary(Cmatrix)
#summary(ABintersect)
#summary(BCintersect)
#summary(ACintersect)
summary(ABCintersect)

#intersection of 75,000 edge networks is consensus network
cons_dups <- get.data.frame(ABCintersect)
write.csv(cons_dups[,c(1,2)], "consensus_edges_duplicates.csv", row.names = F, col.names = T, quote = F)
noselfgraph <- simplify(ABCintersect, remove.loops = T, remove.multiple = T) 
#take out duplicate edges
nodupgraph <- as.undirected(noselfgraph, mode=c("collapse"), edge.attr.comb = "mean")
#convert back to edge list from graph
cons_edges <- get.data.frame(nodupgraph)
write.csv(cons_edges[,c(1,2)], "consensus_edges.csv", row.names = F, col.names = T, quote = F)


#checked this in cytoscape to see if my removing on self loops and duplicates worked
#there were no self loops, and half the edges were duplicates in both my cdoe and cytoscape


#first value is total set of nodes in both graphs
#second value is number of same edges

# divide intersection by 2 because we doubled the edges
