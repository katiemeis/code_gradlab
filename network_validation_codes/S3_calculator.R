library(igraph)
setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")

############################### edge overlap for networks ##########################
#enter adjacency list file names here
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

#first value is total set of nodes in both graphs
#second value is number of same edges

# divide intersection by 2 because we doubled the edges


##################### plot degree distribution and get r-squared ####################

#from https://chengjunwang.com/web_data_analysis/demo2_simulate_networks/


fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  # sequence from 1 to max degree
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # lm prob follows degree (log both)
  reg = lm(log(probability) ~ log(degree))
  # grab the coefficients
  cozf = coef(reg)
  # fit
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  #alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  #print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot data
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Frequency (log)", 
       col = 1, main = "Degree Distribution")
  # add fit line
  curve(power.law.fit, col = "red", add = T, n = length(d))
}


A <- read.csv("consensus_edges.csv", as.is = T)
colnames(A) <- c("V1", "V2", "V3")
sym_edgesA <- data.frame(A$V2, A$V1, A$V3)
A2 <- as.data.frame(rbind(as.matrix(A), as.matrix(sym_edgesA)))

B <- read.csv("regulon_150000.csv", as.is = T)
colnames(B) <- c("V1", "V2", "V3")
sym_edgesB <- data.frame(B$V2, B$V1, B$V3)
B2 <- as.data.frame(rbind(as.matrix(B), as.matrix(sym_edgesB)))

C <- read.csv("rf_150000.csv", as.is = T)
colnames(C) <- c("V1", "V2", "V3")
sym_edgesC <- data.frame(C$V2, C$V1, C$V3)
C2 <- as.data.frame(rbind(as.matrix(C), as.matrix(sym_edgesC)))

Amatrix <- graph.data.frame(A2)
Bmatrix <- graph.data.frame(B2)
Cmatrix <- graph.data.frame(C2)

fit_power_law(Amatrix)
fit_power_law(Bmatrix)
fit_power_law(Cmatrix)



################### edge overlap for networks w/ PPI #########################
setwd("~/ferdig_rotation/regulon_validation/series_thresholds/PPI/")

A <- read.csv("../rf_183932.csv", as.is = T)
colnames(A) <- c("V1", "V2", "V3")
sym_edgesA <- data.frame(A$V2, A$V1, A$V3)
A2 <- as.data.frame(rbind(as.matrix(A), as.matrix(sym_edgesA)))

PPI <- read.table("nat_2_nn_nodup.txt", header=T, as.is=T, sep = " ")
colnames(PPI) <- c("V1", "V2", "V3")
sym_edgesPPI <- data.frame(PPI$V2, PPI$V1, PPI$V3)
PPI2 <- as.data.frame(rbind(as.matrix(PPI), as.matrix(sym_edgesPPI)))


Amatrix <- graph.data.frame(A2)
PPImatrix <- graph.data.frame(PPI2)
ABintersect <- intersection(Amatrix, PPImatrix)
summary(ABintersect)


################# same as regulon rewire analysis ####################
setwd("~/ferdig_rotation/regulon_validation/series_thresholds/PPI/")

hc_genes <- read.csv("HC_subnet_newnames.txt", as.is=T, header=F)
#read in network and find AP2 neighbors
net <- read.csv("../rf_183932.csv", as.is = T)
#net <- read.csv("../../sizematch_nodups_nets/fullreg_size_names_order_nodups.csv", as.is=T, header=T)
colnames(net) <- c("V1", "V2", "V3")
neighbors1 <- net[which(net$V1 == "PF3D7_1007700"),2]
neighbors2 <- net[which(net$V2 == "PF3D7_1007700"),1]
nei_total <- c(neighbors1, neighbors2)
intersect(nei_total, hc_genes$V1)

