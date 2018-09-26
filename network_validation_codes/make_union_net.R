
setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")
library(igraph)

lasso <- read.csv("../final_nets/full_adalasso_renamed.csv", header=T, as.is=T)
lasso5000 <- lasso[1:5000,]

reg5000 <- read.csv("regulon_7500.csv", header=T, as.is=T)

pear5000 <- read.csv("pearson_7500.csv", header=T, as.is=T)

rf5000 <- read.csv("rf_7500.csv", header=T, as.is=T) 

combined_el <- rbind(lasso5000,reg5000,pear5000,rf5000)

combined_graph <- graph_from_edgelist(as.matrix(combined_el[,1:2]), directed = F)
simple_g <- simplify(combined_graph)
summary(combined_graph)
summary(simple_g)

union_el <- get.edgelist(simple_g)
union_el_mat <- as.matrix(union_el)
weights <- rep(1, 13775)

final_el <- cbind(union_el_mat, rep(1,13775))
write.csv(final_el, "union_network.csv", quote=F, row.names = F, col.names = F)
write.table(final_el, "union_network.txt", sep="\t", quote=F, row.names = F, col.names = F)


########################################################################################

#use the size of adaptive lasso as the cutoff - top 6439 edges of each

lasso <- read.csv("../final_nets/full_adalasso_renamed.csv", header=T, as.is=T)

reg <- read.csv("regulon_7500.csv", header=T, as.is=T)
regsize <- reg[1:6439,]

pear <- read.csv("pearson_7500.csv", header=T, as.is=T) 
pearsize <- pear[1:6439,]

rf <- read.csv("rf_7500.csv", header=T, as.is=T) 
rfsize <- rf[1:6439,]

combined_el <- rbind(lasso,regsize,pearsize,rfsize)

combined_graph <- graph_from_edgelist(as.matrix(combined_el[,1:2]), directed = F)
simple_g <- simplify(combined_graph)
summary(combined_graph)
summary(simple_g)

union_el <- get.edgelist(simple_g)
union_el_mat <- as.matrix(union_el)

final_el <- cbind(union_el_mat, rep(1,17661))
write.csv(final_el, "../union_network/union_network_alsize.csv", quote=F, row.names = F, col.names = F)
write.table(final_el, "../union_network/union_network_alsize.txt", sep="\t", quote=F, row.names = F, col.names = F)
