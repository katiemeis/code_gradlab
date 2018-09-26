library(igraph)
#load new data - Gabe's curation
setwd("~/ferdig_rotation/regulon_validation/new_curation/")
new_map_data <- read.csv("GSE19468_Final.csv")
row.names(new_map_data) <- new_map_data[,1]
dat <- new_map_data[,-1]
dat2 <- t(dat)
#build pearson network
pear_mat <- cor(dat2, use = "pairwise.complete.obs")
#get edge list with no self loops and no duplicates
g2 <- graph.adjacency(pear_mat,weighted=T, mode="upper", diag = F)
edges2 <- get.data.frame(g2)


edges2_ordered <- edges2[order(-edges2$weight),]
write.csv(edges2_ordered, "Pearson_allonall_EL.csv", quote = F, row.names = F, col.names = T)


#need the absolute value edges also 
edges2$weight <- abs(edges2$weight)
abs_ordered <- edges2[order(-edges2$weight),]
write.csv(abs_ordered, "Pearson_absval_allonall_EL.csv", quote = F, row.names = F, col.names = T)
