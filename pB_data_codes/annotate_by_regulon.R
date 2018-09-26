
setwd("~/ferdig_rotation/regulon_validation/")
library(igraph)

regulon_data <- read.table("full_regulon_newnames.txt", header=T)

setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/wgcna/wgcna_hardthreshold/")

reg_g <- graph.data.frame(regulon_data, directed = F)
reg_mat <- as_adjacency_matrix(reg_g, type = "both", sparse = F, names = T, attr="MutualInformation")


top_2 <- which(regulon_data$Regulator == "PF3D7_1308700")
top_2_data <- regulon_data[top_2,]
which(top_2_data$MutualInformation > .5)
write.csv(top_2_data$Target[1:59], "PF3D7_1308700_reg_neighbors.csv")


top_7 <- which(regulon_data$Regulator == "PF3D7_0729100")
top_7_data <- regulon_data[top_7,]
which(top_7_data$MutualInformation > .5)
write.csv(top_7_data$Target[1:64], "PF3D7_0729100_reg_neighbors.csv")
