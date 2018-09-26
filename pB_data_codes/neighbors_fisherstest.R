
setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/wgcna/compare_to_regulon/")

nfk13 <- read.csv("NF_ht_K13_neighbors.txt", header = FALSE)
pbk13 <- read.csv("PB_ht_K13_neighbors.txt", header = FALSE)
regk13 <- read.csv("regulon_k13_neighbors.txt", header = FALSE)

intersect(regk13$V1,pbk13$V1)
intersect(regk13$V1,nfk13$V1)

fisher.test(matrix(c(5173, 154, 258, 6), nrow=2), alternative="greater")
fisher.test(matrix(c(5344, 152, 87, 8), nrow=2), alternative="greater")