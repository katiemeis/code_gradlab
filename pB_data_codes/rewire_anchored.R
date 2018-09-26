
setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/")

#load in DTWMIC matrices
nf_data <- read.csv("DTWMIC_NF54_results.csv", header=TRUE)
row.names(nf_data) <- nf_data[,1]
nf_data <- nf_data[,-1]

pb_data <- read.csv("DTWMIC_PB58_results.csv", header=TRUE)
row.names(pb_data) <- pb_data[,1]
pb_data <- pb_data[,-1]

#pull out the K13 row in each data set
nf_k13 <- nf_data[which(row.names(nf_data) == "PF3D7_1343700"),]
pb_k13 <- pb_data[which(row.names(pb_data) == "PF3D7_1343700"),]

#Fisher transformation then calculate p-value based on standard normal distribution
rewire_p <- vector()
for(i in 1:5540){
  Fnf <- (1/2)*log((1+nf_k13[1,i])/(1-nf_k13[1,i]))
  Fpb <- (1/2)*log((1+pb_k13[1,i])/(1-pb_k13[1,i]))
  score <- abs(Fnf-Fpb)/sqrt((1/(12-3))+(1/(12-3)))
  p <- pnorm(score,lower.tail=FALSE)
  rewire_p <- c(rewire_p, p)
}

#get p-values
fdr_pvals <- p.adjust(rewire_p, method="fdr")
#pull out anything less than 0.05
sig_vect <- which(rewire_p < 0.05)
row.names(nf_data)[sig_vect]



#for full rewire, would need to loop over i and j
#put p-vals into edgelist w/ i,j,pval
#add another column with fdrs
#threshold
#convert thresholded edgelist to graph w/ igraph
#get degree w/ igraph     degree(graph, v=V(graph))






#rewiring of significantly rewired edges to K13
#pull out the K13 row in each data set
#k13_rewired_nb <- read.csv("./wgcna/wgcna_hardthreshold/K13_neighbors.txt", header = F, as.is=T)

NF_k13_rewire <- vector()
for(i in sig_vect){
  NF_k13_rewire <- c(NF_k13_rewire, nf_k13[1,i])
}

PB_k13_rewire <- vector()
for(i in sig_vect){
  PB_k13_rewire <- c(PB_k13_rewire, pb_k13[1,i])
}

rewired_k13_df <- cbind(NF_k13_rewire, PB_k13_rewire)
row.names(rewired_k13_df) <- row.names(nf_data)[sig_vect]
#all high in NF and low in PB except PF3D7_1017300.2
#golgi re-assembly stacking protein 1, golgi re-assembly stacking protein 2


ubi_neighbors <- read.csv("./wgcna/wgcna_hardthreshold/ubiquitin_neighbors.txt", header=F, as.is = T)

index_vec <- vector()
for(i in ubi_neighbors$V1){
  index <- which(row.names(nf_data) == i)
  index_vec <- c(index_vec, index)
}


nf_ubi_vals <- nf_data[which(row.names(nf_data) == "PF3D7_0305700"), index_vec]
row.names(nf_ubi_vals) <- "NF"
pb_ubi_vals <- pb_data[which(row.names(pb_data) == "PF3D7_0305700"), index_vec]
row.names(pb_ubi_vals) <- "PB"
rewired_ubi_df <- t(rbind(nf_ubi_vals, pb_ubi_vals))
#all low in NF54 and high in PB58 


