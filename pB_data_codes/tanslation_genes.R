setwd('~/ferdig_rotation/pB_data/my_results/DTWMIC/wgcna/')
nf54 <- read.csv("wgcna_NF54_p6/clustering_nf54.csv", header=T)
pb58 <- read.csv("wgcna_PB58_p5/clustering_pb58.csv", header=T)

translation_genes <- read.table("old_wgcna_results/pb58_up_dha_down_translation_genes.txt", as.is=T)

genes_vec <- as.vector(translation_genes$V1)

nf_clust_trans = data.frame(genes_vec,rep(0,length(genes_vec)), rep(0,length(genes_vec)))
nf_clust_trans[,1] <- genes_vec

for(ii in 1:length(genes_vec)){
  nf_clust_trans[ii,2] <- nf54[which(nf54$row.names.test_nf==genes_vec[ii]),2]
}
for(ii in 1:length(genes_vec)){
  nf_clust_trans[ii,3] <- pb58[which(pb58$row.names.test_pb==genes_vec[ii]),2]
}

write.csv(nf_clust_trans, file="tranlation_genes_clusters.csv")


ordered_clusters_nf <- nf54[with(nf54, order(moduleLabels)),]
ordered_clusters_pb <- pb58[with(pb58, order(moduleLabels)),]

