
setwd("~/ferdig_rotation/pB_data/Data/")
smoothed_data <- read.csv("3d7_smoothed_katie.csv")

row.names(smoothed_data) <- smoothed_data[,1]
smoothed_data <- smoothed_data[,-1]
gene_index_max <- vector()

for(ii in 1:dim(smoothed_data)[1]){
  gene_index_max[ii] <- which.max(smoothed_data[ii,])
  
}

peaks_df <- data.frame(gene_index_max)
row.names(peaks_df) <- row.names(smoothed_data)

k13_nf <- read.csv('../my_results/DTWMIC/wgcna/wgcna_NF54_p6/K13_cluster_nf.csv')

k13_genes_peaks = vector()

for (ii in k13_nf$row.names.test_nf.) {
    k13_genes_peaks = c(k13_genes_peaks, peaks_df[which(row.names(peaks_df) == ii),1])
}

length(which(k13_genes_peaks < 4))
length(which(k13_genes_peaks > 40))
#out of the 430 genes in NF K13 cluster, 372 have epression profiles in 3D7, of those 158 have peak expression between 40 and 4


k13_pb <- read.csv('../my_results/DTWMIC/wgcna/wgcna_PB58_p5/K13_cluster_pb.csv')
k13_genes_peaks = vector()

for (ii in k13_pb$row.names.test_pb.) {
  k13_genes_peaks = c(k13_genes_peaks, peaks_df[which(row.names(peaks_df) == ii),1])
}

length(which(k13_genes_peaks < 4))
length(which(k13_genes_peaks > 40))

wr <- which(k13_genes_peaks < 29)
length(which(k13_genes_peaks[wr] > 12))
