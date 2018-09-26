
setwd("~/ferdig_rotation/pB_data/my_results/diff_gene_exp/")
data <- read.csv("triplicates_data.csv")

dim(data) # 18 by 5545, first 5 columns are metadata

#6 hr compare
pvals_6hr <- vector()
for(i in 6:5545){
  test_res <- t.test(data[1:3,i],data[10:12,i])
  pval <- test_res$p.value
  pvals_6hr <- c(pvals_6hr, pval)
  
}

#26 hr compare
pvals_26hr <- vector()
for(i in 6:5545){
  test_res <- t.test(data[4:6,i],data[13:15,i])
  pval <- test_res$p.value
  pvals_26hr <- c(pvals_26hr, pval)
  
}

#38 hr compare
pvals_38hr <- vector()
for(i in 6:5545){
  test_res <- t.test(data[7:9,i],data[16:18,i])
  pval <- test_res$p.value
  pvals_38hr <- c(pvals_38hr, pval)
  
}

#fdr correction
fdr_vals_6hr <- p.adjust(pvals_6hr, method = "fdr")
fdr_vals_26hr <- p.adjust(pvals_26hr, method = "fdr")
fdr_vals_38hr <- p.adjust(pvals_38hr, method = "fdr")

difgenes_6hr <- which(pvals_6hr < 0.01)
difgenes_26hr <- which(pvals_26hr < 0.01)
difgenes_38hr <- which(pvals_38hr < 0.01)

for(i in 1:22){
  difgenes_6hr[i] <- difgenes_6hr[i] + 5
}

for(i in 1:length(difgenes_26hr)){
  difgenes_26hr[i] <- difgenes_26hr[i] + 5
}

for(i in 1:length(difgenes_38hr)){
  difgenes_38hr[i] <- difgenes_38hr[i] + 5
}

dif_exp_6hr <- data[c(1:3,10:12),difgenes_6hr]
dif_exp_26hr <- data[c(4:6,13:15),difgenes_26hr]
dif_exp_38hr <- data[c(7:9,16:18),difgenes_38hr]


#find average and plot
avg_6_pb <- colMeans(dif_exp_6hr[1:3,])
avg_6_nf <- colMeans(dif_exp_6hr[4:6,])
diff_6 <- cbind(avg_6_nf, avg_6_pb) 

subtract <- vector()
for(i in 1:22){
  subtract[i] <- diff_6[i,1] - diff_6[i,2]
}

diff_6 <- cbind(diff_6, subtract)
ordered_diff_6 <- diff_6[order(-diff_6[,3]),]

mat1=data.matrix(ordered_diff_6[,-3])
heatmap.2(mat1, Rowv=FALSE, Colv=FALSE, scale="row", dendrogram = "none", distfun=dist, hclustfun=hclust, xlab="Line", ylab="Genes", key=TRUE, keysize=1, col=greenred(2000), trace="none", density.info=c("none"), margins=c(10, 8), cexRow=0.5, sepcolor="white")



avg_26_pb <- colMeans(dif_exp_26hr[1:3,])
avg_26_nf <- colMeans(dif_exp_26hr[4:6,])
diff_26 <- cbind(avg_26_nf, avg_26_pb) 

subtract <- vector()
for(i in 1:126){
  subtract[i] <- diff_26[i,1] - diff_26[i,2]
}

diff_26 <- cbind(diff_26, subtract)
ordered_diff_26 <- diff_26[order(-diff_26[,3]),]

mat2=data.matrix(ordered_diff_26[,-3])
heatmap.2(mat2, Rowv=FALSE, Colv=FALSE, scale = "row", dendrogram="none", distfun=dist, hclustfun=hclust, xlab="Line", ylab="Genes", key=TRUE, keysize=1, col=greenred(2000), trace="none", density.info=c("none"), margins=c(10, 8), cexRow=0.5, sepcolor="white")



avg_38_pb <- colMeans(dif_exp_38hr[1:3,])
avg_38_nf <- colMeans(dif_exp_38hr[4:6,])
diff_38 <- cbind(avg_38_nf, avg_38_pb) 

subtract <- vector()
for(i in 1:97){
  subtract[i] <- diff_38[i,1] - diff_38[i,2]
}

diff_38 <- cbind(diff_38, subtract)
ordered_diff_38 <- diff_38[order(-diff_38[,3]),]

mat3=data.matrix(ordered_diff_38[,-3])
heatmap.2(mat3, Rowv=FALSE, Colv=TRUE, scale = "row",dendrogram="none", distfun=dist, hclustfun=hclust, xlab="Line", ylab="Genes", key=TRUE, keysize=1, col=greenred(2000), trace="none", density.info=c("none"), margins=c(10, 8), cexRow=0.5, sepcolor="white")


#find names of diff exp genes
diff_6_names <- row.names(ordered_diff_6)
diff_26_names <- row.names(ordered_diff_26)
diff_38_names <- row.names(ordered_diff_38)




#plot all together
all_genes <- c(diff_6_names, diff_26_names, diff_38_names)

difgeneindex <- vector()
for(i in all_genes){
  difgeneindex <- c(difgeneindex, which(colnames(data) == i))
}

all_genes_data <- data[,difgeneindex]

PB58_6hr <- colMeans(all_genes_data[1:3,])
NF54_6hr <- colMeans(all_genes_data[10:12,])
PB58_26hr <- colMeans(all_genes_data[4:6,])
NF54_26hr <- colMeans(all_genes_data[13:15,])
PB58_38hr <- colMeans(all_genes_data[7:9,])
NF54_38hr <- colMeans(all_genes_data[16:18,])
diff_all <- rbind(PB58_6hr, NF54_6hr, PB58_26hr, NF54_26hr, PB58_38hr, NF54_38hr) 
diff_all_mat <- data.matrix(diff_all)
heatmap.2(diff_all_mat, Rowv=FALSE, Colv=FALSE, scale = "column",dendrogram="none", distfun=dist, hclustfun=hclust, key=TRUE, keysize=1, col=greenred(2000), trace="none", density.info=c("none"), margins=c(10, 8), cexRow=0.5, sepcolor="white")

diff_t <- data.matrix(t(diff_all))
heatmap.2(diff_t, Rowv=FALSE, Colv=FALSE, scale = "row",dendrogram="none", distfun=dist, hclustfun=hclust, key=TRUE, keysize=1, col=greenred(2000), trace="none", density.info=c("none"), margins=c(10, 8), cexRow=0.5, sepcolor="white")


#compare differentially expressed genes to network to see how much overlap
rw_net <- read.csv("../DTWMIC/wgcna/wgcna_hardthreshold/rewired_edgelist_PB58_KBS.csv")



