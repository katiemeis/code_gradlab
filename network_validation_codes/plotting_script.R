setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")

#node overlap plot
overlap_data <- read.table("../edge_node_overlaps.txt", sep=" ", header=T)
edges_percent <- edges_in_common/threshs
plot(overlap_data[,1], overlap_data[,2], xlab="Number of edges", ylab="# node overlap")
lines(overlap_data[,1], overlap_data[,2])


plot(overlap_data[,1], overlap_data[,3], xlab="Number of edges", ylab="# edge overlap")
lines(overlap_data[,1], overlap_data[,3])
percent_edge <- overlap_data[,3]/overlap_data[,1]
plot(overlap_data[,1], percent_edge, xlab="Number of edges", ylab="% edge overlap", ylim=c(0.26,0.41))
lines(overlap_data[,1], percent_edge)

#scale-free plot
sf_fit <- read.table("../scalefree_fit.txt", sep = " ", header=T)
plot(sf_fit$Threshold, sf_fit$Pearson, xlab="Number of edges", ylab="Power law R2", ylim = c(0,1), type="o", lwd=2)
lines(sf_fit$Threshold, sf_fit$MI, col="red", lwd=2, type="o")
lines(sf_fit$Threshold, sf_fit$RF, col="blue", lwd = 2, type="o")
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#modularity
threshs <- c(5000,7500,10000,25000,50000,75000,100000,150000,183932)
pear_mod <- c(0.09717827,0.09270397,0.09225182,0.06807143,0.04988604,0.0403864,0.03457516,0.02697402,0.0236206)
mi_mod <- c(0.1066846,0.1073326,0.1022454,0.06854573,0.04733121,0.03754072,0.03195851,0.02623645,0.02359826)
rf_mod <- c(0.09804214,0.09440995,0.09015398,0.06815542,0.05280602,0.04237684,0.03747424,0.0301244,0.02716064)
plot(threshs, pear_mod, xlab="Number of edges", ylab="Modularity", ylim = c(0,.12), type="o", lwd=2)
lines(threshs, mi_mod, col="red", lwd=2, type="o")
lines(threshs, rf_mod, col="blue", lwd = 2, type="o")
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#fscore
threshs <- c(5000,7500,10000,25000,50000,75000,100000,150000,183932)
pear_fscore <- c(0.00785829,0.010871351,0.014067421,0.022586355,0.026791221,0.027911355,0.027737796,0.026732008,0.026093966)
mi_fscore <- c(0.008902581,0.013045622,0.016542699,0.023952716,0.024982115,0.024105261,0.023169218,0.023176163,0.027729731)
rf_fscore <- c(0.009163653,0.012691671,0.015366329,0.023642179,0.02840297,0.029180053,0.029870576,0.029304128,0.028731539)
plot(threshs, pear_fscore, xlab="Number of edges", ylab="F-score", ylim = c(0,.035), type="o", lwd=2)
lines(threshs, mi_fscore, col="red", lwd=2, type="o")
lines(threshs, rf_fscore, col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#precision
resu_prec <- read.csv("results_from_wrong_betaCVcalc/precision_results.csv", header=T)
plot(resu_prec$Edges, resu_prec$Pearson, xlab="Number of edges", ylab="Precision", ylim = c(0,.07), type="o", lwd=2)
lines(resu_prec$Edges,resu_prec$MI, col="red", type="o", lwd=2)
lines(resu_prec$Edges,resu_prec$RF, col="blue", type="o", lwd=2)
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#recall
resu_recall <- read.csv("results_from_wrong_betaCVcalc/recall_results.csv", header=T)
plot(resu_recall$Number.of.edges, resu_recall$Pearson, xlab="Number of edges", ylab="Recall", ylim = c(0,.07), type="o", lwd=2)
lines(resu_recall$Number.of.edges,resu_recall$MI, col="red", type="o", lwd=2)
lines(resu_recall$Number.of.edges,resu_recall$RF, col="blue", type="o", lwd=2)
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#precision vs recall
plot(resu_recall$Pearson, resu_prec$Pearson, xlab="Recall", ylab="Precision", ylim=c(0,.08), xlim=c(0,0.08), type="o", lwd=2)
lines(resu_recall$MI, resu_prec$MI, col="red", type="o", lwd=2)
lines(resu_recall$RF, resu_prec$RF, col="blue", type="o", lwd=2)
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#number of nodes
num_nodes <- read.table("../node_counts.txt", sep=" ", header=T)
plot(num_nodes$Threshold, num_nodes$Pearson, xlab="Number of edges", ylab="Number of nodes", ylim=c(0,3700), type="o", lwd=2)
lines(num_nodes$Threshold, num_nodes$MI, col="red", type="o", lwd=2)
lines(num_nodes$Threshold, num_nodes$RF, col="blue", type="o", lwd=2)
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)



#GO histogram
#highest 259
go_counts <- read.csv("../new_go_Qi/newPIDgenecounts.csv", header=F)
hist(as.numeric(go_counts[1,]), breaks=258, xlab="Number of genes", ylab="Number of GO terms")

go_counts_num <- as.numeric(go_counts[1,])
length(which(go_counts_num > 1))

gt_edge_count = 0
for(i in 1:length(go_counts)){
  gt_edge_count = gt_edge_count + choose(go_counts_num[i],2)
  
}











data <- read.delim("results_for_rplot.txt", header = F)

plot(data[1:9,1], data[1:9,2], ylim = c(0,6), xlab = "Number of edges", ylab = "BetaCV ratio", main = "BetaCV Ratio vs Number of Edges", lwd =2)
lines(data[1:9,1], data[1:9,2], lwd = 2)
points(data[10:19,1], data[10:19,2], col = "red", lwd=2)
lines(data[10:19,1], data[10:19,2], col = "red", lwd=2)
points(data[20:28,1], data[20:28,2], col = "blue", lwd =2)
lines(data[20:28,1], data[20:28,2], col = "blue", lwd=2)
legend("topright", c("Random Forest", "Mutual Information", "Pearson Correlation"), col=c(2, 1, "blue"), pch=c(1,1), lwd=2)

plot(data[1:8,1], data[1:8,2], ylim = c(0,7), xlab = "Number of edges", ylab = "BetaCV ratio", main = "BetaCV Ratio vs Number of Edges", lwd =2)
lines(data[1:8,1], data[1:8,2], lwd = 2)
points(data[9:16,1], data[9:16,2], col = "red", lwd=2)
lines(data[9:16,1], data[9:16,2], col = "red", lwd=2)
points(data[17:24,1], data[17:24,2], col = "blue", lwd =2)
lines(data[17:24,1], data[17:24,2], col = "blue", lwd=2)
legend("topright", c("Random Forest", "Mutual Information", "Pearson Correlation"), col=c(2, 1, "blue"), pch=c(1,1), lwd=2)

