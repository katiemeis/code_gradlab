
####### precision/recall between two edge lists
setwd("~/ferdig_rotation/regulon_validation/original_nets/gold_standards/")
library(igraph)

#first load in the gold standard/ground truth
string <- read.csv("string_network/string_renamed_noloops.csv")
plasmoMap <- read.csv("plasmoMAP/plasmoMAP_renamed_noloops.csv")

string_graph <- graph.data.frame(string, directed = F)
pmap_graph <- graph.data.frame(plasmoMap, directed = F)

#lets loop over the input networks and get results for both the gold standards
#compared to each input network

setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")
#my_nets <- c("pearson_5000.csv", "pearson_7500.csv", "pearson_10000.csv", "pearson_25000.csv",
#               "pearson_50000.csv", "pearson_75000.csv", "pearson_100000.csv", "pearson_150000.csv")

#my_nets <- c("regulon_5000.csv", "regulon_7500.csv", "regulon_10000.csv", "regulon_25000.csv",
#             "regulon_50000.csv", "regulon_75000.csv", "regulon_100000.csv", "regulon_150000.csv")

my_nets <- c("rf_5000.csv", "rf_7500.csv", "rf_10000.csv", "rf_25000.csv",
             "rf_50000.csv", "rf_75000.csv", "rf_100000.csv", "rf_150000.csv")

#lets store the values here, one matrix per gold standard, matrix is threshold by precision/recall/fscore
pmap_results <- matrix(, nrow = 8, ncol = 3)
colnames(pmap_results) <- c("Precision", "Recall", "Fscore")
string_results <- matrix(, nrow = 8, ncol = 3)
colnames(string_results) <- c("Precision", "Recall", "Fscore")
row.names(pmap_results) <- c(5000, 7500, 10000, 25000, 50000, 75000, 100000, 150000)
row.names(string_results) <- c(5000, 7500, 10000, 25000, 50000, 75000, 100000, 150000)

for(i in 1:length(my_nets)){
  #then read in the predicted network
  net1 <- read.csv(my_nets[i], header = T, as.is = T)
  colnames(net1) <- c("V1", "V2", "V3")
  net1_graph <- graph.data.frame(net1, directed = F)
  
  #true positive edges are those that appear in both networks
  int_pmap <- intersection(pmap_graph, net1_graph)
  tps_pmap <- dim(get.edgelist(int_pmap))[1]
  
  #precision is tp/(tp+fp) so true positives over total number of edges in predicted network
  precision_pmap <- tps_pmap/nrow(net1)
  pmap_results[i,1] <- precision_pmap
  #recall is tp/(tp+fn) so true positive over total number of edges in ground truth
  recall_pmap <- tps_pmap/nrow(plasmoMap)
  pmap_results[i,2] <- recall_pmap
  pmap_results[i,3] <- 2*((precision_pmap*recall_pmap)/(precision_pmap+recall_pmap))
  
  #do the same for the string network
  int_st <- intersection(string_graph, net1_graph)
  tps_st <- dim(get.edgelist(int_st))[1]
  precision_st <- tps_st/nrow(net1)
  string_results[i,1] <- precision_st
  recall_st <- tps_st/nrow(string)
  string_results[i,2] <- recall_st
  string_results[i,3] <- 2*((precision_st*recall_st)/(precision_st+recall_st))
}


#write.csv(pmap_results, "../gold_standards/plasmoMAP/pearson_to_Pmap.csv", row.names = T, col.names = T, quote = F)
#write.csv(string_results, "../gold_standards/string_network/pearson_to_string.csv", row.names = T, col.names = T, quote = F)

#write.csv(pmap_results, "../gold_standards/plasmoMAP/regulon_to_Pmap.csv", row.names = T, col.names = T, quote = F)
#write.csv(string_results, "../gold_standards/string_network/regulon_to_string.csv", row.names = T, col.names = T, quote = F)

write.csv(pmap_results, "../gold_standards/plasmoMAP/rf_to_Pmap.csv", row.names = T, col.names = T, quote = F)
write.csv(string_results, "../gold_standards/string_network/rf_to_string.csv", row.names = T, col.names = T, quote = F)



#now lets plot the results 
pear_pmap <- read.csv("../gold_standards/plasmoMAP/pearson_to_Pmap.csv")
pear_st <- read.csv("../gold_standards/string_network/pearson_to_string.csv")
reg_pmap <- read.csv("../gold_standards/plasmoMAP/regulon_to_Pmap.csv")
reg_st <- read.csv("../gold_standards/string_network/regulon_to_string.csv")
rf_pmap <- read.csv("../gold_standards/plasmoMAP/rf_to_Pmap.csv")
rf_st <- read.csv("../gold_standards/string_network/rf_to_string.csv")

#pmap precision/recall plot
plot(pear_pmap$Recall, pear_pmap$Precision, xlab="Recall", ylab="Precision", ylim=c(0,.06), type="o", lwd=2)
lines(reg_pmap$Recall, reg_pmap$Precision, col="red", type="o", lwd=2)
lines(rf_pmap$Recall, rf_pmap$Precision, col="blue", type="o", lwd=2)
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#string precision/recall plot
plot(pear_st$Recall, pear_st$Precision, xlab="Recall", ylab="Precision", ylim = c(0.5,1), xlim = c(0,0.045), type="o", lwd=2)
lines(reg_st$Recall, reg_st$Precision, col="red", type="o", lwd=2)
lines(rf_st$Recall, rf_st$Precision, col="blue", type="o", lwd=2)
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#fscore pmap
threshs <- pear_pmap$X
plot(threshs, pear_pmap$Fscore, xlab="Number of edges", ylab="F-score", ylim = c(0,.05), type="o", lwd=2)
lines(threshs, reg_pmap$Fscore, col="red", lwd=2, type="o")
lines(threshs, rf_pmap$Fscore, col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#fscore string
threshs <- pear_st$X
plot(threshs, pear_st$Fscore, xlab="Number of edges", ylab="F-score", ylim = c(0,.09), type="o", lwd=2)
lines(threshs, reg_st$Fscore, col="red", lwd=2, type="o")
lines(threshs, rf_st$Fscore, col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)



################################### test ######################################

#then read in the predicted network
net1 <- read.csv("../test_precision_recall/betaCV_test_net.csv", header=T, as.is=T)
colnames(net1) <- c("V1", "V2", "V3")
net1_graph <- graph.data.frame(net1, directed = F)

plasmoMap <- read.csv("../test_precision_recall/pc_testnet.csv")
pmap_graph <- graph.data.frame(plasmoMap, directed = F)

#true positive edges are those that appear in both networks
int_pmap <- intersection(pmap_graph, net1_graph)
tps_pmap <- dim(get.edgelist(int_pmap))[1]

#precision is tp/(tp+fp) so true positives over total number of edges in predicted network
precision_pmap <- tps_pmap/nrow(net1)
#recall is tp/(tp+fn) so true positive over total number of edges in ground truth
recall_pmap <- tps_pmap/nrow(plasmoMap)
fMeasure_pmap <- 2*((precision_pmap*recall_pmap)/(precision_pmap+recall_pmap))

#do the same for the string network
int_st <- intersection(string_graph, net1_graph)
tps_st <- dim(get.edgelist(int_st))[1]

precision_st <- tps_st/nrow(net1)
recall_st <- tps_st/nrow(string)
fMeasure_st <- 2*((precision_st*recall_st)/(precision_st+recall_st))

test_results <- matrix(, nrow = 8, ncol = 3)
colnames(test_results) <- c("Precision", "Recall", "Fscore")
test_results[1,1] <- precision_pmap
test_results[1,2] <- recall_pmap
test_results[1,3] <- fMeasure_pmap
