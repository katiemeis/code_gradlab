#calculate ratio of false positives to true positives
#calculate ratio of false negatives to true negatives 

# TP: in both network and ground truth
# FP: in network but not in ground truth
# TN: not in either network or ground truth
# FN: not in network but in ground truth

setwd("~/ferdig_rotation/regulon_validation/series_thresholds/")

############################# GROUND TRUTH ###############################

#need to only look at genes associated with GO terms that have 10 or more genes
gocounts <- read.csv("PID_6_10.csv") #PID_6_10
row.names(gocounts) <- gocounts[,1]
gocounts <- gocounts[,-1]

#count how many genes are associated with each term
gosums <- colSums(gocounts)
gosumsvec <- as.integer(gosums)

#only keep go terms with at least ten genes
go_ten_index <- which(gosumsvec > 9) #1 for test
go_ten_table <- gocounts[,go_ten_index]

#delete genes that are all zeros so we only keep nodes associated with an ID n > 9
go_final_ten <- go_ten_table[rowSums(go_ten_table) > 0,]
#node_sums <- data.frame(as.numeric(rowSums(go_final_ten)))
#row.names(node_sums) <- row.names(go_final_ten)

from_vec <- vector()
to_vec <- vector()
#if two genes have the same go term, make an edge
#for each column (GO) make edge for each pair of genes
for(i in colnames(go_final_ten)){
  rows_mapped <- which(go_final_ten[,i] > 0.5)
  for(j in 1:(length(rows_mapped)-1)){
    f = j+1
    for(k in f:length(rows_mapped)){
      from_vec <- c(from_vec, row.names(go_final_ten)[rows_mapped[j]])
      to_vec <- c(to_vec, row.names(go_final_ten)[rows_mapped[k]])
    }
  }
}

#double checked on test and with QI's code - I get the same value
#no duplicates - edges only go one way
gt_edgelist <- cbind(from_vec, to_vec)
gt_adj_mat <- as.matrix(get.adjacency((graph.data.frame(gt_edgelist[,1:2]))))


############################### EXPERIMENTAL NETWORK #######################################
#no duplicates - edges only go one way
graph1 <- read.csv("rf_5000.csv", header=T, as.is=T)
colnames(graph1) <- c("V1", "V2", "V3")
sym_edges <- data.frame(graph1$V2, graph1$V1, graph1$V3)
sym_graph <- as.data.frame(rbind(as.matrix(graph1), as.matrix(sym_edges)))
g1_adj_mat <- as.matrix(get.adjacency((graph.data.frame(sym_graph[,1:2]))))
g1_adj_mat[g1_adj_mat > 2] <- 0

#add genes not in network but in GO to g1_adj_mat
for(i in row.names(gt_adj_mat)){
  if(!(i %in% row.names(g1_adj_mat))){
    new_row = rep(0, length(colnames(g1_adj_mat)))
    g1_adj_mat <- rbind(g1_adj_mat, new_row)
    row.names(g1_adj_mat) <- replace(row.names(g1_adj_mat), length(row.names(g1_adj_mat)), i)
    new_col = rep(0, length(row.names(g1_adj_mat)))
    g1_adj_mat <- cbind(g1_adj_mat, new_col)
    colnames(g1_adj_mat) <- replace(colnames(g1_adj_mat), length(colnames(g1_adj_mat)), i)
  }
}

#add genes in network but not in GO to gt_adj_mat
for(i in row.names(g1_adj_mat)){
  if(!(i %in% row.names(gt_adj_mat))){
    new_row = rep(0, length(colnames(gt_adj_mat)))
    gt_adj_mat <- rbind(gt_adj_mat, new_row)
    row.names(gt_adj_mat) <- replace(row.names(gt_adj_mat), length(row.names(gt_adj_mat)), i)
    new_col = rep(0, length(row.names(gt_adj_mat)))
    gt_adj_mat <- cbind(gt_adj_mat, new_col)
    colnames(gt_adj_mat) <- replace(colnames(gt_adj_mat), length(colnames(gt_adj_mat)), i)
  }
}


#for each cell, compare the matrices, add to proper count
# tp when both 1
# fp when g1 cell is 1 but gt cell is 0
# tn when both 0
# fn when g1 cell is 0 but gt cell is 1
tp_count = 0
fp_count = 0
tn_count = 0
fn_count = 0
for(i in 1:(length(row.names(gt_adj_mat))-1)){
  f = i+1
  for(j in f:length(row.names(gt_adj_mat))){
    x = row.names(gt_adj_mat)[i]
    y = row.names(gt_adj_mat)[j]
    if(g1_adj_mat[x,y]==1 & gt_adj_mat[x,y]==1){
      tp_count = tp_count+1
    }
    if(g1_adj_mat[x,y]==1 & gt_adj_mat[x,y]==0){
      fp_count = fp_count+1
    }
    if(g1_adj_mat[x,y]==0 & gt_adj_mat[x,y]==0){
      tn_count=tn_count+1
    }
    if(g1_adj_mat[x,y]==0 & gt_adj_mat[x,y]==1){
      fn_count=fn_count+1
    }
  }
}

#divide all by 2 because symmetric matrix
