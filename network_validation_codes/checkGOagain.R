setwd("~/ferdig_rotation/regulon_validation/new_go_Qi/")
gocounts <- read.csv("PID_6_10.csv")
row.names(gocounts) <- gocounts[,1]
gocounts <- gocounts[,-1]
gosums <- colSums(gocounts)
gosumsvec <- as.integer(gosums)

#find which column the 1 is in and add to vector
index_sp <- vector()
for(i in 1:30){
  col_num <- which(gocounts[i,] == 1)
  index_sp <- c(index_sp, col_num)
}

index_sp2 <- vector()
for(j in 1928:2062){
  col_num2 <- which(gocounts[j,] == 1)
  index_sp2 <- c(index_sp2, col_num2)
}

gosumsvec[102] #40, GO.0006260, DNA replication
gosumsvec[129] #93, GO.0006351, transcription
gosumsvec[157] #259, go GO.0006412, translation

length(which(index_sp == 157)) #26
length(which(index_sp2 == 157)) #116
#so 142 of the 259 are wrong species
length(which(index_sp == 129)) #4
length(which(index_sp2 == 129)) #18
#so 22 of the 93 are wrong species
length(which(index_sp == 102)) #0
length(which(index_sp2 == 102)) #1
#so 1 of the 40 are wrong species

#remove the rows with wrong species
#can leave the columns alone since all GOs still have at least 1
new_PID <- gocounts[-c(1:30, 1928:2062),]
new_ID <- cbind(GeneID=rownames(new_PID), new_PID)

#write new GO file
write.csv(new_PID, "testnew_PID_6_10.csv")
#write.csv(new_ID, "newquotes_PID_6_10.csv", row.names = F, quote = T)
library(dplyr)
new_ID_char <- new_ID %>% mutate_all(as.character)
write.csv(new_ID_char, "newquotesall_PID_6_10.csv", row.names = F, quote = T)

#get new histogram
new_col_sums <- colSums(new_PID)
hist(new_col_sums, breaks=129)
#get new GO network size, gives 43150 edges
new_gt_edge_count = 0
for(i in 1:length(new_col_sums)){
  new_gt_edge_count = new_gt_edge_count + choose(new_col_sums[i],2)
  
}

length(which(new_col_sums >1)) #308 GO terms


#only nodes in the intersection of largest full RF network and GO file
rnet <- read.csv("rf_size_names_nodups2.csv")
g <- graph_from_edgelist(as.matrix(rnet[,c(1,2)]))
names_vect <- V(g)$name
#find the intersection of the nodes
int_nodes <- intersect(names_vect, row.names(gocounts))
#only keep intersection
int_PID <- gocounts[int_nodes,]
#see if any GO terms are now zero and remove them
int_col_sums <- colSums(int_PID)
length(which(int_col_sums > 0)) #658
int_PID_final <- int_PID[,which(int_col_sums > 0)]
max(colSums(int_PID_final))
hist(colSums(int_PID_final), breaks=108, xlab= "Number of genes")



