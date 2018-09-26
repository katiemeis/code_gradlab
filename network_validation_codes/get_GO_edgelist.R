
#turn Qi's PID file into an edgelist

#read in the PID file
gocounts <- read.csv("~/ferdig_rotation/regulon_validation/original_nets/GO_file/PID_6_10_NEW.csv", as.is = T)
row.names(gocounts) <- gocounts[,1]
gocounts <- gocounts[,-1]

#loop over each column - get all-on-all edges from each column
all_edges <- vector()
for(i in 1:length(colnames(gocounts))){
  genes_sameGO <- row.names(gocounts)[which(gocounts[,i] > 0)]
  if(length(genes_sameGO) > 1){
    edges_GO <- t(combn(genes_sameGO,2))
    all_edges <- rbind(all_edges, edges_GO)  
  }
}

