#calculate modularity for a network - clusters defined by GO terms
#input is edgelist with no duplicates, no self loops
#output is modularity score

calc_modularity <- function(file_name){
  #setwd("~/ferdig_rotation/regulon_validation/series_thresholds/")
  setwd("~/ferdig_rotation/regulon_validation/original_nets/series_thresholds/")
  library(igraph)
  
  #need to only look at genes associated with GO terms that have 10 or more genes
  gocounts <- read.csv("../GO_file/PID_6_10_NEW.csv") #PID_6_10
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
  node_sums <- data.frame(as.numeric(rowSums(go_final_ten)))
  row.names(node_sums) <- row.names(go_final_ten)
  
  #get the list of nodes included in GOs w/ at least 10 genes in them
  #these are the only genes we want to consider
  gonodes <- row.names(go_final_ten)
  
  #read in the graph
  graph1 <- read.csv(file_name, header=T, as.is=T)
  colnames(graph1) <- c("V1", "V2", "V3")
  #since graph doesn't have duplicate edges, make symmetric
  sym_edges <- data.frame(graph1$V2, graph1$V1, graph1$V3)
  sym_graph <- as.data.frame(rbind(as.matrix(graph1), as.matrix(sym_edges)))
  
  #find rows of graph1 are edges between GO terms with at least 10 genes
  #only keep edges between gonodes
  ind_vec <- which(sym_graph$V1 %in% gonodes & sym_graph$V2 %in% gonodes)
  final_edges <- sym_graph[ind_vec,]
  #make adjacency matrix for edges left from graph
  adj_mat <- as.matrix(get.adjacency((graph.data.frame(final_edges[,1:2]))))
  network_nodes <- row.names(adj_mat)
  #adj_mat_df <- as.data.frame(as.matrix(get.adjacency((graph.data.frame(final_edges[,1:2])))))
  #degrees of each node
  graph_degree <- degree(graph_from_edgelist(as.matrix(final_edges[,1:2])))
  
  m=dim(final_edges)[1] * 0.5 #only edges between GO terms nodes in network, 0.5 because symmetric
  #m=dim(sym_graph)[1] * 0.5 #all edges
  
  #modulatiy
  # Q = (1/2m) sumc sumij [Aij - (kjki)/2m](1/OiOj)
  # m is number of edges in graph
  # Aij is i,j element of adjacency
  # K is degree of node
  # O is number of clusters i belongs to
  
  #loop over cluster
  #for each node pair in each cluster calculate inside, only if both nodes appear in the network
  clust_sum = 0
  for(c in 1:dim(go_final_ten)[2]){
    clust_ID <- colnames(go_final_ten[c]) #which GO term for this cluster
    clust_nodes <- row.names(go_final_ten)[which(go_final_ten[,c] > 0.5)] #which nodes are associated with this term
    #check to see if all clust_nodes are in go_nodes THIS IS FINE
    #print(length(clust_nodes))
    #print(length(intersect(gonodes, clust_nodes)))
    for(i in 1:length(clust_nodes)){
      #j = i+1
      for(j in i+1:length(clust_nodes)){
        if(clust_nodes[i] %in% network_nodes & clust_nodes[j] %in% network_nodes){
          clust_sum <- clust_sum + ((adj_mat[clust_nodes[i],clust_nodes[j]] - ((graph_degree[[clust_nodes[i]]]*0.5*0.5*graph_degree[[clust_nodes[j]]])/(2*m))) * (1/(node_sums[clust_nodes[i],1]*node_sums[clust_nodes[j],1])))

        }
      }
    }
  }
  
  Q_val <- (1/(2*m))*clust_sum
  return(Q_val)
}


calc_modularity("rf_5000.csv")
