

go.map <- read.table("~/network_science/project/GeneByLocusTag_Summary.txt", header= TRUE, sep= "\t") #LOAD
colnames(go.map) <- c("gene","id","GO","x")
new <- read.table("/Users/emilyherring/Desktop/new.txt")
new.communities<- load("/Users/emilyherring/Desktop/new_communities.R")

#network totals
new.nodes <- row.names(new) #LOAD
new.nodes <- as.data.frame(matrix(unlist(new.nodes)))
colnames(new.nodes) <- "gene"

library(dplyr)
network <- left_join(new.nodes,go.map)

for (j in 1:100){
for (i in 1:dim(network)[1]){
  if(as.character(network[i,3])== "N/A" || is.na(network[i,3])){
    network <- network[-i,]
  }
}}


GOterms.network <- as.data.frame(unique(as.character(network[,3])))
GOcounts.network <- network %>% count(GO)

new_df <- as.data.frame(new.communities)
sig_communities <- array()
for(ii in 1:dim(new_df)[1]){
  if(length(unlist(new_df[ii,])) > 9){
    sig_communities[ii] = new_df[ii,]
    }
}



#####RUN FROM HERE#######


for(ii in 1:length(sig_communities)){
    test.cluster <- (sig_communities[[ii]]) #LOAD
    
    #cluster totals
    test.cluster<- as.data.frame(matrix(unlist(test.cluster)))
    colnames(test.cluster) <- "gene"
    mapped <- left_join(test.cluster,go.map)
    
    for (j in 1:100){
    for (i in 1:dim(mapped)[1]){
      if(as.character(mapped[i,3])== "N/A" || is.na(mapped[i,3])){
        mapped <- mapped[-i,]
      }
    }
    }
    
    GOterms.cluster <- as.data.frame(unique(as.character(mapped[,3])))
    GOcounts.cluster <- mapped %>% count(GO)
    
    #cycle through hiding genes and creating p-values
    GO_sig <- data_frame(ncol = 2, nrow= dim(mapped[1]))
    total.terms <- sum(GOcounts.network$n) #col1+col2
    
    for (i in 1:dim(mapped)[1]){ #cycle through genes in cluster
      count = 0
      cross.check <- mapped[-i,] #create list without gene
      GOcounts.cross <- cross.check %>% count(GO) #count GO terms
      r <- sum(GOcounts.cross$n) #sample size
      pvalue <- array()
      for (j in 1:dim(GOcounts.cross)[1]){ #determine significance of GO terms
        #success is gene with specific GO Term
        #failure is gene without specific GO Term
        a <- as.numeric(GOcounts.cross[j,2]) # success in sample
        rownum <- which(GOcounts.network$GO == GOcounts.cross[j,1]$GO)
        col1 <- as.numeric(GOcounts.network[rownum,2])#success in background
        col2 <- as.numeric(total.terms-col1)#failure in sample
        pvalue <- c(pvalue,phyper(a,col1,col2,r,lower.tail = TRUE))
        } #hypergeometric test
      c<-as.data.frame(cbind(GOcounts.cross[,1],pvalue[-1]))#put pvalues with assoc GO terms
      ind <- which(c$GO == mapped[i,3])
      if((length(ind))!= 0){
         if(c[ind,2] > 0.95){ #if the hidden gene's GO term is in the rest of the cluster
          #if the GO term is significant in the rest of the cluster
          count <- 1}}
     GO_sig[i,] <- c(mapped[i,1],count)
    }
    write.table(GO_sig, file=paste("~/network_science/project/clustering/cluster_new_", ii, ".txt"))
}

