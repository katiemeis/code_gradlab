###   5/23/17 Script to Compare GCN clusters of Mok networks (correlation based)

###   Set working directory
setwd("/Users/kbuttons/Documents/Research/ART rewiring/full network/")


### Network rewiring following Zhou 2013
troph_RN <- read.delim("troph_resistant_network.csv",sep=",",header=TRUE,as.is=TRUE)
troph_SN <- read.delim("troph_sensitive_network.csv",sep=",",header=TRUE,as.is=TRUE)

all_RN <- read.delim("AllStages_resistant_network.csv",sep=",",header=TRUE,as.is=TRUE)
all_SN <- read.delim("AllStages_sensitive_network.csv",sep=",",header=TRUE,as.is=TRUE)

ring_RN <- read.delim("ring_resistant_network.csv",sep=",",header=TRUE,as.is=TRUE)
ring_SN <- read.delim("ring_sensitive_network.csv",sep=",",header=TRUE,as.is=TRUE)

rewiring_score <- function(i) {
  cor_r <- troph_RN[i,"cor"]
  cor_s <- troph_SN[i,"cor"]
  Fr <- (1/2)*log((1+cor_r)/(1-cor_r))
  Fs <- (1/2)*log((1+cor_s)/(1-cor_s))
  score <- abs(Fr-Fs)/sqrt((1/(49-3))+(1/(49-3)))
  p <- pnorm(score,lower.tail=FALSE)
  if(i%%100000==0){
    write.csv(i,paste("out",i,".csv"))
  }
  return(c(troph_RN[i,"i"],troph_RN[i,"j"],cor_r,cor_s,score,p))
}

require(parallel)
cores=6
setwd("/Users/kbuttons/Documents/Research/ART rewiring/full network/ring_rewiring/")
troph_rewiring <- matrix(unlist(mclapply(1:14723451, FUN=rewiring_score, mc.cores=cores)),ncol=6,byrow=TRUE)
colnames(troph_rewiring) <- c("Gene_i","Gene_j","cor_i","cor_j","score","p")
adj.p <- p.adjust(troph_rewiring[,"p"],method="fdr")

write.csv(cbind(troph_rewiring,adj.p),"all_rewiring_sig.csv")
troph_rewiring <- cbind(troph_rewiring,adj.p)
sig_troph_rewired <- all_rewiring[which(troph_rewiring[,"adj.p"]<0.05),]
write.csv(sig_troph_rewired,"troph_rewiring_fdr0.001.csv")
