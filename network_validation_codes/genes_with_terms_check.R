
setwd("~/ferdig_rotation/regulon_validation/new_go_Qi/testing_code/")
gene0 <- read.csv("GeneGos.csv", as.is = T)
gene0 <- gene0[,c(1,2,3)]
names(gene0) <- c("GeneID","CompProcessID", "CurProcessID")

#find rows that have all N/A
indexvec <- which(gene0$CompProcessID == "N/A" & gene0$CurProcessID == "N/A")

#dataframe with only genes that have at least one go term
genesgo <- gene0[-indexvec,]

smallnet <- read.csv("../../series_thresholds/rf_5000.csv")
largenet <- read.csv("../../series_thresholds/rf_150000.csv")

library(igraph)
smat <- as.matrix(get.adjacency((graph.data.frame(smallnet))))
dim(smat)
lmat <- as.matrix(get.adjacency((graph.data.frame(largenet))))
dim(lmat)
snode <- row.names(smat)
lnode <- row.names(lmat)

golist <- genesgo[,1]
length(intersect(golist, snode))
length(intersect(golist, lnode))

#need to only look at genes associated with GO terms that have 10 or more genes
gocounts <- read.csv("../../series_thresholds/PID_6_10.csv")
row.names(gocounts) <- gocounts[,1]
gocounts <- gocounts[,-1]
gosums <- colSums(gocounts)
gosumsvec <- as.integer(gosums)
#only keep go terms with at least ten genes
go_ten_index <- which(gosumsvec > 9)
go_ten_table <- gocounts[,go_ten_index]
#delete genes that are all zeros 
go_final_ten <- go_ten_table[rowSums(go_ten_table) > 0,]

gonodes <- row.names(go_final_ten)
length(intersect(gonodes, snode))
length(intersect(gonodes, lnode))


#plot histogram of GO term gene counts
dat_hist <- read.csv("../new_go_Qi/newPIDgenecounts.csv", header=F)
y <- as.numeric(dat_hist[1,])
barplot(y, xlab = "GO terms", ylab = "Gene Counts")
abline(h=10, col="red", lwd=2, lty=2)
