#setwd("/Users/qili/Dropbox/Genetic Networks/GO_IDs");
setwd("~/ferdig_rotation/regulon_validation/new_go_Qi/testing_code/")
gene0 <- read.csv("GeneGos.csv")
gene0 <- gene0[,c(1,2,3)]
names(gene0) <- c("GeneID","CompProcessID", "CurProcessID")

##remove duplication, by means all three columns are exactly same
gene <- unique( gene0[ , 1:3 ] )

#### Substract ProcessID from computed Computed ProcessID
gene$PID1 <- substr(gene$CompProcessID,1, 10)
gene$PID2 <- substr(gene$CompProcessID,13, 22)
gene$PID3 <- substr(gene$CompProcessID,25, 34)
gene$PID4 <- substr(gene$CompProcessID,37, 46)
gene$PID5 <- substr(gene$CompProcessID,49, 58)
gene$PID6 <- substr(gene$CompProcessID,61, 70)
gene$PID7 <- substr(gene$CompProcessID,73, 82)
gene$PID8 <- substr(gene$CompProcessID,85, 94) ## empty, remove later, use table function to check whether it is empty or not ##

gene1 <- gene[,c(1,3,4,5,6,7,8,9,10)] # remove CompProcessID and PID8 since it is empty #

## substract Process ID from Curated Process ID ##
gene1$PID8 <- substr(gene$CurProcessID, 1, 10)
gene1$PID9 <- substr(gene$CurProcessID, 13, 22)
gene1$PID10 <- substr(gene$CurProcessID, 25, 34)
gene1$PID11 <- substr(gene$CurProcessID, 37, 46)
gene1$PID12 <- substr(gene$CurProcessID, 49, 58)
gene1$PID13 <- substr(gene$CurProcessID, 61, 70) ## empty, remove later ##

#remove curProcessID and PID13 ##
gene2 <- gene1[,-c(2,15)]

### check PID8 - PID13, whether they are same as PID1 - PID7 ##
for (i in 1:nrow(gene2)){
  for (j in 9:13){
    #print (j)
    if (gene2[i,j] == gene2[i,'PID1'] | gene2[i,j] == gene2[i,'PID2'] | gene2[i,j] == gene2[i,'PID3'] |
    gene2[i,j] == gene2[i,'PID4'] | gene2[i,j] == gene2[i,'PID5'] | gene2[i,j] == gene2[i,'PID6'] | 
    gene2[i,j] == gene2[i,'PID7']){
      gene2[i,j] <- ""
    }
    else {
      gene2[i,j] = gene2[i,j]
    }
  }
}

## merge all the PIDs into one column then remove duplicate ##
d <- data.frame(GeneID=character(), PID= character(),stringsAsFactors=FALSE)
for (i in 1:12){
  d_i <- gene2[,c(1, i+1)]
  names(d_i) <- c("GeneID", "PID")
  d<- rbind(d, d_i)
}
dd <- d[(d$PID !="" & d$PID !="N/A"),]

##remove duplication
dd2 <- unique( dd[ , 1:2 ] )

library(reshape2)
ddd <- dcast(dd2, GeneID~PID)
ddd[is.na(ddd)] <- 0

## replace char to 1 ##
for (i in (1:nrow(ddd))){
  for (j in (2:ncol(ddd))){
    if (ddd[i,j] != 0){ddd[i,j] = 1}
  }
}
write.csv(ddd, "PID_6_10.csv",row.names = F)

