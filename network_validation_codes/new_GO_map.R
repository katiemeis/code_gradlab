
##### cut down Qi's GO file to keep only genes in final map and only 

setwd("~/ferdig_rotation/regulon_validation/original_nets/GO_file/")

#################################### compare to GO output if only search these genes #################

newGO_search <- read.table("new_GO_search_only1to1.txt", sep = "\t", header=T, as.is = T)

#if cols 5 and 6 are both NA meaning it has no computed or curated process, remove
remove_ind <- vector()
for(i in 1:nrow(newGO_search)){
  if(newGO_search[i,5] == "N/A" & newGO_search[i,6] == "N/A"){
    remove_ind <- c(remove_ind, i)
  }
}

test_GO_noNA <- newGO_search[-remove_ind,]
dim(test_GO_noNA)
#1900 genes left

#if .2, remove them, don't want isoform copies
library(stringr)
new_cols_go <- str_split_fixed(test_GO_noNA$X.source_id., fixed("."), 2)
test_GO_noNA_cols <- cbind(test_GO_noNA, new_cols_go)
test_GO_noNA_cols[,9] <- as.numeric(test_GO_noNA_cols[,9])
first_copy <- which(test_GO_noNA_cols[,9] == 1)
final_GOlist <- test_GO_noNA_cols[first_copy,]

dim(final_GOlist) #1878 genes with at least 1 process ID
write.table(final_GOlist, "final_GO_dataframe.txt",quote=F, sep="\t", col.names = T, row.names = F)
#write.csv(final_GOlist, "final_GO_dataframe.csv", quote=F)

########################################################## convert to Qi's PID format with her code ########################
gene <- final_GOlist[,c(1,5,6)]
colnames(gene) <- c("GeneID","CompProcessID", "CurProcessID")
#check if any duplicate rows
gene0 <- unique( gene[ , 1:3 ] ) #none, good

#expand this into Qi's format
#most is 8
gene$PID1 <- substr(gene$CompProcessID,1, 10)
gene$PID2 <- substr(gene$CompProcessID,13, 22)
gene$PID3 <- substr(gene$CompProcessID,25, 34)
gene$PID4 <- substr(gene$CompProcessID,37, 46)
gene$PID5 <- substr(gene$CompProcessID,49, 58)
gene$PID6 <- substr(gene$CompProcessID,61, 70)
gene$PID7 <- substr(gene$CompProcessID,73, 82)
gene$PID8 <- substr(gene$CompProcessID,85, 94) 

#table(last PID) if empty, then done

gene1 <- gene[,c(1,3,4,5,6,7,8,9,10,11)] # remove CompProcessID#

## substract Process ID from Curated Process ID ##
#most is 5
gene1$PID9 <- substr(gene$CurProcessID, 1, 10)
gene1$PID10 <- substr(gene$CurProcessID, 13, 22)
gene1$PID11 <- substr(gene$CurProcessID, 25, 34)
gene1$PID12 <- substr(gene$CurProcessID, 37, 46)
gene1$PID13 <- substr(gene$CurProcessID, 49, 58)

#remove curated process ID
gene2 <- gene1[,-2]

#checking for duplicates between computed and curated, get list of unique go terms for each gene
# 
# for (i in 1:nrow(gene2)){
#   for (j in 9:13){
#     #print (j)
#     if (gene2[i,j] == gene2[i,'PID1'] | gene2[i,j] == gene2[i,'PID2'] | gene2[i,j] == gene2[i,'PID3'] |
#         gene2[i,j] == gene2[i,'PID4'] | gene2[i,j] == gene2[i,'PID5'] | gene2[i,j] == gene2[i,'PID6'] | 
#         gene2[i,j] == gene2[i,'PID7']){
#       gene2[i,j] <- ""
#     }
#     else {
#       gene2[i,j] = gene2[i,j]
#     }
#   }
# }
# 
## merge all the PIDs into one column then remove duplicate ##
# gene to PID1, then append gene to PID2, etc.
d <- data.frame(GeneID=character(), PID= character(),stringsAsFactors=FALSE)
for (i in 1:13){
  d_i <- gene2[,c(1, i+1)]
  names(d_i) <- c("GeneID", "PID")
  d<- rbind(d, d_i)
}

#remove anything empty or N/A
dd <- d[(d$PID !="" & d$PID !="N/A"),]

##remove duplication
dd2 <- unique( dd[ , 1:2 ] )

#make adj matrix
library(reshape2)
ddd <- dcast(dd2, GeneID~PID)
ddd[is.na(ddd)] <- 0

## replace char to 1 ##
for (i in (1:nrow(ddd))){
  for (j in (2:ncol(ddd))){
    if (ddd[i,j] != 0){ddd[i,j] = 1}
  }
}

write.csv(ddd, "PID_6_10_NEW.csv",row.names = F)

input <- read.csv("PID_6_10_NEW.csv", as.is=T)
dim(ddd)
ddd_nonames <- input[,-1]
sums <- colSums(ddd_nonames)
which(sums == 0)
