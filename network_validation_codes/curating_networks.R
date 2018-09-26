#curating to final networks
setwd("~/ferdig_rotation/regulon_validation/original_nets/")



##################################### getting list of 1:1 mappings ##########################################

#names_key1 <- read.table("test_renaming_files/test_map_cpplasmo.txt", sep="\t", header = T, as.is=T)
names_key1 <- read.table("plasmoDB_names_hu.txt", sep = "\t", header=T, as.is = T)

#first there are isoforms in second colunm need to get rid of
#separate out the name from period, remove anything > 1
library(stringr)
new_cols <- str_split_fixed(names_key1$X.source_id., fixed("."), 2)
names_key <- cbind(names_key1, new_cols)
names_key[,8] <- as.numeric(names_key[,8])
first_copy <- which(names_key[,8] == 1)
names_key_final <- names_key[first_copy,]
#after removing isoforms have 3659 genes - this means 46 genes have no new name?

#only keep old genes that only appear once
#when the same input goes to different output (same old gene maps to more than one new gene)
n_occur_2 <- data.frame(table(names_key_final$X.Input.ID.))
length(which(n_occur_2[,2] == 1)) #3597 only occur once
to_keep <- as.vector(n_occur_2[which(n_occur_2[,2] == 1),1]) #what are the gene names of non-dups
names_key_noolddups <- names_key_final[names_key_final$X.Input.ID. %in% to_keep,] #keeping only single takes 3705-3597

#check that all inputs are unique, dim should match
length(unique(names_key_noolddups$X.Input.ID.))
dim(names_key_noolddups)

#check again
check_occur <- data.frame(table(names_key_noolddups$X.Input.ID.))
length(which(check_occur$Freq >1)) #returns 0, good

#now only keep new names that only occur once
#in file, if same new maps to two+ old, old are comma sep in input column
#split by comma
new_cols2 <- str_split_fixed(names_key_noolddups$X.Input.ID., fixed(","), 2)
#get index of which genes have empty second column (only one gene maps)
not_multiples <- which(new_cols2[,2] == "")
#only keep these rows of the names file
names_key1to1 <- names_key_noolddups[not_multiples,] #3563 genes out of original 3705, took out 142 genes

#we only want the old and new names, not the rest of the columns
map_1to1 <- names_key1to1[,c(1,4)]
length(unique(map_1to1$X.Input.ID.))

write.csv(map_1to1, "names_map_1to1.csv", quote=F, row.names = F)

#check which ones are removed
#hudata <- read.csv("../Hu_data.csv", header = T, as.is=T)
#full_hu_genes <- hudata[,1]
#length(setdiff(full_hu_genes,map_1to1[,2])) #142, matches above map

#write missing genes to file to search in plasmoDB
#missing_genes <- setdiff(full_hu_genes,map_1to1[,2])
#write.csv(missing_genes, "genes_removed.csv")





############################### Convert original network outputs to final ######################################################

### now we want to rename the edgelists (use map1to1 as the key, new followed by old)
### then remove selfloops and duplicate edges

#net_el <- read.table("test_renaming_files/test_renaming_and_loopsdups.txt", sep = " ", header=T, as.is=T)
#net_el <- read.table("full_regulon.txt", sep=" ", header = T, as.is=T)
#net_el <- read.table("hu_randomForest_edgelist.txt", sep="\t", header=F)
net_el <- read.csv("adalasso_hu_edgelist.csv", header = T, as.is=T)
#net_el <- read.table("full_pearson.txt", sep = " ", header = T, as.is=T)

colnames(net_el) <- c("From", "To", "Weight")

#name old in map_1to1 "From" so we merge from first
colnames(map_1to1) <- c("NewName", "From")
merge1 <- merge(map_1to1, net_el, by="From")
one_merge_el <- merge1[,c(2,3,4)]
colnames(one_merge_el) <- c("From", "To", "Weight")

#now name old in map_1to1 "To" so we merge to also
colnames(map_1to1) <- c("NewName", "To")
merge2 <- merge(map_1to1, one_merge_el, by="To")
two_merge_el <- merge2[,c(2,3,4)]
colnames(two_merge_el) <- c("To", "From", "Weight")

#lets check number of nodes
library(igraph)
g <- graph_from_edgelist(as.matrix(two_merge_el[,c(1,2)]))
V(g) # 3563 nodes, matches total set form above in the mapping so good!

#we want to order from strongest to weakest
net_el_ord <- two_merge_el[order(-two_merge_el$Weight),]

#get rid of self loops and duplicate edges
#make graph with weights
simple_g <- graph.data.frame(net_el_ord)
#take out self loops
simplifygraph <- simplify(simple_g, remove.loops = T, remove.multiple = F) 
#take out duplicate edges
fingraph <- as.undirected(simplifygraph, mode=c("collapse"), edge.attr.comb = "mean")
#convert back to edge list from graph
final_edges <- get.data.frame(fingraph)
#order again to check
fin_el_ord_check <- final_edges[order(-final_edges$Weight),] #compare to net_el_ord
#write.csv(fin_el_ord_check, "full_regulon_renamed.csv", quote = F, row.names = F, col.names = T)
#write.csv(fin_el_ord_check, "full_randonForest_renamed.csv", quote = F, row.names = F, col.names = T)
#write.csv(fin_el_ord_check, "full_adalasso_renamed.csv", quote=F, row.names = F, col.names = T)
#write.csv(fin_el_ord_check, "full_pearson_renamed.csv", quote=F, row.names = F, col.names = T)




######################################################################################



#### EXTRAS #####################

#lets make a test for the merging
nodes_1 <- c("gene1","gene1", "gene2", "gene3", "gene4", "gene1")
nodes_2 <- c("gene2", "gene5","gene1", "gene3", "gene3", "gene3")
mis <- c(.2, .3, .4, .5, .6, .7)
testdf <- data.frame(nodes_1, nodes_2, mis)
testdf

nn <- c("geneA", "geneB", "geneC")
nodes_1 <- c("gene1", "gene2", "gene3")
nk <- data.frame(nn, nodes_1)

mt <- merge(nk, testdf, by="nodes_1")
mt
halfdf <- mt[,c(2,3,4)]
colnames(halfdf)[2] <- "nodes_1"
mt2 <- merge(nk, halfdf, by="nodes_1")
mt2
finaldf <- mt2[,c(2,3,4)]
finaldf

