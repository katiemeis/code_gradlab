

################## need to rename the gold standards ##########################
setwd("~/ferdig_rotation/regulon_validation/original_nets/gold_standards/")

# read in the two graphs and get the union of their nodes lists
# then can search this list on plasmoDB to get the input file to my renaming script
string <- read.table("gold_standards/string_network/string_edgelist.txt", sep=" ", header=T)
string_net <- graph_from_edgelist(as.matrix(string[,c(1,2)]))
string_nodes <- V(string_net)$name
plasmoMap <- read.graph("gold_standards/plasmoMAP/highConfPairs.lgl", format = "lgl")
pm_nodes <- V(plasmoMap)$name
length(intersect(string_nodes, pm_nodes))
length(union(string_nodes, pm_nodes)) #correct length
#union of the two nodes lists is what we need to search in plasmoDB to get names file. 
#write.csv(union(string_nodes, pm_nodes), "full_nodes_forsearch.csv")



names_key1 <- read.table("full_nodes_names_plasmoDB.txt", sep = "\t", header=T, as.is = T)
library(stringr)
new_cols <- str_split_fixed(names_key1$X.source_id., fixed("."), 2)
names_key <- cbind(names_key1, new_cols)
names_key[,8] <- as.numeric(names_key[,8])
first_copy <- which(names_key[,8] == 1)
names_key_final <- names_key[first_copy,]
#after removing isoforms have 5191 genes

#only keep old genes that only appear once
#when the same input goes to different output (same old gene maps to more than one new gene)
n_occur_2 <- data.frame(table(names_key_final$X.Input.ID.))
length(which(n_occur_2[,2] == 1)) #5155 only occur once
to_keep <- as.vector(n_occur_2[which(n_occur_2[,2] == 1),1]) #what are the gene names of non-dups
names_key_noolddups <- names_key_final[names_key_final$X.Input.ID. %in% to_keep,] #keeping only single

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
names_key1to1 <- names_key_noolddups[not_multiples,] #5026 genes left

#we only want the old and new names, not the rest of the columns
map_1to1 <- names_key1to1[,c(1,4)]
length(unique(map_1to1$X.Input.ID.))

#write.csv(map_1to1, "names_map_1to1_goldSt.csv", quote=F, row.names = F)


################################################# rename the edgelists #######################################

#net_el <- read.table("string_network/string_edgelist.txt", sep=" ", header=T)
net_el1 <- read.csv("plasmoMAP/plasmoMAP_el.csv", header = T, as.is=T)

#for plasmomap, no edge weights so lets just add a column of 1s so this runs smoothly
fake_weigths <- rep(1, nrow(net_el1))
net_el <- cbind(net_el1, fake_weigths)

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
V(g) #5008 nodes, matches total set form above in the mapping so good!

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
fin_el_ord_check <- final_edges[order(-final_edges$Weight),]

#write.csv(fin_el_ord_check, "string_network/string_renamed_noloops.csv", quote = F, row.names = F, col.names = T)
write.csv(fin_el_ord_check, "plasmoMAP/plasmoMAP_renamed_noloops.csv", quote = F, row.names = F, col.names = T)
