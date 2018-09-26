
######## if want to use Qi's original file

#load in Qi's original GO file
go_file <- read.csv("PID_6_10.csv", header=T, as.is=T)
row.names(go_file) <- go_file$GeneID
go_file1 <- go_file[,-1]
dim(go_file1)

#load in names map from curating_networks.R
namesmap <- read.csv("names_map_1to1.csv", header = T, as.is = T)
#want new names, column 1
new_keep <- namesmap[,1]
new_keep[c(1:5)]
length(new_keep)

#keep only these rows of the PID
length(which(row.names(go_file1) %in% new_keep))
gene_index_keep <- which(row.names(go_file1) %in% new_keep)
go_cut_genes <- go_file1[gene_index_keep,]
dim(go_cut_genes) #1863 genes, still 665 GOs

#now remove any GO terms that have no genes associated
gene_counts <- colSums(go_cut_genes)
length(gene_counts)
empty_gos <- which(gene_counts == 0) #only 2 drop out
gos_with_1plus <- which(gene_counts > 0) #just check
go_genes_GOs <- go_cut_genes[,-empty_gos]
dim(go_genes_GOs) # 1863 genes, 663 GOs

#check GOsums again to make sure non are zero
final_sum_check <- colSums(go_genes_GOs)
which(final_sum_check == 0) #none, looks ok

write.csv(go_genes_GOs, "reduced_PID_6_10.csv")

