
##### integrity testing for new curation for Hu data #####
setwd("~/ferdig_rotation/regulon_validation/")
old_map_data <- read.csv("Hu_data.csv", as.is=T)
new_names <- read.delim("Newnames_pfal.txt", sep=" ", header=F, as.is=T)
names <-read.csv("Namelists2.csv", header = F, as.is=T)
#my names list
for(i in 1:dim(old_map_data)[1]){
  if(old_map_data[i,1] %in% new_names[,2]){
    index <- which(new_names[,2] == old_map_data[i,1])[1]
    old_map_data[i,1] <- new_names[index,1]
  }
}


#gabes names list
old_map_data2 <- read.csv("Hu_data.csv", as.is=T)
for(i in 1:dim(old_map_data2)[1]){
  if(old_map_data2[i,1] %in% names[,2]){
    index <- which(names[,2] == old_map_data2[i,1])[1]
    old_map_data2[i,1] <- names[index,1]
  }
}

#they change the names the same way
count =0
for(i in 1:dim(old_map_data)[1]){
  if(old_map_data[i,1] == old_map_data2[i,1]){
    count = count + 1
  }
}

#gets rid of 1429 genes if I throw out ones without a new name
old_orig_names <- read.csv("Hu_data.csv", as.is=T)
count =0
for(i in 1:dim(old_map_data)[1]){
  if(old_map_data[i,1] == old_orig_names[i,1]){
    count = count + 1
  }
}





#load in new curation and check correlation of genes at intersection
#load new and subset
new_map_data <- read.csv("GSE19468_Final.csv")
new_sub <- new_map_data[,c(1,9:29)]
#load old and subset
old_sub <- read.csv("Hu_subset_curation_check.csv", header = T, as.is=T)
colnames(old_sub)[1] <- "GeneID"
#load names mapping file and merge
new_names <- read.csv("Namelists2.csv", header = F, as.is=T)
colnames(new_names) <- c("NewGeneID", "GeneID")
old_names <- merge(new_names, old_sub, by = "GeneID")
genes_int <- intersect(new_sub[,1], old_names[,2])

corr_vals <- vector()
for(i in genes_int){
  old_ind <- which(old_names[,2] == i)[1]
  new_ind <- which(new_sub[,1] == i)[1]
  c1 <- cor(as.numeric(old_names[old_ind,3:23]), as.numeric(new_sub[new_ind,2:22]), use = "complete.obs")
  corr_vals <- c(corr_vals, c1)
}

hist(corr_vals)
length(which(corr_vals < 0.75))
95/3509    #0.02707324
length(which(corr_vals < 0.8))
125/3509   #0.03562268

#count missing vals
missing_by_sample = vector()
for(i in 2:ncol(new_map_data)){
  sample_missing <- sum(is.na(new_map_data[,i]))
  missing_by_sample <- c(missing_by_sample, sample_missing)
}

missing_by_gene = vector()
for(i in 1:nrow(new_map_data)){
  gene_missing <- sum(is.na(new_map_data[i,]))
  missing_by_gene <- c(missing_by_gene, gene_missing)
}

hist(missing_by_sample, breaks=247)
hist(missing_by_gene, breaks = 247)
