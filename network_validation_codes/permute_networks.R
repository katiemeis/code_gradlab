
### 1. Code for creating permuted data from the Hu dataset to use for permutation tests ###
### 1. Outputs 1000 files with permuted expression profiles in format for GENIE3 (space delimited, no sample names) - input into network construction methods for null distributions ### 
### 2. Does 10 permutations to give smaller null ditribution - necessary as files too large to save all data ###
### 3. Does permutation tests for Pearson network to set threshold ###
### 4. Does perm tests in a way similar to ARACNE ###
### 5. Same as section 4 but with different sampling to test for rf ###


####################### 1. Permute the Hu dataset and print permuted files #############################
### load in Hu dataset
setwd("~/ferdig_rotation/regulon_validation")
hu_data <- read.csv('Hu_data.csv', header=TRUE)
hu_data_t <- t(hu_data)
colnames(hu_data_t) = hu_data_t[1,]
hu_data_t = hu_data_t[-1,]
class(hu_data_t) <- "numeric"

### making files for permutation tests
for(jj in 1:1000){
  hu_perm_mat <- hu_data_t
  for(ii in 1:ncol(hu_data_t)){
    # for each column, swap all the values
    hu_perm_mat[,ii] <- sample(hu_data_t[,ii])
  }
  write.table(hu_perm_mat, file=paste("~/ferdig_rotation/regulon_validation/permutation_networks/hu_perm_", jj, ".txt", sep=""), row.names=FALSE, quote=FALSE, sep=" ")
}



################# 2. Running permutations for idea of null distribution - Pearson ##################

library(nettools)
library(igraph)
### to read back in for pearson perm testing
#test_space <- read.delim("~/ferdig_rotation/regulon_validation/permutation_networks/hu_perm_1.txt", sep=' ')

#empty list for weights of perm tests
null_dist_values <- list()

#do only 10 to get a good idea of what the distribution looks like
for(ii in 1:10){
  #read in permuted file
  perm_mat <- read.delim(paste("~/ferdig_rotation/regulon_validation/permutation_networks/hu_perm_", ii, ".txt", sep=""), sep=' ')
  #calculate pearson correlations
  hu_pearson_result <- mat2adj(perm_mat, method="cor")
  #matrix to edgelist
  perm_edges <- graph.adjacency(hu_pearson_result, weighted=TRUE, mode=c("upper"))
  perm_edge_df <- get.data.frame(perm_edges)
  # add weights to list
  null_dist_values[[ii]] <- perm_edge_df$weight
}

#convert list to matrix
edges_mat <- do.call(cbind, null_dist_values)

#write.csv(edges_mat, file="~/ferdig_rotation/regulon_validation/perm_weight_distribution_10.csv")

#plot distribution
hist(edges_mat)
sig_vals_cutoff <- unname(quantile(edges_mat, 0.95)) #0.124897
abline(v=sig_vals_cutoff)


############### 3. Run actual permutation test to get distribution ######################
count <- 0
sig_vec <- vector()

#find 5%, if corr val is above this then keep, if below then just count
for(ii in 1:1000){
  #read in permuted file
  perm_mat <- read.delim(paste("~/ferdig_rotation/regulon_validation/permutation_networks/hu_perm_", ii, ".txt", sep=""), sep=' ')
  #calculate pearson correlations
  hu_pearson_result <- mat2adj(perm_mat, method="cor")
  #matrix to edgelist
  perm_edges <- graph.adjacency(hu_pearson_result, weighted=TRUE, mode=c("upper"))
  perm_edge_df <- get.data.frame(perm_edges)
  #keep weights if greater than cutoff, just add to count if not
  for(weight in perm_edge_df$weight){
    if(weight >= 0.124897){
      sig_vec <- c(sig_vec, weight)
    }else{
      count <- count + 1
    }
  } 
}

print(length(sig_vec))
print(count)



################# 4. Perm tests comparable to ARACNE ##############################
#set working directory and load libraries
setwd("~/ferdig_rotation/regulon_validation")
library(nettools)
library(igraph)

#load in data, transpose, make numeric
hu_data <- read.csv('Hu_data.csv', header=TRUE)
hu_data_t <- t(hu_data)
colnames(hu_data_t) = hu_data_t[1,]
hu_data_t = hu_data_t[-1,]
class(hu_data_t) <- "numeric"

#set seed for reproducibility and empty vector for null values
set.seed(0)
perm_null_results <- vector()

for(ii in 1:100000){
  #choose 2 random genes and extract their profiles
  random_genes <- sample(1:3705, 2, replace=F)
  exp_profiles <- hu_data_t[,random_genes]
  #permute profiles separately
  exp_perm <- matrix(0, nrow=nrow(exp_profiles), ncol=ncol(exp_profiles))
  for(jj in 1:ncol(exp_profiles)){
    # for each column, scramble the values
    exp_perm[,jj] <- sample(exp_profiles[,jj])
  }
  #calculate pearson correlations, input rows samples and columns genes
  perm_pearson_result <- mat2adj(exp_perm, method="cor")
  #save single value to vector, can keep top right or bottom left, same
  perm_null_results[ii] <- perm_pearson_result[1,2]
}

#write.csv(perm_null_results, file="sampling_null/sampling_op1_million_seed0.csv")

#graph null distribution and density line
hist(perm_null_results, prob=T)
lines(density(perm_null_results))

#find 10^-5 threshold
test_quant <- unname(quantile(perm_null_results, 0.9999)) #0.248537601774596

#apply threshold
pearson_mat <- read.csv("hu_pearson_network_mat.csv")
rownames(pearson_mat) <- pearson_mat[,1]
pearson_mat <- pearson_mat[,-1]
pearson_mat1 <- as.matrix(pearson_mat)
pear_edges <- graph.adjacency(pearson_mat1, weight=T, mode=c("upper"))
edges_df <- get.data.frame(pear_edges)
length(edges_df$weight[edges_df$weight > .248537601774596]) #gives 3,462,979 edges
ordered_pearson_edges <- edges_df[order(edges_df$weight),]
pearson_thresholded <- ordered_pearson_edges[which(ordered_pearson_edges$weight > .248537601774596),]
write.csv(pearson_thresholded, "pearson_thresholded.csv", row.names = FALSE, col.names = FALSE, quote=FALSE)





################# 5. Perm tests comparable to ARACNE with new sampling ##############################
#### want to take 10,000 from each of 100 networks
#### generate sampling without replacement within networks, but not between

#set working directory and load libraries
setwd("~/ferdig_rotation/regulon_validation")
library(nettools)
library(igraph)

#load in data, transpose, make numeric
hu_data <- read.csv('Hu_data.csv', header=TRUE)
hu_data_t <- t(hu_data)
colnames(hu_data_t) = hu_data_t[1,]
hu_data_t = hu_data_t[-1,]
class(hu_data_t) <- "numeric"

#set seed for reproducibility and empty vector for null values
set.seed(0)
perm_null_results <- vector()

for(ii in 1:54){
  #read in permuted file
  perm_mat <- read.delim(paste("~/ferdig_rotation/regulon_validation/permutation_networks/hu_perm_", ii, ".txt", sep=""), sep=' ')
  #calculate pearson correlations
  hu_pearson_result <- mat2adj(perm_mat, method="cor")
  #matrix to edgelist
  #perm_edges <- graph.adjacency(hu_pearson_result, weighted=TRUE)
  #perm_edge_df <- get.data.frame(perm_edges)
  ############change to choosing 1852 random pairs w/o replacement and saving those values
  random_genes <- sample(1:3705, 3704, replace=F)
  rand_mat <- matrix(random_genes, ncol=2, nrow=1852)
  #take pairs and pull out the weight associated with the pair (pair is row and column)
  for(jj in 1:1852){
    #take each pair, use row and column to find weight and add to list
    rand_row <- rand_mat[jj,1]
    rand_column <- rand_mat[jj,2]
    perm_null_results <- c(perm_null_results, hu_pearson_result[rand_row,rand_column])
  }

}

#write.csv(perm_null_results, file="sampling_null_op2_seed0.csv")




#graph null distribution and density line
hist(perm_null_results, prob=T)
lines(density(perm_null_results))

#find 10^-5 threshold
test_quant <- unname(quantile(perm_null_results, 0.9999))

#qq plot
sample_1 <- read.csv("sampling_null/sampling_null_op1_ns.csv")
sample_2 <- read.csv("sampling_null/sampling_null_op2_ns.csv")
ordered_s1 <- sample_1[order(sample_1$x),]
orderes_s2 <- sample_2[order(sample_2$x),]
ordered_s2 <- orderes_s2[9:100008,]
plot(ordered_s1$x, ordered_s2$x)
lines(x=c(0,0.3),y=c(0,0.3))

pearson_mat <- read.csv("hu_pearson_network_mat.csv")
rownames(pearson_mat) <- pearson_mat[,1]
pearson_mat <- pearson_mat[,-1]
pearson_mat1 <- as.matrix(pearson_mat)
pear_edges <- graph.adjacency(pearson_mat1, weight=T)
edges_df <- get.data.frame(pear_edges)
length(edges_df$weight[edges_df$weight > .248537601774596])
ordered_pearson_edges <- edges_df[order(edges_df$weight),]


###read in rf file and plot
rf_null_df <- read.csv("rf_full_perms_vals_1col.txt", header = FALSE, sep = '\t') 
rf_null_size <- rf_null_df[1:1000000,]
test_quant <- unname(quantile(rf_null_size, 0.9999)) #0.0041739134301528


#make thresholded network
rf_net <- read.csv("../GENIE3_python/hu_randomForest_edgelist.txt", header = FALSE, sep = '\t')
ordered_rf_edges <- rf_net[order(rf_net$V3),]
rf_thresholded <- rf_net[which(rf_net$V3 > 0.0041739134301528),] #gives 7,247,795 edges
write.csv(rf_thresholded, "../rf_thresholded.csv", row.names = FALSE, col.names = FALSE, quote=FALSE)

###make null histograms and fit exponential to extrapolate lower p-vals 
rf_data <- hist(rf_null_size, breaks = 500, main = "Random Forest Null Histogram")
lines(density(rf_null_size))
abline(v=test_quant)

pearson_null <- read.csv("sampling_null/sampling_op1_million_seed0.csv")
pear_test_quant <- unname(quantile(pearson_null$x, 0.9999))
pear_data <- hist(pearson_null$x, breaks = 500,  main = "Pearson Null Histogram")
lines(density(pearson_null$x))
abline(v=pear_test_quant)


plot(pear_data$breaks[1:753], pear_data$counts)
plot(rf_data$breaks[1:913], rf_data$counts)

#fit exponential decay y = e^(a*xvalues)
rf_df <- data.frame(rf_data$breaks[1:913], rf_data$counts)
fit_lines <- nls(counts ~ exp(-b0*breaks), data=rf_df, start=list(b0=1))




