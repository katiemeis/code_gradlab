setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(nettools)

####Load in expression data used to build each network

#data with names parsed
Final_Data <- read.delim("Final_Data_avg.csv",header=TRUE,sep=",",as.is=TRUE) 

#pull out rows for NF54 and order by hpi
NF54_data <- Final_Data[which(Final_Data$pBLine == "NF54"),]
NF54_data <- NF54_data[order(NF54_data$hpi),]
NF54_data2 <- NF54_data[,-1]
rownames(NF54_data2) <- NF54_data[,1]
NF54_data <- NF54_data2
#new matrx without metadata
n <- length(NF54_data)
NF54_data_exp <- NF54_data[, 5:n]
NF54_wgcna <- t(NF54_data_exp)


#pull out rows for PB-58 and order by hpi
PB58_data <- Final_Data[which(Final_Data$pBLine == "PB-58"),]
PB58_data <- PB58_data[order(PB58_data$hpi),]
PB58_data2 <- PB58_data[,-1]
rownames(PB58_data2) <- PB58_data[,1]
PB58_data <- PB58_data2
#new matrx without metadata
n <- length(PB35_data)
PB58_data_exp <- PB58_data[, 5:n] 
PB58_wgcna <- t(PB58_data_exp)

#put into set for wgcna
library(WGCNA)
nSets = 2
setLabels = c("NF54", "PB58")
shortLabels = c("NF54", "PB58")
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(NF54_wgcna)));
names(multiExpr[[1]]$data) = row.names(NF54_wgcna);
rownames(multiExpr[[1]]$data) = names(NF54_wgcna);
multiExpr[[2]] = list(data = as.data.frame(t(PB58_wgcna)));
names(multiExpr[[2]]$data) = row.names(PB58_wgcna);
rownames(multiExpr[[2]]$data) = names(PB58_wgcna);

exprSize = checkSets(multiExpr)
exprSize



######## load in networks made using DTWMIC 
setwd('~/ferdig_rotation/pB_data/my_results/DTWMIC/wgcna/')

PB_mat <- read.csv("DTWMIC_PB58_results.csv", header=T)
test_pb <- PB_mat[,-1]
rownames(test_pb) <- PB_mat[,1]

NF_mat <- read.csv("DTWMIC_NF54_results.csv", header=T)
test_nf <- NF_mat[,-1]
rownames(test_nf) <- NF_mat[,1]

#thresholding
powers = c(seq(4,10,by=1), seq(12,20, by=2))
powerTables = vector(mode = "list", length = nSets)

powerTables[[1]] = list(data = pickSoftThreshold.fromSimilarity(as.matrix(test_nf), powerVector=powers, verbose = 2)[[2]])
powerTables[[2]] = list(data = pickSoftThreshold.fromSimilarity(as.matrix(test_pb), powerVector=powers, verbose = 2)[[2]])
#for (set in 1:nSets)
#  powerTables[[set]] = list(data = pickSoftThreshold.fromSimilarity(multiExpr[[set]]$data, powerVector=powers,
#                                                     verbose = 2)[[2]]);
collectGarbage();

colors = c("black", "red")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}

sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

softPower = 5
nGenes = 5540

#apply sotf threshold to DTWMIC networks
new_nf54_network <- matrix(nrow=nGenes, ncol=nGenes)
for(i in 1:nrow(test_nf)) {
  for(j in 1:ncol(test_nf)) {
    new_nf54_network[i,j] = test_nf[i,j]^6 
  }
}

new_pb58_network <- matrix(nrow=nGenes, ncol=nGenes)
for(i in 1:nrow(test_pb)) {
  for(j in 1:ncol(test_pb)) {
    new_pb58_network[i,j] = test_pb[i,j]^5 
  }
}

#put into array
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
adjacencies[1, , ] = new_nf54_network
adjacencies[2, , ] = new_pb58_network

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])

#########scaling TOM results
# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(1)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
sizeGrWindow(6,6)
#pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))

consensusTOM = pmin(TOM[1, , ], TOM[2, , ])
# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)

sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
######merging similar clusters
# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs
#plot unmerged and merged
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

write.csv(df_nfpb_clust, "wgcna_together/redone_samepower_5/cluster_assignments.csv")
df_nfpb_clust = cbind(row.names(NF54_wgcna), consMEs[[1]]$validColors, moduleColors)

####################################### part 4 Exploring results of network analysis






######################################## part 5
sizeGrWindow(8,10);
#pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(consMEs, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

save(list=ls(), file='wgcna_together/redone_samepower_5/wgcna_redo_925.RData')

