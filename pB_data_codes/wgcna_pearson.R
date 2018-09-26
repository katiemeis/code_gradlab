########################################################################
# Running WGCNA analysis on DTWMIC results 
########################################################################
setwd("~/ferdig_rotation/pB_data/my_results/Pipeline")
library(nettools)

####Load in expression data used to build each network

#data with names parsed
Final_Data <- read.delim("Final_Data_avg.csv",header=TRUE,sep=",",as.is=TRUE) 
options(stringsAsFactors = FALSE)
#pull out rows for NF54 and order by hpi
NF54_data <- Final_Data[which(Final_Data$pBLine == "NF54"),]
NF54_data <- NF54_data[order(NF54_data$hpi),]
NF54_data2 <- NF54_data[,-1]
rownames(NF54_data2) <- NF54_data[,1]
NF54_data <- NF54_data2
#new matrx without metadata
n <- length(NF54_data)
NF54_data_exp <- NF54_data[, 5:n]


#pull out rows for PB-58 and order by hpi
PB58_data <- Final_Data[which(Final_Data$pBLine == "PB-58"),]
PB58_data <- PB58_data[order(PB58_data$hpi),]
PB58_data2 <- PB58_data[,-1]
rownames(PB58_data2) <- PB58_data[,1]
PB58_data <- PB58_data2
#new matrx without metadata
n <- length(PB35_data)
PB58_data_exp <- PB58_data[, 5:n] 



library(WGCNA)
setwd('~/ferdig_rotation/pB_data/my_results/DTWMIC/wgcna/')


###########################  PB58  ######################################

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(PB58_data_exp, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")
 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
-sign(sft$fitIndices[,3])*sft$fitIndices[,2]

softPower = 8
adjacency = adjacency(PB58_data_exp, power = softPower)
TOM = TOMsimilarity(adjacency)

dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
minModuleSize =30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize);

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(PB58_data_exp, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(PB58_data_exp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#dev.off()
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
#save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_individual/PB58_results/wcgna_cluster_pb.RData")

cluster_df <- data.frame(colnames(PB58_data_exp), moduleLabels, moduleColors)

write.csv(cluster_df, file='wgcna_pearson/clustering_pb58.csv', row.names=F)

ordered_clusters_pb <- cluster_df[with(cluster_df, order(moduleLabels)),]
write.csv(ordered_clusters_pb, file='wgcna_pearson/ordered_clusters_pb.csv', row.names=F)

which(ordered_clusters_pb$colnames.PB58_data_exp. =='PF3D7_1343700')
ordered_clusters_pb[4449,] #darkturquoise cluster
K13_cluster <-ordered_clusters_pb[which(ordered_clusters_pb$moduleColors=='darkturquoise'),]
write.csv(K13_cluster, file='wgcna_pearson/K13_cluster_pb.csv', row.names=F)

par(cex=0.9)
plotEigengeneNetworks(MEs, setLabels = c("PB58"), marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

PBmoduleColors = moduleColors
PBMEs = MEs
#save(list=ls(), file='wgcna_pearson/wgcna_pb.RData')
plotTOM = dissTOM^7
diag(plotTOM) = NA
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot PB58, TOM")
#TOMplot(adjacency, geneTree, moduleColors, main = "Network heatmap plot PB58, adjacency")

#conpare clusters between individual PB58 and consensus
load('wgcna_PB58_p5/wgcna_pb_alone.RData')
load("wgcna_together/redone_diff_powers/wgcna_redo_diff_925.RData")
consModuleColors = moduleColors
PBaloneModuleLabels = substring(names(MEs), 3)
consModuleLabels = substring(names(consMEs[[2]]$data), 3)
consModules = labels2colors(as.numeric(consModuleLabels))

nPBMods = length(PBaloneModuleLabels)
nConsMods = length(consModules)
pTable = matrix(0, nrow = nPBMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nPBMods, ncol = nConsMods);


for (fmod in 1:nPBMods)
  for (cmod in 1:nConsMods)
  {
    PBMembers = (PBmoduleColors == PBaloneModuleLabels[fmod]);
    consMembers = (consModuleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(PBMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(PBmoduleColors == PBaloneModuleLabels[fmod] & consModuleColors ==
                                 consModules[cmod])
  }

pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50
PBModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)

sizeGrWindow(10,7 )
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(10, 11, 3, 3)+0.3)

labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", PBaloneModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("PB58 ", PBaloneModuleLabels, ": ", PBModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of PB58 set-specific and NF54-PB58 consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);






###############################  NF54  ##################################################

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(NF54_data_exp, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
-sign(sft$fitIndices[,3])*sft$fitIndices[,2]

softPower = 8
adjacency = adjacency(NF54_data_exp, power = softPower)
TOM = TOMsimilarity(adjacency)

dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
minModuleSize =30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize);

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

MEList = moduleEigengenes(NF54_data_exp, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(NF54_data_exp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
#dev.off()
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
#save(MEs, moduleLabels, moduleColors, geneTree, file = "wcgna_cluster_nf.RData")

cluster_df <- data.frame(colnames(NF54_data_exp), moduleLabels, moduleColors)
write.csv(cluster_df, file='wgcna_pearson/clustering_nf54.csv', row.names=F)

ordered_clusters_pb <- cluster_df[with(cluster_df, order(moduleLabels)),]
write.csv(ordered_clusters_pb, file='wgcna_pearson/ordered_clusters_nf.csv', row.names=F)

which(ordered_clusters_pb$colnames.NF54_data_exp. =='PF3D7_1343700')
ordered_clusters_pb[576,] #black cluster
K13_cluster <-ordered_clusters_pb[which(ordered_clusters_pb$moduleColors=='black'),]
write.csv(K13_cluster, file='wgcna_pearson/K13_cluster_nf.csv', row.names=F)

par(cex=0.9)
plotEigengeneNetworks(MEs, setLabels = c("NF54"), marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

NFmoduleColors = moduleColors
NFMEs = MEs
#save(list=ls(), file='wgcna_pearson/wgcna_nf.RData')
plotTOM = dissTOM^12
diag(plotTOM) = NA
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot NF54, TOM")



#compare NF54 to consensus
load("wgcna_NF54_p6/wgcna_p6_nf.RData")
load("wgcna_together/redone_samepower_5/wgcna_redo_925.RData")
load("wgcna_together/redone_diff_powers/wgcna_redo_diff_925.RData")
consModuleColors = moduleColors
NFaloneModuleLabels = substring(names(NFMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
consModules = labels2colors(as.numeric(consModuleLabels))

nNFMods = length(NFaloneModuleLabels)
nConsMods = length(consModules)
pTable = matrix(0, nrow = nNFMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nNFMods, ncol = nConsMods);


for (fmod in 1:nNFMods)
  for (cmod in 1:nConsMods)
  {
    NFMembers = (NFmoduleColors == NFaloneModuleLabels[fmod]);
    consMembers = (consModuleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(NFMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(NFmoduleColors == NFaloneModuleLabels[fmod] & consModuleColors ==
                                 consModules[cmod])
  }

pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50
NFModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)

sizeGrWindow(10,7 )
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(10, 11, 4, 1)+0.3)

labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", NFaloneModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("NF54 ", NFaloneModuleLabels, ": ", NFModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of NF54 set-specific and NF54-PB58 consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);


#compare NF54 to PB58
load('wgcna_pearson/wgcna_pb.RData')
load("wgcna_NF54_p6/wgcna_nf.RData")
NFaloneModuleLabels = substring(names(NFMEs), 3)
PBaloneModuleLabels = substring(names(PBMEs), 3)

nNFMods = length(NFaloneModuleLabels)
nPBMods = length(PBaloneModuleLabels)
pTable = matrix(0, nrow = nNFMods, ncol = nPBMods);
CountTbl = matrix(0, nrow = nNFMods, ncol = nPBMods);


for (fmod in 1:nNFMods)
  for (cmod in 1:nPBMods)
  {
    NFMembers = (NFmoduleColors == NFaloneModuleLabels[fmod]);
    PBMembers = (PBmoduleColors == PBaloneModuleLabels[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(NFMembers, PBMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(NFmoduleColors == NFaloneModuleLabels[fmod] & PBmoduleColors ==
                                 PBaloneModuleLabels[cmod])
  }

pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50
NFModTotals = apply(CountTbl, 1, sum)
PBModTotals = apply(CountTbl, 2, sum)

sizeGrWindow(12,7 )
#pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1))
par(cex = .8)
par(mar=c(10, 12, 3, .5))

labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", PBaloneModuleLabels),
               yLabels = paste(" ", NFaloneModuleLabels),
               colorLabels = TRUE,
               xSymbols = paste("PB58 ", PBaloneModuleLabels, ": ", PBModTotals, sep=""),
               ySymbols = paste("NF54 ", NFaloneModuleLabels, ": ", NFModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of NF54 and PB58 modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);



