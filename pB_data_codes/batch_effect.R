source("http://bioconductor.org/biocLite.R")

library("sva")
library(limma)


setwd("C://Users//katie//Documents//ferdig_rotation//pB_data//Data//")
BCGene_data <- t(norm_Probe_data)


meta_data <- read.delim("Metadata2_katie.csv",header=TRUE,sep=",",as.is=TRUE,row.names=1)

BCGene_metadata <- merge(BCGene_data,meta_data,by="row.names")
BCGene_data <- BCGene_metadata[,2:62977]
rownames(BCGene_data) <- BCGene_metadata[,1]
rownames(BCGene_metadata) <- BCGene_metadata[,1]
BCGene_pca <- prcomp(BCGene_data,center=TRUE,scale=TRUE)
barplot(BCGene_pca$sdev^2,xlab="Eigenvalues",ylab="variation")
plot(BCGene_pca$x[,1],BCGene_pca$x[,2],main="PCA of BC Gene expression profiles",xlab="PC1",ylab="PC2")
points(BCGene_pca$x[which(BCGene_metadata[,"Time.point"]=="6HR"),1],BCGene_pca$x[which(BCGene_metadata[,"Time.point"]=="6HR"),2],col="red")
points(BCGene_pca$x[which(BCGene_metadata[,"Time.point"]=="24HR"),1],BCGene_pca$x[which(BCGene_metadata[,"Time.point"]=="24HR"),2],col="blue")


evdata <- t(BCGene_metadata[,2:62977])

mod <- model.matrix(~as.factor(pBLine)+as.factor(Replicate)+as.factor(hpi)+as.factor(hybe.batch)+as.factor(cDNA.batch)+as.factor(RNA.batch),data=BCGene_metadata[,62978:62985])
mod0 <- model.matrix(~as.factor(hybe.batch)+as.factor(cDNA.batch)+as.factor(RNA.batch),data=BCGene_metadata[,62978:62985])


n.sv <- num.sv(evdata,mod,method="leek")
svobj <- sva(evdata,mod,mod0)

pValues <- f.pvalue(evdata,mod,mod0)
qValues <- p.adjust(pValues,method="BH")
length(which(qValues < 0.05))

meta_dataSv <- cbind(BCGene_metadata[,62978:62985],svobj$sv)
colnames(meta_dataSv) <- c(colnames(BCGene_metadata[,62978:62985]),paste("sv",seq(1,11),sep="_"))
modSv <- cbind(mod,svobj$sv)
mod0Sv <- cbind(mod0,svobj$sv)
pValuesSv <- f.pvalue(evdata,modSv,mod0Sv)
qValuesSv <- p.adjust(pValuesSv,method="BH")
length(which(qValuesSv < 0.05))

fit <- lmFit(evdata,modSv)

### Batch effects for known batch variables

batch <- meta_dataSv[,6]
modcombat <- model.matrix(~1,data=meta_dataSv)

combat_evdata <- ComBat(dat=evdata,batch=batch,mod=modcombat,par.prior=TRUE,prior.plots=FALSE)

batch <- meta_dataSv[,7]
modcombat <- model.matrix(~1,data=meta_dataSv)

combat1_evdata <- ComBat(dat=combat_evdata,batch=batch,mod=modcombat,par.prior=TRUE,prior.plots=FALSE)

batch <- meta_dataSv[,8]
modcombat <- model.matrix(~1,data=meta_dataSv)

combat2_evdata <- ComBat(dat=combat1_evdata,batch=batch,mod=modcombat,par.prior=TRUE,prior.plots=FALSE)
colnames(combat2_evdata) <- 

pValuesCombat <- f.pvalue(combat_evdata,mod,mod0)
qValuesCombat <- p.adjust(pValuesCombat,method="BH")

length(which(qValuesCombat < 0.05))
combat_corrected_DE <- rownames(evdata)[(which(qValuesCombat < 0.05))]
write.csv(combat_corrected_DE,"combat_corrected_DE.csv")

### PCA of corrected data
Combat_pca <- prcomp(t(combat2_evdata),center=TRUE,scale=TRUE)
barplot(Combat_pca$sdev^2,xlab="Eigenvalues",ylab="variation")

### Plot PC of corrected data by stage
plot(Combat_pca$x[,1],Combat_pca$x[,2],main="PCA of Combat Gene expression profiles outliers removed",xlab="PC1",ylab="PC2")
points(Combat_pca$x[which(meta_dataSv[,5]==2),1],Combat_pca$x[which(meta_dataSv[,5]==2),2],col="deeppink")
points(Combat_pca$x[which(meta_dataSv[,5]==6),1],Combat_pca$x[which(meta_dataSv[,5]==6),2],col="red")
points(Combat_pca$x[which(meta_dataSv[,5]==10),1],Combat_pca$x[which(meta_dataSv[,5]==10),2],col="darkorange3")
points(Combat_pca$x[which(meta_dataSv[,5]==14),1],Combat_pca$x[which(meta_dataSv[,5]==14),2],col="orange")
points(Combat_pca$x[which(meta_dataSv[,5]==18),1],Combat_pca$x[which(meta_dataSv[,5]==18),2],col="yellow")
points(Combat_pca$x[which(meta_dataSv[,5]==22),1],Combat_pca$x[which(meta_dataSv[,5]==22),2],col="gold")
points(Combat_pca$x[which(meta_dataSv[,5]==26),1],Combat_pca$x[which(meta_dataSv[,5]==26),2],col="darkseagreen1")
points(Combat_pca$x[which(meta_dataSv[,5]==30),1],Combat_pca$x[which(meta_dataSv[,5]==30),2],col="green")
points(Combat_pca$x[which(meta_dataSv[,5]==34),1],Combat_pca$x[which(meta_dataSv[,5]==34),2],col="green4")
points(Combat_pca$x[which(meta_dataSv[,5]==38),1],Combat_pca$x[which(meta_dataSv[,5]==38),2],col="darkslategray1")
points(Combat_pca$x[which(meta_dataSv[,5]==42),1],Combat_pca$x[which(meta_dataSv[,5]==42),2],col="blue")
points(Combat_pca$x[which(meta_dataSv[,5]==46),1],Combat_pca$x[which(meta_dataSv[,5]==46),2],col="purple")

plot(Combat_pca$x[,2],Combat_pca$x[,3],main="PCA of Combat Gene expression profiles outliers removed",xlab="PC2",ylab="PC3")
points(Combat_pca$x[which(meta_dataSv[,5]==2),2],Combat_pca$x[which(meta_dataSv[,5]==2),3],col="deeppink")
points(Combat_pca$x[which(meta_dataSv[,5]==6),2],Combat_pca$x[which(meta_dataSv[,5]==6),3],col="red")
points(Combat_pca$x[which(meta_dataSv[,5]==10),2],Combat_pca$x[which(meta_dataSv[,5]==10),3],col="darkorange3")
points(Combat_pca$x[which(meta_dataSv[,5]==14),2],Combat_pca$x[which(meta_dataSv[,5]==14),3],col="orange")
points(Combat_pca$x[which(meta_dataSv[,5]==18),2],Combat_pca$x[which(meta_dataSv[,5]==18),3],col="yellow")
points(Combat_pca$x[which(meta_dataSv[,5]==22),2],Combat_pca$x[which(meta_dataSv[,5]==22),3],col="gold")
points(Combat_pca$x[which(meta_dataSv[,5]==26),2],Combat_pca$x[which(meta_dataSv[,5]==26),3],col="darkseagreen1")
points(Combat_pca$x[which(meta_dataSv[,5]==30),2],Combat_pca$x[which(meta_dataSv[,5]==30),3],col="green")
points(Combat_pca$x[which(meta_dataSv[,5]==34),2],Combat_pca$x[which(meta_dataSv[,5]==34),3],col="green4")
points(Combat_pca$x[which(meta_dataSv[,5]==38),2],Combat_pca$x[which(meta_dataSv[,5]==38),3],col="darkslategray1")
points(Combat_pca$x[which(meta_dataSv[,5]==42),2],Combat_pca$x[which(meta_dataSv[,5]==42),3],col="blue")
points(Combat_pca$x[which(meta_dataSv[,5]==46),2],Combat_pca$x[which(meta_dataSv[,5]==46),3],col="purple")

write.csv(combat2_evdata,"batch_corrected_quantile_normalized_probe_data_katie.csv")

(Combat_pca$sdev[1]^2+Combat_pca$sdev[2]^2+Combat_pca$sdev[3]^2+Combat_pca$sdev[4]^2)/sum(Combat_pca$sdev^2)






###Back to pipeline to get final data, then back here

png("PCA_batcheffects.png", height=10, width=10, units="in", res=220)
par(oma=c(0,0,3,0)) # Will use this space for the "overall" title
layout(matrix(c(1,2,3,4),2,2,byrow=FALSE)) #Setup up the layout of the figure


lA_pca <- prcomp(t(norm_Probe_data[,-1]),center=TRUE,scale=TRUE)
plot(lA_pca$x[,1],lA_pca$x[,2],main="Non batch corrected",xlab="PC1",ylab="PC2")
points(lA_pca$x[which(Filename_data[,6]==2),1],lA_pca$x[which(Filename_data[,6]==2),2],col="deeppink")
points(lA_pca$x[which(Filename_data[,6]==6),1],lA_pca$x[which(Filename_data[,6]==6),2],col="red")
points(lA_pca$x[which(Filename_data[,6]==10),1],lA_pca$x[which(Filename_data[,6]==10),2],col="darkorange3")
points(lA_pca$x[which(Filename_data[,6]==14),1],lA_pca$x[which(Filename_data[,6]==14),2],col="orange")
points(lA_pca$x[which(Filename_data[,6]==18),1],lA_pca$x[which(Filename_data[,6]==18),2],col="yellow")
points(lA_pca$x[which(Filename_data[,6]==22),1],lA_pca$x[which(Filename_data[,6]==22),2],col="gold")
points(lA_pca$x[which(Filename_data[,6]==26),1],lA_pca$x[which(Filename_data[,6]==26),2],col="darkseagreen1")
points(lA_pca$x[which(Filename_data[,6]==30),1],lA_pca$x[which(Filename_data[,6]==30),2],col="green")
points(lA_pca$x[which(Filename_data[,6]==34),1],lA_pca$x[which(Filename_data[,6]==34),2],col="green4")
points(lA_pca$x[which(Filename_data[,6]==38),1],lA_pca$x[which(Filename_data[,6]==38),2],col="darkslategray1")
points(lA_pca$x[which(Filename_data[,6]==42),1],lA_pca$x[which(Filename_data[,6]==42),2],col="blue")
points(lA_pca$x[which(Filename_data[,6]==46),1],lA_pca$x[which(Filename_data[,6]==46),2],col="purple")


plot(lA_pca$x[,3],lA_pca$x[,2],xlab="PC3",ylab="PC2")
points(lA_pca$x[which(Filename_data[,6]==2),3],lA_pca$x[which(Filename_data[,6]==2),2],col="deeppink")
points(lA_pca$x[which(Filename_data[,6]==6),3],lA_pca$x[which(Filename_data[,6]==6),2],col="red")
points(lA_pca$x[which(Filename_data[,6]==10),3],lA_pca$x[which(Filename_data[,6]==10),2],col="darkorange3")
points(lA_pca$x[which(Filename_data[,6]==14),3],lA_pca$x[which(Filename_data[,6]==14),2],col="orange")
points(lA_pca$x[which(Filename_data[,6]==18),3],lA_pca$x[which(Filename_data[,6]==18),2],col="yellow")
points(lA_pca$x[which(Filename_data[,6]==22),3],lA_pca$x[which(Filename_data[,6]==22),2],col="gold")
points(lA_pca$x[which(Filename_data[,6]==26),3],lA_pca$x[which(Filename_data[,6]==26),2],col="darkseagreen1")
points(lA_pca$x[which(Filename_data[,6]==30),3],lA_pca$x[which(Filename_data[,6]==30),2],col="green")
points(lA_pca$x[which(Filename_data[,6]==34),3],lA_pca$x[which(Filename_data[,6]==34),2],col="green4")
points(lA_pca$x[which(Filename_data[,6]==38),3],lA_pca$x[which(Filename_data[,6]==38),2],col="darkslategray1")
points(lA_pca$x[which(Filename_data[,6]==42),3],lA_pca$x[which(Filename_data[,6]==42),2],col="blue")
points(lA_pca$x[which(Filename_data[,6]==46),3],lA_pca$x[which(Filename_data[,6]==46),2],col="purple")

Combat_pca <- prcomp(t(Final_Data[,-1]),center=TRUE,scale=TRUE)
plot(Combat_pca$x[,1],Combat_pca$x[,2],main="Batch Corrected",xlab="PC1",ylab="PC2")
points(Combat_pca$x[which(meta_dataSv[,5]==2),1],Combat_pca$x[which(meta_dataSv[,5]==2),2],col="deeppink")
points(Combat_pca$x[which(meta_dataSv[,5]==6),1],Combat_pca$x[which(meta_dataSv[,5]==6),2],col="red")
points(Combat_pca$x[which(meta_dataSv[,5]==10),1],Combat_pca$x[which(meta_dataSv[,5]==10),2],col="darkorange3")
points(Combat_pca$x[which(meta_dataSv[,5]==14),1],Combat_pca$x[which(meta_dataSv[,5]==14),2],col="orange")
points(Combat_pca$x[which(meta_dataSv[,5]==18),1],Combat_pca$x[which(meta_dataSv[,5]==18),2],col="yellow")
points(Combat_pca$x[which(meta_dataSv[,5]==22),1],Combat_pca$x[which(meta_dataSv[,5]==22),2],col="gold")
points(Combat_pca$x[which(meta_dataSv[,5]==26),1],Combat_pca$x[which(meta_dataSv[,5]==26),2],col="darkseagreen1")
points(Combat_pca$x[which(meta_dataSv[,5]==30),1],Combat_pca$x[which(meta_dataSv[,5]==30),2],col="green")
points(Combat_pca$x[which(meta_dataSv[,5]==34),1],Combat_pca$x[which(meta_dataSv[,5]==34),2],col="green4")
points(Combat_pca$x[which(meta_dataSv[,5]==38),1],Combat_pca$x[which(meta_dataSv[,5]==38),2],col="darkslategray1")
points(Combat_pca$x[which(meta_dataSv[,5]==42),1],Combat_pca$x[which(meta_dataSv[,5]==42),2],col="blue")
points(Combat_pca$x[which(meta_dataSv[,5]==46),1],Combat_pca$x[which(meta_dataSv[,5]==46),2],col="purple")

plot(Combat_pca$x[,2],Combat_pca$x[,3],xlab="PC2",ylab="PC3")
points(Combat_pca$x[which(meta_dataSv[,5]==2),2],Combat_pca$x[which(meta_dataSv[,5]==2),3],col="deeppink")
points(Combat_pca$x[which(meta_dataSv[,5]==6),2],Combat_pca$x[which(meta_dataSv[,5]==6),3],col="red")
points(Combat_pca$x[which(meta_dataSv[,5]==10),2],Combat_pca$x[which(meta_dataSv[,5]==10),3],col="darkorange3")
points(Combat_pca$x[which(meta_dataSv[,5]==14),2],Combat_pca$x[which(meta_dataSv[,5]==14),3],col="orange")
points(Combat_pca$x[which(meta_dataSv[,5]==18),2],Combat_pca$x[which(meta_dataSv[,5]==18),3],col="yellow")
points(Combat_pca$x[which(meta_dataSv[,5]==22),2],Combat_pca$x[which(meta_dataSv[,5]==22),3],col="gold")
points(Combat_pca$x[which(meta_dataSv[,5]==26),2],Combat_pca$x[which(meta_dataSv[,5]==26),3],col="darkseagreen1")
points(Combat_pca$x[which(meta_dataSv[,5]==30),2],Combat_pca$x[which(meta_dataSv[,5]==30),3],col="green")
points(Combat_pca$x[which(meta_dataSv[,5]==34),2],Combat_pca$x[which(meta_dataSv[,5]==34),3],col="green4")
points(Combat_pca$x[which(meta_dataSv[,5]==38),2],Combat_pca$x[which(meta_dataSv[,5]==38),3],col="darkslategray1")
points(Combat_pca$x[which(meta_dataSv[,5]==42),2],Combat_pca$x[which(meta_dataSv[,5]==42),3],col="blue")
points(Combat_pca$x[which(meta_dataSv[,5]==46),2],Combat_pca$x[which(meta_dataSv[,5]==46),3],col="purple")

mtext("Principal Component Analysis and Batch Effect Correction",outer=TRUE,cex=1.4)

dev.off()
