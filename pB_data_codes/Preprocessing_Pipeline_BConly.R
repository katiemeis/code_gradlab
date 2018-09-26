### Set working directory -- START HERE!!!
setwd("C://Users//katie//Documents//ferdig_rotation//pB_data//Data//")

### attempt to load packages, may need to install if not already installed
library(reshape)
library(preprocessCore)

### Get list of txt files in working directory and then remove .txt extension
FileNames <- list.files(pattern = "*.txt")
name <- sub(".txt", "", FileNames)
### Output filenames without .txt to create a file mapping filename to sample - use text to columns in excel and then 
write.csv(FileNames,"Sample_Names2_katie.csv")
### Read back in the file that maps file name to sample ID
Filename_data <- read.delim("Sample_Names2.csv",sep=",",header=FALSE,as.is=TRUE)
name <- Filename_data[,2]

### Read in probe files
# Probe.csv maps ProbeName to gene and exon identifiers
probes <- read.csv("Probe.csv",header=TRUE)
Probe_data <- matrix(probes$ProbeName,ncol=1,nrow=length(probes[,1]))
colnames(Probe_data) <- "ProbeName"

# Read in expression probes (have identifier 0 in design file)
ExpProbes <- read.csv("ExpProbes.csv",header=FALSE)

# Read in all probes on array (not parsed by expression and control probes)
All.probes <- read.csv("AllProbes.csv",header=TRUE)

# Read in control probes, these have ControlType identifier 1 
Control.probes <- read.csv("Control_Probes.csv",header=FALSE)

# Read in all probe indentifiers
All_Probe_data <- matrix(ncol=length(FileNames)+1,nrow=62976)
rownames(All_Probe_data) <- All.probes[,2]
colnames(All_Probe_data) <- c(t(name), "ControlType")


### Read in all data to RMA all probes together ###
# For all txt files in the data folder get gProcessedSignal column and add to All_Probe_data
# All output from Agilent feature extraction files are identical for a given array design, so can just fill in matrix, no need to merge
for(j in 1:length(FileNames)) {
  data <- read.delim(file=FileNames[j], header=TRUE, sep="\t",skip=9)
  All_Probe_data[,j] <- data[,"gProcessedSignal"]
}
# Pull controlType information from last file
All_Probe_data[,length(FileNames)+1] <- data[,"ControlType"]

### Save unnormalized probe summary file
write.csv(All_Probe_data,"unnormalized_raw_transcript_expression_final_katie.csv")

### Perform quantile normalization using BioConductor package preprocessCore
norm_Probe_data <- normalize.quantiles(All_Probe_data[,-(length(FileNames)+1)])
colnames(norm_Probe_data) <- t(name)
rownames(norm_Probe_data) <- rownames(All_Probe_data)
# Save quantile normalized probe file
write.csv(norm_Probe_data,"quantilenormalized_raw_transcript_expression_final_katie.csv")




### Batch Effect Section - Need to add this
#didn't work for me, missing Metadata2.csv file




### Create a new data.frame to hold quantile normalized results
Probe_Data_Full <- data.frame(combat2_evdata)
colnames(Probe_Data_Full) <- colnames(combat2_evdata)
# Add a column for ProbeName
Probe_Data_Full$ProbeName <- All.probes[,1]


# Merge ControlType info into Probe_Data_Full
Probe_Data_Full <- merge(Probe_Data_Full,All_Probe_data[,"ControlType"],by="row.names") #ControlTypes is row names
rownames(Probe_Data_Full) <-Probe_Data_Full[,"Row.names"]
Probe_Data_Full <- Probe_Data_Full[,-1]
colnames(Probe_Data_Full)[length(Probe_Data_Full[1,])] <- "ControlType"

### Find control probes ##
Control_Probes <- Probe_Data_Full[Probe_Data_Full$ControlType==1,]

### Subset out expression probes ##
Expression_Probes <- Probe_Data_Full[Probe_Data_Full$ControlType==0,]
Expression_Probes$ProbeID <- rownames(Expression_Probes)

### Need Gene and Exon info in the Probe file
Probe_data <- merge(Expression_Probes,probes,by="ProbeName")

### Count up probes for each Gene_Exon identifier and save as Replicate column in Probe_data
list <- unique(Probe_data[,"Gene_Exon"])
for (i in 1:length(list)){
  ind_gene <- which(Probe_data[,"Gene_Exon"]==list[i])
  Probe_data[ind_gene,"Replicate"] <- seq(1:length(ind_gene))
}

### Create Final Data Matrix to merge data into
Final_Data <- data.frame(unique(Probe_data[,"GeneID"]))
colnames(Final_Data) <- c("GeneID")  

###  Loop through samples and summarize data for exon and gene
for(j in 1:length(FileNames)) {
  # Create a new dataset with ProbeName, j+1th column, ProbeID, GeneID, ExonID and Replicate
  Probe_data_sample <- Probe_data[,c("ProbeName",colnames(Probe_data)[j+1],"ProbeID","GeneID","ExonID","Replicate","Gene_Exon")]
  # Use recast function to get orientation from long format to wide format so we can use rowMeans
  ProbeData_WF <-recast(Probe_data_sample,GeneID+ExonID~Replicate,measure.var=colnames(Probe_data)[j+1],id.var=c("GeneID","ExonID","Replicate"))
  
  # Create an Exons dataset to fill in probe counts, Exon averages and SD
  Exons <- ProbeData_WF[,1:5]
  for (i in 1:length(Exons[,3])){
    Exons[i,3] <- length(which(is.na(ProbeData_WF[i,3:54])==FALSE))
    temp <- t(ProbeData_WF[i,3:(3+Exons[i,3]-1)])
    if(Exons[i,3]>1) {
      Exons[i,5] <- apply(temp,2,sd,na.rm=TRUE)
    }
    else {
      Exons[i,5] <- 0
    }
  }
  Exons[,4] <- rowMeans(ProbeData_WF[,3:54],na.rm=TRUE)
  colnames(Exons) <- c("GeneID","ExonID","Num_Probes","Mean","SD")
  
  # If we are on first sample, don't merge just create Exon_Data, otherwise for subsequent samples merge into Exon_Data
  if(j==1){
    Exon_Data <- Exons[,c(1,2,4)]
  }
  else{
    Exon_Data <- merge(Exon_Data,Exons[,c(1,2,4)],by=c("GeneID","ExonID"))
  }
  
  # Now summarize at gene level
  Test <- as.data.frame(Exons[,c("GeneID","ExonID","Mean")])  #gets rid of extra information stored during first recast
  # Get unique GeneID values from Exons file
  list_gene <- unique(Test[,"GeneID"])
  # Get number of Exons per GeneID
  for (i in 1:length(list_gene)){
    ind_gene <- which(Test[,"GeneID"]==list_gene[i])
    Test[ind_gene,"ExonID"] <- seq(1:length(ind_gene))
  }
  #Test[which(is.na(Test[,"Mean"])==TRUE),"Mean"] <- 0
  
  # Recast data in wide format to take rowMeans
  Exon_WF <- recast(Test,GeneID~ExonID,measure.var="Mean",id.var=c("GeneID","ExonID"))
  
  # Create Genes file and calculate gene summaries
  Genes <- Exon_WF[,1:4]
  for (i in 1:length(Genes[,2])){
    # Get number of exons
    Genes[i,2] <- length(which(is.na(Exon_WF[i,2:32])==FALSE))
    # store exon values in temp
    temp <- t(Exon_WF[i,2:(2+Genes[i,2]-1)])
    # calculate SD
    if(Genes[i,2]>1) {
      Genes[i,4] <- apply(temp,2,sd,na.rm=TRUE)
    }
    else {
      Genes[i,4] <- 0
    }
  }
  # Calculate rowmeans to get average Gene level summary
  Genes[,3] <- rowMeans(Exon_WF,na.rm=TRUE)
  colnames(Genes) <- c("GeneID","Num_Exons","Mean","SD")
  
  # Output Exon and Gene files for each entry in name
  write.csv(Exons,paste(colnames(Probe_data)[j+1],"Exon_FDR.csv",sep="_"))
  write.csv(Genes,paste(colnames(Probe_data)[j+1],"Gene_FDR.csv",sep="_"))
  
  # merge Gene info into Final_Data file
  Final_Data<- merge(Final_Data,Genes[,c("GeneID","Mean")],by="GeneID")
}

colnames(Final_Data) <- c("GeneID",colnames(Probe_data[,2:147]))
colnames(Exon_Data) <- c("GeneID","ExonID",colnames(Probe_data[,2:147]))

write.csv(Exon_Data,"Exon_Data_BC_all_katie.csv")  ##same as Katie's
write.csv(Final_Data,"Gene_Data_BC_all_katie.csv") ##same as Katie's

BCGene_pca <- prcomp(t(Final_Data[,-1]),center=TRUE,scale=TRUE)
barplot(BCGene_pca$sdev^2,xlab="Eigenvalues",ylab="variation")
plot(BCGene_pca$x[,1],BCGene_pca$x[,2],main="PCA of BC Gene expression profiles outliers removed",xlab="PC1",ylab="PC2")
(BCGene_pca$sdev[1]^2+BCGene_pca$sdev[2]^2+BCGene_pca$sdev[3]^2+BCGene_pca$sdev[4]^2)/sum(BCGene_pca$sdev^2)

