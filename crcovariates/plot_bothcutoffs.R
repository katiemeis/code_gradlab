

############### pvalue plots - with both 0.01 pval cutoff and FDR cutoff ##################

setwd("~/ferdig_rotation/gabes_CRC/")

#Read files, remove non-gene rows

uncorrected <- read.csv("GrpBuncorrectedcorr.csv", stringsAsFactors = F, header = T, row.names = 1)
uncorrected <- uncorrected[-c(1,2),]


corrected <- read.csv("GrpBCorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
corrected <- corrected[-c(1,2,3,4),]

#Sort lists by FDR
uncor_ord <- uncorrected[order(uncorrected$FDR),]
cor_ord <- corrected[order(corrected$FDR),]

############################### plot anything that meets FDR 0.05 in BOTH ##########################

#want to plot anything with FDR 0.05 or better
#which(cor_ord$PValue > 0.01)
uncor_pval <- uncor_ord[which(uncor_ord$PValue < 0.01),]
cor_pval <- cor_ord[which(cor_ord$PValue < 0.01),]

uncor_fdr <- uncor_pval[which(uncor_pval$FDR < 0.05),]
cor_fdr <- cor_pval[which(cor_pval$FDR < 0.05),] 

####### plot anything that meets FDR 0.05 in BOTH
both <- merge(uncor_fdr, cor_fdr, by = "row.names")
#x in uncor, y is cor

plot(both$FDR.y, both$FDR.x, xlab = "Corrected", ylab = "Unorrected")


####### plot anything that meets FDR 0.05 in at least 1
#need union of 0.05 cutoff lists by rowname
genes_in_either <- union(row.names(uncor_fdr), row.names(cor_fdr))

#for each gene, pull it's uncorrected fdr value
genes_uncor <- uncor_ord[row.names(uncor_ord) %in% genes_in_either,]
max(genes_uncor$FDR)

#for each gene, pull it's corrected fdr value
genes_cor <- cor_ord[row.names(cor_ord) %in% genes_in_either,]
max(genes_cor$FDR)

genes_both <- merge(genes_uncor, genes_cor, by="row.names")

plot(genes_both$FDR.y, genes_both$FDR.x, xlab = "Corrected", ylab = "Uncorrected")
abline(h=0.05, col="red")
abline(v=0.05, col="red")


#lets add the middle section
full_gunc <- uncor_ord[which(uncor_ord$FDR < max(genes_uncor$FDR)),]
full_gc <- cor_fdr <- cor_ord[which(cor_ord$FDR < max(genes_cor$FDR)),]

full_both <- merge(full_gunc, full_gc, by="row.names")

plot(full_both$FDR.y, full_both$FDR.x, xlab = "Corrected", ylab = "Uncorrected")
abline(h=0.05, col="red")
abline(v=0.05, col="red")









#For this to work, you have to split by + and -

posuncorrect <- uncorrected[which(uncorrected$Slope > 0),]
neguncorrect <- uncorrected[which(uncorrected$Slope < 0),]

poscorrect <- corrected[which(corrected$Slope > 0),]
negcorrect <- corrected[which(corrected$Slope < 0),]

#now sort those, positive from low to high, negative from high to low

posuncorrect <- posuncorrect[order(posuncorrect$FDR),]
neguncorrect <- neguncorrect[order(-neguncorrect$FDR),]

poscorrect <- poscorrect[order(poscorrect$FDR),]
negcorrect <- negcorrect[order(-negcorrect$FDR),]

#stitch them back together

uncorrectedsorted <- rbind(posuncorrect, neguncorrect)
correctedsorted <- rbind(poscorrect, negcorrect)

#Add rank
for (num in 1:nrow(uncorrectedsorted)){
  uncorrectedsorted[num,5] <- num
}

for (num in 1:nrow(correctedsorted)){
  correctedsorted[num,5] <- num
}
#Rename rank columns
colnames(uncorrectedsorted)[5] <- "UncorRank"

colnames(correctedsorted)[5] <- "CorRank"

uncorrectedsorted <- uncorrectedsorted[which(uncorrectedsorted$FDR < 0.05),]

correctedsorted <- correctedsorted[which(correctedsorted$FDR < 0.05),]

#merge them for figures
combine <- merge(uncorrectedsorted, correctedsorted, by = "row.names", sort=FALSE)

######################################################################
#Figure Code Goes Here
######################################################################

plot(combine$UncorRank, combine$CorRank, xlab = "Uncorrected Rank", ylab = "Corrected Rank")

plot(combine$FDR.x, combine$FDR.y, xlab = "Uncorrected FDR", ylab = "Corrected FDR")

######################################################################

#Pull up/down reg genes
uncorsignificant <- uncorrected[which(uncorrected$PValue < 0.01),]
uncorsignificant <- uncorsignificant[which(uncorsignificant$FDR < 0.05),]

uncorupreg <- uncorsignificant[which(uncorsignificant$Correlation>0),]

uncordownreg <- uncorsignificant[which(uncorsignificant$Correlation<0),]

corsignificant <- corrected[which(corrected$PValue < 0.01),]
corsignificant <- corsignificant[which(corsignificant$FDR < 0.05),]

corupreg <- corsignificant[which(corsignificant$Correlation>0),]

cordownreg <- corsignificant[which(corsignificant$Correlation<0),]

