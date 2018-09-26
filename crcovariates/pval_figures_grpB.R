setwd("~/ferdig_rotation/gabes_CRC/")

#Read files, remove non-gene rows

uncorrected <- read.csv("GrpACorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
uncorrected <- uncorrected[-c(1,2),]

corrected <- read.csv("GrpAandBCorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
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

plot(both$FDR.y, both$FDR.x, xlab = "AB Corrected", ylab = "A Corrected")
#plot(both$FDR.x, both$FDR.y, xlim = c(0.05, 0), ylim = c(0.05,0), xlab = "Uncorrected", ylab = "Corrected")
plot(both$FDR.y, both$FDR.x, xlim = c(0.05, 0), ylim = c(0.05,0), xlab = "Corrected", ylab = "Uncorrected")


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

plot(genes_both$FDR.y, genes_both$FDR.x, xlab = "AB Corrected", ylab = "A Corrected")
abline(h=0.05, col="red")
abline(v=0.05, col="red")
#plot(genes_both$FDR.x, genes_both$FDR.y, xlim = c(0.05, 0), ylim = c(0.05,0), xlab = "Uncorrected", ylab = "Corrected")
plot(genes_both$FDR.y, genes_both$FDR.x, xlim = c(0.05, 0), ylim = c(0.05,0), xlab = "Corrected", ylab = "Uncorrected")


#lets add the middle section
full_gunc <- uncor_ord[which(uncor_ord$FDR < max(genes_uncor$FDR)),]
full_gc <- cor_fdr <- cor_ord[which(cor_ord$FDR < max(genes_cor$FDR)),]

full_both <- merge(full_gunc, full_gc, by="row.names")

plot(full_both$FDR.y, full_both$FDR.x, xlab = "AB Corrected", ylab = "A Corrected")
abline(h=0.05, col="red")
abline(v=0.05, col="red")

plot(full_both$FDR.y, full_both$FDR.x, xlim = c(max(genes_cor$FDR),0), ylim=c(max(genes_uncor$FDR),0), xlab = "Corrected", ylab = "Uncorrected")
abline(h=0.05, col="red")
abline(v=0.05, col="red")


#color the points by intersect, out, in
full_set <- union()


