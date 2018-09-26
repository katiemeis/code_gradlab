
########################## adalasso, MI, consensus #######################################

#plot the precision/recall
setwd("~/ferdig_rotation/regulon_validation/original_nets/CRC_MCL_clustering/")
pr_res <- read.csv("best_LOO_precision_recall.csv", header = F, as.is = T)

cols_vec <- c("darkcyan", "cyan3", "lightseagreen",
              "darkolivegreen", "darkolivegreen4", "darkolivegreen3",
              "orangered3", "orangered", "orange2", "orange")


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

n_go_terms <- c(41,41,31,14,30,7,13,28,28,12)/10
n_genes <- c(615,436,151,3563,2745,436,2211,910,663,96)/200

plot(pr_res$V2, pr_res$V3, xlab="Precision", ylab="Recall")
points(pr_res$V2, pr_res$V3, col = cols_vec, pch=20, cex=n_genes)
legend("topright", inset=c(-0.2,0), pr_res$V1, pch=20, col = cols_vec)


#plot fscores of top best
fscore <- 2*(pr_res$V2*pr_res$V3)/(pr_res$V2+pr_res$V3)
barplot(fscore, col = cols_vec, ylab = "LOO Fscore", names.arg = pr_res$V1, las=2)

#AUC
library(pracma)
loo_al <- read.table("al_LOOprecrec.txt", sep="\t", header=T)
loo_reg <- read.table("reg_LOOprecrec.txt", sep="\t", header=T)
loo_con <- read.table("con_LOOprecrec.txt", sep="\t", header=T)
loo_uni <- read.table("union_LOOprecrec.txt", sep="\t", header=T)

plot(loo_al$LOO.prec, loo_al$LOO.recall, type="l")
plot(loo_reg$LOO.prec, loo_reg$LOO.recall, type="l")
plot(loo_con$LOO.prec, loo_con$LOO.recall, type="l")

ord_al <- loo_al[order(loo_al$LOO.prec),]
ord_reg <- loo_reg[order(loo_reg$LOO.prec),]
ord_con <- loo_con[order(loo_con$LOO.prec),]
ord_uni <- loo_uni[order(loo_uni$LOO.prec),]

plot(ord_al$LOO.prec, ord_al$LOO.recall, type="l")
plot(ord_reg$LOO.prec, ord_reg$LOO.recall, type="l")
plot(ord_con$LOO.prec, ord_con$LOO.recall, type="l")

plot(ord_al$LOO.prec, ord_al$LOO.recall, type="l", col = "cyan3", lwd = 3, 
     xlim=c(0,0.76), ylim=c(0,0.12), xlab="Precision", ylab="Recall")
lines(ord_reg$LOO.prec, ord_reg$LOO.recall, type="l", col = "darkolivegreen", lwd=3)
lines(ord_con$LOO.prec, ord_con$LOO.recall, type="l", col = "orangered", lwd=3)
#lines(ord_uni$LOO.prec, ord_uni$LOO.recall, type="l", col = "darkorchid3", lwd=3)
legend("topright", c("Adaptive Lasso", "Mutual Information", "Consensus"), lty=c(1,1,1), lwd=c(3,3,3),
       col=c("cyan3", "darkolivegreen", "orangered"))

trapz(ord_al$LOO.prec, ord_al$LOO.recall)
trapz(ord_reg$LOO.prec, ord_reg$LOO.recall)
trapz(ord_con$LOO.prec, ord_con$LOO.recall)
#trapz(ord_uni$LOO.prec, ord_uni$LOO.recall)

#plot fscore distributions
al13_p <- read.csv("lasso_clustering_results/I13_res/cons_sig_precision.txt", header=F)
al13_r <- read.table("lasso_clustering_results/I13_res/cons_sig_recall.txt", sep= ",", header=F)
al13_f <- as.numeric(as.vector(2*(al13_p*al13_r)/(al13_p+al13_r)))
al14_p <- read.csv("lasso_clustering_results/I14_res/cons_sig_precision.txt", header=F)
al14_r <- read.table("lasso_clustering_results/I14_res/cons_sig_recall.txt", sep= ",", header=F)
al14_f <- as.numeric(as.vector(2*(al14_p*al14_r)/(al14_p+al14_r)))
al26_p <- read.csv("lasso_clustering_results/I26_res/cons_sig_precision.txt", header=F)
al26_r <- read.table("lasso_clustering_results/I26_res/cons_sig_recall.txt", sep= ",", header=F)
al26_f <- as.numeric(as.vector(2*(al26_p*al26_r)/(al26_p+al26_r)))

mi17_p <- read.csv("reg_clustering_results/I17_res/cons_sig_precision.txt", header=F)
mi17_r <- read.csv("reg_clustering_results/I17_res/cons_sig_recall.txt", header=F)
mi17_f <- as.numeric(as.vector(2*(mi17_p*mi17_r)/(mi17_p+mi17_r)))
mi22_p <- read.csv("reg_clustering_results/I22_res/cons_sig_precision.txt", header=F)
mi22_r <- read.csv("reg_clustering_results/I22_res/cons_sig_recall.txt", header=F)
mi22_f <- as.numeric(as.vector(2*(mi22_p*mi22_r)/(mi22_p+mi22_r)))
mi33_p <- read.csv("reg_clustering_results/I33_res/cons_sig_precision.txt", header=F)
mi33_r <- read.csv("reg_clustering_results/I33_res/cons_sig_recall.txt", header=F)
mi33_f <- as.numeric(as.vector(2*(mi33_p*mi33_r)/(mi33_p+mi33_r)))

con13_p <- read.csv("consensus_clustering_results/I13_res/cons_sig_precision.txt", header=F)
con13_r <- read.csv("consensus_clustering_results/I13_res/cons_sig_recall.txt", header=F)
con13_f <- as.numeric(as.vector(2*(con13_p*con13_r)/(con13_p+con13_r)))
con2_p <- read.csv("consensus_clustering_results/I2_res/cons_sig_precision.txt", header=F)
con2_r <- read.csv("consensus_clustering_results/I2_res/cons_sig_recall.txt", header=F)
con2_f <- as.numeric(as.vector(2*(con2_p*con2_r)/(con2_p+con2_r)))
con22_p <- read.csv("consensus_clustering_results/I22_res/cons_sig_precision.txt", header=F)
con22_r <- read.csv("consensus_clustering_results/I22_res/cons_sig_recall.txt", header=F)
con22_f <- as.numeric(as.vector(2*(con22_p*con22_r)/(con22_p+con22_r)))
con33_p <- read.csv("consensus_clustering_results/I33_res/cons_sig_precision.txt", header=F)
con33_r <- read.csv("consensus_clustering_results/I33_res/cons_sig_recall.txt", header=F)
con33_f <- as.numeric(as.vector(2*(con33_p*con33_r)/(con33_p+con33_r)))


boxplot(al13_f, al14_f, al26_f, mi17_f, mi22_f, mi33_f, con13_f,con2_f,con22_f,con33_f,
        col = cols_vec, names = pr_res$V1, las=2, ylab = "Fscore (per enriched cluster)")







########################## UNION_ALSIZE, adalasso, MI, consensus #######################################

#plot the precision/recall
setwd("~/ferdig_rotation/regulon_validation/original_nets/CRC_MCL_clustering/")
pr_res <- read.csv("best_LOO_precision_recall_union.csv", header = F, as.is = T)

cols_vec <- c("darkcyan", "cyan3", "lightseagreen",
              "darkolivegreen", "darkolivegreen4", "darkolivegreen3",
              "orangered3", "orangered", "orange2", "orange",
              "darkorchid4", "darkorchid3", "darkorchid2", "darkorchid1")


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

plot(pr_res$V2, pr_res$V3, xlab="Precision", ylab="Recall")
points(pr_res$V2, pr_res$V3, col = cols_vec, pch=20, cex=4)
legend("topright", inset=c(-0.2,0), pr_res$V1, pch=20, col = cols_vec)

n_go_terms <- c(41,41,31,14,30,7,13,28,28,12,31,51,20,7)/10
n_genes <- c(615,436,151,3563,2745,436,2211,910,663,96,2469,955,317,203)/300

plot(pr_res$V2, pr_res$V3, xlab="Precision", ylab="Recall")
points(pr_res$V2, pr_res$V3, col = cols_vec, pch=20, cex=n_genes)
legend("topright", inset=c(-0.2,0), pr_res$V1, pch=20, col = cols_vec)


#plot fscores of top best
fscore <- 2*(pr_res$V2*pr_res$V3)/(pr_res$V2+pr_res$V3)
barplot(fscore[c(11,12,13,14,1,2,3,4,5,6,7,8,9,10)], col = cols_vec[c(11,12,13,14,1,2,3,4,5,6,7,8,9,10)], 
        ylab = "LOO Fscore", names.arg = pr_res$V1[c(11,12,13,14,1,2,3,4,5,6,7,8,9,10)], las=2)



#plot fscore distributions
al13_p <- read.csv("lasso_clustering_results/I13_res/cons_sig_precision.txt", header=F)
al13_r <- read.table("lasso_clustering_results/I13_res/cons_sig_recall.txt", sep= ",", header=F)
al13_f <- as.numeric(as.vector(2*(al13_p*al13_r)/(al13_p+al13_r)))
al14_p <- read.csv("lasso_clustering_results/I14_res/cons_sig_precision.txt", header=F)
al14_r <- read.table("lasso_clustering_results/I14_res/cons_sig_recall.txt", sep= ",", header=F)
al14_f <- as.numeric(as.vector(2*(al14_p*al14_r)/(al14_p+al14_r)))
al26_p <- read.csv("lasso_clustering_results/I26_res/cons_sig_precision.txt", header=F)
al26_r <- read.table("lasso_clustering_results/I26_res/cons_sig_recall.txt", sep= ",", header=F)
al26_f <- as.numeric(as.vector(2*(al26_p*al26_r)/(al26_p+al26_r)))

mi17_p <- read.csv("reg_clustering_results/I17_res/cons_sig_precision.txt", header=F)
mi17_r <- read.csv("reg_clustering_results/I17_res/cons_sig_recall.txt", header=F)
mi17_f <- as.numeric(as.vector(2*(mi17_p*mi17_r)/(mi17_p+mi17_r)))
mi22_p <- read.csv("reg_clustering_results/I22_res/cons_sig_precision.txt", header=F)
mi22_r <- read.csv("reg_clustering_results/I22_res/cons_sig_recall.txt", header=F)
mi22_f <- as.numeric(as.vector(2*(mi22_p*mi22_r)/(mi22_p+mi22_r)))
mi33_p <- read.csv("reg_clustering_results/I33_res/cons_sig_precision.txt", header=F)
mi33_r <- read.csv("reg_clustering_results/I33_res/cons_sig_recall.txt", header=F)
mi33_f <- as.numeric(as.vector(2*(mi33_p*mi33_r)/(mi33_p+mi33_r)))

con13_p <- read.csv("consensus_clustering_results/I13_res/cons_sig_precision.txt", header=F)
con13_r <- read.csv("consensus_clustering_results/I13_res/cons_sig_recall.txt", header=F)
con13_f <- as.numeric(as.vector(2*(con13_p*con13_r)/(con13_p+con13_r)))
con2_p <- read.csv("consensus_clustering_results/I2_res/cons_sig_precision.txt", header=F)
con2_r <- read.csv("consensus_clustering_results/I2_res/cons_sig_recall.txt", header=F)
con2_f <- as.numeric(as.vector(2*(con2_p*con2_r)/(con2_p+con2_r)))
con22_p <- read.csv("consensus_clustering_results/I22_res/cons_sig_precision.txt", header=F)
con22_r <- read.csv("consensus_clustering_results/I22_res/cons_sig_recall.txt", header=F)
con22_f <- as.numeric(as.vector(2*(con22_p*con22_r)/(con22_p+con22_r)))
con33_p <- read.csv("consensus_clustering_results/I33_res/cons_sig_precision.txt", header=F)
con33_r <- read.csv("consensus_clustering_results/I33_res/cons_sig_recall.txt", header=F)
con33_f <- as.numeric(as.vector(2*(con33_p*con33_r)/(con33_p+con33_r)))

union13_p <- read.csv("union_clustering_results/I13_res/cons_sig_precision.txt", header=F)
union13_r <- read.csv("union_clustering_results/I13_res/cons_sig_recall.txt", header=F)
union13_f <- as.numeric(as.vector(2*(union13_p*union13_r)/(union13_p+union13_r)))
union16_p <- read.csv("union_clustering_results/I16_res/cons_sig_precision.txt", header=F)
union16_r <- read.csv("union_clustering_results/I16_res/cons_sig_recall.txt", header=F)
union16_f <- as.numeric(as.vector(2*(union16_p*union16_r)/(union16_p+union16_r)))
union25_p <- read.csv("union_clustering_results/I25_res/cons_sig_precision.txt", header=F)
union25_r <- read.csv("union_clustering_results/I25_res/cons_sig_recall.txt", header=F)
union25_f <- as.numeric(as.vector(2*(union25_p*union25_r)/(union25_p+union25_r)))
union33_p <- read.csv("union_clustering_results/I33_res/cons_sig_precision.txt", header=F)
union33_r <- read.csv("union_clustering_results/I33_res/cons_sig_recall.txt", header=F)
union33_f <- as.numeric(as.vector(2*(union33_p*union33_r)/(union33_p+union33_r)))


boxplot(union13_f, union16_f, union25_f, union33_f, al13_f, al14_f, al26_f, mi17_f, mi22_f, mi33_f, con13_f,con2_f,con22_f,con33_f,
        col = cols_vec[c(11,12,13,14,1,2,3,4,5,6,7,8,9,10)], 
        names = pr_res$V1[c(11,12,13,14,1,2,3,4,5,6,7,8,9,10)], las=2, 
        ylab = "Fscore (per enriched cluster)")
