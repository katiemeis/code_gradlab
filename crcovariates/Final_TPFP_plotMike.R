setwd("~/ferdig_rotation/gabes_CRC/")

######################## Group A ###############################

######## TP in uncorrected out of total TP
######## FP in uncorrected out of total predicted
uncorrected <- read.csv("GrpAuncorrectedcorr.csv", stringsAsFactors = F, header = T, row.names = 1)
uncorrected <- uncorrected[-c(1,2),]

corrected <- read.csv("GrpACorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
corrected <- corrected[-c(1,2,3,4),]

#want to plot anything with Pvaluse 0.01 AND FDR 0.05 or better
corrected_filt <- corrected[which(corrected$PValue < 0.01 & corrected$FDR < 0.05),]
uncorrected_filt <- uncorrected[which(uncorrected$PValue < 0.01 & uncorrected$FDR < 0.05),]

# this is total TP
ground_truth <- row.names(corrected_filt)

my_seq <- seq(0, 0.01, 0.000001)
tp_vec <- vector()
fp_vec <- vector()
tp_vec_corrected <- vector()

tp_count_vec <- vector()
fp_count_vec <- vector()
tp_count_corrected <-vector()

for(i in my_seq){
  #threshold the uncorrected
  row_ind <- which(uncorrected_filt$PValue < i)
  all_pos <- row.names(uncorrected_filt)[row_ind]
  #threshold the corrected
  row_ind_cor <- which(corrected_filt$PValue < i)
  all_pos_cor <- row.names(corrected_filt)[row_ind_cor]
  #counts
  tp_count <- length(intersect(ground_truth, all_pos))
  fp_count <- length(all_pos) - tp_count
  tp_count_cor <- length(all_pos_cor)
  tp_count_vec <- c(tp_count_vec, tp_count)
  fp_count_vec <- c(fp_count_vec, fp_count)
  tp_count_corrected <- c(tp_count_corrected, tp_count_cor)
  #percents
  tp_percent <- tp_count/length(ground_truth)
  fp_percent <- fp_count/length(all_pos)
  tp_percent_cor <- tp_count_cor/length(ground_truth)
  tp_vec <- c(tp_vec, tp_percent)
  fp_vec <- c(fp_vec, fp_percent)
  tp_vec_corrected <- c(tp_vec_corrected, tp_percent_cor)
}

#plot the counts on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_count_vec[-1], xlab = "P-value", ylab = "True Positives (count)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_count_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, 50, 100, 150, 200, 250))
mtext("False Positives (count)", side=4, line=3)


#plot the percents on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, .05, 0.1, 0.15))
mtext("False Positives (%)", side=4, line=3)
#legend("topright",legend = c("TP", "FP"), col = c("blue", "green"))


#counts but added the orange line
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_count_vec[-1], xlab = "P-value", ylab = "True Positives (count)", type = "o", col = "blue", lwd = 3, ylim=c(500,1400)) # first plot
lines(my_seq[-1], tp_count_corrected[-1], type="o", col = "orange", lwd=3)
par(new = TRUE)
plot(my_seq[-1], fp_count_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, 50, 100, 150, 200, 250))
mtext("False Positives (count)", side=4, line=3)


#percents but added the orange line
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3, ylim=c(0.37,1)) # first plot
lines(my_seq[-1], tp_vec_corrected[-1], type="o", col = "orange", lwd=3)
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, .05, 0.1, 0.15))
mtext("False Positives (%)", side=4, line=3)
legend("bottomright", legend = c("TP uncorrected", "TP corrected", "FP uncorrected"), col=c("blue", "orange", "green"), 
       pch=c(19, 19, 19), lwd = c(3,3,3))



################################### Group B ################################
setwd("~/ferdig_rotation/gabes_CRC/")
######## TP in uncorrected out of total TP
######## FP in uncorrected out of total predicted
uncorrected <- read.csv("GrpBuncorrectedcorr.csv", stringsAsFactors = F, header = T, row.names = 1)
uncorrected <- uncorrected[-c(1,2),]

corrected <- read.csv("GrpBCorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
corrected <- corrected[-c(1,2,3,4),]

#want to plot anything with Pvaluse 0.01 AND FDR 0.05 or better
corrected_filt <- corrected[which(corrected$PValue < 0.01 & corrected$FDR < 0.05),]
uncorrected_filt <- uncorrected[which(uncorrected$PValue < 0.01 & uncorrected$FDR < 0.05),]

# this is total TP
ground_truth <- row.names(corrected_filt)

my_seq <- seq(0, 0.01, 0.000001)
tp_vec <- vector()
fp_vec <- vector()
tp_vec_corrected <- vector()

tp_count_vec <- vector()
fp_count_vec <- vector()
tp_count_corrected <-vector()

for(i in my_seq){
  #threshold the uncorrected
  row_ind <- which(uncorrected_filt$PValue < i)
  all_pos <- row.names(uncorrected_filt)[row_ind]
  #threshold the corrected
  row_ind_cor <- which(corrected_filt$PValue < i)
  all_pos_cor <- row.names(corrected_filt)[row_ind_cor]
  #counts
  tp_count <- length(intersect(ground_truth, all_pos))
  fp_count <- length(all_pos) - tp_count
  tp_count_cor <- length(all_pos_cor)
  tp_count_vec <- c(tp_count_vec, tp_count)
  fp_count_vec <- c(fp_count_vec, fp_count)
  tp_count_corrected <- c(tp_count_corrected, tp_count_cor)
  #percents
  tp_percent <- tp_count/length(ground_truth)
  fp_percent <- fp_count/length(all_pos)
  tp_percent_cor <- tp_count_cor/length(ground_truth)
  tp_vec <- c(tp_vec, tp_percent)
  fp_vec <- c(fp_vec, fp_percent)
  tp_vec_corrected <- c(tp_vec_corrected, tp_percent_cor)
}

#plot the percents on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_count_vec[-1], xlab = "P-value", ylab = "True Positives (count)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_count_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, 50, 100, 150, 200, 250))
mtext("False Positives (count)", side=4, line=3)


#plot the percents on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, .05, 0.1, 0.15))
mtext("False Positives (%)", side=4, line=3)
#legend("topright",legend = c("TP", "FP"), col = c("blue", "green"))

par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3, ylim=c(0.17,1)) # first plot
lines(my_seq[-1], tp_vec_corrected[-1], type="o", col = "orange", lwd=3)
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0.50, 0.55, 0.6, 0.65, 0.70))
mtext("False Positives (%)", side=4, line=3)
legend("bottomright", legend = c("TP uncorrected", "TP corrected", "FP uncorrected"), col=c("blue", "orange", "green"), 
       pch=c(19, 19, 19), lwd = c(3,3,3))









############################# Group AB ###############################
setwd("~/ferdig_rotation/gabes_CRC/")
######## TP in uncorrected out of total TP
######## FP in uncorrected out of total predicted
uncorrected <- read.csv("GrpAandBuncorrectedcorr.csv", stringsAsFactors = F, header = T, row.names = 1)
uncorrected <- uncorrected[-c(1,2),]

corrected <- read.csv("GrpAandBCorrectedCorr.csv", stringsAsFactors = F, header = T, row.names = 1)
corrected <- corrected[-c(1,2,3,4),]

#want to plot anything with Pvaluse 0.01 AND FDR 0.05 or better
corrected_filt <- corrected[which(corrected$PValue < 0.01 & corrected$FDR < 0.05),]
uncorrected_filt <- uncorrected[which(uncorrected$PValue < 0.01 & uncorrected$FDR < 0.05),]

# this is total TP
ground_truth <- row.names(corrected_filt)

my_seq <- seq(0, 0.01, 0.000001)
tp_vec <- vector()
fp_vec <- vector()
tp_vec_corrected <- vector()

tp_count_vec <- vector()
fp_count_vec <- vector()
tp_count_corrected <-vector()

for(i in my_seq){
  #threshold the uncorrected
  row_ind <- which(uncorrected_filt$PValue < i)
  all_pos <- row.names(uncorrected_filt)[row_ind]
  #threshold the corrected
  row_ind_cor <- which(corrected_filt$PValue < i)
  all_pos_cor <- row.names(corrected_filt)[row_ind_cor]
  #counts
  tp_count <- length(intersect(ground_truth, all_pos))
  fp_count <- length(all_pos) - tp_count
  tp_count_cor <- length(all_pos_cor)
  tp_count_vec <- c(tp_count_vec, tp_count)
  fp_count_vec <- c(fp_count_vec, fp_count)
  tp_count_corrected <- c(tp_count_corrected, tp_count_cor)
  #percents
  tp_percent <- tp_count/length(ground_truth)
  fp_percent <- fp_count/length(all_pos)
  tp_percent_cor <- tp_count_cor/length(ground_truth)
  tp_vec <- c(tp_vec, tp_percent)
  fp_vec <- c(fp_vec, fp_percent)
  tp_vec_corrected <- c(tp_vec_corrected, tp_percent_cor)
}

#plot the percents on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_count_vec[-1], xlab = "P-value", ylab = "True Positives (count)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_count_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, 50, 100, 150, 200, 250))
mtext("False Positives (count)", side=4, line=3)


#plot the percents on the same graph, TP is blue FP is green
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3) # first plot
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0, .05, 0.1, 0.15))
mtext("False Positives (%)", side=4, line=3)
#legend("topright",legend = c("TP", "FP"), col = c("blue", "green"))

par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(my_seq[-1], tp_vec[-1], xlab = "P-value", ylab = "True Positives (%)", type = "o", col = "blue", lwd = 3, ylim=c(0.40,1)) # first plot
lines(my_seq[-1], tp_vec_corrected[-1], type="o", col = "orange", lwd=3)
par(new = TRUE)
plot(my_seq[-1], fp_vec[-1], type = "o", col="green", axes = FALSE, bty = "n", xlab = "", ylab = "", lwd = 3)
axis(side=4, at = c(0.3, 0.35, 0.4, 0.45, 0.50))
mtext("False Positives (%)", side=4, line=3)
legend("bottomright", legend = c("TP uncorrected", "TP corrected", "FP uncorrected"), col=c("blue", "orange", "green"), 
       pch=c(19, 19, 19), lwd = c(3,3,3))

