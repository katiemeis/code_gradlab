
setwd("~/ferdig_rotation/regulon_validation/original_nets/precision_recall/")
p_res <- read.csv("pearson_results.csv")
reg_res <- read.csv("regulon_results.csv")
rf_res <- read.csv("rf_results.csv")

thresholds <- c(5000, 7500, 10000, 25000, 50000, 75000, 100000, 150000)

#plot precision
plot(thresholds, p_res[4, 2:9], xlab="Number of edges", ylab="Precision", ylim = c(0,.1), type="o", lwd=2)
lines(thresholds, reg_res[4, 2:9], col="red", lwd=2, type="o")
lines(thresholds, rf_res[4, 2:9], col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#plot recall
plot(thresholds, p_res[5, 2:9], xlab="Number of edges", ylab="Recall", ylim = c(0,.1), type="o", lwd=2)
lines(thresholds, reg_res[5, 2:9], col="red", lwd=2, type="o")
lines(thresholds, rf_res[5, 2:9], col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)


#plot fscore
plot(thresholds, p_res[6, 2:9], xlab="Number of edges", ylab="Fscore", ylim = c(0,.1), type="o", lwd=2)
lines(thresholds, reg_res[6, 2:9], col="red", lwd=2, type="o")
lines(thresholds, rf_res[6, 2:9], col="blue", lwd = 2, type="o")
legend("bottomright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

#plot precision vs recall
plot(as.numeric(as.vector(p_res[4, 2:9])), as.numeric(as.vector(p_res[5, 2:9])), xlab="Precision", ylab="Recall", lwd=2, xlim=c(0,0.1), ylim=c(0,0.1))
points(as.numeric(as.vector(reg_res[4, 2:9])), as.numeric(as.vector(reg_res[5, 2:9])), col="red", lwd=2)
points(as.numeric(as.vector(rf_res[4, 2:9])), as.numeric(as.vector(rf_res[5, 2:9])), col="blue", lwd = 2)
legend("topleft", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)
