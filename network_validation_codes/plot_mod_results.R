
setwd("~/ferdig_rotation/regulon_validation/original_nets/modularity/")
#plot modularity results

mod_res <- read.csv("mod_results_all.csv")
thresholds <- mod_res[,1]

plot(thresholds, mod_res$Pearson, xlab="Number of edges", ylab="Modularity", ylim = c(0,.1), type="o", lwd=2)
lines(thresholds, mod_res$MI, col="red", lwd=2, type="o")
lines(thresholds, mod_res$RF, col="blue", lwd = 2, type="o")
legend("topright", c("Pearson Correlation", "Mutual Information", "Random Forest"), col=c("black", "red", "blue"), pch=c(1,1), lwd=2)

