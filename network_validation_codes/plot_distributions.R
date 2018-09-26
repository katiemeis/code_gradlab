#make histograms of null distributions for modularity
setwd("~/ferdig_rotation/regulon_validation/series_thresholds/modularity_null_distributions/")

#percent difference
(6-3)/3 #1, 6 is 100% more than 3
(6+3)/3 #3, 6 300% more than -3


##### load the results files, lists for random to loop over

RF_mod_res <- read.table("../modularity_res_RF.txt", header=T, sep = "\t")
reg_mod_res <- read.table("../modularity_res_Reg.txt", header=T, sep = "\t")
pear_mod_res <- read.table("../modularity_res_Pear.txt", header=T, sep = "\t")

#plot modularity
plot(reg_mod_res$Edges, reg_mod_res$Modularity, type = "l", ylim = c(-0.01,0.1), lwd=2, xlab = "Number of edges", ylab = "Modularity")
lines(RF_mod_res$Edges, RF_mod_res$Modularity, lwd=2, col = "red")
lines(pear_mod_res$Edges, pear_mod_res$Modularity, lwd=2, col="blue")
legend("topright", c("Random Forest", "Mutual Information", "Pearson Correlation"), col=c("red", "black", "blue"), pch=c(1,1), lwd=2)

rf_rand_files <- c("mod_perm_rf_5000.csv", "mod_perm_rf_7500.csv", "mod_perm_rf_10000.csv", 
              "mod_perm_rf_25000.csv", "mod_perm_rf_50000.csv", "mod_perm_rf_75000.csv", 
              "mod_perm_rf_100000.csv", "mod_perm_rf_150000.csv")

reg_rand_files <- c("mod_perm_reg_5000.csv", "mod_perm_reg_7500.csv", "mod_perm_reg_10000.csv", 
                   "mod_perm_reg_25000.csv", "mod_perm_reg_50000.csv", "mod_perm_reg_75000.csv", 
                   "mod_perm_reg_100000.csv", "mod_perm_reg_150000.csv")

pear_rand_files <- c("mod_perm_pear_5000.csv", "mod_perm_pear_7500.csv", "mod_perm_pear_10000.csv", 
                   "mod_perm_pear_25000.csv", "mod_perm_pear_50000.csv", "mod_perm_pear_75000.csv", 
                   "mod_perm_pear_100000.csv", "mod_perm_pear_150000.csv")


############ random forest ###########
rf_pvals <- vector()
rf_perc_inc <- vector()
rf_which <- vector()
for(i in 1:8){
  rf_rand <- read.csv(rf_rand_files[i], header = F, as.is = T)
  mod_score <- RF_mod_res[i,2]
  #get p-val
  p <- pnorm(mod_score, mean = mean(rf_rand$V1), sd = sd(rf_rand$V1), lower.tail = F)
  rf_pvals <- c(rf_pvals, p)
  inc <- (mod_score-mean(rf_rand$V1))/abs(mean(rf_rand$V1))
  rf_perc_inc <- c(rf_perc_inc, inc)
  whichp <- length(which(rf_rand$V1 > mod_score))
  rf_which <- c(rf_which, whichp)
}
rf_pvals
rf_perc_inc
rf_which


############ regulon/mutual information ##########

reg_pvals <- vector()
reg_perc_inc <- vector()
reg_which <- vector()
for(i in 1:8){
  reg_rand <- read.csv(reg_rand_files[i], header = F, as.is = T)
  mod_score <- reg_mod_res[i,2]
  #get p-val
  p <- pnorm(mod_score, mean = mean(reg_rand$V1), sd = sd(reg_rand$V1), lower.tail = F)
  reg_pvals <- c(reg_pvals, p)
  inc <- (mod_score-mean(reg_rand$V1))/abs(mean(reg_rand$V1))
  reg_perc_inc <- c(reg_perc_inc, inc)
  whichp <- length(which(reg_rand$V1 > mod_score))
  reg_which <- c(reg_which, whichp)
}
reg_pvals
reg_perc_inc
reg_which

################ Pearosn correlation ##################

pear_pvals <- vector()
pear_perc_inc <- vector()
pear_which <- vector()
for(i in 1:8){
  pear_rand <- read.csv(pear_rand_files[i], header = F, as.is = T)
  mod_score <- pear_mod_res[i,2]
  #get p-val
  p <- pnorm(mod_score, mean = mean(pear_rand$V1), sd = sd(pear_rand$V1), lower.tail = F)
  pear_pvals <- c(pear_pvals, p)
  inc <- (mod_score-mean(pear_rand$V1))/abs(mean(pear_rand$V1))
  pear_perc_inc <- c(pear_perc_inc, inc)
  whichp <- length(which(pear_rand$V1 > mod_score))
  pear_which <- c(pear_which, whichp)
}
pear_pvals
pear_perc_inc
pear_which

#plot percent change for all
plot(reg_mod_res$Edges, reg_perc_inc, type = "l", lwd = 2, ylim = c(-70, 1100))
lines(reg_mod_res$Edges, rf_perc_inc, lwd = 2, col = "red")
lines(reg_mod_res$Edges, pear_perc_inc, lwd = 2, col = "blue")
legend("topright", c("Random Forest", "Mutual Information", "Pearson Correlation"), col=c(2, 1, "blue"), pch=c(1,1), lwd=2)








######### random forest networks ##########
#all p-values return 0
rf_5k <- read.csv("mod_perm_rf_5000.csv", header = F, as.is = T)
hist(rf_5k$V1)


rf_7k <- read.csv("mod_perm_rf_7500.csv", header = F, as.is = T)
hist(rf_7k$V1)


rf_10k <- read.csv("mod_perm_rf_10000.csv", header = F, as.is = T)
hist(rf_10k$V1)


rf_25k <- read.csv("mod_perm_rf_25000.csv", header = F, as.is = T)
hist(rf_25k$V1)


rf_50k <- read.csv("mod_perm_rf_50000.csv", header = F, as.is = T)
hist(rf_50k$V1)


rf_75k <- read.csv("mod_perm_rf_75000.csv", header = F, as.is = T)
hist(rf_75k$V1)


rf_100k <- read.csv("mod_perm_rf_100000.csv", header = F, as.is = T)
hist(rf_100k$V1)


rf_150k <- read.csv("mod_perm_rf_150000.csv", header = F, as.is = T)
hist(rf_150k$V1)



######## regulon/mutual informaiton netwroks ##########
reg_5k <- read.csv("mod_perm_reg_5000.csv", header = F, as.is = T)
hist(reg_5k$V1)


reg_7k <- read.csv("mod_perm_reg_7500.csv", header = F, as.is = T)
hist(reg_7k$V1)


reg_10k <- read.csv("mod_perm_reg_10000.csv", header = F, as.is = T)
hist(reg_10k$V1)


reg_25k <- read.csv("mod_perm_reg_25000.csv", header = F, as.is = T)
hist(reg_25k$V1)


reg_50k <- read.csv("mod_perm_reg_50000.csv", header = F, as.is = T)
hist(reg_50k$V1)


reg_75k <- read.csv("mod_perm_reg_75000.csv", header = F, as.is = T)
hist(reg_75k$V1)


reg_100k <- read.csv("mod_perm_reg_100000.csv", header = F, as.is = T)
hist(reg_100k$V1)


reg_150k <- read.csv("mod_perm_reg_150000.csv", header = F, as.is = T)
hist(reg_150k$V1)


######## pearson networks ########
pear_5k <- read.csv("mod_perm_pear_5000.csv", header = F, as.is = T)
hist(pear_5k$V1)


pear_7k <- read.csv("mod_perm_pear_7500.csv", header = F, as.is = T)
hist(pear_7k$V1)


pear_50k <- read.csv("mod_perm_pear_50000.csv", header = F, as.is = T)
hist(pear_50k$V1)

