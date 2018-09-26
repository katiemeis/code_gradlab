#testing DTWMIC

library(nettools)
require(reshape2)
#mat2adj takes in matrix or dataframe with N sample rows and P gens columns

geneA <- c(5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3, 4)
geneB <- c(4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3) #move A by one time step, want to be higher
geneC <- c(5, 12, 28, 40, 28, 12, 5, -2, -20, -33, -20, -2) #amplify A, want to be higher
geneD <- c(2, 12, 20, 14, 30, 5, 22, 28, 14, 14, 4, 25) #random, shouldnt do well with repect to anything
geneE <- c(-33, -20, -2, 5, 12, 28, 40, 28, 12, 5, -2, -20) #moved C by 3 time points, want to be lower
geneF <- c(5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3, 4) #repeat of A
plot(c(1:12), geneC)
lines(geneA)
lines(geneB)
lines(geneD)
lines(geneE)

test <- data.frame(geneA, geneB, geneC, geneD, geneE, geneF)
DTWMIC_test <- mat2adj(test, method='DTWMIC')

pearson_test <- mat2adj(test, method="cor")
