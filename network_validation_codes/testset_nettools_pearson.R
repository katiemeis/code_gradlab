

geneA <- c(5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3, 4)
geneB <- c(4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3) #move A by one time step, want to be higher
geneC <- c(5, 12, 28, 40, 28, 12, 5, -2, -20, -33, -20, -2) #amplify A, want to be higher, nonlinear
geneD <- c(2, 12, 20, 14, 30, 5, 22, 28, 14, 14, 4, 25) #random
geneE <- c(1,2,3,4,5,6,7,8,9,10,11,12)
geneF <- c(5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 3, 4) #repeat of A
geneG <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24) # multiply E by 2
geneH <- geneD*5 #multiply D by 5
geneI <- c(5, 4, 3, 2, 3, 4, 5, 6, 7, 8, 7, 6) # negative correlation with A 
plot(c(1:12), geneC)
lines(geneA)
lines(geneB)
lines(geneD)
lines(geneE)

test <- data.frame(geneA, geneB, geneC, geneD, geneE, geneF, geneG, geneH, geneI)
#matrix or dataframe, samples by genes
#does absolute value pearson
library(nettools)
hu_pearson_result <- mat2adj(test, method="cor")
cor(geneA, geneB)
cor(geneA, geneC)
cor(geneA, geneD)
cor(geneA, geneI)

#checks out