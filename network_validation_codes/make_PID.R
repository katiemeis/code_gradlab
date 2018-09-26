#to make a PID file

goslim <- read.table("goslim2.dat", header = T)
library(reshape)
newmap <- cast(goslim, GENE_SOURCE_ID ~ GO_ID)
newmap[newmap > 1] <- 1
write.csv(newmap, "PID_function1.csv", quote = F, row.names = F)