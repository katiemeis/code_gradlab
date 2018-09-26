setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/top5/")

NF54_degree_data <- read.csv("NF54_degree_sorted.txt", sep='\t', header=FALSE)
x = as.vector(NF54_degree_data[,2])
x2 = as.numeric(x)
NF_points = x2[1:5523]


#plot y=#, then y=%
library(igraph)
plot(tabulate(NF_points), log="xy", xlab='Degree', ylab='Number of nodes', main="NF54 Degree Distribution")

percent_vec_NF = tabulate(NF_points)/length(NF_points)
plot(percent_vec_NF, log="xy", xlab='Degree', ylab='Number of nodes (%)', main="NF54 Degree Distribution")





PB58_degree_data <- read.csv("PB58_degree_sorted.txt", sep='\t', header=FALSE)
x = as.vector(PB58_degree_data[,2])
x2 = as.numeric(x)
PB_points = x2[1:5531]

#plot y=#, then y=%
library(igraph)
plot(tabulate(PB_points), log="xy", xlab='Degree', ylab='Number of nodes', main="PB58 Degree Distribution")

percent_vec_PB = tabulate(PB_points)/length(PB_points)
plot(percent_vec_PB, log="xy", xlab='Degree', ylab='Number of nodes (%)', main="PB58 Degree Distribution")
