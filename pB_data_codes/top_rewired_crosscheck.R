#check the lists of top rewired genes and see which overlap for all the networks

setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/rewired_nets/")

pb58_table <- read.csv("rewired_PB58/largest_cc_degree1.csv") 
pb54_table <- read.csv("PB54_nodetable_bydegree.csv")
pb55_table <- read.csv("PB55_nodetable_bydegree.csv")
pb57_table <- read.csv("PB57_nodetable_bydegree.csv")

p58_100 <- pb58_table[1:100,c("name", "Degree")]
p54_100 <- pb54_table[1:100,c("name", "Degree")]
p55_100 <- pb55_table[1:100,c("name", "Degree")]
p57_100 <- pb57_table[1:100,c("name", "Degree")]

p58_100 <- pb58_table[1:500,"name"]
p54_100 <- pb54_table[1:500,"name"]
p55_100 <- pb55_table[1:500,"name"]
p57_100 <- pb57_table[1:500,"name"]

Reduce(intersect, list(p58_100, p54_100, p55_100, p57_100))
