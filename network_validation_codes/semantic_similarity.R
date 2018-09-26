
################################# Semantic Similarity ########################

setwd("~/ferdig_rotation/regulon_validation/original_nets/")

#source("https://bioconductor.org/biocLite.R")
#biocLite("GOSemSim")
library(GOSemSim)

#load the GO term map to get the list of GO terms in Hu data
#need check.names = F to leave the colons as colons instead of changing to periods
go_mat <- read.csv("GO_file/PID_6_10_NEW-Qi.csv", as.is=T, check.names = F)
go_terms <- colnames(go_mat)[2:length(go_mat)]

#make the semData object
d <- godata(OrgDb = 'org.Pf.plasmo.db', keytype = "SYMBOL", ont="BP", computeIC=TRUE)

semsim_mat <- termSim(go_terms, go_terms, semData = d, method = "Rel")

#need to remove obsolete terms
to_remove <- c(76, 90, 158, 208, 448, 529, 530)
go_terms[to_remove]

go_terms[448] <- "GO:0006325"
go_terms[529] <- "GO:0043087"
go_terms[530] <- "GO:0043087"
to_remove <- c(76, 90, 158, 208)

#"GO:0006118" "GO:0006184" "GO:0006410" "GO:0006512" "GO:0016568" "GO:0032312" "GO:0032313"
#obs, obs, obs, obs, GO:0006325, GO:0043087, GO:0043087
#obsolete or alternate IDs - should I combine 529 and 530 into one column in the GO file
#should I get rid of all of these, only obsolete, none of them?

go_terms1 <- go_terms[-to_remove]


semsim_mat <- termSim(go_terms1, go_terms1, semData = d, method = "Rel")

#this does the same as above but rounds to 3 decimals
#testcomb_mat <- mgoSim(go_terms1, go_terms1, semData = d, measure = "Rel", combine = NULL)

#need to convert to a distance matrix - ask Katie
#can do 1- val, 1/val, or euclidean dist - 1-val makes most sense to me
dis_mat <- 1-semsim_mat
test <- as.dist(dis_mat, diag = T)

#heirarchical clustering of dissimilarity matrix hclust
tclust <- hclust(test, method = "average")
final_clust <- cutree(tclust, h=0.7)

#randomly choose a cluster representative






#found for euclidean I think
#http://onertipaday.blogspot.com/2007/04/from-similarity-to-distance-matrix.html
# This function returns an object of class "dist"
#sim2dist <- function(mx) as.dist(sqrt(outer(diag(mx), diag(mx), "+") - 2*mx)) 
#from similarity to distance matrix 
#s.mx = as.matrix(s.mx) 
#.mx = sim2dist(s.mx)

#heirarchical clustering of dissimilarity matrix hclust
