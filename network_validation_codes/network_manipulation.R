#get rid of self loops in rf net
rf_net <- read.csv("../sizematch_nodups_nets/rf_size_names_nodups.csv", as.is = T)

idex_vec <- vector()
for(i in 1:dim(rf_net)[1]){
  if(rf_net[i,1] == rf_net[i,2]){
    idex_vec <- c(idex_vec, i)
  }
}

rf_net2 <- rf_net[-idex_vec,]
write.csv(rf_net2,file="../sizematch_nodups_nets/rf_size_names_nodups2.csv",row.names=F, quote=F)

#order reg net
full_reg <- read.csv("~/ferdig_rotation/regulon_validation/sizematch_nodups_nets/fullreg_size_names_nodups.csv", as.is=T)
full_rego <- full_reg[order(-full_reg$MutualInformation),]
write.csv(full_rego, "~/ferdig_rotation/regulon_validation/sizematch_nodups_nets/fullreg_size_names_order_nodups.csv", row.names=F, quote=F)

#order pearson edgelist
pear<-read.csv("~/ferdig_rotation/regulon_validation/hu_pearson_network_mat.csv")
row.names(pear) <- pear$X
pear_n <- pear[,-1]
edges <- flattenSquareMatrix(data.matrix(pear_n))
edges2 <- edges[order(-edges$cor),]
write.csv(edges2, "~/ferdig_rotation/regulon_validation/sizematch_nodups_nets/full_pearson.csv", row.names=F, quote=F)
#printed, changed names, removed duplicate edges with remove_hem_dups.py in code folder
#remove self loops pearson
nodup_pear <- read.table("~/ferdig_rotation/regulon_validation/sizematch_nodups_nets/full_pear_nodup_names.txt", sep=" ", as.is = T)

idex_vec <- vector()
for(i in 1:dim(nodup_pear)[1]){
  if(nodup_pear[i,1] == nodup_pear[i,2]){
    idex_vec <- c(idex_vec, i)
  }
}

pear_net2 <- nodup_pear[-idex_vec,]
write.csv(pear_net2[-1,], "~/ferdig_rotation/regulon_validation/sizematch_nodups_nets/full_pear_nodup_self_names.csv", row.names=F, quote=F)
