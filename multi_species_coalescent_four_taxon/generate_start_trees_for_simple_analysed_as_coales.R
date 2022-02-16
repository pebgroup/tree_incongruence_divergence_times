library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)

trees <- vector("list", 400)

for (i in 1:400){
trees[[i]] <- read.tree(paste("../simple_four_taxon/data_incongruence_0.5/", i, ".tre", sep=""))
}

for (i in 1:length(trees)){
trees[[i]]$edge.length[[1]] <- trees[[i]]$edge.length[[1]]+1
trees[[i]]$edge.length[[2]] <- trees[[i]]$edge.length[[2]]+3
trees[[i]]$edge.length[[3]] <- trees[[i]]$edge.length[[3]]+3
trees[[i]]$edge.length[[4]] <- trees[[i]]$edge.length[[4]]+0.5
trees[[i]]$edge.length[[5]] <- trees[[i]]$edge.length[[5]]+3.5
trees[[i]]$edge.length[[6]] <- trees[[i]]$edge.length[[6]]+3.5
}

dir.create("data/data_incongruence_0.5")
for (i in 1:length(trees)){
write.tree(trees[[i]], paste("data/data_incongruence_0.5/extended_start_", i, ".tre", sep=""))
}





