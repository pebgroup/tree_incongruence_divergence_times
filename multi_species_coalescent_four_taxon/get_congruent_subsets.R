library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)

species_tree <- read.tree("entire_tree_balanced.tre")
n_gene_trees <- 200

gene_trees <- vector("list", n_gene_trees)
for (i in 1:n_gene_trees){
gene_trees[[i]] <- read.tree(paste("data/gene_trees/", i, ".nexus", sep=""))
}

species_tree_clades <- vector("list", 2)
species_tree_clades[[1]] <- sort(extract.clade(species_tree, 6)$tip.label)
species_tree_clades[[2]] <- sort(extract.clade(species_tree, 7)$tip.label)

contain <- vector(mode="numeric", length=0)
for (i in 1:length(gene_trees)){
gene_tree_clades <- vector("list", 2)
gene_tree_clades[[1]] <- sort(extract.clade(gene_trees[[i]], 6)$tip.label)
gene_tree_clades[[2]] <- sort(extract.clade(gene_trees[[i]], 7)$tip.label)
if (((identical(species_tree_clades[[1]], gene_tree_clades[[1]]) == TRUE) || (identical(species_tree_clades[[1]], gene_tree_clades[[2]]) == TRUE))
&
((identical(species_tree_clades[[2]], gene_tree_clades[[1]]) == TRUE) || (identical(species_tree_clades[[2]], gene_tree_clades[[2]]) == TRUE))){
contain <- append(contain, i)
}
}

contain_loci <- vector("list", length(contain))
for (i in 1:length(contain)){
contain_loci[[i]] <- read.phyDat(paste("data/gene_sequences/", contain[[i]], ".nexus", sep = ""), format = "nexus", type ="DNA")
} 

concat <- contain_loci[[1]]
for (i in 2:length(contain_loci)){
concat <- append(concat, contain_loci[[i]])
}

write.phyDat(concat, file = "data/concat_sequence/congruent.nexus", format = "nexus")
