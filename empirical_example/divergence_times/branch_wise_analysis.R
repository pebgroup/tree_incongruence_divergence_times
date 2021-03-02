library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(coda)
library(ggplot2)
library(Rmisc)

###############
###GET_TREES###
###############

`%notin%` <- Negate(`%in%`)
n_gene_trees <- 351
gene_trees <- vector("list", n_gene_trees)
outgroup_names <- vector("list", n_gene_trees)
setwd("../gene_shopping/gene_trees")
for (a in 1:n_gene_trees){
gene_trees[[a]] <- read.tree(list.files(pattern = "FBP.rt.tre")[[a]])
}

overall_species_tree <- read.tree("../species_rt_tree.tre")

######################
###DO_CORE_ANALYSIS###
######################

species_tree_branch_estimate <- vector("list", nrow(overall_species_tree[[1]]))

for (a in 1:length(gene_trees)){

ignore <- which(gene_trees[[a]][[1]][,1] == length(gene_trees[[a]]$tip.label)+1)
use <- seq(1, nrow(gene_trees[[a]][[1]]), 1)[-ignore]

for (b in use){

if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
gene_trees_descendant <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,2][[b]])
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
} else {
gene_trees_descendant <- gene_trees[[a]]$tip.label[[gene_trees[[a]][[1]][,2][[b]]]]
gene_trees_ancestral <- extract.clade(gene_trees[[a]], gene_trees[[a]][[1]][,1][[b]])
}

if (gene_trees[[a]][[1]][,2][[b]] > length(gene_trees[[a]]$tip.label)){
species_tree_descendant <- findMRCA(overall_species_tree, gene_trees_descendant$tip.label, "node")
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
} else {
species_tree_descendant <- which(overall_species_tree$tip.label == gene_trees_descendant)
species_tree_ancestral <- findMRCA(overall_species_tree, gene_trees_ancestral$tip.label, "node")
}

for (c in 1:nrow(overall_species_tree[[1]])){
if (overall_species_tree[[1]][,2][[c]] == species_tree_descendant & overall_species_tree[[1]][,1][[c]] == species_tree_ancestral){
species_tree_branch_estimate[[c]] <- append(species_tree_branch_estimate[[c]], gene_trees[[a]]$edge.length[[b]])
}
}

}
}

########################
###BRANCH_LENGTH_TREE###
########################

branch_length_tree <- overall_species_tree

for (i in 1:length(species_tree_branch_estimate)){
species_tree_branch_estimate[[i]] <- mean(species_tree_branch_estimate[[i]])
}

for (i in 1:length(species_tree_branch_estimate)){
if (is.na(species_tree_branch_estimate[[i]]) == TRUE){
species_tree_branch_estimate[[i]] <- 0
}
}


branch_length_tree$edge.length <- unlist(species_tree_branch_estimate)