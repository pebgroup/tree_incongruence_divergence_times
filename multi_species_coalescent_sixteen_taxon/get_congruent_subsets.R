library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)

n_reps <- 50
n_loci <- 400
incongruent_incorporated <- 0
congruence_genes_per_species_tree <- vector(mode="numeric", length=0)

for (a in 1:n_reps){
species_tree <- read.tree("entire_tree_unbalanced.tre")

species_tree_clades <- vector("list", length(species_tree$tip.label) - 2)
for (i in 1:length(species_tree_clades)){
species_tree_clades[[i]] <- extract.clade(species_tree, length(species_tree$tip.label)+1+i)
}

gene_trees <- vector("list", n_loci)
for (i in 1:n_loci){
gene_trees[[i]] <- read.tree(paste("unbalanced/", a, "_gene_trees/", i, ".nexus", sep=""))
}

species_tree_clade_support <- vector("list", length(species_tree_clades))
for (j in 1:length(species_tree_clades)){
for (i in 1:length(gene_trees)){
for (x in 1:(length(gene_trees[[i]]$tip.label)-2)){
if (setequal(extract.clade(gene_trees[[i]], length(gene_trees[[i]]$tip.label)+1+x)$tip.label, species_tree_clades[[j]]$tip.label) == TRUE){
species_tree_clade_support[[j]] <- append(species_tree_clade_support[[j]], i)
}
}
}
}

for (j in 1:n_loci){
counter <- vector(mode="numeric", length=0)
for (i in 1:length(species_tree_clade_support)){
if (j %in% species_tree_clade_support[[i]]){
counter <- append(counter, 1)
}
}

if (length(counter) >= ((length(species_tree$tip.label)-2)-incongruent_incorporated)){
congruence_genes_per_species_tree <- append(congruence_genes_per_species_tree, j)
}
}

if ((length(congruence_genes_per_species_tree)) > 0){
high_support_loci <- vector("list", length(congruence_genes_per_species_tree))
for (i in 1:length(congruence_genes_per_species_tree)){
high_support_loci[[i]] <- read.phyDat(paste("unbalanced/", a, "_gene_sequences/", congruence_genes_per_species_tree[[i]], ".nexus", sep=""), format = "nexus")  
}
} 

if (length(high_support_loci) > 0){
concat <- high_support_loci[[1]]
if (length(high_support_loci) > 1){
for (i in 2:length(high_support_loci)){
concat <- append(concat, high_support_loci[[i]])
}
}
write.phyDat(concat, file = paste("unbalanced/", a, "_concat_sequence/congruent.nexus", sep = ""), format="nexus")
}
}
