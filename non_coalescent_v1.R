library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)

n_randomisations <- 2 #for each gene tree the number of topology swaps that will be attempted
proceed_probability <- 0.5 #for each gene tree the probability that attempted topology swaps are successful
n_tips <- 30 #number of tips in species trees and gene trees
n_gene_trees <- 30 #number of gene trees to simulate for each species tree
locus_size <- 1500 #length in bp for each locus simulated along each gene tree
consistent_incongruence <- TRUE #if TRUE, topology swaps will ititiate on the same branch
n_reps <- 10 #total number of replicate experiments - number of species trees
consistent_clade_size <- 2  

##########

lambda <- 3
mu <- 0

###########

file_names <- paste(seq(1, n_gene_trees, 1), ".nexus", sep="")

##########

tip_vector <- seq(1, n_tips, 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

##########

for (a in 1:n_reps){

entire_tree <- sim.bd.taxa(n_tips, 1, lambda, mu, frac = 1, complete = TRUE, stochsampling = FALSE)
entire_tree <- entire_tree[[1]]

#########

gene_trees  <- vector("list", n_gene_trees)
for (i in 1:length(gene_trees)){
gene_trees[[i]] <- entire_tree
}

########

for(j in 1:n_randomisations){

if (consistent_incongruence == TRUE){
clade <- sample(entire_tree$tip.label, consistent_clade_size)

for(i in 1:length(gene_trees)){
proceed <- sample(seq(0, proceed_probability, 0.0001), 1)

if (proceed < proceed_probability){
node <- findMRCA(gene_trees[[i]], clade, type="node")
if (node != length(gene_trees[[i]]$tip.label) +1){
exclusion <- c(node, getDescendants(gene_trees[[i]], node), ancestors(phylo4(gene_trees[[i]]), node, "all"))
swap_node <- sample(gene_trees[[i]][[1]][,2], 1)

while(swap_node %in% exclusion){
swap_node <- sample(gene_trees[[i]][[1]][,2], 1)
}
gene_trees[[i]][[1]][,2][c(which(gene_trees[[i]][[1]][,2] == node), which(gene_trees[[i]][[1]][,2] == swap_node))] <- c(swap_node, node) 
write.tree(gene_trees[[i]], "tree.tre")
gene_trees[[i]] <- read.tree("tree.tre")
}
}
}

for (i in 1:length(gene_trees)){
for (x in 1:length(gene_trees[[i]]$tip.label)){
if (node.depth.edgelength(gene_trees[[i]])[[x]] < max(node.depth.edgelength(gene_trees[[i]]))){
gene_trees[[i]]$edge.length[[which(gene_trees[[i]][[1]][,2] == x)]]<- gene_trees[[i]]$edge.length[[which(gene_trees[[i]][[1]][,2] == x)]] + (max(node.depth.edgelength(gene_trees[[i]])) -node.depth.edgelength(gene_trees[[i]])[[x]])
}
}
}
}

if (consistent_incongruence == FALSE){

for(i in 1:length(gene_trees)){
proceed <- sample(seq(0, proceed_probability, 0.0001), 1)

if (proceed < proceed_probability){
node <- sample(gene_trees[[i]][[1]][,2], 1)
exclusion <- c(node, getDescendants(gene_trees[[i]], node), ancestors(phylo4(gene_trees[[i]]), node, "all"))
swap_node <- sample(gene_trees[[i]][[1]][,2], 1)

while(swap_node %in% exclusion){
swap_node <- sample(gene_trees[[i]][[1]][,2], 1)
}
gene_trees[[i]][[1]][,2][c(which(gene_trees[[i]][[1]][,2] == node), which(gene_trees[[i]][[1]][,2] == swap_node))] <- c(swap_node, node) 
write.tree(gene_trees[[i]], "tree.tre")
gene_trees[[i]] <- read.tree("tree.tre")
}
}

for (i in 1:length(gene_trees)){
for (x in 1:length(gene_trees[[i]]$tip.label)){
if (node.depth.edgelength(gene_trees[[i]])[[x]] < max(node.depth.edgelength(gene_trees[[i]]))){
gene_trees[[i]]$edge.length[[which(gene_trees[[i]][[1]][,2] == x)]]<- gene_trees[[i]]$edge.length[[which(gene_trees[[i]][[1]][,2] == x)]] + (max(node.depth.edgelength(gene_trees[[i]])) -node.depth.edgelength(gene_trees[[i]])[[x]])
}
}
}
}

}

#########

gene_clock_trees <- gene_trees
for (i in 1:length(gene_clock_trees)){
gene_clock_trees[[i]]$edge.length <- gene_clock_trees[[i]]$edge.length*0.05
}

sequences <- vector("list", n_gene_trees)
for (i in 1:length(sequences)){
sequences[[i]] <- simSeq(gene_trees[[i]], l = locus_size, type = "DNA", rate = 1)
}

#########

dir.create(paste(a, "gene_sequences", sep=""))
for (i in 1:length(file_names)){
write.phyDat(sequences[[i]], file = paste(a, paste("gene_sequences/", file_names[[i]], sep=""), sep=""), format = "nexus")
}

########

concat <- sequences[[1]]
for (i in 2:length(gene_trees)){
concat <- append(concat, sequences[[i]])
}

dir.create(paste(a, "concat_sequence", sep=""))
write.phyDat(concat, file = paste(a, paste("concat_sequence/", "concat.nexus", sep=""), sep=""), format = "nexus")

#########

dir.create(paste(a, "species_tree", sep=""))
write.tree(entire_tree, file = paste(a, paste("species_tree/", "species_tree.tre", sep=""), sep=""))

#########

species_tree_clades <- vector("list", length(entire_tree$tip.label) - 2)
for (i in 1:length(species_tree_clades)){
species_tree_clades[[i]] <- extract.clade(entire_tree, length(entire_tree$tip.label)+1+i)
}


species_tree_clade_support <- vector("list", length(species_tree_clades))
for (j in 1:length(species_tree_clades)){
for (i in 1:length(gene_trees)){
for (x in 1:(length(gene_trees[[i]]$tip.label)-2)){
if (setequal(extract.clade(gene_trees[[i]], length(gene_trees[[i]]$tip.label)+1+x)$tip.label, species_tree_clades[[j]]$tip.label) == TRUE){
species_tree_clade_support[[j]] <- append(species_tree_clade_support[[j]], 1)
}
}
}
}

for (i in 1:length(species_tree_clade_support)){
if (length(species_tree_clade_support[[i]]) > 0){
species_tree_clade_support[[i]] <- length(species_tree_clade_support[[i]])/length(gene_trees)
}	
if (length(species_tree_clade_support[[i]]) == 0){
species_tree_clade_support[[i]] <- 0
}
}

dir.create(paste(a, "clade_support", sep=""))
sink(paste(a, "clade_support/clade_support.txt", sep=""))
print(unlist(species_tree_clade_support))
sink()

#######

}

