library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(TreePar)

n_gene_trees <- 400 
locus_size <- 800 
n_reps <- 50 
effective_population_size_approx <- 120 
n_tips <- 16

###########

file_names <- paste(seq(1, n_gene_trees, 1), ".nexus", sep="")

##########

tip_vector <- seq(1, n_tips, 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

##########SIMULATE_SPECIES_TREE

for (z in 1:n_reps){

entire_tree <- read.tree("entire_tree_balanced.tre")

##########

coalescent_populations <- vector("list", (length(tip_vector)-1))
for (i in 1:length(coalescent_populations)){
coalescent_populations[[i]] <- extract.clade(entire_tree, node_numbers[[i]])$tip.label
}

coalescent_populations <- coalescent_populations[order(sapply(coalescent_populations,length),decreasing=F)]

##########MINIMUM_AND_MAXIMUM_AGES_OF_COALESCENT_POPULATION_AS_DETERMINED_BY_SPECIES_TREE

coalescent_populations_minimum_ages <- vector("list", length(coalescent_populations))
coalescent_populations_maximum_ages <- vector("list", length(coalescent_populations))
for (i in 1:length(coalescent_populations_minimum_ages)){
coalescent_populations_minimum_ages[[i]] <- max(node.depth.edgelength(extract.clade(entire_tree, findMRCA(entire_tree, coalescent_populations[[i]], "node"))))
if(length(coalescent_populations[[i]]) < length(entire_tree$tip.label)){
coalescent_populations_maximum_ages[[i]] <- max(node.depth.edgelength(extract.clade(entire_tree, entire_tree[[1]][,1][[which(entire_tree[[1]][,2] == findMRCA(entire_tree, coalescent_populations[[i]], "node"))]])))
}
}

#############

gene_trees <- vector("list", n_gene_trees)
for (y in 1:length(gene_trees)){

pool <- vector("list", length(entire_tree$tip.label))
for (i in 1:length(pool)){
pool[[i]] <- list(edge=matrix(c(2,1),1,2), tip.label=entire_tree$tip.label[[i]], edge.length=1.0, Nnode=1)
class(pool[[i]]) <- "phylo"
}

#############

coalescent_ages <- vector("list", length(coalescent_populations_minimum_ages))
for (i in 1:length(coalescent_ages)){

age_generation <- vector(mode="numeric", length=0)
for (a in 1:length(pool)){
counter <- vector(mode="numeric", length=0)
for (b in 1:length(pool[[a]]$tip.label)){
if (pool[[a]]$tip.label[[b]] %in% coalescent_populations[[i]]){
counter <- append(counter, 1)
}
}
if (length(counter) == length(pool[[a]]$tip.label)){
age_generation <- append(age_generation, 1)
}
}

coalescent_ages[[i]] <- coalescent_populations_minimum_ages[[i]] + rexp((length(age_generation)-1), 1/(effective_population_size_approx/1000))
coalescent_ages[[i]] <- sort(coalescent_ages[[i]])
remove <- vector(mode="numeric", length=0)
for (a in 1:length(coalescent_ages[[i]])){
if(is.numeric(coalescent_populations_maximum_ages[[i]]) == TRUE)
if(coalescent_ages[[i]][[a]] > coalescent_populations_maximum_ages[[i]]){
remove <- append(remove, a)
}
}

if(length(remove) > 0){
coalescent_ages[[i]] <- coalescent_ages[[i]][-remove]
}

if (length(coalescent_ages[[i]]) > 0){
for (a in 1:length(coalescent_ages[[i]])){
a
correct <- 0
if (length(pool) >= 2){
while (correct == 0){
counter <- vector(mode="numeric", length=0)
to_join_index <- sample(seq(1, length(pool), 1), 2, replace = FALSE)
to_join <- pool[to_join_index]
for (b in 1:length(c(to_join[[1]]$tip.label, to_join[[2]]$tip.label))){
if (c(to_join[[1]]$tip.label, to_join[[2]]$tip.label)[[b]] %in% coalescent_populations[[i]]){
counter <- append(counter, 1)
}
}
if (length(counter) == length(c(to_join[[1]]$tip.label, to_join[[2]]$tip.label))){
correct <- 1
}
}
if(correct == 1){
if ((length(to_join[[1]]$tip.label) == 1) & (length(to_join[[2]]$tip.label) == 1)){
to_join[[1]]$edge.length <- coalescent_ages[[i]][[a]]
to_join[[2]]$edge.length <- coalescent_ages[[i]][[a]]
}
if ((length(to_join[[1]]$tip.label) > 1) & (length(to_join[[2]]$tip.label) == 1)){
to_join[[1]] <- addroot(to_join[[1]], coalescent_ages[[i]][[a]] - max(node.depth.edgelength(to_join[[1]])))
to_join[[2]]$edge.length <- coalescent_ages[[i]][[a]]
}
if ((length(to_join[[1]]$tip.label) == 1) & (length(to_join[[2]]$tip.label) > 1)){
to_join[[1]]$edge.length <- coalescent_ages[[i]][[a]]
to_join[[2]] <- addroot(to_join[[2]], coalescent_ages[[i]][[a]] - max(node.depth.edgelength(to_join[[2]])))
}
if ((length(to_join[[1]]$tip.label) > 1) & (length(to_join[[2]]$tip.label) > 1)){
to_join[[1]] <- addroot(to_join[[1]], coalescent_ages[[i]][[a]] - max(node.depth.edgelength(to_join[[1]])))
to_join[[2]] <- addroot(to_join[[2]], coalescent_ages[[i]][[a]] - max(node.depth.edgelength(to_join[[2]])))
}
new <- bind.tree(to_join[[1]], to_join[[2]])
pool <- pool[-c(to_join_index[[1]], to_join_index[[2]])]
pool[[length(pool)+1]] <- new
}
}
}
}
}

gene_trees[[y]] <- pool[[1]]

}

#########

gene_clock_trees <- gene_trees
for (i in 1:length(gene_clock_trees)){
gene_clock_trees[[i]]$edge.length <- gene_clock_trees[[i]]$edge.length*0.05
}

sequences <- vector("list", n_gene_trees)
for (i in 1:length(sequences)){
sequences[[i]] <- simSeq(gene_clock_trees[[i]], l = locus_size, type = "DNA", rate = 1)
}

#########

dir.create("balanced")

dir.create(paste("balanced/", paste(z, "gene_sequences_medium", sep=""), sep = ""))
for (i in 1:length(file_names)){
write.phyDat(sequences[[i]], file = paste("balanced/", paste(z, paste("gene_sequences_medium/", file_names[[i]], sep=""), sep=""), sep = ""), format = "nexus")
}

dir.create(paste("balanced/", paste(z, "gene_trees_medium", sep=""), sep = ""))
for (i in 1:length(file_names)){
write.tree(gene_trees[[i]], file = paste("balanced/", paste(z, paste("gene_trees_medium/", file_names[[i]], sep=""), sep=""), sep = ""))
}

#########

concat <- sequences[[1]]
for (i in 2:length(gene_trees)){
concat <- append(concat, sequences[[i]])
}

dir.create(paste("balanced/", paste(z, "concat_sequence_medium", sep=""), sep = ""))
write.phyDat(concat, file = paste("balanced/", paste(z, paste("concat_sequence_medium/", "concat.nexus", sep=""), sep=""), sep = ""), format = "nexus")

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

dir.create(paste("balanced/", paste(z, "clade_support_medium", sep=""), sep = ""))
sink(paste("balanced/", paste(z, "clade_support_medium/clade_support.txt", sep=""), sep=""))
print(unlist(species_tree_clade_support))
sink()

#######

}
