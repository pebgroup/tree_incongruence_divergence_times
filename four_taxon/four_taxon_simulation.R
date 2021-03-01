library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)

tree <- read.tree("four_taxon_tree.tre")
n_gene_trees <- 400
locus_size <- 800
probs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

###########

file_names <- paste(seq(1, n_gene_trees, 1), ".nexus", sep="")
tree_file_names <- paste(seq(1, n_gene_trees, 1), ".tre", sep="")

###########

for (a in 1:length(probs)){
	
gene_trees  <- vector("list", n_gene_trees)
for (i in 1:length(gene_trees)){
gene_trees[[i]] <- tree
}

##########

trees_to_distort <- sample(n_gene_trees, n_gene_trees*probs[[a]])

#########

alternative_labels <- vector("list", 2)
alternative_labels[[1]] <- c("A", "C", "B", "D")
alternative_labels[[2]] <- c("A", "D", "B", "C")

if (length(trees_to_distort) > 0){
for (i in trees_to_distort){
gene_trees[[i]]$tip.label <- unlist(sample(alternative_labels, 1))
}
}

#########

gene_clock_trees <- gene_trees
for (i in 1:length(gene_clock_trees)){
gene_clock_trees[[i]]$edge.length <- gene_clock_trees[[i]]$edge.length*0.05
}

##########

sequences <- vector("list", n_gene_trees)
for (i in 1:length(sequences)){
sequences[[i]] <- simSeq(gene_clock_trees[[i]], l = locus_size, type = "DNA", rate = 1)
}

#########

dir.create(paste("data_incongruence_", probs[[a]], sep=""))
for (i in 1:length(file_names)){
write.phyDat(sequences[[i]], file = paste(paste("data_incongruence_", probs[[a]], sep=""), paste("/", file_names[[i]], sep=""), sep=""), format = "nexus")
write.tree(gene_trees[[i]], file = paste(paste("data_incongruence_", probs[[a]], sep = ""), paste("/", tree_file_names[[i]], sep=""), sep=""))
}

########

concat_incongruence <- sequences[[1]]
for (i in 2:length(gene_trees)){
concat_incongruence <- append(concat_incongruence, sequences[[i]])
}

write.phyDat(concat_incongruence, file = paste(paste("data_incongruence_", probs[[a]], sep = ""), paste("/", "concat.nexus", sep=""), sep=""), format = "nexus")

#########

write.tree(tree, file = paste(paste("data_incongruence_", probs[[a]], sep = ""), paste("/", "species_tree.tre", sep=""), sep=""))

#########

}