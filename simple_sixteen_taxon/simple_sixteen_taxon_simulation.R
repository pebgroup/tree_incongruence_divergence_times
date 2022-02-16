library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)

congruence_level <- "unbalanced_sixteen_cong"
tree <- read.tree("unbalanced_sixteen_cong.tre")
incong_tree <- read.tree("unbalanced_sixteen_cong.tre")
n_gene_trees <- 400
locus_size <- 800

####################
###INITIAL_PARAMS###
####################

file_names <- paste(seq(1, n_gene_trees, 1), ".nexus", sep="")
tree_file_names <- paste(seq(1, n_gene_trees, 1), ".tre", sep="")

#######################
###DEFINE_GENE_TREES###
#######################
	
gene_trees  <- vector("list", n_gene_trees)
for (i in 1:200){
gene_trees[[i]] <- tree
}

for (i in 201:400){
gene_trees[[i]] <- incong_tree
}

#############################
###DEFINE_GENE_CLOCK_TREES###
#############################

gene_clock_trees <- gene_trees
for (i in 1:length(gene_clock_trees)){
gene_clock_trees[[i]]$edge.length <- gene_clock_trees[[i]]$edge.length*0.05
}

########################
###SIMULATE_SEQUENCES###
########################

sequences <- vector("list", n_gene_trees)
for (i in 1:length(sequences)){
sequences[[i]] <- simSeq(gene_clock_trees[[i]], l = locus_size, type = "DNA", rate = 1)
}

###################
###WRITE_TO_FILE###
###################

dir.create("data_incongruence/")
dir.create(paste("data_incongruence/", congruence_level, sep=""))

for (i in 1:length(file_names)){
write.phyDat(sequences[[i]], file = paste("data_incongruence/", congruence_level, "/", file_names[[i]], sep=""), format = "nexus")
write.tree(gene_trees[[i]], file = paste("data_incongruence/", congruence_level, "/", tree_file_names[[i]], sep=""))
}

###########################
###CONCATENATE_SEQUENCES###
###########################

concat_incongruence <- sequences[[1]]
for (i in 2:length(gene_trees)){
concat_incongruence <- append(concat_incongruence, sequences[[i]])
}

write.phyDat(concat_incongruence, file = paste("data_incongruence/", congruence_level, "/concat.nexus", sep = ""), format = "nexus")

#########

write.tree(tree, file = paste("data_incongruence/", congruence_level, "/species_tree.tre", sep=""))