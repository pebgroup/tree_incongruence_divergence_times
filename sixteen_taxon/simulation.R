library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)

congruence_level <- "sixteen_four_incong"
tree <- read.tree("../sixteen_cong.tre")
incong_tree <- read.tree("../sixteen_one_incong.tre")
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
gene_trees[[i]] <- read.tree(paste(paste("../", congruence_level, sep = ""), ".tre", sep=""))
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

dir.create(congruence_level)

for (i in 1:length(file_names)){
write.phyDat(sequences[[i]], file = paste(congruence_level, paste("/", file_names[[i]], sep=""), sep=""), format = "nexus")
write.tree(gene_trees[[i]], file = paste(congruence_level, paste("/", tree_file_names[[i]], sep=""), sep=""))
}

###########################
###CONCATENATE_SEQUENCES###
###########################

concat_incongruence <- sequences[[1]]
for (i in 2:length(gene_trees)){
concat_incongruence <- append(concat_incongruence, sequences[[i]])
}

write.phyDat(concat_incongruence, file = paste(congruence_level, paste("/", "concat.nexus", sep=""), sep=""), format = "nexus")

#########

write.tree(tree, file = paste(congruence_level, paste("/", "species_tree.tre", sep=""), sep=""))