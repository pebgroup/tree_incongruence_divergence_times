library(phangorn)
library(devtools)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(coda)
library(ggplot2)

burnin <- 0.25
n_reps <- 10
n_tips <- 30

########

tip_vector <- seq(1, n_tips, 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

#########

congruence_storage <- vector(mode="numeric", length=0)
error_storage <- vector(mode="numeric", length=0)

#########

for (a in 1:n_reps){

correct_tree <- read.tree(paste(a, "species_tree/species_tree.tre", sep=""))
correct_branches <- correct_tree$edge.length
mpe_storage <- vector("list", length(correct_branches))

############

trees_table <- read.table(paste0("output/", paste0(a,"output.trees")))
trees <- vector("list", nrow(trees_table)-round(burnin*nrow(trees_table), digits = 0))
for (i in seq(round(burnin*nrow(trees_table), digits = 0), nrow(trees_table), 1)){
write(as.character(trees_table[,5][[i]]), "tree.tre")
trees[[i - (round(burnin*nrow(trees_table), digits = 0)-1)]] <- read.tree("tree.tre")
}

########

branch_posteriors <- vector("list", nrow(trees[[1]][[1]]))
for (j in 1:length(trees)){
for (i in 1:length(trees[[1]]$edge.length)){
branch_posteriors[[i]] <- append(branch_posteriors[[i]], trees[[j]]$edge.length[[i]])
}
}

#########

if (identical(trees[[1]][[1]], correct_tree[[1]]) == FALSE){
print("STOP_ANALYSIS")
print(a)
}

#########

mpe <- vector("list", length(branch_posteriors))
hpd <- vector("list", length(branch_posteriors))

for (i in 1:length(branch_posteriors)){
mpe[[i]] <- mean(branch_posteriors[[i]])
hpd[[i]] <- HPDinterval(as.mcmc(branch_posteriors[[i]]), prob = 0.95) 
}


##########

source(paste(a, "clade_support/clade_support_out.txt", sep=""))
clade_support <- c(rep(1, length(correct_tree$tip.label)), clade_support)
tree_table <- data.frame(correct_tree[[1]], correct_branches, unlist(mpe))
tree_table <- tree_table[order(tree_table[,2]),]
tree_table <- cbind(tree_table, clade_support)
error <- (tree_table[,4] - tree_table[,3])/tree_table[,3]
tree_table <- cbind(tree_table, error)

congruence_storage <- append(congruence_storage, tree_table[,5])
error_storage <- append(error_storage, tree_table[,6])

}


data_frame <- data.frame(congruence_storage, error_storage)

##########

result_plot <- ggplot(data = data_frame, aes(x = data_frame[,1], y =data_frame[,2])) +
geom_point(aes(y=data_frame[,2]), size = 1, shape = 19, fill = "grey") +
scale_y_continuous(breaks = seq(-1, 4, 0.1), labels = seq(-1, 4, 0.1)) +
scale_x_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) +
coord_cartesian(xlim = c(0,1), ylim = c(-1,4), FALSE)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(size=1, colour = "black"), text = element_text(size = 32, colour = "black"), axis.ticks = element_line(size=1, colour="black"))