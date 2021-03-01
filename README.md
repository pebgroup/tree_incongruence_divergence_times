# The implications of topological incongruence between gene trees and the species tree for divergence time estimation

The scripts presented here are relevent to the article: The implications of incongruence between gene trees and the species tree for divergence time estimation 
The branch specific approach discussed in the article will be available to download as an r package shortly.

## Simple four taxon simulations

### To simulate the data
Run "four_taxon/four_taxon_simulation.R" to perform the four taxon simulations discussed in the article.

#### Arguments: 
~ n_gene_trees refers to the number of gene trees to simulate for a given species tree
~ locus size refers to the length in bp of each locus
~ probs refers to the probavbility that gene trees are topologically incongruent with the species tree

### To analyse the data as a concatenated alignment
Run "four_taxon/rev_shell.Rev" in RevBayes. 

### To estimate individual gene trees, as a basis for using the branch specific method 
Run "four_taxon/rev_shell_gene_trees.Rev" in RevBayes.

## Simple sixteen taxon simulations
Run "sixteen_taxon/
