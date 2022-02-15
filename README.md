# The implications of incongruence between gene tree and species tree topologies for divergence time estimation

The scripts presented here are relevent to the article: Carruthers et al. The implications of incongruence between gene tree and species tree topologies for divergence time estimation. Syst Biol.  

## Generating simulated data
### Simple four taxon
Run "four_taxon/four_taxon_simulation.R" to perform the four taxon simulations discussed in the article.
#### Arguments: 
n_gene_trees: refers to the number of gene trees to simulate for a given species tree\
locus_size: refers to the length in bp of each locus\
probs: refers to the probavbility that gene trees are topologically incongruent with the species tree\






### To analyse the data as a concatenated alignment
Run "four_taxon/rev_shell.Rev" in RevBayes. 

### To estimate individual gene trees, as a basis for using the branch specific method 
Run "four_taxon/rev_shell_gene_trees.Rev" in RevBayes.


## Simple sixteen taxon simulations

### To simulate the data
Run "sixteen_taxon/simulation.R" to perform the simple sixteen taxon simulations discussed in the article.
#### Arguments:
tree: the species tree topology, and the topology of the topologically congruent gene trees\
incong_tree: the topology of the topologically incongruent gene trees\
n_gene_trees: the number of gene trees to simulate for a given species tree\
locus_size: the length in bp of each locus\

### To analyse the data as a concatenated alignment
Run "sixteen_taxon/rev_shell.Rev" in RevBayes. Note: set the root times to the correct age. Depends on whether the balanced tree or imbalanced tree is being analysed.


## Multispecies coalescent simulations

### To simulate the data
Run "multispecies_coalescent/coalescent.R". 
#### Arguments:
entire_tree: the species tree topology

### To analyse the simulated data as a concatenated alignment
Run "multispecies_coalescent/rev_shell.Rev"

### To estimate individual gene trees as a basis for the branch specific method
Run "rev_shell_by_gene.Rev" in RevBayes


## Empirical example

Tree inference are linked [here](https://github.com/pebgroup/Seed_Plant_BackBone).  

### To perform gene shopping
Run "empirical_example/gene_shopping.R" to identify loci with gene trees that are topologically congruent with the species tree. Ensure that script is searching correct directory for the gene trees.

### To estimate species tree branch lengths with congruent branches method
Run "empirical_example/branch_wise_analysis.R". Ensure that script is searching correct directory for gene trees. 

### To perform treePL analysis
Run relevant ".cppr8s" file. 

























