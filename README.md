# The implications of incongruence between gene tree and species tree topologies for divergence time estimation
The scripts presented here are relevent to the article: Carruthers et al. The implications of incongruence between gene tree and species tree topologies for divergence time estimation. Syst Biol.  

--- detailed instructions coming shortly ---

## Simple four taxon
#### To generate simulated data [run four_taxon_simulation.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/simple_four_taxon_simulation.R).
**n_gene_trees:** number of gene trees to simulate.\
**locus_size:** length in base pairs of each simulated locus.\
**probs:** the probability that gene tree topologies are incongruent with the species tree topology in each simulated dataset.\

#### To analyse simulated data [run overall.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/overall.Rev) in RevBayes. 
This script will perform analyses of:
1) The entire dataset\ 
........using [simple_clock_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/simple_clock_script.Rev) to estimate _t_, [more_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/more_relaxed_script_fixed.Rev) and [less_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/less_relaxed_script_fixed.Rev) to estimate _r_, and [simple_branch_length_estimation.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/simple_branch_length_estimation.Rev) to estimate _n_.
3) Concatenated loci with gene trees that are topologically congruent with the species tree
4) Individual loci, such that parameters can be estimated from congruent branches in gene trees

## Simple sixteen taxon
Use for generating simple sixteen taxon simulated data, and analysing simulated data. 

## Multi-species coalescent sixteen taxon
Use for generating data in a multispecies coalescent framework with 16 taxon species tree, and analysing simulated data.

## Multi-species coalescent sixteen taxon
Use for generating data in a multispecies coalescent framework with 4 taxon species tree, and analysing simulated data.

## Empirical Example
Scripts for empirical example presented in the study.






