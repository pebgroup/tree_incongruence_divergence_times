# The implications of incongruence between gene tree and species tree topologies for divergence time estimation
The scripts presented here are relevent to the article: Carruthers et al. The implications of incongruence between gene tree and species tree topologies for divergence time estimation. Syst Biol.  

--- detailed instructions coming shortly ---

## Simple four taxon
### To generate simulated data run [simlpe_four_taxon_simulation.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/simple_four_taxon_simulation.R).
##### _Arguments:_
**n_gene_trees:** number of gene trees to simulate.\
**locus_size:** length in base pairs of each simulated locus.\
**probs:** the probability that gene tree topologies are incongruent with the species tree topology in each simulated dataset.

### To analyse simulated data [run overall.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/overall.Rev) in RevBayes. 
This script will perform analyses of:
1) **The entire dataset (concatenated)** - using [simple_clock_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/simple_clock_script.Rev) to estimate _t_, [more_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/more_relaxed_script_fixed.Rev) and [less_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/less_relaxed_script_fixed.Rev) to estimate _r_, and [simple_branch_length_estimation.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/simple_branch_length_estimation.Rev) to estimate _n_.
2) **Concatenated loci with gene trees that are topologically congruent with the species tree** - using same scripts as for 1
3) **Individual loci, such that parameters can be estimated from congruent branches in gene trees** - using [simple_clock_scripts_gene_trees.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_four_taxon/analyse/simple_clock_script_gene_trees.Rev)

## Simple sixteen taxon
### To generate simulated data run [simple_sixteen_taxon_simulation.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/simple_sixteen_taxon_simulation.R)
Simulation needs to be repeated for imbalanced and imbalanced tree, and different levels of topological incongruence.
##### _Arguments:_
**congruence_level:** output_file (needs to refer to whether tree is balanced, and level of incongruence)\
**tree:** species tree topology. Either [balanced_sixteen_taxon_cong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/balanced_sixteen_cong.tre) or [unbalanced_sixteen_taxon_cong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/unbalanced_sixteen_cong.tre)\
**incong_tree:** incongurent gene tree topology. Either [balanced_sixteen_one_incong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/balanced_sixteen_one_incong.tre), [balanced_sixteen_two_incong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/balanced_sixteen_two_incong.tre), [balanced_sixteen_three_incong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/balanced_sixteen_three_incong.tre), [balanced_sixteen_four_incong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/balanced_sixteen_four_incong.tre), [unbalanced_sixteen_incong.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/unbalanced_sixteen_incong.tre)\
**n_gene_trees:** number of gene trees to simulate\
**locus_size:** number of base pairs per locus

### To analyse simulated data run [overall.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/analyse/overall.Rev) in RevBayes.
Thie script will perform analyses of **the entire dataset (concatenated)** - using [simple_clock_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/analyse/simple_clock_script.Rev) to estimate _t_; [simple_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/analyse/simple_relaxed_script.Rev) and [simple_less_relaxed_script.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/analyse/simple_less_relaxed_script.Rev) to estimate _r_; and [simple_branch_length_estimation.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/simple_sixteen_taxon/analyse/simple_branch_length_estimation.Rev) to estimate _n_.

## Multi-species coalescent sixteen taxon
### To generate simulated data use [multispecies_coalescent_sixteen_taxon_simulation.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/multispecies_coalescent_sixteen_taxon_simulation.R)
##### _Arguments:_
**n_gene_trees:** number of gene trees to simulate\
**locus_size:** size in base pairs of each gene tree\
**n_reps:** number of times to repeat entire simulation\
**effective_population_size_approx:** _Ne_\
**n_tips:** number of taxa in species tree\
**species_tree:** species tree, either [entire_tree_unbalanced.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/entire_tree_unbalanced.tre), or [entire_tree_balanced.tre](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/entire_tree_balanced.tre)\
**type:** name of output file
### To generate concatenated alignment of simulated data based only on loci with gene trees that are topologically congruent with the species tree use [get_congruent_subsets.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/get_congruent_subsets.R)
### To analyse simulated data run [rev_shell.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/analyse/rev_shell.Rev)
This script will estimate _t_ in the balanced species tree with [entire_script_balanced.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/analyse/simple_clock_script_balanced.Rev); _t_ in the unbalanced species tree with [entire_script_unbalanced.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/analyse/simple_clock_script_unbalanced.Rev)
### also use [rev_shell_balanced_by_gene.Rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_sixteen_taxon/analyse/rev_shell_balanced_by_gene.Rev) to estimate _t_ in individual gene trees. 

## Multi-species coalescent four taxon
### To generate simulated data use [simulation.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_four_taxon/simulation.R)
##### _Arguments:_ 
**n_gene_trees:** number of gene trees to simulate\
**locus_size:** size in base pairs of each gene tree\
**n_reps:** number of times to repeat entire simulation\
**effective_population_size_approx:** _Ne_\
**n_tips:** number of taxa in species tree
### To generate concatenated alignment of simulated data based only on loci with gene trees that are topologically congruent with the species tree use [get_congruent_subsets.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_four_taxon/get_congruent_subsets.R)
### To generate start trees from initial four taxon simulation for analysis in multi species coalescent framework use [generate_start_trees_for_simple_analysed_as_coales.R](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_four_taxon/generate_start_trees_for_simple_analysed_as_coales.R)
### To analyse simulated data use [overall.rev](https://github.com/pebgroup/tree_incongruence_divergence_times/blob/master/multi_species_coalescent_four_taxon/analyse/overall.rev) 

## Empirical Example
Scripts for empirical example presented in the study.






