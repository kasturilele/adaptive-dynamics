# adaptive-dynamics
This repository contains all the code used in the numerical solutions and calculations, population genetics simulations, and analysis and figure generation for the manuscript "Fitness landscapes for species interactions: when do population genetics and adaptive dynamics agree?".

Scripts used: 

Part A: Python and Mathematica scripts used for calculations and numerical solutions-
- eqn_solver.py - this script contains the code for trial tradeoff functions and pairwise invasibility plots, numerical solutions to selection gradients (derivative of invasion fitness), and predictions of coexistence over a range of constants for tradeoff functions. These results are discussed in Results section 1.
- AD_single.nb - this script contains Mathematica code to double-check the calculations from Appendix Section 3 (single-species adaptive dynamics predictions)
- AD_paired.nb - this script contains Mathematica code to double-check the calculations from Appendix Section 4 (two-species adaptive dynamics predictions)

Part B: SLiM scripts + bash scripts used for simulations-
- script9_exp_3.slim - this script runs single-species simulations over fixed time scales, for a determined range of mutation supply parameters (takes as input single_params_3.csv and mutation_kernel_extra.csv). compiled outputs from this file generate log_combined_single3.csv (data for figure 3)
- run_reps_sing3.sh and combine_file_sing.sh - bash scripts for running the previous SLiM script
- script10_exp_2.slim - this script runs two-species simulations over fixed time scales, for a determined range of mutation supply parameters (takes as input paired_params_2.csv and mutation_kernel_extra.csv). compiled outputs from this file generate log_combined_pair2new.csv (data for Figure 4)
- script10_3.slim - this script runs two-species simulations over a larger set of initial parameters, for a single set of mutation supply parameters (takes as input pairparams_new.csv and pairparams_mk_new.csv). compiled outputs from this file generate log_comb_init_00.csv (data for Figure 5)
- script10_2_expanded.slim - this script runs two-species simulations that terminate for predetermined end conditions, mutation rates sampled from distribution (takes as input paired_params_2.csv). compiled outputs from this file generate log_comb_exp_15_ext2.csv and log_comb_exp_16_ext2.csv (data for Figure 6)
- run_reps_pair_2.sh, run_reps_pair_expanded.sh and combine_file_sing.sh - bash scripts for running the previous SLiM scripts

Part C: R scripts used for analysis- 
The script AD_all_analysis.R contains all the R scripts used for statistical analysis and making figures. Some of the file locations need to be changed for the script to work correctly. 


Simulation data (only including data files for main figures):

 - tradeoff_fun_ESS1.csv and tradeoff_fun_ESS2.csv - files containing predictions of coexistence for different tradeoff function constants generated from the script eqn.solver.py
 - mutation_kernel_extra.csv - file containing all the mutation supply parameters (for most simulations)
 - mutation_kernel_expanded.csv - file with a smaller set of mutation kernel parameters, for script10_2_expanded.slim
 - pairparams_mk_new.csv - file with specific mutation parameters, for script10_3.slim
 - mutkern_15_ext2.csv and mutkern_16_ext2.csv - mutation supply parameters for simulations (for figure 6)
 - single_params_3.csv - file containing initial values for growth parameters (ri, aii) for single species simulations
 - paired_params_2.csv - file containing initial values for growth parameters (ri, aii) for two species simulations
 - pairparams_new.csv - file containing initial values for growth parameters (ri, aii) for extended two species simulations (figure 5)
 - log_combined_single3.csv - simulation outputs for single-species simulations (population size and mean values of ri and aii)
 - fixed_combined_single3.csv - information about fixed mutations for simulations in previous file
 - log_combined_pair2new.csv - simulation outputs for first set of two-species simulations (population size and mean values of ri and aii)
 - log_combined_extinct.csv - simulation outputs for some parameter combinations not covered in this previous file
 - fixed_combined_pair2new.csv - information about fixed mutations for simulations in previous file
 - log_comb_exp_15_ext2.csv and log_comb_exp_16_ext2.csv - simulation outputs for extended two-species simulations (for figure 6, population size and mean values of ri and aii)
 - log_comb_init_00.csv - simulation outputs for extended two-species simulations (for figure 5, population size and mean values of ri and aii)
