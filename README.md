# adaptive-dynamics
This repository contains all the code used in the numerical solutions and calculations, population genetics simulations, and analysis and figure generation for the manuscript "Fitness landscapes for species interactions: when do population genetics and adaptive dynamics agree?". 

Scripts used: Part A: Python and Mathematica scripts used for calculations and numerical solutions-
eqn_solver.py - this script contains the code for trial tradeoff functions and pairwise invasibility plots, numerical solutions to selection gradients (derivative of invasion fitness), and predictions of coexistence over a range of constants for tradeoff functions. These results are discussed in Results section 1.
AD_single.nb - this script contains Mathematica code to double-check the calculations from Appendix Section 3 (single-species adaptive dynamics predictions)
AD_paired.nb - this script contains Mathematica code to double-check the calculations from Appendix Section 4 (two-species adaptive dynamics predictions)

Part B: SLiM scripts + bash scripts used for simulations

Part C: R scripts used for analysis- 
The script AD_all_analysis.R contains all the R scripts used for statistical analysis and making figures. Some of the file locations need to be changed for the script to work correctly. 
