# NCTSimulations
**Monte Carlo simulations to determine power and actual rejection size of the network comparison test for the difference in edge weights**

This repository contains scripts and functions that implement the Monte Carlo simulation to assess the power and actual rejection size of the permutation test that was implemented in the NetworkComparisonTest package of Van Borkulo.
It also allows to reproduce the Monte Carlo simulation described in the paper.

Reference:
de Rooij, B., et al. (2021). Symptom clusters in survivors of 7 seven cancer types from the PROFILES registry: A Network analysis. Manuscript submitted for publication. 

To reproduce the analysis of the influenza data discussed in the paper, take the following steps.

## 1. CREATE DATA

Either,
-Run the R script *R/SCRIPT_READDATA.R*. This script reads and processes the original SPSS data that were obtained from the PROFILES registry. The data are not publicly available but access can be requested. 
or,
-Run the R script *R/SCRIPT_MAKEDATA.R* to obtain artificial data with the same structure as the empirical data. These artificial data are mainly provided for convenience, to allow trying out the provided functions and scripts.

## 2. MONTE CARLO SIMULATIONS

1. Adapted Network Comparison Test for mixed graphical models, in case of both binary and continuous variables
  *R/Functions_NCTforMGM.R*
2. Script for performing the Network Comparison Test and the Monte Carlo simulation
  *R/Script_MCsimulation_size_power.R*
