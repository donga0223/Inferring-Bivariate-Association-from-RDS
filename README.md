# Inferring-Bivariate-Association-from-RDS

This repository refers to the paper Kim D., et al.(2020); Inferring Bivariate Association from Respondent-Driven Sampling Data. In this paper, we proposed a method to semi-parametrically estimate the null distributions of standard test statistics in the presence of sampling dependence, allowing for more valid statistical testing for dependence between pairs of variables within the sample. 

[data](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/tree/master/data) contains simulated network data for Section 4.3.3.


[R_functions](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/tree/master/R_functions) contains all of the R code for all simulations in this paper. 

[SPRTBA_function.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/SPRTBA_function.R) is a main function for SPRTBA. 

[Generate.data.1stMC.function.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/Generate.data.1stMC.function.R) contains a function for generating data based on 1st order Markov dependence, then code to do a SPRTBA test. (This code is for Sections 4.3.1 and 4.3.2)

[simulation.study.1stMC.noassociation.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/simulation.study.1stMC.noassociation.R) code for the simulation study for Section 4.3.1 and 4.3.2.

[simulation.study.1stMC.association.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/simulation.study.1stMC.association.R) code for the simulation study for Section 4.3.3.

[latent_space_ftn.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/latent_space_ftn.R) R functions for generating latent space model data and sampling using Respondent Driven sampling, then doing SPRTBA.

[simulation_latent_space_dep_unif01.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/simulation_latent_space_dep_unif01.R) obtaining simulation results.
