# Inferring-Bivariate-Association-from-RDS

This repository refers to the paper Dongah Kim, Krista J. Gile (2020). Inferring Bivariate Association from Respondent-Driven Sampling Data. In this paper, we proposed a method to semi-parametrically estimate the null distributions of standard test statistics in the presence of sampling dependence, allowing for more valid statistical testing for dependence between pairs of variables within the sample. 

[R_functions](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/tree/master/R_functions) contains all of the R code for all simulations we run on this paper. 
[SPRTBA_function.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/SPRTBA_function.R) is a main function for SPRTBA. 
[Generate.data.1stMC.function.R](https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/Generate.data.1stMC.function.R) contains a function for generating data based on 1st order Markov dependence, then do a SPRTBA test. (This code is for Chapter 4.3.1 and 4.3.2)
[simulation.study.1stMC.noassociation.R]https://github.com/donga0223/Inferring-Bivariate-Association-from-RDS/blob/master/R_functions/simulation.study.1stMC.noassociation.R contains the simulation study code.
