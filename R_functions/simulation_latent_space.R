## This function shows how we get the result values for Section 4.4.
## Setting maxunif = 1 generates weakly clustered network data.
## Setting maxunif = 15 generates stronly clustered network data.
## beta requires two values for beta1 and beta2.
### As the beta approaches zero, it implies weak homophily. If beta >>1, it implies strong homophily.  
## a = 0 : two variables do not have association.
## a != 0, two variables have association.
## As ``a'' value increases, the two variables have a strong association.
## Sigma == NULL : using Uniform distribution

source("latent_space_ftn.R")

library(RDS)
library(MASS)
library(network)

##### This is one example for getting result for Section 4.4. 
## You can change maxunif, alpha, beta, a, b values for getting other results.
res1 <- latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-0.5, beta = c(0.1,0.5), a = 0.1, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)

