source("latent_space_ftn.R")

library(foreach)
library(doParallel)

myCluster <- makeCluster(20, # number of cores to use
                         type = "PSOCK")
registerDoParallel(myCluster)
print(detectCores())


iter <- 200
sim_res1 <- sim_res2 <- sim_res3 <- sim_res4 <- sim_res5 <- sim_res6 <- list()
sim_res7 <- sim_res8 <- sim_res9 <- sim_res10 <- sim_res11 <- sim_res12 <- list()


sim_res1 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-0.5, beta = c(0.1,0.5), a = 0.1, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res2 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(0.1,3), a = 0.1, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res3 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(1,3), a = 0.1, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}


sim_res4 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-0.5, beta = c(0.1,0.5), a = 0.3, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res5 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(0.1,3), a = 0.3, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res6 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(1,3), a = 0.3, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res7 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-0.5, beta = c(0.1,0.5), a = 0.5, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res8 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(0.1,3), a = 0.5, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res9 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(1,3), a = 0.5, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res10 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-0.5, beta = c(0.1,0.5), a = 0.2, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res11 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(0.1,3), a = 0.2, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

sim_res12 <- foreach(i = 1:iter, .combine='rbind', .packages=c('network', 'MASS', 'RDS')) %dopar% {
  latentspace_SPRTBA(gen_latent_space(Sigma = NULL, maxunif = 1, sample.size = 2000, myseed = i, alpha = 0-3, beta = c(1,3), a = 0.2, b = 1)
                   , RDS.size = 500, number.of.seeds = 10, per.iter = 1000, correct = T)
}

save.image(file="latent_space_SPRTBA_dep_unif01.RData")

plot(sort(sim_res1$chisq.class.group), type = "l")
lines(sort(sim_res1$RDS.chisq), type = "l", col = 2)
lines(sort(sim_res1$RDS1per), type = "l", col = 3)
lines(sort(sim_res1$RDS2per), type = "l", col = 4)
abline(h = 0.05, lty = 2)

mean(sim_res1$chisq.class.group <= 0.05)
mean(sim_res1$RDS.chisq <= 0.05)
mean(sim_res1$RDS1per <= 0.05)
mean(sim_res1$RDS2per <= 0.05)

