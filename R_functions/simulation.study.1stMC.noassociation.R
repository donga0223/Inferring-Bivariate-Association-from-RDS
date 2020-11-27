## This is R code for getting results for our paper Section 4.3.1. and 4.3.2. 

##########################################################################
####### Simple Two variables; Generate the data
##########################################################################
### Independent Case

library(RDS)
library(network)

source("SPRTBA_function.R")
source("Generate.data.1stMC.function.R")

AA <- c(0.91, 0.75, 0.85, 0.74, 0.88, 0.82, 0.95, 0.85, 0.74, 0.63, 0.89, 0.79, 0.85)
BB <- c(0.56, 0.9, 0.82, 0.81, 0.82, 0.61, 0.63, 0.54, 0.81, 0.89, 0.85, 0.82, 0.74)
x00 <- c(0.48, 0.48, 0.64, 0.46, 0.5, 0.55, 0.6, 0.59, 0.62, 0.49, 0.62, 0.52, 0.75)
x11 <- c(0.65, 0.63, 0.46, 0.63, 0.64, 0.73, 0.52, 0.48, 0.6, 0.66, 0.5, 0.53, 0.43)

name <- as.character(1:13)

for(i in 1:length(name)){
  iter <- 1000
  number.of.seeds <- 10
  sample.size <- 500
  
  
  assign(paste("sim.study.1stMC", name[i], "bothstrong", sep=".")
         , Generate.data(iter, number.of.seeds, sample.size, AA[i], BB[i], AA[i], BB[i], correct = FALSE, per.iter = 1000))
  
  save.image(file=paste("sim.study.1stMC", name[i], "bothstrong", "RData", sep = "."))
  
}


for(i in 1:length(name)){
  iter <- 1000
  number.of.seeds <- 10
  sample.size <- 500
  
  
  assign(paste("sim.study.1stMC", name[i], "bothweak", sep=".")
         , Generate.data(iter, number.of.seeds, sample.size, x00[i], x11[i], x00[i], x11[i], correct = FALSE, per.iter = 1000))
  
  save.image(file=paste("sim.study.1stMC", name[i], "bothweak", "RData", sep = "."))
  
}

for(i in 1:length(name)){
  iter <- 1000
  number.of.seeds <- 10
  sample.size <- 500
  
  
  assign(paste("sim.study.1stMC", name[i], "onestrongoneweak", sep=".")
         , Generate.data(iter, number.of.seeds, sample.size, AA[i], BB[i], x00[i], x11[i], correct = F, per.iter = 1000))
  
  save.image(file=paste("sim.study.1stMC", name[i], "onestrongoneweak", "RData", sep = "."))
  
}

