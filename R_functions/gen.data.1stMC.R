
##########################################################################
####### Simple Two variables; Generate the data
##########################################################################
### Independent Case

library(RDS)
library(network)

source("sim.permutation.test.R")
source("Generate.data.function.R")


AA <- c(.63, .91, .89, .95, .75, .85, .88, .74, .74, .79, .85, .85, .82)
BB <- c(.89, .56, .85, .63, .9, .54, .82, .81, .81, .82, .82, .74, .61)
x00 <- c(.49, .48, .62, .6, .48, .59, .5, .62, .46, .52, .64, .75, .55)
x11 <- c(.66, .65, .5, .52, .63, .48, .64, .60, .63, .53, .46, .43, .73)

name <- c("01", "03", "04", "06", "08", "09", "10", "12", "13", "14", "15", "16", "17")

for(i in 1:length(name)){
  iter <- 1000
  number.of.seeds <- 10
  sample.size <- 500
  
  
  assign(paste("gen.data.1stMC", name[i], "bothstrong", sep=".")
         , Generate.data(iter, number.of.seeds, sample.size, AA[i], BB[i], AA[i], BB[i], correct = F, per.iter = 1000))
  
  save.image(file=paste("gen.data.1stMC", name[i], "bothstrong", "RData", sep = "."))
  
}



##############################################################################
## To test for other cases, change the assign and save.image parts as follows:
##############################################################################

## for one weak and one strong ###############################################
#assign(paste("gen.data.1stMC", name[i], "onestrongoneweak", sep=".")
#       , Generate.data(iter, number.of.seeds, sample.size, AA[i], BB[i], x00[i], x11[i], correct = F, per.iter = 1000))

#save.image(file=paste("gen.data.1stMC", name[i], "onestrongoneweak", "RData", sep = "."))


## for both weak ###############################################

#assign(paste("gen.data.1stMC", name[i], "onestrongoneweak", sep=".")
#       , Generate.data(iter, number.of.seeds, sample.size, x00[i], x11[i], x00[i], x11[i], correct = F, per.iter = 1000))

#save.image(file=paste("gen.data.1stMC", name[i], "bothweak", "RData", sep = "."))

