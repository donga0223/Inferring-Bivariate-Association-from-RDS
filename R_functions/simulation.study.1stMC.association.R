## This is R code for getting results for our paper Section 4.3.3. 

source("getRDSsample.R")
source("rdssamplecode.krista.R")
source("as.network.uncompressed.R")
source("SPRTBA_function.R")

library(RDS)
library(network)


asso.SPRTBA <- function(data, iter, sample.size = 500, number.of.seeds = 10, per.iter = 100, correct = FALSE){
  
  net <- data
  cg <- get.vertex.attribute(net, "cg")
  class <- get.vertex.attribute(net, "class")
  group <- get.vertex.attribute(net, "group")
  id <- get.vertex.attribute(net, "vertex.names")
  
  deg <-sapply(net$iel,length)+sapply(net$oel,length)
  head(deg)
  
  undirect.edge2 <- undirect.edge <- as.edgelist(net)
  colnames(undirect.edge2) <- colnames(undirect.edge) <- c("id1", "id2")
  undirect.edge2[,1] <- undirect.edge[,2]
  undirect.edge2[,2] <- undirect.edge[,1]
  undirect.edge <- undirect.edge3 <- rbind(undirect.edge, undirect.edge2)
  undirect.edge <- undirect.edge[order(undirect.edge[,1],undirect.edge[,2]),]
  
  net.data <- data.frame(id, class, group, cg, deg)
  
  net.data.att <- merge(undirect.edge, net.data, by.x = "id2", by.y = "id", all.x = TRUE)
  net.data.att <- merge(net.data.att, net.data, by.x = "id1", by.y = "id",all.x = TRUE, suffixes = c(".2",".1"))
  #head(net.data.att)
  
  ####################
  chisq.test.att.class <- chisq.test(net.data.att$class.1, net.data.att$class.2, correct = correct)
  chisq.test.att.group <- chisq.test(net.data.att$group.1, net.data.att$group.2, correct = correct)
  
  chisq.test.class.group <- chisq.test(net.data$class, net.data$group, correct = correct)
  
  
  chisq.res.net <- rep(NA, iter)
  permute.res1.mc.net <- rep(NA, iter)
  permute.res2.mc.net <- rep(NA, iter)
  
  res1.permute.var.11 <- rep(NA, iter)
  res1.permute.var.22 <- rep(NA, iter)
  res1.fix.var.11 <- rep(NA, iter)
  res1.fix.var.22 <- rep(NA, iter)
  res2.permute.var.11 <- rep(NA, iter)
  res2.permute.var.22 <- rep(NA, iter)
  res2.fix.var.11 <- rep(NA, iter)
  res2.fix.var.22 <- rep(NA, iter)
  
  for (i in 1:iter){
    set.seed(i)
    net.RDS <- getRDSsample(net, "group", sample.size = sample.size, number.of.seeds = number.of.seeds)
    N <- get.population.size(net.RDS)
    net.RDS <- as.data.frame(net.RDS)
    merge.net.data <- merge(net.RDS, net.data, by.x = "id", by.y = "id", all.x = TRUE)
    
    net.RDS.data <- as.rds.data.frame(merge.net.data, id = "id"
                                      , recruiter.id = "recruiter.id"
                                      , network.size = "degree" 
                                      , population.size = N)
    chisq.res.net[i] <- chisq.test(net.RDS.data$group, net.RDS.data$class, correct = correct)$p.value
    
    permute.res1.net <- SPRTBA(net.RDS.data, "group", "class", chisq.test, per.iter, correct = correct)
    permute.res1.mc.net[i] <- permute.res1.net$p.value.mc
    res1.permute.var.11[i] <- permute.res1.net$permute.var.prop[1,1]
    res1.permute.var.22[i] <- permute.res1.net$permute.var.prop[2,2]
    res1.fix.var.11[i] <- permute.res1.net$fix.var.prop[1,1]
    res1.fix.var.22[i] <- permute.res1.net$fix.var.prop[2,2]
    
    permute.res2.net <- SPRTBA(net.RDS.data, "class", "group", chisq.test, per.iter, correct = correct)
    permute.res2.mc.net[i] <- permute.res2.net$p.value.mc
    res2.permute.var.11[i] <- permute.res2.net$permute.var.prop[1,1]
    res2.permute.var.22[i] <- permute.res2.net$permute.var.prop[2,2]
    res2.fix.var.11[i] <- permute.res2.net$fix.var.prop[1,1]
    res2.fix.var.22[i] <- permute.res2.net$fix.var.prop[2,2]
    
    rm(net.RDS, merge.net.data, net.RDS.data)
    print(i)
  }
  
  res.net <- data.frame(chisq.res.net, permute.res1.mc.net, permute.res2.mc.net
                        ,res1.permute.var.11 
                        ,res1.permute.var.22 
                        ,res1.fix.var.11 
                        ,res1.fix.var.22 
                        ,res2.permute.var.11 
                        ,res2.permute.var.22 
                        ,res2.fix.var.11 
                        ,res2.fix.var.22 )
  
  names(res.net) <- c("chisq.res", "permute.res1.mc", "permute.res2.mc"
                      ,"res1.per.11" ,"res1.per.22", "res1.fix.11", "res1.fix.22"
                      ,"res2.per.11", "res2.per.22", "res2.fix.11", "res2.fix.22")
  
  
  
  return(list(chisq.test.att.class =  chisq.test.att.class
              ,chisq.test.att.group = chisq.test.att.group
              ,chisq.test.class.group = chisq.test.class.group
              ,res.net = res.net))
}


name <- as.character(5:13)

for(i in 1:length(name)){
  load(paste("pops",name[i],"RData",sep="."))
  data <- get(paste("pops", name[i], sep="."))
  
  iter <- 1000
  sample.size <- 500
  number.of.seeds <- 10
  per.iter <- 100
  
  assign(paste("sim.association", name[i], sep=".")
         , asso.SPRTBA(data, iter, sample.size = sample.size, correct = FALSE))
  
  save.image(file=paste("sim.association", name[i], "RData", sep = "."))
}

