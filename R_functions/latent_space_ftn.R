library(network)
library(MASS)
if(2==3){
  Sigma <- matrix(c(100,0,0,100),2,2)
  Sigma = NULL
  sample.size = 100
  myseed = 1234
  alpha = -min(mybeta[c(1,2)])
  beta = c(1,2)
  a = 0
  b = 1
  net <- gen_latent_space(Sigma = NULL, maxunif =10, sample.size = 2000, myseed = 1234, alpha = 2, beta = c(1,2), a = 0, b = 1)
    
}


#a = 0 independent
#a != 0 dependent
# Sigma = NULL -> two independent uniform distribution
# Sigma not NULL and 2 by 2 matrix -> multivariate normal
# beta : can control the homophily effects. 
## When beta closes to 0, there is no homophily. WHen beta >> 1 , there is a strong homophily.
gen_latent_space <- function(Sigma = NULL, maxunif =10, sample.size = 2000, myseed = myseed, alpha = 2, beta = c(1,2), a = 0, b = 1){
  set.seed(myseed)
  if(is.null(Sigma)==FALSE){
    myspace <- mvrnorm(n = sample.size, rep(0, 2), Sigma)
  }else if(is.null(Sigma)){
    myspace <- matrix(runif(sample.size*2, 0, maxunif), ncol = 2)
  }
  myspace <- as.data.frame(myspace)
  colnames(myspace) <- c("nodeX1", "nodeX2")
  head(myspace)
  plot(myspace)
  summary(myspace)
  
  myspace$nodeX22 <- a*myspace$nodeX1 + b*myspace$nodeX2
  myspace$X1p <- (myspace$nodeX1-min(myspace$nodeX1))/(max(myspace$nodeX1)-min(myspace$nodeX1))
  myspace$X2p <- (myspace$nodeX22-min(myspace$nodeX22))/(max(myspace$nodeX22)-min(myspace$nodeX22))
  #summary(myspace)
  
  myspace$X1 <- rbinom(sample.size,1,myspace$X1p)
  myspace$X2 <- rbinom(sample.size,1,myspace$X2p)
  head(myspace)
  chisq.test(myspace$X1, myspace$X2)
  
  X1ij <- X2ij <- distf <- etaij <- Pyij <- yij <- matrix(0, nrow = sample.size, ncol = sample.size)
  for(i in 1:sample.size){
    if(i < sample.size){
      for(j in (i+1):sample.size){
        X1ij[i,j] <- X1ij[j,i] <- ifelse(myspace$X1[i] == myspace$X1[j], 1, 0)
        X2ij[i,j] <- X2ij[j,i] <- ifelse(myspace$X2[i] == myspace$X2[j], 1, 0)
        distf[i,j] <- distf[j,i] <- sqrt((myspace$nodeX1[i]-myspace$nodeX1[j])^2 + (myspace$nodeX2[i]-myspace$nodeX2[j])^2)
        etaij[i,j] <- etaij[j,i] <- alpha + beta[1]*X1ij[i,j] + beta[2]*X2ij[i,j] - distf[i,j]
        Pyij[i,j] <- Pyij[j,i] <- exp(etaij[i,j] - log(1+exp(etaij[i,j])))
        yij[i,j] <- yij[j,i] <- rbinom(1,1,Pyij[i,j])
      } 
    }else{}
  }
  
  #Pyij[1:20,1:20]
  #yij[1:20,1:20]
  #table(apply(yij, 1, sum))
  yij <- as.network(yij, matrix.type = "adjacency", directed = FALSE)
  set.vertex.attribute(yij, "group", myspace$X1)
  set.vertex.attribute(yij, "class", myspace$X2)
  
  return(yij)
}


source("getRDSsample.R")
source("rdssamplecode.krista.R")
source("as.network.uncompressed.R")
source("SPRTBA_function.R")

library(RDS)

if(2 == 3){
  #net <- mydata
  RDS.size <- 500
  number.of.seeds <- 10
  per.iter <- 100
  correct = T
  net <- gen_latent_space(Sigma = matrix(c(100,0,0,100),2,2), sample.size = 500, myseed = i, alpha = 2, beta = c(1,2), a = 0, b = 1)
  net <- gen_latent_space(Sigma = NULL, maxunif =1, sample.size = 500, myseed = i, alpha = -5, beta = c(0,1), a = 0, b = 1)
  net <- gen_latent_space(Sigma = NULL, maxunif =15, sample.size = 1000, myseed = i, alpha = -min(c(3,5)), beta = c(3,5), a = 0, b = 1)
  net
  plot(net)
  latentspace_SPRTBA(net, RDS.size = 100, number.of.seeds = 10, per.iter = 100, correct = T)
}


##net : latent space model data generated from the ``gen_latent_space'' function.
## RDS.size : number of sample size sampling RDS method.

latentspace_SPRTBA <- function(net, RDS.size = 500, number.of.seeds = 10, per.iter = 100, correct = T){
  group <- get.vertex.attribute(net, "group")
  class <- get.vertex.attribute(net, "class")
  id <- get.vertex.attribute(net, "vertex.names")
  
  deg <-sapply(net$iel,length)+sapply(net$oel,length)
  head(deg)
  
  undirect.edge2 <- undirect.edge <- as.edgelist(net)
  colnames(undirect.edge2) <- colnames(undirect.edge) <- c("id1", "id2")
  undirect.edge2[,1] <- undirect.edge[,2]
  undirect.edge2[,2] <- undirect.edge[,1]
  undirect.edge <- undirect.edge3 <- rbind(undirect.edge, undirect.edge2)
  undirect.edge <- undirect.edge[order(undirect.edge[,1],undirect.edge[,2]),]
  
  net.data <- data.frame(id, group, class, deg)
  
  net.data.att <- merge(undirect.edge, net.data, by.x = "id2", by.y = "id", all.x = TRUE)
  net.data.att <- merge(net.data.att, net.data, by.x = "id1", by.y = "id",all.x = TRUE, suffixes = c(".2",".1"))
  #head(net.data.att)
  
  ####################
  chisq.class <- chisq.test(net.data.att$class.1, net.data.att$class.2, correct = correct)
  chisq.group <- chisq.test(net.data.att$group.1, net.data.att$group.2, correct = correct)
  
  chisq.class.group <- chisq.test(net.data$class, net.data$group, correct = correct)
  
  net.RDS <- getRDSsample(net, "group", sample.size = RDS.size, number.of.seeds = number.of.seeds)
  net.RDS <- as.data.frame(net.RDS)
  merge.net.data <- merge(net.RDS, net.data, by.x = "id", by.y = "id", all.x = TRUE)
  
  net.RDS.data <- as.rds.data.frame(merge.net.data, id = "id"
                                    , recruiter.id = "recruiter.id"
                                    , network.size = "degree" )
  chisq.res.net <- chisq.test(net.RDS.data$group, net.RDS.data$class, correct = correct)$p.value
  
  permute.res1 <- SPRTBA(net.RDS.data, "group", "class", chisq.test, per.iter)
  permute.res2 <- SPRTBA(net.RDS.data, "class", "group", chisq.test, per.iter)
  
  permute.res.both <- SPRTBA.both(net.RDS.data, "group", "class", chisq.test, per.iter)
  
  res <- data.frame(mean_deg = mean(deg)
                    , chisq.class= chisq.class$p.value
                    , chisq.group = chisq.group$p.value
                    , chisq.class.group = chisq.class.group$p.value
                    , RDS.chisq = chisq.res.net
                    , RDS.Pclass = mean(merge.net.data$class)
                    , RDS.Pgroup = mean(merge.net.data$group)
                    , RDS1pervar = permute.res1$permute.variable
                    , RDS1per_00 = permute.res1$permute.var.prop[1,1]
                    , RDS1per_11 = permute.res1$permute.var.prop[2,2]
                    , RDS1per_pvalue = permute.res1$p.value.mc
                    , RDS2pervar = permute.res2$permute.variable
                    , RDS2per_00 = permute.res2$permute.var.prop[1,1]
                    , RDS2per_11 = permute.res2$permute.var.prop[2,2]
                    , RDS2per_pvalue = permute.res2$p.value.mc
                    , RDSbothper_pvalue = permute.res.both$p.value.mc)
  return(res)
}


