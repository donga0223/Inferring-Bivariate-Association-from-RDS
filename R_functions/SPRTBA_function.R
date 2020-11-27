if(2==3){
  net <- dat
  fix.variable <- "char1"
  permute.variable <- "char2"
  test <- chisq.test
  n <- 10
  permute.level <- c("0","1")
}


if(2==3){
  net <- net.RDS.data
  fix.variable <- "group"
  permute.variable <- "class"
  test <- chisq.test
  n <- 10
  permute.level <- c("A", "B")
}

################################
## This is the SPRTBA function that permutes one variable and fixes the other.
## net : RDS data
## fix.variable : variable name for not permuting
## permute.variable : variable name for permuting
## n : number of permutations
################################

SPRTBA <- function(net, fix.variable, permute.variable, test = chisq.test, n, correct = TRUE)
{
  if(is.rds.data.frame(net)==FALSE){
    stop(sprintf("net must be an rds.data.frame"))
  }
  net[] <- lapply(net, factor)
  
  ID <- get.id(net)
  rec.id <- get.rid(net)
  degree <- get.net.size(net)
  fix.var <- net[[fix.variable]]
  permute.var <- net[[permute.variable]]
  wave <- get.wave(net)
  seed.id <- get.seed.id(net)
  
  rec.fix.variable <- paste("rec",fix.variable,sep=".")
  rec.permute.variable <- paste("rec",permute.variable,sep=".")
  new.permute.variable <- paste("new",permute.variable,sep=".")
  
  ## merge data : including recruitment information
  a1 <- data.frame(ID, rec.id, degree, fix.var, permute.var, wave, seed.id)
  colnames(a1) <- c("ID", "rec.id", "degree", fix.variable, permute.variable, "wave", "seed.id")
  a2 <- data.frame(ID, fix.var, permute.var)
  colnames(a2) <- c("rec.id", rec.fix.variable, rec.permute.variable)
  
  net.with.rec <- merge(a1, a2, by = "rec.id", all.x = TRUE)
  order.net.with.rec <- net.with.rec[order(net.with.rec[["seed.id"]], net.with.rec[["wave"]]),]
  
  ## 
  permute.var.table <- table(order.net.with.rec[[rec.permute.variable]], order.net.with.rec[[permute.variable]])
  permute.var.prop <- prop.table(permute.var.table, 1)
  permute.var.prop[which(is.na(permute.var.prop)==TRUE)] <- 1/length(which(is.na(permute.var.prop)==TRUE))
  
  ##
  fix.var.table <- table(order.net.with.rec[[rec.fix.variable]], order.net.with.rec[[fix.variable]])
  fix.var.prop <- prop.table(fix.var.table, 1)
  fix.var.prop[which(is.na(fix.var.prop)==TRUE)] <- 1/length(which(is.na(fix.var.prop)==TRUE))
  

  ## permutation test
  obs.test.res <- test(order.net.with.rec[[permute.variable]], order.net.with.rec[[fix.variable]], correct = correct)
  
  order.net.with.rec[[new.permute.variable]] <- rep(NA, dim(order.net.with.rec)[1])
  #per.table.newper.fix <- matrix(NA, nrow = n, ncol = )
  per.test.res <- rep(NA, n)
  
  for(j in 1:n){
    for(i in 1:dim(net)[1]){
      if(is.na(order.net.with.rec[["rec.id"]][i])){
        order.net.with.rec[[new.permute.variable]][i] <- NA
      }else if(order.net.with.rec[["rec.id"]][i] == 0){
        ccc <- order.net.with.rec[[permute.variable]][i]
        order.net.with.rec[[new.permute.variable]][i] <- levels(order.net.with.rec[[permute.variable]])[ccc]
        rm(ccc)
      }else if(order.net.with.rec[["rec.id"]][i] != 0){
        aaa <- order.net.with.rec[[new.permute.variable]][match(order.net.with.rec[["rec.id"]][i], order.net.with.rec[["ID"]][1:i])]
        bbb <- rmultinom(1, 1, prob = permute.var.prop[aaa,])==1      
        order.net.with.rec[[new.permute.variable]][i] <- rownames(bbb)[which(bbb)]
        rm(aaa, bbb)
      }
    }
    if(length(table(order.net.with.rec[[new.permute.variable]])) == 1){
      per.test.res[j] <- NA
    }else{
      per.test.res[j] <- test(order.net.with.rec[[new.permute.variable]], order.net.with.rec[[fix.variable]], correct = correct)$statistic
    }
  }
    
  test.res <- na.omit(per.test.res)
  b <- sum(ifelse(test.res >= obs.test.res$statistic, 1, 0))
  p.value <- b/length(test.res)
  p.value.mc <- (b+1)/(length(test.res)+1)
  return(list(fix.variable = fix.variable
              ,permute.variable = permute.variable
              ,n.iteration = n
              ,per.test.res = per.test.res
              ,observed.test = obs.test.res
              #,p.value = p.value
              ,p.value.mc = p.value.mc
              ,permute.var.prop = permute.var.prop
              ,fix.var.prop = fix.var.prop
              ,permute.var.table = permute.var.table
              ,fix.var.table = fix.var.table))

}


################################
## This is the SPRTBA function that permutes two variables together.
## net : RDS data
## fix.variable : variable 1 for bivariate association test
## permute.variable : variable 2 for bivariate association test
## In this test, we permute both variables, so you don't need to consider which one is fixed and which one is permuted
## n : number of permutations
################################

SPRTBA.both <- function(net, fix.variable, permute.variable, test = chisq.test, n, correct = TRUE)
{
  if(is.rds.data.frame(net)==FALSE){
    stop(sprintf("net must be an rds.data.frame"))
  }
  net[] <- lapply(net, factor)
  
  ID <- get.id(net)
  rec.id <- get.rid(net)
  degree <- get.net.size(net)
  fix.var <- net[[fix.variable]]
  permute.var <- net[[permute.variable]]
  wave <- get.wave(net)
  seed.id <- get.seed.id(net)
  
  rec.fix.variable <- paste("rec",fix.variable,sep=".")
  rec.permute.variable <- paste("rec",permute.variable,sep=".")
  new.permute.variable <- paste("new",permute.variable,sep=".")
  new.fix.variable <- paste("new",fix.variable,sep=".")
  
  ## merge data : including recruitment information
  a1 <- data.frame(ID, rec.id, degree, fix.var, permute.var, wave, seed.id)
  colnames(a1) <- c("ID", "rec.id", "degree", fix.variable, permute.variable, "wave", "seed.id")
  a2 <- data.frame(ID, fix.var, permute.var)
  colnames(a2) <- c("rec.id", rec.fix.variable, rec.permute.variable)
  
  net.with.rec <- merge(a1, a2, by = "rec.id", all.x = TRUE)
  order.net.with.rec <- net.with.rec[order(net.with.rec[["seed.id"]], net.with.rec[["wave"]]),]
  
  ## 
  permute.var.table <- table(order.net.with.rec[[rec.permute.variable]], order.net.with.rec[[permute.variable]])
  permute.var.prop <- prop.table(permute.var.table, 1)
  permute.var.prop[which(is.na(permute.var.prop)==TRUE)] <- 1/length(which(is.na(permute.var.prop)==TRUE))
  
  ##
  fix.var.table <- table(order.net.with.rec[[rec.fix.variable]], order.net.with.rec[[fix.variable]])
  fix.var.prop <- prop.table(fix.var.table, 1)
  fix.var.prop[which(is.na(fix.var.prop)==TRUE)] <- 1/length(which(is.na(fix.var.prop)==TRUE))
  
  
  ## permutation test
  obs.test.res <- test(order.net.with.rec[[permute.variable]], order.net.with.rec[[fix.variable]], correct = correct)
  
  order.net.with.rec[[new.permute.variable]] <- rep(NA, dim(order.net.with.rec)[1])
  order.net.with.rec[[new.fix.variable]] <- rep(NA, dim(order.net.with.rec)[1])
  
  per.test.res <- rep(NA, n)
  
  for(j in 1:n){
    for(i in 1:dim(net)[1]){
      if(is.na(order.net.with.rec[["rec.id"]][i])){
        order.net.with.rec[[new.permute.variable]][i] <- NA
      }else if(order.net.with.rec[["rec.id"]][i] == 0){
        ccc <- order.net.with.rec[[permute.variable]][i]
        order.net.with.rec[[new.permute.variable]][i] <- levels(order.net.with.rec[[permute.variable]])[ccc]
        rm(ccc)
      }else if(order.net.with.rec[["rec.id"]][i] != 0){
        aaa <- order.net.with.rec[[new.permute.variable]][match(order.net.with.rec[["rec.id"]][i], order.net.with.rec[["ID"]][1:i])]
        bbb <- rmultinom(1, 1, prob = permute.var.prop[aaa,])==1      
        order.net.with.rec[[new.permute.variable]][i] <- rownames(bbb)[which(bbb)]
        rm(aaa, bbb)
      }
    }
    
    for(i in 1:dim(net)[1]){
      if(is.na(order.net.with.rec[["rec.id"]][i])){
        order.net.with.rec[[new.fix.variable]][i] <- NA
      }else if(order.net.with.rec[["rec.id"]][i] == 0){
        ccc <- order.net.with.rec[[fix.variable]][i]
        order.net.with.rec[[new.fix.variable]][i] <- levels(order.net.with.rec[[fix.variable]])[ccc]
        rm(ccc)
      }else if(order.net.with.rec[["rec.id"]][i] != 0){
        aaa <- order.net.with.rec[[new.fix.variable]][match(order.net.with.rec[["rec.id"]][i], order.net.with.rec[["ID"]][1:i])]
        bbb <- rmultinom(1, 1, prob = fix.var.prop[aaa,])==1      
        order.net.with.rec[[new.fix.variable]][i] <- rownames(bbb)[which(bbb)]
        rm(aaa, bbb)
      }
    }
    if(length(table(order.net.with.rec[[new.permute.variable]])) == 1){
      per.test.res[j] <- NA
    }else{
      per.test.res[j] <- test(order.net.with.rec[[new.permute.variable]], order.net.with.rec[[fix.variable]], correct = correct)$statistic
    }
    
  }
  
  test.res <- na.omit(per.test.res)
  b <- sum(ifelse(test.res >= obs.test.res$statistic, 1, 0))
  p.value <- b/length(test.res)
  p.value.mc <- (b+1)/(length(test.res)+1)
  
  return(list(fix.variable = fix.variable
              ,permute.variable = permute.variable
              ,iteration = n
              ,per.test.res = per.test.res
              ,observed.test = obs.test.res
              #,p.value = p.value
              ,p.value.mc = p.value.mc
              ,permute.var.prop = permute.var.prop
              ,fix.var.prop = fix.var.prop
              ,permute.var.table = permute.var.table
              ,fix.var.table = fix.var.table))
}


