
##########################################
## iter : number of repeat the whole process
## number.of.seeds : this is RDS data, so we need the number of seeds size for sampling
## sample.size : desired sample size
## AA, BB : probabilites of transition from X1 vales "A" to "A" or B" to "B"
## x00, x11 : probabilites of transition from X2 vales "0" to "0" or 1" to "1"

## per.iter : number of iteration for the permutation test
##########################################

Generate.data <- function(iter, number.of.seeds, sample.size, AA, BB, x00, x11, correct = TRUE, per.iter = 1000){

chisq.res <- rep(NA, iter)
permute.res1 <- rep(NA, iter)
permute.res2 <- rep(NA, iter)
permute.res1.mc <- rep(NA, iter)
permute.res2.mc <- rep(NA, iter)

for(j in 1:iter){
  ## for plot set.seed(18374)
  set.seed(j)
  
  id <- 1:sample.size
  
  rec.id <- rep(0, number.of.seeds)
  c <- length(rec.id)
  initial <- 1
  finish <- number.of.seeds
  while(length(rec.id) <= sample.size){
    a <- sample(1:3, c, replace = T)
    b <- sample(rep(initial:finish, a))
    c <- length(b)
    rec.id <- c(rec.id, b)
    initial <- max(b)+1
    finish <- max(b)+c
  }
  
  rec.id <- rec.id[1:sample.size]
  degree <- sample(1:30, sample.size, replace=T)
  dat <- data.frame(id, rec.id, degree)
  
  
  ## cori.03.1: c(0.91, 0.09, 0.44, 0.56)
  ## cori.15.1: c(0.85, 0.15, 0.18, 0.82)
  myprob.char1 <- matrix(c(AA, 1-AA, 1-BB, BB),nrow=2)
  
  char1 <- rep(NA,sample.size)
  for(i in 1:sample.size){
    if(rec.id[i]==0){char1[i] <- sample(1:2,1)}
    else{
      a <- char1[match(dat$rec.id[i],dat$id[1:i])]
      b <- rmultinom(1, 1, prob = myprob.char1[,a])==1
      char1[i] <- which(b)
      rm(a,b)
    }
  }
  
  char1 <- as.factor(char1)
  levels(char1) <- c("A", "B")
  
  
  ## cori.03.1: c(0.48, 0.52, .35, .65)
  ## cori.15.1: c(0.64, 0.36, .54, .46)
  myprob.char2 <- matrix(c(x00, 1-x00, 1-x11, x11), nrow=2)
  
  char2 <- rep(NA,sample.size)
  for(i in 1:sample.size){
    if(rec.id[i]==0){char2[i] <- sample(1:2,1)}
    else{
      a <- char2[match(dat$rec.id[i],dat$id[1:i])]
      b <- rmultinom(1, 1, prob = myprob.char2[,a])==1
      char2[i] <- which(b)
      rm(a,b)
    }
  }
  
  char2 <- as.factor(char2)
  levels(char2) <- c("0", "1")
  
  
  
  dat <- cbind(dat, char1, char2)
  dat <-as.rds.data.frame(df=dat,id="id",
                          recruiter.id="rec.id",
                          network.size="degree",
                          population.size=15000)
  
  chisq.res[j] <- chisq.test(char1, char2, correct = correct)$p.value
  
  permute.test1 <- permutation.test.II(dat, "char1", "char2", chisq.test, per.iter, c("0", "1"))
  permute.res1[j] <- permute.test1$p.value
  permute.res1.mc[j] <- permute.test1$p.value.mc
   
  permute.test2 <- permutation.test.II(dat, "char2", "char1",  chisq.test, per.iter, c("A", "B"))
  permute.res2[j] <- permute.test2$p.value
  permute.res2.mc[j] <- permute.test2$p.value.mc
  
  print(j)
  rm(permute.test1, permute.test2)
  
}

res <- data.frame(chisq.res, permute.res1, permute.res2, permute.res1.mc, permute.res2.mc)
return(res = res)

}

