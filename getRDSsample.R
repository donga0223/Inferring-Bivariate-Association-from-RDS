if(1==2){
  net <- dixon
  trait.variable <- "sex"
  sample.size <- 200
  seed.selection="degree"
  number.of.seeds=10
  number.of.coupons=2
  max.degree=NULL
  use.C.sampler=FALSE
  recruitment.rates=NULL
  recruitment.rates.by="disease"
  sample.with.replacement=FALSE
  verbose=TRUE
}



getRDSsample <- function(net,trait.variable="disease",sample.size=500, 
 seed.selection="degree", number.of.seeds=10, number.of.coupons=2,
 max.degree=NULL, use.C.sampler=FALSE, recruitment.rates=NULL,
recruitment.rates.by="disease",sample.with.replacement=FALSE,
 verbose=TRUE){
#
 N <- network.size(net)
 disease <- as.numeric(net %v% trait.variable)
 if(is.null(disease) || all(is.na(disease))){
  stop(sprintf("No variable called %s appears in the data.", trait.variable))
 }
#
 deg<-sapply(net$iel,length)+sapply(net$oel,length)
 if(is.null(max.degree)){max.degree <- max(deg, na.rm=TRUE)}
 deg[deg>max.degree] <- max.degree
 if(length(seed.selection)!=network.size(net)){
  seed.distribution <- switch(as.character(seed.selection[1]),
   "random"={
       rep(1,length=network.size(net))
     },
   "degree"={
      deg
     },
   "allwithtrait"={
      disease
     },
   "allwithtraitdegree"={
      disease*deg
     },
      deg)
  fixinitial <- switch(seed.selection,
     "random"=-3, "allwithtrait"=1, "degree"=3,
     "allwithtraitdegree"=2, 1)
  }else{
    seed.distribution <- seed.selection
    fixinitial <- 3
  }
     
# dummy<- sort(10*rep((0:max.degree),2)+ rep(c(0,1) ,each= (max.degree + 1)))
# seed.indices<-1:number.of.seeds
# indices<-toindexref.RDS(deg,disease,dummy)
# u <- sort(unique(seed.indices))
# for(i in u){
#  seed.distribution[indices==i] <- 1000
# }

# si <- rep(0,length=number.of.seeds)
# for(i in u){
#  if(any(indices==i)){
#   if(sum(seed.indices==i) <= sum(indices==i)){
#    si[seed.indices==i] <- sample(seq(along=indices)[indices==i],size=sum(seed.indices==i),replace=FALSE)
#   }else{
#    si[seed.indices==i] <- sample(seq(along=indices)[abs(indices-i)<2],size=sum(seed.indices==i),replace=FALSE)
#   }
#  }
# }
# sd <- seed.distribution
# sd[si] <- 100000
  s<-rdssample.krista(net,nsamp0=number.of.seeds, fixinitial=fixinitial,
        nsamp=sample.size, replace=sample.with.replacement, coupons=number.of.coupons,
        select=NULL, disall=TRUE,
        bias=NULL,seed.distribution=seed.distribution, trait.variable=trait.variable, 
        recruitment.rates=recruitment.rates, recruitment.rates.by=recruitment.rates.by)
  s$degsample[s$degsample>max.degree] <- max.degree
  s=as.rds.data.frame(
    data.frame(id=s$nsample,recruiter.id=s$nominators,degree=s$degsample,disease=s$dissample,todisall=s$todisall,tonodisall=s$tonodisall),
    population.size=N, network.size="degree", max.coupons=number.of.coupons,
    check.valid=!sample.with.replacement)
  s
}
