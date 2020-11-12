if(1==2){
  nnodes <- network.size(net)
  nsamp0 <- number.of.seeds
  fixinitial=fixinitial
  nsamp <- sample.size
  replace <- sample.with.replacement
  coupons <- number.of.coupons
  select=NULL
  disall=TRUE
  bias=NULL
  seed.distribution=seed.distribution
  trait.variable=trait.variable
  recruitment.rates=recruitment.rates
  recruitment.rates.by=recruitment.rates.by
  rds.samp=NULL
  }

### working stuff

rdssample.krista<-function(net,nnodes=network.size(net), nsamp0,fixinitial,nsamp, replace=FALSE,
                    coupons,select=NULL,bias=NULL,rds.samp=NULL,seed.distribution=NULL, 	
                    disall=FALSE,trait.variable="disease", recruitment.rates=NULL, 
                    recruitment.rates.by="disease"){ 
	
net<-as.network.uncompressed.rds(net)

#smtx<-as.sociomatrix(net)
disease <- network::get.vertex.attribute(net,trait.variable)
if(is.null(disease) || all(is.na(disease))){
  stop(sprintf("No variable called %s appears in the data.", trait.variable))
}

# Note:  recruitment.rates is a matrix (later maybe an array).  number of columns should be 1 plus max coupons, corresponding to the frequencies of returing 0, 1, ... coupons coupons.  number of rows should be the number of classes of recruitment.rates.by.  Row labels should correspond to levels of recruitment.rates.by                    
if(!is.null(recruitment.rates)){
	rec.var<-network::get.vertex.attribute(net,recruitment.rates.by)
	if(is.null(rec.var) || all(is.na(rec.var))){
  		stop(sprintf("No variable called %s appears in the data.", recruitment.rates.by))
	}
	if(any(sort(unique(rec.var)) != sort(dimnames(recruitment.rates)[[1]]))){
  		stop(sprintf("Levels of %s not the same as row names of recruitment.rates.", recruitment.rates.by))
	}
	if(dim(recruitment.rates)[2] != coupons+1){
  		stop(sprintf("recruitment.rates must have columns equal to the number of coupons plus 1, corresponding to 0 through coupons recruits"))
  	}  		
}


degs <- sapply(net$iel,length)+sapply(net$oel,length)
# Choose Initial Sample
#set.seed(seeed)
nsample<-NULL # vector of indices of sampled nodes
wsample<-NULL # vector of waves of each sampled node
degsample<-NULL #vector of degrees of sampled nodes
dissample<-NULL # vector of disease of sampled nodes
todis<-rep(0,nsamp) # vector of numbers of referrals to diseased nodes
tonodis<-rep(0,nsamp) # vector of number of referrans to undiseased
nominators<-NULL # recruiter of each sample

if(is.null(seed.distribution)){
#ps<-apply(smtx,1,sum)
 ps<-sapply(net$iel,length)+sapply(net$oel,length)
}else{
 ps<-seed.distribution
}
if(any(ps<0)){
 stop("Negative probabilities of seed selection are not allowed.")
}
if(all(ps==0)){
 stop("At least one node should have a positive probability of seed selection.\n Are you sure you are selecting some nodes with positive degree?.")
}


if(fixinitial=="FALSE"){fixinitial<--1} #line added 0527
if(!(fixinitial %in% c(1,2))){
# If fixinitial is not 1 or 2 choose seeds proportional to seed.distribution
  if(sum(ps>0) < nsamp0){
    nsample<-sample.int(nnodes,size=nsamp0,prob=ps+0.000000001)
  }else{
    nsample<-sample.int(nnodes,size=nsamp0,prob=ps)
  }
}else{
# If fixinitial is 1 or 2 choose seeds from those with disease==1 (and proportional to seed.distribution)
  if(sum(ps[which(disease==1)]) < nsamp0){
    nsample<-sample(which(disease==1),size=nsamp0,prob=ps[which(disease==1)]+0.000000001)
  }else{
    nsample<-sample(which(disease==1),size=nsamp0,prob=ps[which(disease==1)])
  }
}

		
wsample<-rep(0,nsamp0)
nominators<-rep(0,nsamp0) #who nominated each sample

# Subsequent Sampling
# Note:  for now, I am not keeping track of when limited numbers of unsampled alters reduces the sampled counts below desired.

refnode<-1
while(length(nsample)<nsamp){
        aaa <- .Call("getNeighborhood_R", net, nsample[refnode], "combined", TRUE, PACKAGE = "network")
   	if(!replace){  # choose group from which to choose next one
   	  aaa<-setdiff(aaa,nsample)
        }
	# determine number to sample
	if(!is.null(bias)){
		prob<-disease[aaa]*bias
		prob[which(prob==0)]<-1}else{prob<-rep(1,length(aaa))}
	if(!is.null(recruitment.rates)){
		target_num_recruits<-sample(c(0:coupons), size=1, 
			prob=recruitment.rates[disease[nsample[refnode]]+1,]) 
		}else{target_num_recruits<-coupons}
	if((nsamp-length(nsample))>=target_num_recruits){ntosample<-target_num_recruits
		}else{ntosample<-nsamp-length(nsample)}
	if(ntosample>length(aaa)){ntosample<-length(aaa)}### backhere
	if(length(aaa)>1){thissamp<-sample(aaa,size=ntosample,prob=prob)}else{thissamp<-NULL}
	if(length(aaa)==1){thissamp<-aaa}
	if(ntosample>0){
	 nsample<-c(nsample,thissamp)
	 wsample<-c(wsample,rep((wsample[refnode]+1),ntosample))
	 nominators<-c(nominators,rep(nsample[refnode],ntosample))
	 for(nn in 1:ntosample){
		if(disease[thissamp[nn]]==1){todis[refnode]<-todis[refnode]+1}
		if(disease[thissamp[nn]]==0){tonodis[refnode]<-tonodis[refnode]+1}
		}
    }
   	refnode<-refnode+1
   	if(refnode>length(nsample)){#here, trapped -need new seed!
		cat(sprintf("trapped. fixinitial=%d, sampled so far=%d.\n",fixinitial,length(nsample)))
   		if(fixinitial==-1){
# If fixinitial is -1 choose a reseed proportional to degree
		  www<-setdiff(c(1:nnodes),nsample)
		  nsample<-c(nsample,sample(www,size=1,prob=degs[www]))	    }else{
# If fixinitial is not -1 choose a reseed from those with disease==fixinitial (and proportional to degree)
		  zzz<-setdiff(which(disease==fixinitial),nsample)
		  if(length(zzz)>0){
		    nsample<-c(nsample,sample(zzz,size=1,prob=degs[zzz]))
		  }else{
# If fixinitial is not -1 and there are none with disease==fixinitial, choose a reseed proportional to degree
# This is the typical case!!!
# Note that it is not proportional to the seed.distribution but the nodal degree
		    cat(sprintf("Re-seeding proportional to degree.\n"))
		    www<-setdiff(c(1:nnodes),nsample)
		    nsample<-c(nsample,sample(www,size=1,prob=degs[www]))
		    cat(sprintf("Re-seed %d has degree %d.\n",nsample[length(nsample)],
			degs[nsample[length(nsample)]]))
			}
		}
wsample<-c(wsample,0)
nominators<-c(nominators,0)
   		}	
   		}	

degsample<-degs[nsample]
dissample<-disease[nsample]

    if(disall){
     temp <- matrix(0,ncol=2,nrow=nnodes)
     for (i in 1:nnodes)  {
       temp[i,] <- tabulate(
          as.numeric(disease[.Call("getNeighborhood_R", net, i, "combined", TRUE,
                        PACKAGE = "network")])+1,nbin=2)
     }
    } 

    if(disall){
     list(nsample=nsample,wsample=wsample, nominators=nominators, degsample=degsample,dissample=dissample,todis=todis,tonodis=tonodis,
       todisall=temp[nsample,2],tonodisall=temp[nsample,1])
    }else{
     list(nsample=nsample,wsample=wsample, nominators=nominators, degsample=degsample,dissample=dissample,todis=todis,tonodis=tonodis)
    }
}
