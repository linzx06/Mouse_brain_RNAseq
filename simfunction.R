getStatesSim <- function(paraMRF, nbr=3, ng=500, nt=5, numsex=2, iter=3, iprob = 0.5)
{
  g <- paraMRF[1]; bsl <- paraMRF[2]; bsi <- paraMRF[3]; 
  bli <- paraMRF[4]; bt <- paraMRF[5]; bsex <- paraMRF[6]
  
  states <- array(sample(2, nbr*ng*nt*numsex, replace=T, prob=c(1-iprob, iprob)) - 1, dim=c(nbr, ng, nt, numsex))
  for (it in 1:iter){
    br <- 1
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
      sex <- 2
      nsex <- 1
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsl*states[2,,t,sex] + bsi*states[3,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[2,,t,sex]) + bsi*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
    }
    br <- 2
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsl*states[1,,t,sex] + bli*states[3,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsl*(1-states[1,,t,sex]) + bli*(1-states[3,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
    }
    br <- 3
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*states[br,,t+1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*states[br,,t-1,sex] + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t-1,sex]) + bsex*(1-states[br,,t,nsex])
      } else {
        t1 <- g + bsi*states[1,,t,sex] + bli*states[2,,t,sex] + bt*(states[br,,t-1,sex]+states[br,,t+1,sex]) + bsex*states[br,,t,nsex]
        t0 <- bsi*(1-states[1,,t,sex]) + bli*(1-states[2,,t,sex]) + bt*(1-states[br,,t-1,sex]+1-states[br,,t+1,sex]) + bsex*(1-states[br,,t,nsex])
      }
      states[br,,t,sex] <- ((exp(t1)/(exp(t1) + exp(t0))) > runif(ng)) + 0
    }   
  }
  return(states)
}

getStatesSimhmm <- function(nbr=3, ng=500, nt=5, numsex=2, nDE, tch, perc){
  states <- array(0, dim=c(nbr, ng, nt, numsex))
  stmp <- matrix(0,nrow=ng, ncol=nt)
  tmp <- sample(ng, round(ng*nDE), replace=F)
  stmp[tmp, 1] <- 1
  for (t in 2:nt){
    e1 <- which(stmp[, t-1]==1)
    e0 <- which(stmp[, t-1]==0)
    
    ne1 <- sum(stmp[, t-1])
    ct <- round(ne1*tch)
    
    e1n <- c(sample(e1, ne1-ct, replace=F), sample(e0, ct, replace=F))
    stmp[e1n, t] <- 1    
  }
  for (b in 1:nbr){
    for (s in 1:numsex){
      states[b,,,s] <- stmp
    }
  }
  l1 <- which(states==1)
  l0 <- which(states==0)
  c10 <- sample(l1, round(length(l1)*perc), replace=F)
  c01 <- sample(l0, round(length(l1)*perc), replace=F)
  states[c(c10,c01)] <- 1 - states[c(c10,c01)]
  return(states)
}

get_z <- function(states,f1_mean,f1_sd)
{
  z <- states*0
  z[states==0] <- rnorm(sum(states==0),mean=0,sd=1)
  z[states==1] <- sample(c(rnorm(floor(sum(states==1)/2),mean=f1_mean,sd=f1_sd), rnorm(ceiling(sum(states==1)/2),mean=-f1_mean,sd=f1_sd)))
  return(z)
}

####simulate the expression values
getsimSeq <- function(b, nbr=3, ng=500, nt=5, numsex=2,  
                      rep = 3, librarySize=NULL, lamda=NULL, shape = 0.85, scale = 0.5,
                      option=c("mrf", "hmm"), paraMRF=NULL, iter=3, iprob = 0.5,
                      nDE=NULL, tch=NULL, perc=NULL)
{
  ###simulate the latent states
  if (option=="mrf"){
    if (is.null(paraMRF)){
      stop(stop("Please enter a valid paraMRF"))
    }
    states <- getStatesSim(paraMRF=paraMRF, nbr=nbr, ng=ng, nt=nt, numsex=numsex, iter=iter, iprob=iprob)
  } else {
    if (is.null(nDE) | is.null(tch) | is.null(perc)){
      stop(stop("Please enter a valid nDE, tch, perc"))
    }
    states <- getStatesSimhmm(nbr=nbr, ng=ng, nt=nt, numsex=numsex, nDE=nDE, tch=tch, perc=perc)
  }
  
  if (is.null(lamda)){
    lamda <- array(1/ng, dim=c(nbr, ng, numsex))
  }
  if (is.null(librarySize)){
    ###median of read count in the brain data around 200
    librarySize <- ng * runif(nbr * (nt + 1) * numsex * rep, min = 100, max = 300)
  }
  lall <- array(dim=c(nbr, ng, nt + 1, numsex))
  lall[,,1,] <- lamda
  ###the threshold for lamda
  lthres <- 1/ng/20
  getl <- function(lamda, states, b, lthres){
    ###don't want lamda to be too small to cause small counts
    tmp <- states[states==1] * (sample(2,length(states[states==1]),replace=T) - 1.5) * 2
    ln <- lamda
    ln[states == 1] <- lamda[states == 1] * b ^ tmp
    ln[ln <= lthres & states == 1] <- lamda[ln <= lthres & states == 1] * b
    return(ln)
  }
  for (t in 1:nt){
    lall[,,t + 1,] <- getl(lall[,,t,], states[,,t,], b, lthres)
  }
  #disp <- array(rgamma(length(lall), shape = 0.85, scale = 0.5), dim = dim(lall))
  ###assume each gene have one dispersion
  disp <- rgamma(ng, shape = 0.85, scale = 0.5)
  data <- c()
  breg <- c()
  tp <- c()
  sex <- c()
  group <- c()
  gid <- 1
  count <- 0
  for (b in 1:nbr){
    for (t in 1:(nt + 1)){
      for (s in 1:numsex){
        for (r in 1:rep){
          count <- count + 1
          lib <- librarySize[count]
          data <- cbind(data, matrix(rnbinom(n = ng, mu = lall[b,,t,s] * lib, 
                                             size = 1/disp),byrow = F, ncol = 1))
          group <- c(group, gid)
        } 
        breg <- c(breg, rep(b, rep))
        tp <- c(tp, rep(t, rep))
        sex <- c(sex, rep(s, rep))
        if (rep!=1){
          gid <- gid + 1
        }
      }
      if (rep==1){
        gid <- gid + 1
      }
    }
  }
  #print(count)
  return(list(breg = breg, tp = tp, sex = sex, data = data, states = states, lamdaall = lall, 
              librarySize = librarySize, dispersion = disp, group = group))
}

getspesen <- function(fdr, states, alpha=0.05){
  tmp <- sort(as.numeric(fdr))
  if (tmp[1]>=alpha){
    sen <- 0; spe <- 1; fdr <- 0; cutoff <- alpha  
  } else {
    tot <- length(tmp)
    lab <- tot
    low <- 1; up <- tot
    FDR <- 0.2
    nsat <- TRUE
    while(nsat){
      if (FDR-alpha > 0){
        up <- lab
        labpre <- lab
        lab <- floor((low + lab)/2)
        if (labpre==lab){
          nsat <- FALSE
        }
        FDR <- sum(tmp[1:lab])/lab
      } else if (FDR-alpha<= -0.001) {
        low <- lab
        labpre <- lab
        lab <- floor((up + lab)/2)
        if (labpre==lab){
          nsat <- FALSE
        }
        FDR <- sum(tmp[1:lab])/lab
      } else {
        nsat <- FALSE
      }
      FDRpre <- FDR
    }
    cutoff <- tmp[lab]
    ss <- (fdr<=cutoff)+0
    spe<- sum((ss+states)==0)/sum(states==0)
    sen <- sum((ss+states)==2)/sum(states==1)
    fdr <- sum((2*ss+states)==2)/sum(ss==1)
  }
  return(list(spe=spe, sen=sen, fdr=fdr, cutoff=cutoff))
}

getspesenBH <- function(fdr, states, alpha=0.05){
  tmp <- sort(as.numeric(fdr))
  lab <- which(((tmp<=((1:length(tmp))/length(tmp)*alpha))+0)==0)[1] - 1
  cutoff <- tmp[lab]
  ss <- (fdr<=cutoff)+0
  spe<- sum((ss+states)==0)/sum(states==0)
  sen <- sum((ss+states)==2)/sum(states==1)
  fdr <- sum((2*ss+states)==2)/sum(ss==1)
  return(list(spe=spe, sen=sen, fdr=fdr, cutoff=cutoff))
}

getspesenA <- function(fdr, states, breaks=200){
  tmp <- sort(as.numeric(fdr))
  totp <- sum(states)*2
  labs <- ceiling(seq(0,totp,totp/breaks))
  labs[1] <- 1
  p <- c(); tp <- c()
  for (cutoff in tmp[labs]){
    ss <- (fdr<=cutoff)+0
    p <- c(p, sum(ss))
    tp <- c(tp, sum((ss+states)==2))
  }
  return(list(p=p, tp=tp))
}
