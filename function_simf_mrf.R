calmf1I  <- function(meanf, fy1, fy0, paraMRF, option = c("mf", "simf"), iter=100, br=50)
{
  g <- paraMRF[1]; bsl <- paraMRF[2]; bsi <- paraMRF[3]; bli <- paraMRF[4]; bt <- paraMRF[5]; bsex <- paraMRF[6]
  tmp <- dim(fy1)
  nbr <- tmp[1]; ng <- tmp[2]; nt <- tmp[3]; nsex <- tmp[4];
  mfsum <- meanf*0
  for (it in 1:iter){
    br <- 1
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
      sex <- 2
      nsex <- 1
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
    }  
    br <- 2
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
    }
    br <- 3
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
      } else {
        t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
        t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      if (option == "mf"){
        meanf[br,,t,sex] <- prob
      } else if (option == "simf"){
        meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      } else {
        print("please enter a valid option: mf or simf")
      }
    }
    if (it >= br){
      mfsum <- mfsum + meanf
    } 
  }
  if (option == "mf"){
    return(meanf)
  } else {
    return(mfsum/(iter-br+1))
  }
}

calw <- function(fy1, fy0, paraMRF, meanf){
  g <- paraMRF[1]; bsl <- paraMRF[2]; bsi <- paraMRF[3]; bli <- paraMRF[4]; bt <- paraMRF[5]; bsex <- paraMRF[6]
  tmp <- dim(fy1)
  nbr <- tmp[1]; ng <- tmp[2]; nt <- tmp[3]; nsex <- tmp[4];
  w1 <- meanf * 0
  
  br <- 1
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      ##the term when z=1
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      ##the term when z=0
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
    sex <- 2
    nsex <- 1
    if (t==1){
      ##the term when z=1
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      ##the term when z=0
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
  }  
  br <- 2
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
    sex <- 2
    nsex <- 1
    if (t==1){
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
  }
  br <- 3
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
    sex <- 2
    nsex <- 1
    if (t==1){
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    w1[br,,t,sex] <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex] + exp(t0)*fy0[br,,t,sex])  
  }
  return(w1)
}


optimNR <- function(paraIni = rep(0, 6), meanf, w1, alpha = 10^(-6), maxiter = 500)
{
  tmp <- dim(fy1)
  nbr <- tmp[1]; ng <- tmp[2]; nt <- tmp[3]; nsex <- tmp[4];
  w0 <- 1 - w1
  iter <- 1
  para <- rep(1, 6)
  parapre <- paraIni
  while (max(abs(para - parapre)) >= alpha & iter <= maxiter){
    #print(para)
    if (iter != 1){
      parapre <- para
    }
    g <- parapre[1]; bsl <- parapre[2]; bsi <- parapre[3]; bli <- parapre[4]; bt <- parapre[5]; bsex <- parapre[6]
    hes <- matrix(0, nrow = 6, ncol = 6)
    del <- rep(0, 6)  
    br <- 1
    for (t in 1:nt){
      for (sex in 1:2){
        nsex <- 3 - sex
        if (t==1){
          ##the term when z=1
          t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
          ##the term when z=0
          t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t+1,sex] 
        } else if(t==nt) {
          t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
          t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] 
        } else {
          t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
          t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] + meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] + 1 - meanf[br,,t+1,sex] 
        }
        denom <- exp(t1) + exp(t0)
        zl1 <- meanf[2,,t,sex]
        zl0 <- 1 - meanf[2,,t,sex]
        zi1 <- meanf[3,,t,sex]
        zi0 <- 1 - meanf[3,,t,sex]
        zsex1 <- meanf[br,,t,nsex]
        zsex0 <- 1 - meanf[br,,t,nsex]
        ###g
        numerg <- exp(t1)
        del[1] <- del[1] + sum(w1[br,,t,sex]) - sum(numerg/denom)
        ###bsl
        numersl <- zl1*exp(t1) + zl0*exp(t0)
        del[2] <- del[2] + sum(w1[br,,t,sex]*zl1) + sum(w0[br,,t,sex]*zl0) - sum(numersl/denom)
        ###bsi
        numersi <- zi1*exp(t1) + zi0*exp(t0)
        del[3] <- del[3] + sum(w1[br,,t,sex]*zi1) + sum(w0[br,,t,sex]*zi0) - sum(numersi/denom)
        ###bli
        ###bt
        numert <- zt1*exp(t1) + zt0*exp(t0)
        del[5] <- del[5] + sum(w1[br,,t,sex]*zt1) + sum(w0[br,,t,sex]*zt0) - sum(numert/denom)
        ###bsex
        numersex <- zsex1*exp(t1) + zsex0*exp(t0)
        del[6] <- del[6] + sum(w1[br,,t,sex]*zsex1) + sum(w0[br,,t,sex]*zsex0) - sum(numersex/denom)
        #########diagonal of the hessian matrix########
        ###g, g
        numer <- exp(t1) * exp(t0)
        hes[1, 1] <- hes[1, 1] + sum(numer/denom^2)
        ###bsl, bsl
        term1 <- ( zl1^2*exp(t1) + zl0^2*exp(t0) ) / denom
        term2 <- - numersl^2 / denom^2
        hes[2, 2] <- hes[2, 2] + sum(term1 + term2)
        ###bsi, bsi
        term1 <- ( zi1^2*exp(t1) + zi0^2*exp(t0) ) / denom
        term2 <- - numersi^2 / denom^2
        hes[3, 3] <- hes[3, 3] + sum(term1 + term2)
        ###bli, bli
        ###bt, bt
        term1 <- ( zt1^2*exp(t1) + zt0^2*exp(t0) ) / denom
        term2 <- - numert^2 / denom^2
        hes[5, 5] <- hes[5, 5] + sum(term1 + term2)
        ###bsex, bsex
        term1 <- ( zsex1^2*exp(t1) + zsex0^2*exp(t0) ) / denom
        term2 <- - numersex^2 / denom^2
        hes[6, 6] <- hes[6, 6] + sum(term1 + term2)
        #########off-diagonal of the hessian matrix########
        ###g, bsl
        term1 <- zl1*exp(t1) / denom
        term2 <- - numersl * numerg / denom^2
        hes[1, 2] <- hes[1, 2] + sum(term1 + term2)
        ###g, bsi
        term1 <- zi1*exp(t1) / denom
        term2 <- - numersi * numerg / denom^2
        hes[1, 3] <- hes[1, 3] + sum(term1 + term2)
        ###g, bli
        ###g, bt
        term1 <- zt1*exp(t1) / denom
        term2 <- - numert * numerg / denom^2
        hes[1, 5] <- hes[1, 5] + sum(term1 + term2)
        ###g, bsex
        term1 <- zsex1*exp(t1) / denom
        term2 <- - numersex * numerg / denom^2
        hes[1, 6] <- hes[1, 6] + sum(term1 + term2)  
        ###bsl, bsi
        term1 <- ( zl1*zi1*exp(t1) + zl0*zi0*exp(t0) ) / denom
        term2 <- - numersl * numersi / denom^2
        hes[2, 3] <- hes[2, 3] + sum(term1 + term2)
        ###bsl, bli
        ###bsl, bt
        term1 <- ( zl1*zt1*exp(t1) + zl0*zt0*exp(t0) ) / denom
        term2 <- - numersl * numert / denom^2
        hes[2, 5] <- hes[2, 5] + sum(term1 + term2)
        ###bsl, bsex
        term1 <- ( zl1*zsex1*exp(t1) + zl0*zsex0*exp(t0) ) / denom
        term2 <- - numersl * numersex / denom^2
        hes[2, 6] <- hes[2, 6] + sum(term1 + term2)
        ###bsi, bli
        ###bsi, bt
        term1 <- ( zi1*zt1*exp(t1) + zi0*zt0*exp(t0) ) / denom
        term2 <- - numersi * numert / denom^2
        hes[3, 5] <- hes[3, 5] + sum(term1 + term2)
        ###bsi, bsex
        term1 <- ( zi1*zsex1*exp(t1) + zi0*zsex0*exp(t0) ) / denom
        term2 <- - numersi * numersex / denom^2
        hes[3, 6] <- hes[3, 6] + sum(term1 + term2)
        ###bli, bt
        ###bli, bsex
        ###bt, bsex
        term1 <- ( zt1*zsex1*exp(t1) + zt0*zsex0*exp(t0) ) / denom
        term2 <- - numert * numersex / denom^2
        hes[5, 6] <- hes[5, 6] + sum(term1 + term2)       
      }
    }
    br <- 2
    for (t in 1:nt){
      for (sex in 1:2){
        nsex <- 3 - sex
        if (t==1){
          t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
          t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t+1,sex] 
        } else if(t==nt) {
          t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
          t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] 
        } else {
          t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
          t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] + meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] + 1 - meanf[br,,t+1,sex] 
        }
        denom <- exp(t1) + exp(t0)
        zs1 <- meanf[1,,t,sex]
        zs0 <- 1 - meanf[1,,t,sex]
        zi1 <- meanf[3,,t,sex]
        zi0 <- 1 - meanf[3,,t,sex]
        zsex1 <- meanf[br,,t,nsex]
        zsex0 <- 1 - meanf[br,,t,nsex]
        ###g
        numerg <- exp(t1)
        del[1] <- del[1] + sum(w1[br,,t,sex]) - sum(numerg/denom)
        ###bsl
        numersl <- zs1*exp(t1) + zs0*exp(t0)
        del[2] <- del[2] + sum(w1[br,,t,sex]*zs1) + sum(w0[br,,t,sex]*zs0) - sum(numersl/denom)
        ###bsi
        ###bli
        numerli <- zi1*exp(t1) + zi0*exp(t0)
        del[4] <- del[4] + sum(w1[br,,t,sex]*zi1) + sum(w0[br,,t,sex]*zi0) - sum(numerli/denom)
        ###bt
        numert <- zt1*exp(t1) + zt0*exp(t0)
        del[5] <- del[5] + sum(w1[br,,t,sex]*zt1) + sum(w0[br,,t,sex]*zt0) - sum(numert/denom)
        ###bsex
        numersex <- zsex1*exp(t1) + zsex0*exp(t0)
        del[6] <- del[6] + sum(w1[br,,t,sex]*zsex1) + sum(w0[br,,t,sex]*zsex0) - sum(numersex/denom)
        #########diagonal of the hessian matrix########
        ###g, g
        numer <- exp(t1) * exp(t0)
        hes[1, 1] <- hes[1, 1] + sum(numer/denom^2)
        ###bsl, bsl
        term1 <- ( zs1^2*exp(t1) + zs0^2*exp(t0) ) / denom
        term2 <- - numersl^2 / denom^2
        hes[2, 2] <- hes[2, 2] + sum(term1 + term2)
        ###bsi, bsi
        ###bli, bli
        term1 <- ( zi1^2*exp(t1) + zi0^2*exp(t0) ) / denom
        term2 <- - numerli^2 / denom^2
        hes[4, 4] <- hes[4, 4] + sum(term1 + term2)
        ###bt, bt
        term1 <- ( zt1^2*exp(t1) + zt0^2*exp(t0) ) / denom
        term2 <- - numert^2 / denom^2
        hes[5, 5] <- hes[5, 5] + sum(term1 + term2)
        ###bsex, bsex
        term1 <- ( zsex1^2*exp(t1) + zsex0^2*exp(t0) ) / denom
        term2 <- - numersex^2 / denom^2
        hes[6, 6] <- hes[6, 6] + sum(term1 + term2)
        #########off-diagonal of the hessian matrix########
        ###g, bsl
        term1 <- zs1*exp(t1) / denom
        term2 <- - numersl * numerg / denom^2
        hes[1, 2] <- hes[1, 2] + sum(term1 + term2)
        ###g, bsi
        ###g, bli
        term1 <- zi1*exp(t1) / denom
        term2 <- - numerli * numerg / denom^2
        hes[1, 4] <- hes[1, 4] + sum(term1 + term2)
        ###g, bt
        term1 <- zt1*exp(t1) / denom
        term2 <- - numert * numerg / denom^2
        hes[1, 5] <- hes[1, 5] + sum(term1 + term2)
        ###g, bsex
        term1 <- zsex1*exp(t1) / denom
        term2 <- - numersex * numerg / denom^2
        hes[1, 6] <- hes[1, 6] + sum(term1 + term2)  
        ###bsl, bsi
        ###bsl, bli
        term1 <- ( zs1*zi1*exp(t1) + zs0*zi0*exp(t0) ) / denom
        term2 <- - numersl * numerli / denom^2
        hes[2, 4] <- hes[2, 4] + sum(term1 + term2)
        ###bsl, bt
        term1 <- ( zs1*zt1*exp(t1) + zs0*zt0*exp(t0) ) / denom
        term2 <- - numersl * numert / denom^2
        hes[2, 5] <- hes[2, 5] + sum(term1 + term2)
        ###bsl, bsex
        term1 <- ( zs1*zsex1*exp(t1) + zs0*zsex0*exp(t0) ) / denom
        term2 <- - numersl * numersex / denom^2
        hes[2, 6] <- hes[2, 6] + sum(term1 + term2)
        ###bsi, bli
        ###bsi, bt 
        ###bsi, bsex
        ###bli, bt
        term1 <- ( zi1*zt1*exp(t1) + zi0*zt0*exp(t0) ) / denom
        term2 <- - numerli * numert / denom^2
        hes[4, 5] <- hes[4, 5] + sum(term1 + term2)
        ###bli, bsex
        term1 <- ( zi1*zsex1*exp(t1) + zi0*zsex0*exp(t0) ) / denom
        term2 <- - numerli * numersex / denom^2
        hes[4, 6] <- hes[4, 6] + sum(term1 + term2)
        ###bt, bsex
        term1 <- ( zt1*zsex1*exp(t1) + zt0*zsex0*exp(t0) ) / denom
        term2 <- - numert * numersex / denom^2
        hes[5, 6] <- hes[5, 6] + sum(term1 + term2)       
      }
    }  
    br <- 3
    for (t in 1:nt){
      for (sex in 1:2){
        nsex <- 3 - sex
        if (t==1){
          t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
          t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t+1,sex] 
        } else if(t==nt) {
          t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
          t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] 
        } else {
          t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
          t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
          zt1 <- meanf[br,,t-1,sex] + meanf[br,,t+1,sex] 
          zt0 <- 1 - meanf[br,,t-1,sex] + 1 - meanf[br,,t+1,sex] 
        }
        denom <- exp(t1) + exp(t0)
        zs1 <- meanf[1,,t,sex]
        zs0 <- 1 - meanf[1,,t,sex]
        zl1 <- meanf[2,,t,sex]
        zl0 <- 1 - meanf[2,,t,sex]
        zsex1 <- meanf[br,,t,nsex]
        zsex0 <- 1 - meanf[br,,t,nsex]
        ###g
        numerg <- exp(t1)
        del[1] <- del[1] + sum(w1[br,,t,sex]) - sum(numerg/denom)
        ###bsl
        ###bsi
        numersi <- zs1*exp(t1) + zs0*exp(t0)
        del[3] <- del[3] + sum(w1[br,,t,sex]*zs1) + sum(w0[br,,t,sex]*zs0) - sum(numersi/denom)
        ###bli
        numerli <- zl1*exp(t1) + zl0*exp(t0)
        del[4] <- del[4] + sum(w1[br,,t,sex]*zl1) + sum(w0[br,,t,sex]*zl0) - sum(numerli/denom)
        ###bt
        numert <- zt1*exp(t1) + zt0*exp(t0)
        del[5] <- del[5] + sum(w1[br,,t,sex]*zt1) + sum(w0[br,,t,sex]*zt0) - sum(numert/denom)
        ###bsex
        numersex <- zsex1*exp(t1) + zsex0*exp(t0)
        del[6] <- del[6] + sum(w1[br,,t,sex]*zsex1) + sum(w0[br,,t,sex]*zsex0) - sum(numersex/denom)
        #########diagonal of the hessian matrix########
        ###g, g
        numer <- exp(t1) * exp(t0)
        hes[1, 1] <- hes[1, 1] + sum(numer/denom^2)
        ###bsl, bsl
        ###bsi, bsi
        term1 <- ( zs1^2*exp(t1) + zs0^2*exp(t0) ) / denom
        term2 <- - numersi^2 / denom^2
        hes[3, 3] <- hes[3, 3] + sum(term1 + term2)
        ###bli, bli
        term1 <- ( zl1^2*exp(t1) + zl0^2*exp(t0) ) / denom
        term2 <- - numerli^2 / denom^2
        hes[4, 4] <- hes[4, 4] + sum(term1 + term2)
        ###bt, bt
        term1 <- ( zt1^2*exp(t1) + zt0^2*exp(t0) ) / denom
        term2 <- - numert^2 / denom^2
        hes[5, 5] <- hes[5, 5] + sum(term1 + term2)
        ###bsex, bsex
        term1 <- ( zsex1^2*exp(t1) + zsex0^2*exp(t0) ) / denom
        term2 <- - numersex^2 / denom^2
        hes[6, 6] <- hes[6, 6] + sum(term1 + term2)
        #########off-diagonal of the hessian matrix########
        ###g, bsl
        ###g, bsi
        term1 <- zs1*exp(t1) / denom
        term2 <- - numersi * numerg / denom^2
        hes[1, 3] <- hes[1, 3] + sum(term1 + term2)
        ###g, bli
        term1 <- zl1*exp(t1) / denom
        term2 <- - numerli * numerg / denom^2
        hes[1, 4] <- hes[1, 4] + sum(term1 + term2)
        ###g, bt
        term1 <- zt1*exp(t1) / denom
        term2 <- - numert * numerg / denom^2
        hes[1, 5] <- hes[1, 5] + sum(term1 + term2)
        ###g, bsex
        term1 <- zsex1*exp(t1) / denom
        term2 <- - numersex * numerg / denom^2
        hes[1, 6] <- hes[1, 6] + sum(term1 + term2)  
        ###bsl, bsi
        ###bsl, bli
        ###bsl, bt
        ###bsl, bsex
        ###bsi, bli
        term1 <- ( zs1*zl1*exp(t1) + zs0*zl0*exp(t0) ) / denom
        term2 <- - numersi * numerli / denom^2
        hes[3, 4] <- hes[3, 4] + sum(term1 + term2)
        ###bsi, bt 
        term1 <- ( zs1*zt1*exp(t1) + zs0*zt0*exp(t0) ) / denom
        term2 <- - numersi * numert / denom^2
        hes[3, 5] <- hes[3, 5] + sum(term1 + term2)
        ###bsi, bsex
        term1 <- ( zs1*zsex1*exp(t1) + zs0*zsex0*exp(t0) ) / denom
        term2 <- - numersi * numersex / denom^2
        hes[3, 6] <- hes[3, 6] + sum(term1 + term2)
        ###bli, bt
        term1 <- ( zl1*zt1*exp(t1) + zl0*zt0*exp(t0) ) / denom
        term2 <- - numerli * numert / denom^2
        hes[4, 5] <- hes[4, 5] + sum(term1 + term2)
        ###bli, bsex
        term1 <- ( zl1*zsex1*exp(t1) + zl0*zsex0*exp(t0) ) / denom
        term2 <- - numerli * numersex / denom^2
        hes[4, 6] <- hes[4, 6] + sum(term1 + term2)
        ###bt, bsex
        term1 <- ( zt1*zsex1*exp(t1) + zt0*zsex0*exp(t0) ) / denom
        term2 <- - numert * numersex / denom^2
        hes[5, 6] <- hes[5, 6] + sum(term1 + term2)       
      }
    }
    hes <- -hes
    hes <- t(hes) + hes
    diag(hes) <- diag(hes)/2
    if (is.na(det(hes))){
      return("NA")
    }
    if (abs(det(hes))<=1/10^6){
      return("NA")
    }
    para <- parapre - solve(hes) %*% del
    iter <- iter + 1
    if (sum(is.na(para))){
      return("NA")
    }
  }
  converge <- 1
  if (iter == (maxiter + 1)){
    print("Not converged, please increase alpha or maxiter")
    converge <- 0
  }
  return(list(para = para, converge = converge))
}

calmf1  <- function(meanf, fy1, fy0, paraMRF, option = c("mf", "simf"))
{
  g <- paraMRF[1]; bsl <- paraMRF[2]; bsi <- paraMRF[3]; bli <- paraMRF[4]; bt <- paraMRF[5]; bsex <- paraMRF[6]
  tmp <- dim(fy1)
  nbr <- tmp[1]; ng <- tmp[2]; nt <- tmp[3]; nsex <- tmp[4];
  br <- 1
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      ##the term when z=1
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      ##the term when z=0
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
    sex <- 2
    nsex <- 1
    if (t==1){
      ##the term when z=1
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      ##the term when z=0
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[2,,t,sex] + bsi*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[2,,t,sex]) + bsi*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
  }  
  br <- 2
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
    sex <- 2
    nsex <- 1
    if (t==1){
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsl*meanf[1,,t,sex] + bli*meanf[3,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsl*(1-meanf[1,,t,sex]) + bli*(1-meanf[3,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
  }
  br <- 3
  for (t in 1:nt){
    sex <- 1
    nsex <- 2
    if (t==1){
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
    sex <- 2
    nsex <- 1
    if (t==1){
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t+1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else if(t==nt) {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*meanf[br,,t-1,sex] + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]) + bsex*(1-meanf[br,,t,nsex])
    } else {
      t1 <- g + bsi*meanf[1,,t,sex] + bli*meanf[2,,t,sex] + bt*(meanf[br,,t-1,sex]+meanf[br,,t+1,sex]) + bsex*meanf[br,,t,nsex]
      t0 <- bsi*(1-meanf[1,,t,sex]) + bli*(1-meanf[2,,t,sex]) + bt*(1-meanf[br,,t-1,sex]+1-meanf[br,,t+1,sex]) + bsex*(1-meanf[br,,t,nsex])
    }
    prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
    if (option == "mf"){
      meanf[br,,t,sex] <- prob
    } else if (option == "simf"){
      meanf[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    } else {
      print("please enter a valid option: mf or simf")
    }
  }
  return(meanf)
}

getMRFDE <- function(fy1, fy0, p1=NULL, mf1=NULL, option="simf", paraMRFIni=NULL, iterEM=200, iterGibbsPost=10500,brPost=500, skip=5){ 
  if(is.null(paraMRFIni)) {
    paraMRFIni <- rep(0,6)
  }
  if(is.null(mf1)) {
    mf1 <- ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  }
  paraMRF <- paraMRFIni
  paraMRFTrace <- c()
  converge <- c()
  cat("Estimating MRF parameters,iteration:\n")  
  for (j in 1:iterEM) {
    write(j, file="iter_EM.txt")
    mf1 <- calmf1(meanf=mf1, fy1=fy1, fy0=fy0, paraMRF=paraMRF, option = option)
    w1 <- calw(fy1, fy0, paraMRF, mf1)
    tmp <- optimNR(paraIni = rep(0, 6), mf1, w1, alpha = 10^(-6), maxiter = 500)
    if ((tmp=="NA")[1]){
      return("NA")
    }
    paraMRF <- as.numeric(tmp$para)
    paraMRFTrace <- rbind(paraMRFTrace, paraMRF)
    converge <- c(converge, tmp$converged)
    cat(paste(j, " ",sep="")); flush.console()
  }
  statesI <- ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  pfdr1 <- get_states_nopara(fy1, fy0, paraMRF, statesI, iterGibbsPost, brPost, skip)
  list(pfdr1=pfdr1,paraMRF=paraMRF,paraMRFTrace=paraMRFTrace)
}

get_states_nopara <- function(fy1, fy0, paraMRF, statesI, iterGibbs, br, skip){
  g <- paraMRF[1]; bsl <- paraMRF[2]; bsi <- paraMRF[3]; bli <- paraMRF[4]; bt <- paraMRF[5]; bsex <- paraMRF[6]
  tmp <- dim(fy1)
  nbr <- tmp[1]; ng <- tmp[2]; nt <- tmp[3]; nsex <- tmp[4];
  count <- 0
  statessum <- statesI*0
  for (iter in 1:iterGibbs){
    write(iter, file="iter_post.txt")
    br <- 1
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      
      sex <- 2
      nsex <- 1
      if (t==1){
        ##the term when z=1
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        ##the term when z=0
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsl*statesI[2,,t,sex] + bsi*statesI[3,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[2,,t,sex]) + bsi*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    }  
    br <- 2
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsl*statesI[1,,t,sex] + bli*statesI[3,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsl*(1-statesI[1,,t,sex]) + bli*(1-statesI[3,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    }
    br <- 3
    for (t in 1:nt){
      sex <- 1
      nsex <- 2
      if (t==1){
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
      
      sex <- 2
      nsex <- 1
      if (t==1){
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*statesI[br,,t+1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else if(t==nt) {
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*statesI[br,,t-1,sex] + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t-1,sex]) + bsex*(1-statesI[br,,t,nsex])
      } else {
        t1 <- g + bsi*statesI[1,,t,sex] + bli*statesI[2,,t,sex] + bt*(statesI[br,,t-1,sex]+statesI[br,,t+1,sex]) + bsex*statesI[br,,t,nsex]
        t0 <- bsi*(1-statesI[1,,t,sex]) + bli*(1-statesI[2,,t,sex]) + bt*(1-statesI[br,,t-1,sex]+1-statesI[br,,t+1,sex]) + bsex*(1-statesI[br,,t,nsex])
      }
      prob <- exp(t1)*fy1[br,,t,sex]/(exp(t1)*fy1[br,,t,sex]+ exp(t0)*fy0[br,,t,sex])
      statesI[br,,t,sex] <- (prob >= runif(length(prob))) + 0
    }
    if (iter >= br && (iter%%skip)==0)
    {
      count <- count+1
      statessum <- statessum + statesI
    }    
  }
  pfdr1 <- statessum/count
  return(pfdr1) 
}

getMRFDEALL <- function(fy1, fy0, p1=NULL, mf1=NULL, option="simf", paraMRFIni=NULL, iterEM=200, iterGibbsPost=10500,brPost=500, skip=5){ 
  if(is.null(paraMRFIni)) {
    paraMRFIni <- rep(0,6)
  }
  if(is.null(mf1)) {
    mf1 <- ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  }
  paraMRF <- paraMRFIni
  paraMRFTrace <- c()
  converge <- c()
  cat("Estimating MRF parameters,iteration:\n")  
  for (j in 1:iterEM) {
    write(j, file="iter_EM.txt")
    mf1 <- calmf1(meanf=mf1, fy1=fy1, fy0=fy0, paraMRF=paraMRF, option = option)
    w1 <- calw(fy1, fy0, paraMRF, mf1)
    tmp <- optimNR(paraIni = rep(0, 6), mf1, w1, alpha = 10^(-6), maxiter = 500)
    if (tmp=="NA"){
      return("NA")
    }
    paraMRF <- as.numeric(tmp$para)
    paraMRFTrace <- rbind(paraMRFTrace, paraMRF)
    converge <- c(converge, tmp$converged)
    cat(paste(j, " ",sep="")); flush.console()
  }
  statesI <- ((p1*fy1/(p1*fy1 + (1-p1)*fy0)) > runif(length(fy0))) + 0
  para_200 <- paraMRF
  pfdr1_200 <- get_states_nopara(fy1, fy0, para_200, statesI, iterGibbsPost, brPost, skip)
  para_100 <- paraMRFTrace[100,]
  pfdr1_100 <- get_states_nopara(fy1, fy0, para_100, statesI, iterGibbsPost, brPost, skip)
  para_mean <- apply(paraMRFTrace[100:200,], 2, mean)
  pfdr1_mean <- get_states_nopara(fy1, fy0, para_mean, statesI, iterGibbsPost, brPost, skip)
  para_median <- apply(paraMRFTrace[100:200,], 2, median)
  pfdr1_median <- get_states_nopara(fy1, fy0, para_median, statesI, iterGibbsPost, brPost, skip)
  list(pfdr1_200=pfdr1_200, pfdr1_100=pfdr1_100, pfdr1_mean=pfdr1_mean, pfdr1_median=pfdr1_median, paraMRF=paraMRF,paraMRFTrace=paraMRFTrace)
}

calZ <- function(data, reg, tp, sex, group, method = "TMM")
{
  y <- DGEList(counts = data, group = group)
  y <- calcNormFactors(y, method = method)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  dispersion <- y$tagwise.dispersion
  pseuy <- y$pseudo.counts
  
  nbr <- length(unique(reg))
  nt <- length(unique(tp)) - 1
  ng <- dim(data)[1]
  z <- array(dim=c(nbr, ng, nt,2))
  for (b in 1:nbr){
    for (t in 1:nt){
      for (s in 1:2){
        l1 <- reg == b & tp == t & sex == s
        l2 <- reg == b & tp == t+1 & sex == s
        #l1 <- reg == b & tp == t 
        #l2 <- reg == b & tp == t+1 
        #c <- c(rep(1, sum(l1)), rep(2, sum(l2))); c <- as.factor(c)
        #y <- DGEList(counts = cbind(data[,l1], data[,l2]), group = c)
        #y <- calcNormFactors(y, method = method)
        #y <- estimateCommonDisp(y)
        #y <- estimateTagwiseDisp(y)
        #dispersion <- y$tagwise.dispersion
        #pseuy <- y$pseudo.counts   
        d1 <- as.matrix(pseuy[, l1])
        d2 <- as.matrix(pseuy[, l2])
        #d <- cbind(d1, d2)
        p <- myexactTestDoubleTail(d1, d2, dispersion=dispersion)
        z[b, ,t,s] <- qnorm(p)
      }
    }
  }
  z[z == Inf] <- max(z[is.finite(z)])
  z[z == -Inf] <- min(z[is.finite(z)])
  return(z)
}

####the modified version of edgeR
myexactTestDoubleTail <- function (y1, y2, dispersion = 0, big.count = 900) 
{
  ntags <- NROW(y1)
  n1 <- NCOL(y1)
  n2 <- NCOL(y2)
  if (n1 > 1) 
    s1 <- round(rowSums(y1))
  else s1 <- round(y1)
  if (n2 > 1) 
    s2 <- round(rowSums(y2))
  else s2 <- round(y2)
  if (length(dispersion) == 1) 
    dispersion <- rep(dispersion, ntags)
  s <- s1 + s2
  mu <- s/(n1 + n2)
  mu1 <- n1 * mu
  mu2 <- n2 * mu
  q <- rep(1, ntags)
  pois <- dispersion <= 0
  if (any(pois)) {
    q[pois] <- mybinomTest(s1[pois], s2[pois], p = n1/(n1 + n2))
  }
  big <- s1 > big.count & s2 > big.count
  if (any(big)) {
    y1 <- as.matrix(y1)
    y2 <- as.matrix(y2)
    q[big] <- myexactTestBetaApprox(y1[big, , drop = FALSE], 
                                    y2[big, , drop = FALSE], dispersion[big])
  }
  p.bot <- size1 <- size2 <- rep(0, ntags)
  left <- s1 <= mu1 & !pois & !big
  if (any(left)) {
    p.bot[left] <- dnbinom(s[left], size = (n1 + n2)/dispersion[left], 
                           mu = s[left])
    size1[left] <- n1/dispersion[left]
    size2[left] <- n2/dispersion[left]
    for (g in which(left)) {
      ##in the case when all counts are 0, will give pval=1/2. 
      ##For the quantile calculation for discrete variables, use 1/2
      if (s1[g] == 0){
        q[g] <-  1/2 * dnbinom(0, size = size1[g], mu = mu1[g]) * 
          dnbinom(s[g] - 0, size = size2[g], mu = mu2[g])
      } else {
        x <- 0:(s1[g] - 1)
        p.top <- dnbinom(x, size = size1[g], mu = mu1[g]) * 
          dnbinom(s[g] - x, size = size2[g], mu = mu2[g])
        ##
        q[g] <-  sum(p.top) + 1/2 * dnbinom(s1[g], size = size1[g], mu = mu1[g]) * 
          dnbinom(s2[g], size = size2[g], mu = mu2[g])
      }
    }
    q[left] <- q[left]/p.bot[left]
  }
  right <- s1 > mu1 & !pois & !big
  if (any(right)) {
    p.bot[right] <- dnbinom(s[right], size = (n1 + n2)/dispersion[right], 
                            mu = s[right])
    size1[right] <- n1/dispersion[right]
    size2[right] <- n2/dispersion[right]
    for (g in which(right)) {
      if (s2[g] == 0){
        q[g] <-  1/2 * dnbinom(s1[g], size = size1[g], mu = mu1[g]) * 
          dnbinom(0, size = size2[g], mu = mu2[g])
      } else {
        x <- (s1[g] + 1):s[g]
        p.top <- dnbinom(x, size = size1[g], mu = mu1[g]) * 
          dnbinom(s[g] - x, size = size2[g], mu = mu2[g])
        ##
        q[g] <- sum(p.top) + 1/2 * dnbinom(s1[g], size = size1[g], mu = mu1[g]) * 
          dnbinom(s2[g], size = size2[g], mu = mu2[g]) 
      }
    }
    ##
    q[right] <- 1-q[right]/p.bot[right]
  }
  ##
  q
}

myexactTestBetaApprox <- function (y1, y2, dispersion = 0)
{
  ntags <- NROW(y1)
  n1 <- NCOL(y1)
  n2 <- NCOL(y2)
  if (n1 > 1) 
    y1 <- rowSums(y1)
  if (n2 > 1) 
    y2 <- rowSums(y2)
  if (length(dispersion) == 1) 
    dispersion <- rep(dispersion, ntags)
  y <- y1 + y2
  mu <- y/(n1 + n2)
  alpha1 <- n1 * mu/(1 + dispersion * mu)
  alpha2 <- n2/n1 * alpha1
  return(pbeta(y1/y, alpha1, alpha2))
}

mybinomTest <- function (y1, y2, p)
{
  if (length(y1) != length(y2)) 
    stop("y1 and y2 must have same length")
  if (any(is.na(y1)) || any(is.na(y2))) 
    stop("missing values not allowed")
  y1 <- round(y1)
  y2 <- round(y2)
  if (any(y1 < 0) || any(y2 < 0)) 
    stop("y1 and y2 must be non-negative")
  if (p <= 0 || p >= 1) 
    stop("p must be between 0 and 1")
  size <- y1 + y2
  q <- pbinom(y1 - 1, size = size, prob = p) + 
    1/2 * dbinom(y1, size = size, prob = p)
  return(q)
}

calf <- function(data, reg, tp, sex, group, method = "TMM", nulltype=1, df=15){
  z <- calZ(data, reg, tp, sex, group, method = method)
  tmp <- locfdr(z, nulltype=nulltype, df=df)
  if (nulltype==1){
    p0 <- (tmp$fp0)[3, 3]
  } else {
    p0 <- (tmp$fp0)[1, 3]
  }
  p1 <- 1 - p0  
  fdr <- array(tmp$fdr, dim=dim(z))
  f0 <- fdr/p0
  f1 <- (1-fdr)/p1 
  return(list(fy1=f1, fy0=f0, fdr=fdr, p1=p1))
}

getpedgeRori <- function(data, reg, tp, sex, method = "TMM")
{
  br <- length(unique(reg))
  ng <- dim(data)[1]
  nt <- length(unique(tp))-1
  ns <- length(unique(sex))
  pv <- array(dim=c(br, ng, nt, ns))
  ###potential issues of different data format of z and zper
  for (b in 1:br){
    for (t in 1:nt){
      for (s in 1:ns){
        l1 <- reg == b & tp == t & sex==s
        l2 <- reg == b & tp == t+1 & sex==s
        if (any(l1) & any(l2)){
          d1 <- data[, l1]
          d2 <- data[, l2]
          d<- cbind(d1, d2)
          c <- c(rep(1, dim(d1)[2]), rep(2, dim(d2)[2])); c <- as.factor(c)       
          y <- DGEList(counts = d, group = c)
          y <- calcNormFactors(y, method = method)
          y <- estimateCommonDisp(y)
          y <- estimateTagwiseDisp(y)
          p <- exactTest(y)
          pv[b, ,t, s] <- (p$table)[,3]
        } else {
          
        }
      }
    }
  } 
  return(pv)
}

calZo <- function(data, reg, tp, sex, group, method = "TMM")
{
  br <- length(unique(reg))
  ng <- dim(data)[1]
  nt <- length(unique(tp))-1
  ns <- length(unique(sex))
  z <- array(dim=c(br, ng, nt, ns))
  ###potential issues of different data format of z and zper
  for (b in 1:br){
    for (t in 1:nt){
      for (s in 1:ns){
        l1 <- reg == b & tp == t & sex==s
        l2 <- reg == b & tp == t+1 & sex==s
        if (any(l1) & any(l2)){
          d1 <- data[, l1]
          d2 <- data[, l2]
          d<- cbind(d1, d2)
          c <- c(rep(1, dim(d1)[2]), rep(2, dim(d2)[2])); c <- as.factor(c)       
          y <- DGEList(counts = d, group = c)
          y <- calcNormFactors(y, method = method)
          y <- estimateCommonDisp(y)
          y <- estimateTagwiseDisp(y)
          dispersion <- y$tagwise.dispersion
          pseuy <- y$pseudo.counts
          d1 <- as.matrix(pseuy[, c==1])
          d2 <- as.matrix(pseuy[, c==2])
          #d <- cbind(d1, d2)
          p <- myexactTestDoubleTail(d1, d2, dispersion=dispersion)
          z[b, ,t,s] <- qnorm(p)
        } else {
          
        }
      }
    }
  } 
  z[z == Inf] <- max(z[is.finite(z)])
  z[z == -Inf] <- min(z[is.finite(z)])
  return(z)
}


getpedgeR <- function(data, reg, tp, sex, group, method = "TMM")
{
  y <- DGEList(counts = data, group = group)
  y <- calcNormFactors(y, method = method)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  dispersion <- y$tagwise.dispersion
  pseuy <- y$pseudo.counts
  
  nbr <- length(unique(reg))
  nt <- length(unique(tp)) - 1
  ng <- dim(data)[1]
  pv <- array(dim=c(nbr, ng, nt,2))
  fdrBH <- pv
  for (b in 1:nbr){
    for (t in 1:nt){
      for (s in 1:2){
        l1 <- unique(group[reg == b & tp == t & sex == s])
        l2 <- unique(group[reg == b & tp == t+1 & sex == s])
        p <- exactTest(y, pair=c(l1, l2))
        pv[b, ,t, s] <- (p$table)[,3]
        tmp <- (topTags(p, n=ng))$table
        tmp <- cbind(as.numeric(row.names(tmp)), tmp[,4])
        tmp <- tmp[order(tmp[,1]),]
        fdrBH[b, , t, s] <- tmp[,2]
      }
    }
  }
  return(list(pv=pv, fdrBH=fdrBH))
}

getpedgeR_data <- function(data, reg, tp, sex, group, method = "TMM")
{
  y <- DGEList(counts = data, group = group)
  y <- calcNormFactors(y, method = method)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  dispersion <- y$tagwise.dispersion
  pseuy <- y$pseudo.counts
  
  nbr <- length(unique(reg))
  nt <- length(unique(tp)) - 1
  ng <- dim(data)[1]
  pv <- array(dim=c(nbr, ng, nt,2))
  fdrBH <- pv
  for (b in 1:nbr){
    for (t in 1:nt){
      for (s in 1:2){
        p <- exactTestDoubleTail(pseuy[,reg == b & tp == t & sex == s], pseuy[,reg == b & tp == t+1 & sex == s], dispersion=dispersion)
        pv[b, ,t, s] <- p
        #tmp <- (topTags(p, n=ng))$table
        #tmp <- cbind(as.numeric(row.names(tmp)), tmp[,4])
        #tmp <- tmp[order(tmp[,1]),]
        #fdrBH[b, , t, s] <- tmp[,2]
      }
    }
  }
  return(pv)
}
