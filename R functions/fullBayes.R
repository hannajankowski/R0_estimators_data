

#rm(list = ls())

########### Fullbayes SIR #####################
fullbayes_SIR <- function(data, sbeta, rbeta, sgamma, rgamma, theta, N){
  set.seed(1)
  
  priorI <- function(theta,m,w,j){
    if(j==1){
      if(m==1){
        I0 <- -rexp(theta)  
        Itime <- I0
      }else{
        I0 <- -rexp(theta)
        I <- runif(m-1,w[j],w[j+1])
        Itime <- c(I0,I)
      }
      Itime <- sort(Itime)  
    }else{
      Itime <- runif(m,w[j],w[j+1])
      Itime <- sort(Itime)  
    }
    return(Itime)
  }
  
  priorR <- function(n,w,j){
    Rtime <- runif(n,w[j],w[j+1])
    Rtime <- sort(Rtime)
    return(Rtime)
  }
  
  genI <- function(Itime,beta,gamma,theta,w,j){
    newI0 <- -rexp(beta+gamma+theta)
    IStime <- sample(Itime[Itime > 0],1)  
    genItime <- runif(1,w[j],w[j+1])
    Itime[which(Itime == IStime)] <- genItime
    Itime <- sort(Itime)  
    return(Itime)
  }
  
  genR <- function(Rtime,w,j){
    RStime <- sample(Rtime,1)  
    genRtime <- runif(1,w[j],w[j+1])
    Rtime[which(Rtime == RStime)] <- genRtime
    Rtime <- sort(Rtime)   
    return(Rtime)
  }
  
  likj1 <- function(Itime,N,beta){
    l <- length(Itime)
    a <- cbind(Itime,rep(1,l))
    numS <- N-cumsum(a[,2])
    numI <- cumsum(a[,2])
    nr <- rep(0,nrow(a))
    numdat <- cbind(a,nr,numS,numI)
    tt <- numdat[,1]
    ss <- numdat[,4]
    ii <- numdat[,5]
    si1 <- ss*ii
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    f1 <- prod(beta*(ss/N)*ii)
    f2 <- exp(-(beta*intSI))
    lik <- f1*f2
    return(c(intSI,intI,lik))
  }
  
  likotherj <- function(Itime,Rtime,cm,cn,N,beta,gamma){
    a <- cbind(Itime,rep(1,sum(cm)))
    b <- cbind(Rtime,rep(-1,sum(cn)))
    c <- cbind(Rtime,rep(0,sum(cn))) 
    dt1 <- rbind(a,b)
    dt2 <- rbind(a,c)
    dt1 <- dt1[order(dt1[,1]),]
    dt2 <- dt2[order(dt2[,1]),]        
    dt <- as.data.frame(cbind(dt1,dt2[,2]))
    numS <- N-cumsum(dt[,3])
    numI <- cumsum(dt[,2])
    numdat <- cbind(dt,numS,numI)
    
    tt <- numdat[,1]
    ss <- numdat[,4]
    ii <- numdat[,5]
    
    si1 <- ss*ii
    l <- length(si1)
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    f1 <- prod(gamma*ii)*prod(beta*(ss/N)*ii)
    f2 <- exp(-(beta*intSI+gamma*intI))
    lik <- f1*f2
    return(c(intSI,intI,lik))
  }
  
  
  sbeta <- sbeta	#shape of beta   
  rbeta <- rbeta 	#rate of beta   
  sgamma <- sgamma	#shape of gamma   
  rgamma <- rgamma	#rate of gamma   
  theta <- theta	#rate of I0   
  N <- N  #the population size
  burn <- 100
  iter <- burn+1000  #burn in 100 rounds
  nc <- length(data)   
  w <- seq(0,nc*7,7)
  
  cm <- c()
  cn <- c()
  resbeta <- rep(0,nc)
  resgamma <- rep(0,nc)
  resR0 <- rep(0,nc)
  SI <- rep(0,nc)
  cItime <- c()
  cRtime <- c()
  
  
  for (j in 1:nc){
    m <-  data[j]    
    cm <- c(cm,m)
    if(j==1){
      n <- 1
      cn <- c(cn,n)
    }else if(j==2){
      n <- data[j-1]-1
      cn <- c(cn,n)
    }else{
      n <- data[j-1]
      cn <- c(cn,n)
    }
    
    if(m==0){
      Itimew <- 0
    }else{
      Itimew <- priorI(theta,m,w,j)  
    }
    if(n==0){
      Rtimew <- 0
    }else{
      Rtimew <- priorR(n,w,j)
    }
    
    ## prior for beta, gamma
    beta <- rgamma(1, sbeta, rbeta)
    gamma <- rgamma(1, sgamma, rgamma)
    
    ##################### Gibbs within Metropolis
    if(m==1 && j==1){
      cItime <- Itimew
      Rtimew <- Rtimew[which(Rtimew!=0)]
      cRtime <- Rtimew  
      
      resbeta[j] <- 0
      resgamma[j] <- 0
      resR0[j] <- 0
      SI <- "Cannot find since no new infectious"
    }else{
      p <- 1
      prob <- rep(0,iter)
      hatbeta <- rep(0,iter)
      hatgamma <- rep(0,iter)
      hatR0 <- rep(0,iter)
      while (p <= iter){
        if(m!=0){
          oldItimew <- Itimew
          Itimew <- genI(Itimew,beta,gamma,theta,w,j)  
          newItimew <- Itimew
        }else{
          oldItimew <- Itimew
          newItimew <- 0
        }
        if(n!=0){
          oldRtimew <- Rtimew
          Rtimew <- genR(Rtimew,w,j)
          newRtimew <- Rtimew
        }else{
          newRtimew <- 0
          oldRtimew <- Rtimew
        }
        
        Itime <- c(cItime,newItimew)
        Itime <- Itime[which(Itime!=0)]
        oldItime <- c(cItime,oldItimew)
        oldItime <- oldItime[which(oldItime!=0)]
        Rtime <- c(cRtime,newRtimew)
        Rtime <- Rtime[which(Rtime!=0)]
        oldRtime <- c(cRtime,oldRtimew)
        oldRtime <- oldRtime[which(oldRtime!=0)]
        
        #likelihood for the whole data
        if(j==1){
          newlik <- likj1(Itime,N,beta)[3]
          oldlik <- likj1(oldItime,N,beta)[3]
        }else{
          newlik <- likotherj(Itime,Rtime,cm,cn,N,beta,gamma)[3]
          oldlik <- likotherj(oldItime,oldRtime,cm,cn,N,beta,gamma)[3]
        }
        
        ##acceptance probability
        likratio <- newlik/oldlik 
        cri <- c(1,likratio)
        prob[p] <- cri[which.min(cri)]
        u <- runif(1)
        if(u <= prob[p]){
          Itimew <- newItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          Rtimew <- newRtimew 
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew) 
        }else{
          Itimew <- oldItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          Rtimew <- oldRtimew
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew) 
        }
        
        ##update beta and gamma
        sh_beta <- sbeta+sum(cm)-1
        sh_gamma <- sgamma+sum(cn)
        if(j==1){
          ra_beta <- rbeta+likj1(Itime,N,beta)[1]
          ra_gamma <- rgamma+likj1(Itime,N,beta)[2]
        }else{
          ra_beta <- rbeta+likotherj(Itime,Rtime,cm,cn,N,beta,gamma)[1]
          ra_gamma <- rgamma+likotherj(Itime,Rtime,cm,cn,N,beta,gamma)[2]
        }
        beta <- rgamma(1,sh_beta,ra_beta)
        gamma <- rgamma(1,sh_gamma,ra_gamma)
        
        ##Estimate beta and gamma 
        hatbeta[p] <- (sh_beta-1)/ra_beta
        hatgamma[p] <- (sh_gamma-1)/ra_gamma
        if(hatgamma[p]==0){
          hatR0[p] <- 0
        }else{
          hatR0[p] <- hatbeta[p]/hatgamma[p]      
        }
        p <- p+1
        
      }
      cItime <- Itime
      cRtime <- Rtime  
      
      #average from 1,000 rounds (100 burn-in)
      resbeta[j] <- mean(hatbeta[burn:iter])
      resgamma[j] <- mean(hatgamma[burn:iter])
      resR0[j] <- mean(hatR0[burn:iter])
      SI[j] <- 1/resgamma[j]
    }
  }
  return(list(R0=resR0, SI=SI))
}






########### Fullbayes SEIR #####################

fullbayes_SEIR <- function(data, sbeta, rbeta, sgamma, rgamma,ssigma, rsigma, theta, N){
  set.seed(1)
  
  priorI <- function(theta,m,w,j){
    if(j==1){
      if(m==1){
        I0 <- -rexp(theta)  
        Itime <- I0
      }else{
        I0 <- -rexp(theta)
        I <- runif(m-1,w[j],w[j+1])
        Itime <- c(I0,I)
      }
      Itime <- sort(Itime)  
    }else{
      Itime <- runif(m,w[j],w[j+1])
      Itime <- sort(Itime)  
    }
    return(Itime)
  }
  
  priorR <- function(n,w,j){
    Rtime <- runif(n,w[j],w[j+1])
    Rtime <- sort(Rtime)
    return(Rtime)
  }
  
  priorE <- function(e,w,j){
    Etime <- runif(e,w[j],w[j+1])
    Etime <- sort(Etime)
    return(Etime)
  }
  
  genI <- function(Itime,beta,gamma,theta,w,j){
    newI0 <- -rexp(beta+gamma+theta)
    IStime <- sample(Itime[Itime > 0],1)  
    genItime <- runif(1,w[j],w[j+1])
    Itime[which(Itime == IStime)] <- genItime
    Itime <- sort(Itime)   
    return(Itime)
  }
  
  genR <- function(Rtime,w,j){
    RStime <- sample(Rtime,1) 
    genRtime <- runif(1,w[j],w[j+1])
    Rtime[which(Rtime == RStime)] <- genRtime
    Rtime <- sort(Rtime)  
    return(Rtime)
  }
  
  genE <- function(Etime,w,j){
    EStime <- sample(Etime,1) 
    genEtime <- runif(1,w[j],w[j+1])
    Etime[which(Etime == EStime)] <- genEtime
    Etime <- sort(Etime)   
    return(Etime)
  }
  
  likj1 <- function(Itime,Etime,cm,ce,N,beta,sigma){
    a <- cbind(Itime,rep(1,sum(cm)))
    b <- cbind(Etime,rep(0,sum(ce)))
    dt1 <- rbind(a,b)
    
    d <- cbind(Itime,rep(0,sum(cm)))
    e <- cbind(Etime,rep(1,sum(ce)))
    dt2 <- rbind(d,e)
    
    dt1 <- dt1[order(dt1[,1]),]
    dt2 <- dt2[order(dt2[,1]),]        
    dt <- as.data.frame(cbind(dt1,dt2[,2]))
    
    numI <- cumsum(dt[,2])
    numE <- cumsum(dt[,3])
    numS <- N-numI-numE
    
    numdat <- cbind(dt,numS,numE,numI)
    
    tt <- numdat[,1]
    ss <- numdat[,4]
    ee <- numdat[,5]
    ii <- numdat[,6]
    si1 <- ss*ii
    l <- length(si1)
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    ee1 <- ee[-1]
    ee2 <- ee[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intE <- sum(((t1-t2)*(ee1+ee2))/2)
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    
    f1a <-  prod(beta*(ss[which(ii!=0)]/N)*ii[which(ii!=0)])
    f1b <- prod(sigma*ee[which(ee!=0)])
    
    if(f1a!=0|f1b!=0){
      f1 <- f1a*f1b
    }else{
      f1 <- 0
    }
    
    f2 <- exp(-(beta*intSI+sigma*intE))
    if(f2==0){
      lik <- 0 
    }else if(f1==Inf){
      lik <- 999999
    }else{
      lik <- f1*f2
    }
    return(c(intSI,intE,intI,lik))
  }
  
  likotherj <- function(Itime,Etime,Rtime,cm,ce,cn,N,beta,sigma,gamma){
    a <- cbind(Itime,rep(1,sum(cm)))
    b <- cbind(Rtime,rep(-1,sum(cn)))
    c <- cbind(Etime,rep(0,sum(ce)))
    dt1 <- rbind(a,b,c)
    
    d <- cbind(Itime,rep(0,sum(cm)))
    e <- cbind(Rtime,rep(0,sum(cn)))
    f <- cbind(Etime,rep(1,sum(ce)))
    dt2 <- rbind(d,e,f)
    
    dt1 <- dt1[order(dt1[,1]),]
    dt2 <- dt2[order(dt2[,1]),]        
    dt <- as.data.frame(cbind(dt1,dt2[,2]))
    
    numI <- cumsum(dt[,2])
    numE <- cumsum(dt[,3])
    numS <- N-numI-numE
    
    numdat <- cbind(dt,numS,numE,numI)
    
    tt <- numdat[,1]
    ss <- numdat[,4]
    ee <- numdat[,5]
    ii <- numdat[,6]
    
    si1 <- ss*ii
    l <- length(si1)
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    ee1 <- ee[-1]
    ee2 <- ee[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intE <- sum(((t1-t2)*(ee1+ee2))/2)
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    
    f1a <- prod(beta*(ss[which(ii!=0)]/N)*ii[which(ii!=0)])
    f1b <- prod(sigma*ee[which(ee!=0)])
    f1c <- prod(gamma*ii[which(ii!=0)])
    if(f1a!=0|f1b!=0|f1c!=0){
      f1 <- f1a*f1b*f1c
    }else{
      f1 <- 0
    }
    f2 <- exp(-(beta*intSI+sigma*intE+gamma*intI))
    if(f2==0){
      lik <- 0 
    }else if(f1==Inf){
      lik <- 999999
    }else{
      lik <- f1*f2
    }
    return(c(intSI,intE,intI,lik))
  }
  
  
  sbeta <- sbeta	#shape of beta   
  rbeta <- rbeta	#rate of beta   
  sgamma <- sgamma	#shape of gamma   
  rgamma <- rgamma	#rate of gamma  
  ssigma <- ssigma   #shape of sigma
  rsigma <- rsigma	 #rate of sigma
  theta <- theta	#rate of I0   
  N <- N  #the population size
  burn <- 100
  iter <- burn+1000  #burn in 100 rounds
  nc <- length(data)
  w <- seq(0,nc*7,7)
  
  cm <- c()
  cn <- c()
  ce <- c()
  resbeta <- rep(0,nc)
  resgamma <- rep(0,nc)
  ressigma <- rep(0,nc)
  resR0 <- rep(0,nc)
  SI <- rep(0,nc)
  cEtime <- c()
  cItime <- c()
  cRtime <- c()
  
  
  for (j in 1:nc){
    #number of I
    m <-  data[j] 
    cm <- c(cm,m)
    
    #number of R
    if(j==1){
      n <- 1
      cn <- c(cn,n)
    }else if(j==2){
      n <- data[j-1]-1
      cn <- c(cn,n)
    }else{
      n <- data[j-1]
      cn <- c(cn,n)
    }
    
    #number of E
    e <- m
    ce <- c(ce,e)
    
    if(m==0 && n==0){
      Itimew <- 0
      Etimew <- 0
      Rtimew <- 0
    }else if(m==0 && n !=0){
      Itimew <- 0
      Etimew <- 0  
      Rtimew <- priorR(n,w,j)
    }else if(m !=0 && n==0){
      Itimew <- priorI(theta,m,w,j)  
      Etimew <- priorE(e,w,j)
      Rtimew <- 0
    }else{
      Itimew <- priorI(theta,m,w,j)  
      Etimew <- priorE(e,w,j)
      Rtimew <- priorR(n,w,j)
    }
    
    ## prior for sigma, beta, gamma
    beta <- rgamma(1, sbeta, rbeta)
    gamma <- rgamma(1, sgamma, rgamma)
    sigma <- rgamma(1, ssigma, rsigma)
    
    ##################### Gibbs within Metropolis
    if(m==1 && j==1){
      Itimew <- Itimew[which(Itimew!=0)]
      cItime <- Itimew
      Rtimew <- Rtimew[which(Rtimew!=0)]
      cRtime <- Rtimew  
      Etimew <- Etimew[which(Etimew!=0)]
      cEtime <- Etimew  
      
      resbeta[j] <- 0
      resgamma[j] <- 0
      ressigma[j] <- 0
      resR0[j] <- 0
      SI <- "Cannot find since no new infectious"
    }else{
      p <- 1
      prob <- rep(0,iter)
      hatbeta <- rep(0,iter)
      hatgamma <- rep(0,iter)
      hatsigma <- rep(0,iter)
      hatR0 <- rep(0,iter)
      while (p <= iter){
        oldItimew <- Itimew
        oldRtimew <- Rtimew
        oldEtimew <- Etimew
        
        if(m > 1){
          Itimew <- genI(Itimew,beta,gamma,theta,w,j)
        }else{
          if(m==0){
            Itimew <- 0
          }else{
            Itimew <- oldItimew
          }
        }
        
        if(e > 1){
          Etimew <- genE(Etimew,w,j)
        }else{
          if(e==0){
            Etimew <- 0
          }else{
            Etimew <- oldEtimew
          }
        }
        
        if(n > 1){
          Rtimew <- genR(Rtimew,w,j)
        }else{
          if(n==0){
            Rtimew <- 0
          }else{
            Rtimew <- oldRtimew
          }
        }
        
        newItimew <- Itimew
        newEtimew <- Etimew
        newRtimew <- Rtimew
        
        
        Itime <- c(cItime,newItimew)
        Itime <- Itime[which(Itime!=0)]
        oldItime <- c(cItime,oldItimew)
        oldItime <- oldItime[which(oldItime!=0)]
        
        Rtime <- c(cRtime,newRtimew)
        Rtime <- Rtime[which(Rtime!=0)]
        oldRtime <- c(cRtime,oldRtimew)
        oldRtime <- oldRtime[which(oldRtime!=0)]
        
        Etime <- c(cEtime,newEtimew)
        Etime <- Etime[which(Etime!=0)]
        oldEtime <- c(cEtime,oldEtimew)
        oldEtime <- oldEtime[which(oldEtime!=0)]     
        
        if(j==1){
          newlik <- likj1(Itime,Etime,cm,ce,N,beta,sigma)[4]
          oldlik <- likj1(oldItime,Etime,cm,ce,N,beta,sigma)[4]
        }else{
          newlik <- likotherj(Itime,Etime,Rtime,cm,ce,cn,N,beta,sigma,gamma)[4]
          oldlik <- likotherj(oldItime,oldEtime,oldRtime,cm,ce,cn,N,beta,sigma,gamma)[4]
        }
        
        ##acceptance probability
        likratio <- newlik/oldlik
        if(oldlik==0){likratio <- 1}
        cri <- c(1,likratio)
        prob[p] <- cri[which.min(cri)]
        u <- runif(1)
        if(u <= prob[p]){
          Itimew <- newItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          
          Etimew <- newEtimew 
          Etimew <- Etimew[which(Etimew!=0)]
          Etime <- c(cEtime,Etimew)
          
          Rtimew <- newRtimew 
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew) 
        }else{
          Itimew <- oldItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          
          Etimew <- oldEtimew
          Etimew <- Etimew[which(Etimew!=0)]
          Etime <- c(cEtime,Etimew)
          
          Rtimew <- oldRtimew
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew)
        }
        
        ##update beta, gamma and sigma
        sh_beta <- sbeta+sum(cm)
        sh_sigma <- ssigma+sum(ce)
        sh_gamma <- sgamma+sum(cn)
        
        if(j==1){
          ra_beta <- rbeta+likj1(Itime,Etime,cm,ce,N,beta,sigma)[1]
          ra_sigma <- rsigma+likj1(Itime,Etime,cm,ce,N,beta,sigma)[2]
          ra_gamma <- rgamma+likj1(Itime,Etime,cm,ce,N,beta,sigma)[3]
        }else{
          ra_beta <- rbeta+likotherj(Itime,Etime,Rtime,cm,ce,cn,N,beta,sigma,gamma)[1]
          ra_sigma <- rsigma+likotherj(Itime,Etime,Rtime,cm,ce,cn,N,beta,sigma,gamma)[2]
          ra_gamma <- rgamma+likotherj(Itime,Etime,Rtime,cm,ce,cn,N,beta,sigma,gamma)[3]
        }
        
        beta <- rgamma(1,shape=sh_beta,rate=ra_beta)
        sigma <- rgamma(1,shape=sh_sigma,rate=ra_sigma)
        gamma <- rgamma(1,shape=sh_gamma,rate=ra_gamma)
        
        ##Estimate beta, sigma, gamma
        hatbeta[p] <- (sh_beta-1)/ra_beta
        hatsigma[p] <- (sh_sigma-1)/ra_sigma
        hatgamma[p] <- (sh_gamma-1)/ra_gamma
        if(hatgamma[p]==0){
          hatR0[p] <- 0
        }else{
          hatR0[p] <- hatbeta[p]/hatgamma[p]      
        }
        p <- p+1
      }
      
      cItime <- Itime
      cEtime <- Etime
      cRtime <- Rtime  
      
      #average from 1,000 rounds (100 burn-in)
      resbeta[j] <- mean(hatbeta[burn:iter])
      ressigma[j] <- mean(hatsigma[burn:iter])
      resgamma[j] <- mean(hatgamma[burn:iter])
      resR0[j] <- mean(hatR0[burn:iter])
      SI[j] <- (1/resgamma[j])+(1/ressigma[j])
    }
  }
  return(list(R0=resR0, SI=SI))
}




########### Fullbayes SEAIR #####################
fullbayes_SEAIR <- function(data, sbeta, rbeta, sgamma, rgamma,ssigma, rsigma, srho, rrho, theta, N){
  set.seed(1)
  
  priorI0 <- function(theta){
    I0 <- -rexp(theta)  
    I0time <- I0
    return(I0time)
  }
  
  priorEAIw1 <- function(q,w,j){
    time <- runif((3*q-1),w[j],w[j+1])
    time <- sort(time)
    if(q==1){
      Etime <- time[1]
      Atime <- time[2]
      Itime <- 0
    }else{
      Etime <- time[1:q]
      Atime <- time[(q+1):(2*q)]
      Itime <- time[(2*q+1):(3*q-1)] 
    }
    return(list(Etime=Etime, Atime=Atime, Itime=Itime))
  }
  
  priorEAI <- function(q,n,w,j){
    time <- runif((3*q+n),w[j],w[j+1])
    time <- sort(time)
    if(q!=0){
      Etime <- time[1:q]
      Atime <- time[(q+1):(2*q)]
      Itime <- time[(2*q+1):(3*q)] 
      et <- c(w[j],time[q])
      at <- c(time[q],time[2*q])
      it <- c(time[2*q],time[3*q])
    }else{
      Etime <- 0
      Atime <- 0
      Itime <- 0
      et <- 0
      at <- 0
      it <- 0
    }
    if(n==0){
      Rtime <-  0
      rt <- 0
    }else{
      Rtime <- time[(3*q+1):(3*q+n)] 
      rt <- c(w[j],time[3*q+n])
    }
    return(list(Etime=Etime, Atime=Atime, Itime=Itime, Rtime=Rtime,et=et,at=at,it=it,rt=rt))
  }
  
  genI <- function(Itime,beta,gamma,theta,it,j){
    newI0 <- -rexp(beta+gamma+theta)
    if(length(Itime)==1){IStime <- Itime}else{
      IStime <- sample(Itime[Itime > 0],1)  #not sample the initial inf.  
    }
    genItime <- runif(1,it[1],it[2])
    Itime[which(Itime == IStime)] <- genItime
    Itime <- sort(Itime)   #new time
    return(list(Itime=Itime,IStime=IStime,genItime=genItime))
  }
  
  genR <- function(Rtime,rt,j){
    if(length(Rtime)==1){RStime <- Rtime}else{
      RStime <- sample(Rtime,1)  
    }
    genRtime <- runif(1,rt[1],rt[2])
    Rtime[which(Rtime == RStime)] <- genRtime
    Rtime <- sort(Rtime)   #new time
    return(list(Rtime=Rtime,RStime=RStime,genRtime=genRtime))
  }
  
  genE <- function(Etime,et,j){
    if(length(Etime)==1){EStime <- Etime}else{
      EStime <- sample(Etime,1)  
    }
    genEtime <- runif(1,et[1],et[2])
    Etime[which(Etime == EStime)] <- genEtime
    Etime <- sort(Etime)   #new time
    return(list(Etime=Etime,EStime=EStime,genEtime=genEtime))
  }
  
  genA <- function(Atime,at,j){
    if(length(Atime)==1){AStime <- Atime}else{
      AStime <- sample(Atime,1)  
    }
    genAtime <- runif(1,at[1],at[2])
    Atime[which(Atime == AStime)] <- genAtime
    Atime <- sort(Atime)   #new time
    return(list(Atime=Atime,AStime=AStime,genAtime=genAtime))
  }
  
  likj1 <- function(Itime,Etime,Atime,cm,ce,ca,N,beta,sigma,rho){
    #for E
    a <- cbind(Atime,rep(-1,sum(ca)))
    b <- cbind(Itime,rep(0,sum(cm)))
    c <- cbind(Etime,rep(1,sum(ce)))
    dt1 <- rbind(a,b,c)
    
    #for A
    d <- cbind(Etime,rep(0,sum(ce)))
    e <- cbind(Itime,rep(-1,sum(cm)))
    e[which(e[,1] < 0),2] <- 0
    f <- cbind(Atime,rep(1,sum(ca)))
    dt2 <- rbind(d,e,f)
    
    #for I
    g <- cbind(Etime,rep(0,sum(ce)))
    h <- cbind(Atime,rep(0,sum(ca)))
    hh <- cbind(Itime,rep(1,sum(cm)))
    dt3 <- rbind(g,h,hh)
    
    dt1 <- dt1[order(dt1[,1]),]
    dt2 <- dt2[order(dt2[,1]),]
    dt3 <- dt3[order(dt3[,1]),]        
    dt <- as.data.frame(cbind(dt1,dt2[,2],dt3[,2]))
    
    numE <- cumsum(dt[,2])
    numA <- cumsum(dt[,3])
    numI <- cumsum(dt[,4])
    numS <- N-numI-numE-numA
    
    numdat <- cbind(dt,numS,numE,numA,numI)
    
    tt <- numdat[,1]
    ss <- numdat[,5]
    ee <- numdat[,6]
    aa <- numdat[,7]
    ii <- numdat[,8]
    
    si1 <- ss*(ii+aa)
    l <- length(si1)
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    ee1 <- ee[-1]
    ee2 <- ee[-l]
    aa1 <- aa[-1]
    aa2 <- aa[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intE <- sum(((t1-t2)*(ee1+ee2))/2)
    intA <- sum(((t1-t2)*(aa1+aa2))/2)
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    
    f1 <- prod(beta*(ss/N)*(ii+aa))*prod(sigma*ee)*prod(rho*aa)
    f2 <- exp(-(beta*intSI+sigma*intE+rho*intA))
    lik <- f1*f2
    return(list(intSI=intSI,intE=intE,intI=intI,intA=intA,lik=lik))
  }
  
  likotherj <- function(Itime,Etime,Atime,Rtime,genItime,genEtime,genAtime,cm,ce,cn,ca,N,beta,sigma,gamma,rho){
    #for E
    a <- cbind(Atime,rep(-1,sum(ca)))
    b <- cbind(Itime,rep(0,sum(cm)))
    c <- cbind(Etime,rep(1,sum(ce)))
    d <- cbind(Rtime,rep(0,sum(cn)))
    dt1 <- rbind(a,b,c,d)
    
    #for A
    e <- cbind(Etime,rep(0,sum(ce)))
    f <- cbind(Itime,rep(-1,sum(cm)))
    f[which(f[,1] < 0),2] <- 0
    g <- cbind(Atime,rep(1,sum(ca)))
    h <- cbind(Rtime,rep(0,sum(cn)))
    dt2 <- rbind(e,f,g,h)
    
    #for I
    h1 <- cbind(Etime,rep(0,sum(ce)))
    h2 <- cbind(Atime,rep(0,sum(ca)))
    h3 <- cbind(Itime,rep(1,sum(cm)))
    h4 <- cbind(Rtime,rep(-1,sum(cn)))
    dt3 <- rbind(h1,h2,h3,h4)
    
    #for R
    r1 <- cbind(Etime,rep(0,sum(ce)))
    r2 <- cbind(Atime,rep(0,sum(ca)))
    r3 <- cbind(Itime,rep(0,sum(cm)))
    r4 <- cbind(Rtime,rep(1,sum(cn)))
    dt4 <- rbind(r1,r2,r3,r4)
    
    dt1 <- dt1[order(dt1[,1]),]
    dt2 <- dt2[order(dt2[,1]),]
    dt3 <- dt3[order(dt3[,1]),] 
    dt4 <- dt4[order(dt4[,1]),] 
    dt <- as.data.frame(cbind(dt1,dt2[,2],dt3[,2],dt4[,2]))
    
    numE <- cumsum(dt[,2])
    numA <- cumsum(dt[,3])
    numI <- cumsum(dt[,4])
    numR <- cumsum(dt[,5])
    numS <- N-numI-numE-numA-numR
    
    numdat <- cbind(dt,numS,numE,numA,numI)
    
    tt <- numdat[,1]
    ss <- numdat[,6]
    ee <- numdat[,7]
    aa <- numdat[,8]
    ii <- numdat[,9]
    
    si1 <- ss*(ii+aa)
    l <- length(si1)
    si2 <- si1[-1]
    si3 <- si1[-l]
    ii1 <- ii[-1]
    ii2 <- ii[-l]
    ee1 <- ee[-1]
    ee2 <- ee[-l]
    aa1 <- aa[-1]
    aa2 <- aa[-l]
    t1 <- tt[-1]
    t2 <- tt[-l]
    intSI <- sum(((t1-t2)*(si2+si3))/(2*N))
    intE <- sum(((t1-t2)*(ee1+ee2))/2)
    intA <- sum(((t1-t2)*(aa1+aa2))/2)
    intI <- sum(((t1-t2)*(ii1+ii2))/2)
    
    if(genItime != 0){
      ssu <- ss[which(numdat[,1]==genItime)]
      iiu <- ii[which(numdat[,1]==genItime)]
      eeu <- ee[which(numdat[,1]==genEtime)]
      aau <- aa[which(numdat[,1]==genAtime)]
      
      f1 <- prod(beta*(ssu/N)*(iiu+aau))*prod(sigma*eeu)*prod(rho*aau)*prod(gamma*iiu)
      f2 <- exp(-(beta*intSI+sigma*intE+rho*intA+gamma*intI))
      lik <- f1*f2
    }else{
      lik <- 0
    }
    
    return(list(intSI=intSI,intE=intE,intI=intI,intA=intA,lik=lik))
  }
  
  sbeta <- sbeta	#shape of beta   
  rbeta <- rbeta	#rate of beta   
  sgamma <- sgamma	#shape of gamma   
  rgamma <- rgamma	#rate of gamma  
  ssigma <- ssigma	 #shape of sigma 
  rsigma <- rsigma   #rate of sigma 
  srho <- srho	   #shape of rho
  rrho <- rrho     #rate of rho
  theta <- theta	#rate of I0   
  N <- N  #the population size
  burn <- 100
  iter <- burn+1000  #burn in 100 rounds
  nc <- length(data)
  w <- seq(0,nc*7,7)
  
  cm <- c()
  cn <- c()
  ce <- c()
  ca <- c()
  resbeta <- rep(0,nc)
  resgamma <- rep(0,nc)
  ressigma <- rep(0,nc)
  resrho <- rep(0,nc)
  resR0 <- rep(0,nc)
  SI <- rep(0,nc)
  cAtime <- c()
  cEtime <- c()
  cItime <- c()
  cRtime <- c()
  
  for (j in 1:nc){
    #number of I
    m <-  data[j]  #New this week
    cm <- c(cm,m)
    
    #number of R
    if(j==1){
      n <- 1
      cn <- c(cn,n)
    }else if(j==2){
      n <- data[j-1]-1
      cn <- c(cn,n)
    }else{
      n <- data[j-1]
      cn <- c(cn,n)
    }
    
    #number of E
    e <- m
    ce <- c(ce,e)
    
    #number of A
    a <- m
    ca <- c(ca,a)
    
    if(m==0 && n==0){
      Itimew <- 0
      Etimew <- 0
      Atimew <- 0
      Rtimew <- 0
    }else if(m==0 && n !=0){
      Itimew <- 0
      Etimew <- 0
      Atimew <- 0
      timewj <- priorEAI(m,n,w,j)
      Rtimew <- timewj$Rtime
    }else if(m !=0 && n==0){
      timewj <- priorEAI(m,n,w,j)
      Itimew <- timewj$Itime
      Etimew <- timewj$Etime
      Atimew <- timewj$Atime
      Rtimew <- 0
    }else{
      timewj <- priorEAI(m,n,w,j)
      Itimew <- timewj$Itime
      Etimew <- timewj$Etime
      Atimew <- timewj$Atime
      Rtimew <- timewj$Rtime
    }
    
    if(m > 1){
      et <- timewj$et
      at <- timewj$at
      it <- timewj$it
    }
    if(n > 1){
      rt <- timewj$rt
    }
    
    ## prior for sigma, beta, gamma, rho
    beta <- rgamma(1, sbeta, rbeta)
    gamma <- rgamma(1, sgamma, rgamma)
    sigma <- rgamma(1, ssigma, rsigma)
    rho <- rgamma(1, srho, rrho)
    
    ##################### Gibbs within Metropolis
    if(m==1 && j==1){
      Itimew <- Itimew[which(Itimew!=0)]
      cItime <- Itimew
      Rtimew <- Rtimew[which(Rtimew!=0)]
      cRtime <- Rtimew  
      Etimew <- Etimew[which(Etimew!=0)]
      cEtime <- Etimew  
      Atimew <- Atimew[which(Atimew!=0)]
      cAtime <- Atimew
      
      resbeta[j] <- 0
      resgamma[j] <- 0
      ressigma[j] <- 0
      resrho[j] <- 0
      resR0[j] <- 0
      SI <- "Cannot find since no new infectious"
    }else{
      
      p <- 1
      prob <- rep(0,iter)
      hatbeta <- rep(0,iter)
      hatgamma <- rep(0,iter)
      hatsigma <- rep(0,iter)
      hatrho <- rep(0,iter)
      hatR0 <- rep(0,iter)
      while (p <= iter){  ## loop MCMC
        oldItimew <- Itimew
        oldRtimew <- Rtimew
        oldEtimew <- Etimew
        oldAtimew <- Atimew
        
        if(m > 1){
          genIt <- genI(oldItimew,beta,gamma,theta,it,j)
          Itimew <- genIt$Itime  
          IStime <- genIt$IStime
          genItime <- genIt$genItime
        }else{
          if(m==0){
            Itimew <- 0
            genItime <- 0
            IStime <- 0
          }else{
            IStime <- oldItimew
            genItime <- oldItimew
          }
        }
        
        if(e > 1){
          genEt <- genE(oldEtimew,et,j)
          Etimew <- genEt$Etime 
          EStime <- genEt$EStime
          genEtime <- genEt$genEtime
        }else{
          if(e==0){
            Etimew <- 0
            genEtime <- 0
            EStime <- 0
          }else{
            EStime <- oldEtimew
            genEtime <- oldEtimew
          }
        }
        
        if(a > 1){
          genAt <- genA(oldAtimew,at,j)
          Atimew <- genAt$Atime
          AStime <- genAt$AStime
          genAtime <- genAt$genAtime
        }else{
          if(a==0){
            Atimew <- 0
            genAtime <- 0
            AStime <- 0
          }else{
            AStime <- oldAtimew
            genAtime <- oldAtimew
          }
        }
        
        if(n > 1){
          genRt <- genR(oldRtimew,rt,j)
          Rtimew <- genRt$Rtime 
        }else{
          if(n==0){
            Rtimew <- 0
          }else{
            RStime <- oldRtimew
            genRtime <- oldRtimew
          }
        }
        
        newItimew <- Itimew
        newEtimew <- Etimew
        newAtimew <- Atimew
        newRtimew <- Rtimew
        
        Itime <- c(cItime,newItimew)
        Itime <- Itime[which(Itime!=0)]
        oldItime <- c(cItime,oldItimew)
        oldItime <- oldItime[which(oldItime!=0)]
        
        Rtime <- c(cRtime,newRtimew)
        Rtime <- Rtime[which(Rtime!=0)]
        oldRtime <- c(cRtime,oldRtimew)
        oldRtime <- oldRtime[which(oldRtime!=0)]
        
        Etime <- c(cEtime,newEtimew)
        Etime <- Etime[which(Etime!=0)]
        oldEtime <- c(cEtime,oldEtimew)
        oldEtime <- oldEtime[which(oldEtime!=0)]
        
        Atime <- c(cAtime,newAtimew)
        Atime <- Atime[which(Atime!=0)]
        oldAtime <- c(cAtime,oldAtimew)
        oldAtime <- oldAtime[which(oldAtime!=0)]
        
        
        #likelihood for the whole data
        if(j==1){
          newlik <- likj1(Itime,Etime,Atime,cm,ce,ca,N,beta,sigma,rho)$lik
          oldlik <- likj1(oldItime,oldEtime,oldAtime,cm,ce,ca,N,beta,sigma,rho)$lik
        }else{
          newlik <- likotherj(Itime,Etime,Atime,Rtime,genItime,genEtime,genAtime,cm,ce,cn,ca,N,beta,sigma,gamma,rho)$lik
          oldlik <- likotherj(oldItime,oldEtime,oldAtime,oldRtime,IStime,EStime,AStime,cm,ce,cn,ca,N,beta,sigma,gamma,rho)$lik
        }
        
        ##acceptance probability
        likratio <- newlik/oldlik
        if(oldlik==0){likratio <- 1}
        cri <- c(1,likratio)
        prob[p] <- cri[which.min(cri)]
        u <- runif(1)
        if(u <= prob[p]){
          Itimew <- newItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          
          Etimew <- newEtimew 
          Etimew <- Etimew[which(Etimew!=0)]
          Etime <- c(cEtime,Etimew)
          
          Atimew <- newAtimew 
          Atimew <- Atimew[which(Atimew!=0)]
          Atime <- c(cAtime,Atimew) 
          
          Rtimew <- newRtimew 
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew) 
        }else{
          Itimew <- oldItimew
          Itimew <- Itimew[which(Itimew!=0)]
          Itime <- c(cItime,Itimew)
          
          Etimew <- oldEtimew
          Etimew <- Etimew[which(Etimew!=0)]
          Etime <- c(cEtime,Etimew)
          
          Atimew <- oldAtimew
          Atimew <- Atimew[which(Atimew!=0)]
          Atime <- c(cAtime,Atimew)
          
          Rtimew <- oldRtimew
          Rtimew <- Rtimew[which(Rtimew!=0)]
          Rtime <- c(cRtime,Rtimew) 
        }
        
        
        ##update beta, gamma sigma and rho
        sh_beta <- sbeta+sum(cm)
        sh_sigma <- ssigma+sum(ce)
        sh_gamma <- sgamma+sum(cn)
        sh_rho <- srho+sum(ca)
        
        if(j==1){
          lik <- likj1(Itime,Etime,Atime,cm,ce,ca,N,beta,sigma,rho)
          ra_beta <- rbeta+lik$intSI
          ra_sigma <- rsigma+lik$intE
          ra_gamma <- rgamma+lik$intI
          ra_rho <- rrho+lik$intA
        }else{
          lik <- likotherj(Itime,Etime,Atime,Rtime,genItime,genEtime,genAtime,cm,ce,cn,ca,N,beta,sigma,gamma,rho)
          ra_beta <- rbeta+lik$intSI
          ra_sigma <- rsigma+lik$intE
          ra_gamma <- rgamma+lik$intI
          ra_rho <- rrho+lik$intA
        }
        beta <- rgamma(1,sh_beta,ra_beta)
        sigma <- rgamma(1,sh_sigma,ra_sigma)
        gamma <- rgamma(1,sh_gamma,ra_gamma)
        rho <- rgamma(1,sh_rho,ra_rho)
        
        ##Estimate beta, sigma, gamma, rho 
        hatbeta[p] <- (sh_beta-1)/ra_beta
        hatsigma[p] <- (sh_sigma-1)/ra_sigma
        hatgamma[p] <- (sh_gamma-1)/ra_gamma
        hatrho[p] <- (sh_rho-1)/ra_rho
        if(hatgamma[p]==0 | hatrho[p]==0){
          hatR0[p] <- 0
        }else{
          hatR0[p] <- (hatbeta[p]/hatrho[p])+(hatbeta[p]/hatgamma[p])
        }
        p <- p+1
      }
      
      cItime <- Itime
      cEtime <- Etime
      cRtime <- Rtime  
      cAtime <- Atime  
      
      #average from 1,000 rounds (100 burn-in)
      resbeta[j] <- mean(hatbeta[burn:iter])
      ressigma[j] <- mean(hatsigma[burn:iter])
      resgamma[j] <- mean(hatgamma[burn:iter])
      resrho[j] <- mean(hatrho[burn:iter])
      resR0[j] <- mean(hatR0[burn:iter])
      SI[j] <- (1/resgamma[j])+(1/ressigma[j])
    }
  }
  return(list(R0=resR0, SI=SI))
}


