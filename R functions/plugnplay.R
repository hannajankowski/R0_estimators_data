
library(pomp)

#rm(list = ls())

############################ Plug-n-Play SIR #############################
plugnplay_SIR <- function(data, sbeta, rbeta, sgamma, rgamma, N){
  set.seed(12345)
  
  nc <- length(data)
  w <- nc
  
  rho <- 1   ## all cases are reported
  difft <- 1/7  ## do filtering for every day
  
  resbeta <- c()
  resgamma <- c()
  R0e <- c()
  SI <- c()
  
  for (j in 1:nc){
    data <- as.numeric(c(1,data[1:j]))
    date <- as.numeric(seq(0,j,1))
    infdat <- cbind(date,data)
    infdat <- as.data.frame(infdat)
    colnames(infdat) <- c("week","cases")
    
    sir.proc.sim <- function (S, I, R, H, beta, gamma, delta.t, ...) {
      N <- sum(S,I,R)
      foi <- beta*I/N
      trans <- c(reulermultinom(n=1,size=S,rate=foi,dt=delta.t),
                 reulermultinom(n=1,size=I,rate=gamma,dt=delta.t))
      S = S-trans[1]
      I = I+trans[1]-trans[2]
      R = R+trans[2]
      H = H+trans[2]
      c(S=S, I=I, R=R, H=H)
    }
    
    f <- function(t,S, I, R, beta, gamma){
      N <- sum(S, I, R)
      foi <- beta*I/N
      terms <- c(
        S*foi,
        I*gamma
      )
      terms <- unname(terms)
      c(
        S=-terms[1],
        I=terms[1]-terms[2],
        R=terms[2],
        H=terms[2]                 
      )
    }
    
    sir_dmeas <- function (cases, H, rho, log, ...) {
      dbinom(x=cases, size=H, prob=rho, log=log)
    }
    sir_rmeas <- function (H, rho, ...) {
      cases=rbinom(n=1, size=H, prob=rho)
    }
    
    a <- sbeta
    b <- rbeta
    c <- sgamma
    d <- rgamma
    
    flu.sir <- pomp(data= infdat,
                    times="week",
                    t0=0,
                    params=c(rho=1,gamma=mean(rgamma(n=10,sgamma,rgamma)),beta=mean(rgamma(n=10,sbeta,rbeta)),
                             S.0=N-1,R.0=0,I.0=1,H.0=0),
                    rmeasure=sir_rmeas,
                    dmeasure=sir_dmeas,
                    rprocess=euler(sir.proc.sim, delta.t=difft)
    )
    
    para <- coef(flu.sir)
    simpar <- c("beta","gamma")
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    fit <- mif2(flu.sir,Nmif=5, 
                rw.sd=rw.sd(beta=0.0001,gamma=0.0005),
                cooling.fraction.50=0.01,
                Np=1000)
    
    options(warn = oldw)
    
    resbeta[j] <- as.numeric(coef(fit)[3])
    resgamma[j] <- as.numeric(coef(fit)[2])  
    IP <- difft/(1-exp(-difft*resgamma[j]))
    R0e[j] <- resbeta[j]*IP
    SI[j] <- 1/resgamma[j]
  }
  return(list(R0=R0e, SI=SI))
}





############################ Plug-n-Play SEIR #############################
plugnplay_SEIR <- function(data, sbeta, rbeta, sgamma, rgamma, ssigma, rsigma, N){
  set.seed(12345)
  
  nc <- length(data)
  w <- nc
  
  rho <- 1   ## all cases are reported
  difft <- 1/7  ## do filtering for every day
  
  resbeta <- c()
  resgamma <- c()
  ressigma <- c()
  R0e <- c()
  SI <- c()
  
  for (j in 1:nc){
    data <- as.numeric(c(1,data[1:j]))
    date <- as.numeric(seq(0,j,1))
    infdat <- cbind(date,data)
    infdat <- as.data.frame(infdat)
    colnames(infdat) <- c("week","cases")
    
    sir.proc.sim <- function (S,E, I, R, H, beta,sigma, gamma, delta.t, ...) {
      N <- sum(S,E,I,R)
      foi <- beta*I/N
      trans <- c(reulermultinom(n=1,size=S,rate=foi,dt=delta.t),
                 reulermultinom(n=1,size=E,rate=sigma,dt=delta.t),
                 reulermultinom(n=1,size=I,rate=gamma,dt=delta.t))
      S = S-trans[1]
      E = E+trans[1]-trans[2]
      I = I+trans[2]-trans[3]
      R = R+trans[3]
      H = H+trans[3]
      c(S=S,E=E,I=I, R=R, H=H)
    }
    
    f <- function(t,S,E,I, R, beta,sigma, gamma){
      N <- sum(S,E, I, R)
      foi <- beta*I/N
      terms <- c(
        S*foi,
        E*sigma,
        I*gamma
      )
      terms <- unname(terms)
      c(
        S = -trans[1],
        E = trans[1]-trans[2],
        I = trans[2]-trans[3],
        R = trans[3],
        H = trans[3]                 
      )
    }
    
    sir_dmeas <- function (cases, H, rho, log, ...) {
      dbinom(x=cases, size=H, prob=rho, log=log)
    }
    sir_rmeas <- function (H, rho, ...) {
      cases=rbinom(n=1, size=H, prob=rho)
    }
    
    s.beta <- sbeta
    r.beta <- rbeta
    s.sigma <- ssigma
    r.sigma <- rsigma
    s.gamma <- sgamma
    r.gamma <- rgamma
    
    flu.sir <- pomp(data= infdat,
                    times="week",
                    t0=0,
                    params=c(rho=1,gamma=mean(rgamma(10,s.gamma,r.gamma)),beta=mean(rgamma(10,s.beta,r.beta)),
                             sigma=mean(rgamma(10,s.sigma,r.sigma)),
                             S.0=N-1,E.0=0, R.0=0,I.0=1,H.0=0),
                    rmeasure=sir_rmeas,
                    dmeasure=sir_dmeas,
                    rprocess=euler(sir.proc.sim, delta.t=difft)
    )
    
    para <- coef(flu.sir)
    simpar <- c("beta","sigma","gamma")
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    fit <- mif2(flu.sir,Nmif=5, 
                rw.sd=rw.sd(beta=0.0001,gamma=0.0005,sigma=0.0001),
                cooling.fraction.50=0.01,
                Np=1000)
    
    options(warn = oldw)
    
    resbeta[j] <- as.numeric(coef(fit)[3])
    resgamma[j] <- as.numeric(coef(fit)[2]) 
    ressigma[j] <- as.numeric(coef(fit)[4]) 
    IP <- difft/(1-exp(-difft*resgamma[j]))
    R0e[j] <- resbeta[j]*IP
    SI[j] <- (1/resgamma[j])+(1/ressigma[j])
  }
  return(list(R0=R0e, SI=SI))
}





############################ Plug-n-Play SEAIR #############################

plugnplay_SEAIR <- function(data, sbeta, rbeta, sgamma, rgamma, ssigma, rsigma, srhoA, rrhoA, N){
  set.seed(12345)
  
  nc <- length(data)
  w <- nc
  
  rho <- 1   ## all cases are reported
  difft <- 1/7  ## do filtering for every day
  
  resbeta <- c()
  resgamma <- c()
  ressigma <- c()
  resrhoA <- c()
  R0e <- c()
  SI <- c()
  
  for (j in 1:nc){
    data <- as.numeric(c(1,data[1:j]))
    date <- as.numeric(seq(0,j,1))
    infdat <- cbind(date,data)
    infdat <- as.data.frame(infdat)
    colnames(infdat) <- c("week","cases")
    
    sir.proc.sim <- function (S,E,A,I,R,H, beta,sigma,rhoA,gamma, delta.t, ...) {
      N <- sum(S,E,A,I,R)
      foi <- beta*I/N
      trans <- c(reulermultinom(n=1,size=S,rate=foi,dt=delta.t),
                 reulermultinom(n=1,size=E,rate=sigma,dt=delta.t),
                 reulermultinom(n=1,size=A,rate=rhoA,dt=delta.t),
                 reulermultinom(n=1,size=I,rate=gamma,dt=delta.t))
      S = S-trans[1]
      E = E+trans[1]-trans[2]
      A = A+trans[2]-trans[3]
      I = I+trans[3]-trans[4]
      R = R+trans[4]
      H = H+trans[4]
      c(S=S,E=E,A=A,I=I, R=R, H=H)
    }
    
    f <- function(t,S,E,A,I, R, beta,sigma,rhoA, gamma){
      N <- sum(S,E,A, I, R)
      foi <- beta*(I+A)/N
      terms <- c(
        S*foi,
        E*sigma,
        A*rhoA,
        I*gamma
      )
      terms <- unname(terms)
      c(
        S = -trans[1],
        E = trans[1]-trans[2],
        A = trans[2]-trans[3],
        I = trans[3]-trans[4],
        R = trans[4],
        H = trans[4]                 
      )
    }
    
    sir_dmeas <- function (cases, H, rho, log, ...) {
      dbinom(x=cases, size=H, prob=rho, log=log)
    }
    sir_rmeas <- function (H, rho, ...) {
      cases=rbinom(n=1, size=H, prob=rho)
    }
    
    s.beta <- sbeta
    r.beta <- rbeta
    s.sigma <- ssigma
    r.sigma <- rsigma
    s.gamma <- sgamma
    r.gamma <- rgamma
    s.rhoA <- srhoA
    r.rhoA <- rrhoA
    
    flu.sir <- pomp(data= infdat,
                    times="week",
                    t0=0,
                    params=c(rho=1,gamma=mean(rgamma(10,s.gamma,r.gamma)),beta=mean(rgamma(10,s.beta,r.beta)),
                             sigma=mean(rgamma(10,s.sigma,r.sigma)),rhoA=mean(rgamma(10,s.rhoA,r.rhoA)),
                             S.0=N-1,E.0=0,A.0=0, R.0=0,I.0=1,H.0=0),
                    rmeasure=sir_rmeas,
                    dmeasure=sir_dmeas,
                    rprocess=euler(sir.proc.sim, delta.t=difft)
    )
    
    para <- coef(flu.sir)
    simpar <- c("beta","sigma","rhoA","gamma")
    
    oldw <- getOption("warn")
    options(warn = -1)
    
    fit <- mif2(flu.sir,Nmif=5, 
                rw.sd=rw.sd(beta=0.0001,gamma=0.0005,sigma=0.0001,rhoA=0.0001),
                cooling.fraction.50=0.01,
                Np=1000)
    
    options(warn = oldw)
    
    resbeta[j] <- as.numeric(coef(fit)[3])
    resgamma[j] <- as.numeric(coef(fit)[2]) 
    ressigma[j] <- as.numeric(coef(fit)[4]) 
    resrhoA[j] <- as.numeric(coef(fit)[5])  
    IP <- difft/(1-exp(-difft*resgamma[j]))
    R0e[j] <- (resbeta[j]/resrhoA[j])+(resbeta[j]*IP)
    SI[j] <- (1/resgamma[j])+(1/ressigma[j])
  }
  return(list(R0=R0e, SI=SI))
}


