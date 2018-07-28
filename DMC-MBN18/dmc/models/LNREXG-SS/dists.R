####################### LNR choice, ExG stop signal

####################### EXG single accumulator ----


  # Modified from gamlss.dist to make cdf in nu > 0.05 * sigma case robust,
  # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
  # robust against small sigma cases. NOTE LATTER NOT VECTORIZED ASSUMES ALL mu
  # sigma AND tau THE SAME
  pexGAUS <- function (q, mu = 5, sigma = 1, nu = 1, lower.tail = TRUE, log.p = FALSE) 
  {
    if (sigma[1] <= 0) return(rep(NA,length(q)))
    if (nu[1] <= 0) return(rep(NA,length(q)))
    
    # if (sigma < 0.05*nu) 
    if (sigma[1] < 1e-4) 
      return(pexp(q-mu,1/nu,log=log.p,lower.tail=lower.tail)) # shfited exponential
    
    ly <- length(q)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    nu <- rep(nu, length = ly)
    index <- seq(along = q)
    z <- q - mu - ((sigma^2)/nu)
    cdf <- ifelse(is.finite(q), 
                  ifelse(nu > 0.05 * sigma, 
                         pnorm((q - mu)/sigma) - 
                           exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/nu))^2 - (mu^2) - 
                                                        2 * q * ((sigma^2)/nu))/(2 * sigma^2)), 
                         pnorm(q, mean = mu, sd = sigma)),
                  ifelse(q<0,0,1)   
    )
    if (lower.tail == TRUE) 
      cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE) 
      cdf <- cdf
    else cdf <- log(cdf)
    cdf
  }
  
  # gamlss.dist function, but returns NA for bad sigma or tau, and
  # robust against small sigma cases. NOTE LATTER NOT VECTORIZED ASSUMES ALL mu
  # sigma AND tau THE SAME
  dexGAUS <- function (x, mu = 5, sigma = 1, nu = 1, log = FALSE) 
  {
    if (sigma[1] <= 0) return(rep(NA,length(x)))
    if (nu[1] <= 0) return(rep(NA,length(x)))
    
    # if (sigma < 0.05*nu) 
    if (sigma[1] < 1e-4) 
      return(dexp(x-mu,1/nu,log=log)) # shfited exponential
    
    ly <- length(x)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    nu <- rep(nu, length = ly)
    z <- x - mu - ((sigma^2)/nu)
    logfy <- ifelse(nu > 0.05 * sigma, -log(nu) - (z + (sigma^2/(2 * 
                                                                   nu)))/nu + log(pnorm(z/sigma)), dnorm(x, mean = mu, sd = sigma, 
                                                                                                         log = TRUE))
    if (log == FALSE) 
      fy <- exp(logfy)
    else fy <- logfy
    fy
  }
  

####################### LNR n-choice ----

  rlnr <- function (n, meanlog, sdlog, t0, st0 = 0) 
    # Race among n_acc accumulators, mealnlog and sdlog can be n_acc length
    # vectors or n_acc x n matrices. t0 can be
    # a) a scalar, b) a vector of length number of accumulators or
    # c) a matrix with 1 row per accumulator, when start times differ on each trial
    # st0, range of non-decison time variability, must be a scalar, as the same
    # variability is assumed in a common encoding/production stage
    
  {
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    dt <- matrix(rlnorm(n = n*n_acc, meanlog = meanlog, sdlog = sdlog), 
                 nrow = n_acc) + t0
    winner <- apply(dt,2,which.min)    
    if (st0[1]==0) data.frame(RT=dt[cbind(winner,1:n)],R=winner) else
      data.frame(RT=dt[cbind(winner,1:n)]+runif(n,0,st0[1]),R=winner)
  }
  
  n1PDFfixedt0.lnr=function(dt,meanlog,sdlog) 
    # Generates defective PDF for responses first among n_acc accumulator at 
    # dt (decison time), a matrix with one row for each accumulator (allowing for
    # different start times per accumulator) 
  {
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
    if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,dim(dt)[2]),nrow=n_acc)
    # winner
    dt[1,] <- dlnorm(dt[1,],meanlog[1,],sdlog[1,])
    # loosers
    if (dim(meanlog)[1]>1) for (i in 2:dim(meanlog)[1])
      dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i,],sdlog[i,],lower.tail=FALSE)
    dt[1,]
  }


####################### Combine LNR choice and EXG stop signal ----
  
  # NB1: no st0
  
  # NB2: TRIALS effect on meanlog for lnr accumulators 
  
  n1PDFfixedt0.lnrexg=function(dt,meanlog,sdlog,tau) 
    # Generates defective PDF for responses first among n_acc accumulator at 
    # dt (decison time), a matrix with one row for each accumulator (allowing for
    # different start times per accumulator). Assumed EXG is n2.  
  {
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    # if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
    # if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,dim(dt)[2]),nrow=n_acc)
    # winner
    dt[1,] <- dlnorm(dt[1,],meanlog[1,],sdlog[1,])
    # loosers
    if ( dim(meanlog)[1]==2 ) dt[1,] <- dt[1,]*
      pexGAUS(q=dt[2,],mu=meanlog[2,],sigma=sdlog[2,],nu=tau,lower.tail=FALSE) else
      {
        dt[1,] <- dt[1,]*pexGAUS(q=dt[2,],meanlog=mu[2,],sigma=sdlog[2,],nu=tau,lower.tail=FALSE)
        if (dim(meanlog)[1]>2) for (i in 3:dim(meanlog)[1])
          dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i,],sdlog[i,],lower.tail=FALSE)
      }
    dt[1,]
  }

  
  n1PDFfixedt0.exglnr=function(dt,meanlog,sdlog,tau) 
    # Generates defective PDF for responses first among n_acc accumulator at 
    # dt (decison time), a matrix with one row for each accumulator (allowing for
    # different start times per accumulator) EXG is n1 
  {
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    # if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(dt)[2]),nrow=n_acc)
    # if (!is.matrix(sdlog))     sdlog <- matrix(rep(  sdlog,dim(dt)[2]),nrow=n_acc)
    # winner
    dt[1,] <- dexGAUS(dt[1,],meanlog[1,],sdlog[1,],tau) 
    # loosers
    if (dim(meanlog)[1]>1) for (i in 2:dim(meanlog)[1])
        dt[1,] <- dt[1,]*plnorm(dt[i,],meanlog[i,],sdlog[i,],lower.tail=FALSE)
    dt[1,]
  }

    
  rlnrexg <- function(n,meanlog,sdlog,t0,tau) 

  # Race among n_acc accumulators, first exg rest lnr.
  # mealnlog and sdlog can be n_acc length
  # vectors or n_acc x n matrices. t0 and tau can be
  # a) a scalar, b) a vector of length number of accumulators or
  # c) a matrix with 1 row per accumulator, when start times differ on each trial

  {
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])-1
    dt <- rbind(
      rnorm(n=length(tau), mean=meanlog[1,],sd=sdlog[1,]) + rexp(n=length(tau),rate=1/tau),
      matrix(rlnorm(n = n*n_acc, meanlog = meanlog[-1,], sdlog = sdlog[-1,]) ,nrow = n_acc)
    ) + t0
    winner <- apply(dt,2,which.min)    
    data.frame(RT=dt[cbind(winner,1:n)],R=winner) 
  }
    
 
 
# n=10; meanlog=c(.3,-1,0); sdlog=c(.05,1,1); tau=(.1); t0=.1; SSD=.1
# tf=0; gf=0; ts = 0; TRIALS = NA; staircase=NA
  rlnrexgss <- function (n, meanlog, sdlog, t0, tau, tf=0, gf=0, ts = 0,
                      SSD=Inf, TRIALS = NA, staircase=NA) 
    # Race among length(meanlog) accumulators (if meanlog a vector), 
    # or dim(meanlog)[1] (if a matrix), first of which is a stop accumulator.
    # Acts the same as rlnr except NA returned for RT when winner = 1. 
    # Optional SSD argument can be used to adjust start time for first
    # accumulator. SSD can be a scalar or vector length n; output has an SSD column 
    # For trials with winning first accumulator RT and R set to NA. 
    # tf = trigger failure probability, gf = go failure probability
    # If any !is.na in staircase runs a staircase
    # t0 is a scalar with the standard interpritaiton for GO accumulators. 
    # Tau is stop accumulator tau (note meanlog = mu and sdlog = sigma for stop accumulator)
  # ts = slope of slowing (speeding if negative) over TRIALS, meanlog - ts*TRIALS
  # This has a linear effect on mean and sd, but no effect on 1st (stop) accumulator
  
  {
    if ( length(SSD)==1 ) SSD <- rep(SSD,n)
    if ( any(is.na(SSD)) || length(SSD) != n )
      stop("SSD cannot have NAs and must be a scalar or same length as n")
    
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    
    if (!is.matrix(t0)) t0 <- matrix(rep(t0,n*n_acc),nrow=n_acc)
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,n),nrow=n_acc)
    if (!is.matrix(sdlog)) sdlog <- matrix(rep(sdlog,n),nrow=n_acc)
    if (!is.matrix(tau)) tau <- matrix(rep(tau,n),nrow=1)

    if ( !any(is.na(TRIALS)) ) {
      if ( length(TRIALS)!=n )
        stop("TRIALS must have length n")
      meanlog[-1,] <- meanlog[-1,] + rep(ts*TRIALS,each=n_acc-1)
    }
    
    if ( gf > 0 ) # Setup for GO failure
      is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
        is.gf <- logical(length(SSD))
    
    if ( all(!is.finite(SSD)) ) {              # ALL GO
      out <- rlnr(n,meanlog[-1,,drop=FALSE],sdlog[-1,,drop=FALSE],t0[-1,,drop=FALSE])
      out$R <- out$R+1
    } else {                                   # SOME STOP
      if ( any(is.na(staircase)) ) {           # STOP fixed SSD
        # add SSD to stop accumulator
        t0[1,] <- t0[1,] + SSD 
        out <- rlnrexg(n,meanlog,sdlog,t0,tau)
        if ( tf>0 ) {
          is.tf <- logical(length(SSD))
          is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE  
          if ( any(is.tf) ) { 
            out[is.tf,] <- rlnr(sum(is.tf),meanlog[-1,is.tf,drop=FALSE],
                                sdlog[-1,is.tf,drop=FALSE],t0=0)
            out[is.tf,"R"] <- out[is.tf,"R"]+1
          }
        }
      } else {                                 # STOP, staircase
        if ( !is.numeric(staircase) | length(staircase)!=1 )
          stop("Staircase must be a numeric vector of length 1 specifying the step.")
        SSDi <- SSD[is.finite(SSD)][1] # begining SSD 
        dt <- rbind(
          rnorm(n=length(tau), mean=meanlog[1,],sd=sdlog[1,]) + rexp(n=length(tau),rate=1/tau),
          matrix(rlnorm(n = n*(n_acc-1), meanlog = meanlog[-1,], sdlog = sdlog[-1,]) ,nrow = n_acc-1)
        ) + t0 
        # Setup
        winner <- numeric(n)
        for ( i in c(1:n) ) {
          if ( !is.finite(SSD[i]) )   # not staircase
            dt[1,i] <- dt[1,i] + SSD[i] else
              dt[1,i] <- dt[1,i] + SSDi # staircase
            if ( runif(1)<tf ) # Trigger failure
              winner[i] <- which.min(dt[2:n_acc,i])+1 else
                winner[i] <- which.min(dt[,i])
              if (is.gf[i]) winner[i] <- 1
              if ( is.finite(SSD[i]) ) { # update staircase
                SSD[i] <- SSDi
                if ( winner[i]==1 ) 
                  SSDi <- SSDi + staircase else
                    SSDi <- SSDi - staircase
                  if (SSDi<1e-10) SSDi <- 0
              }
        }
        out <- data.frame(RT=dt[cbind(winner,1:n)],R=winner)
      }
    }
    out[out$R==1,"RT"] <- NA
    if (gf > 0) {
      out$RT[is.gf] <- NA
      out$R[is.gf] <- 1
    }
    if ( any(is.na(TRIALS)) ) cbind.data.frame(out,SSD=SSD) else
      cbind.data.frame(out,SSD=SSD,TRIALS=TRIALS)
  }
  
  
  
  
  n1PDF.lnrexgss <- function(rt,meanlog,sdlog,t0,tau,tf=0,gf=0,ts=0,
                          SSD=Inf,TRIALS=NA,Si)
    # Same as n1PDF.lnr except SSD is either a scalar or vector of length(rt)
    # stop accumulator must have name "NR". SSD is subtracted stop accumulator time
    # and dt=NA done by integration.
    #
    # tf= probabiliy of trigger failure, where
    # L = trigger fail & respond + trigger and respond + trigger and no-response
    #   = tf*L(N-1)+(1-tf)[L(N)+p(S)],
    # L(N-1) = choice race likelihood (no stop accumulator), 
    # L(N) = full N unit race likelihood given they did respond, 
    # p(S) probability of stop winning
    #
  # gf = probabiliy of go failure. 
  # On go trials:   L = go fail (so no response) + go and L above 
  # L = gf + (1-gf)*[tf*L(N-1)+(1-tf)[L(N)+p(S)]]  or similarly
  # L =    [ p(non-response) ]    +           [ p(response) ] 
  #   = [ gf + (1-gf)(1-tf)p(S) ] + [ (1-gf){(tf*Ln(n-1) + (1-tf)*L(N))} ]
  #
  # NB:rt is NOT decision time, but rather full RT as t0 has to be passed
  #    in order to include properly in cases where RT is NA (i.e., sucessful stop)
 
  
  {
    
    stopfn <- function(t,meanlogj,sdlogj,t0,tau,SSD,Si) 
    {

      t0S <- rep(0,length(meanlogj)) # Set relative to stop accumulator finish time
      t0S[-Si] <- SSD # Finishing time advantage for GO   
      # subtracting min(t0S) keeps all times positive
      dt <- matrix(rep(t,each=length(meanlogj)),nrow=length(meanlogj))+t0S
      i <- c(Si,c(1:length(meanlogj))[-Si])
      n1PDFfixedt0.exglnr(dt[i,,drop=FALSE],meanlogj[i],sdlogj[i],tau)
    }
    
    # NOTE: t0 is not subtracted when making dt but passed to handle RT=NA case
    
    if ( length(SSD)==1 ) SSD <- rep(SSD,length(rt))
    if (length(SSD) != length(rt))
      stop("SSD must be a scalar or same length as rt")
    n_acc <- ifelse(is.null(dim(meanlog)),length(meanlog),dim(meanlog)[1])
    
    rt <- matrix(rep(rt,each=n_acc),nrow=n_acc)
    is.stop <- is.na(rt[1,])  
    
    if (!is.matrix(meanlog)) meanlog <- matrix(rep(meanlog,dim(rt)[2]),nrow=n_acc)
    if (!is.matrix(sdlog)) sdlog <- matrix(rep(sdlog,dim(rt)[2]),nrow=n_acc)
    tau <- rep(tau,length.out=dim(rt)[2])
    
    if ( any(is.na(TRIALS)) | ts == 0 ) { 
      p <- SSD[is.stop]
      pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique SSD
    } else {
      meanlog[-Si,] <- meanlog[-Si,] + ts*TRIALS
      p <- apply(
        rbind(meanlog[,is.stop,drop=FALSE],SSD[is.stop]),
        2,paste,collapse="")
      pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique p and SSD
    }
    
    if ( any(!is.stop) ) 
    {
      rt[Si,!is.stop] <- rt[Si,!is.stop] - SSD[!is.stop]
      rt[-Si,!is.stop] <- rt[-Si,!is.stop]-t0
      if ( tf > 0 ) 
      {
        rt[1,!is.stop] <- (1-gf)*(
          tf*n1PDFfixedt0.lnr(rt[-Si,!is.stop,drop=FALSE],
            meanlog[-Si,!is.stop,drop=FALSE],sdlog[-Si,!is.stop,drop=FALSE]) +
          (1-tf)*n1PDFfixedt0.lnrexg(dt=rt[,!is.stop,drop=FALSE],
            meanlog=meanlog[,!is.stop,drop=FALSE],sdlog=sdlog[,!is.stop,drop=FALSE],
            tau=tau[!is.stop])
        )
      } else 
        rt[1,!is.stop] <- (1-gf)*n1PDFfixedt0.lnrexg(rt[,!is.stop,drop=FALSE],
          meanlog[,!is.stop,drop=FALSE],sdlog[,!is.stop,drop=FALSE],
          tau=tau[!is.stop])
    }
    
    if ( any(is.stop) ) for (j in pj) {
      #     tmp <- ifelse(!is.finite(SSD[j]),0,
      #       try(integrate(f=stopfn,lower=0,upper=Inf,meanlogj=meanlog[,j],
      #         sdlogj=sdlog[,j],t0=t0,t0sg=t0sg,SSD=SSD[j],Si=Si,
      #         NoBallistic = FALSE)$value,silent=TRUE))
      #     if (!is.numeric(tmp)) tmp <- 0
      if ( !is.finite(SSD[j]) ) tmp <- 0 else
        tmp <- my.integrate(f=stopfn,lower=0,meanlogj=meanlog[,j],
                            sdlogj=sdlog[,j],t0=t0,tau=tau[j],SSD=SSD[j],Si=Si)
      rt[1,is.stop][p %in% p[j]] <- gf +(1-gf)*(1-tf)*tmp
    }
    rt[1,]
  }
  
  
#### VERY EXTENSIVE TESTING WITH Two different SSDs ----
  { 
  # ########### TWO ACCUMULATOR CASE

  n=1e5
  meanlog=c(.5,-1); sdlog=c(.05,1); tau=.1; t0=.2
  # SSD = rep(c(1e5,1e6),each=n/2)
  SSD = rep(c(1,10)/10,each=n/2)

  # Run one of the follwing two lines
  do.trials=FALSE
  # do.trials = TRUE # requires very differnet plotting check, can be SLOW!

  ### RUN ONE OF THE FOLLOWING FOUR LINES
  # Without trigger failure or go failure
  tf=0; gf=0
  # With trigger failure, no go failure
  tf=.1;gf=0
  # Without trigger failure, with go failure
  tf=0; gf=.1
  # With trigger failure and go failure
  tf=.1;gf=.1

  if (do.trials) {
    ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
    TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
    # Plot slowing in GO (usually nice and linear, up to smoothing overfitting)
    sim.go <- rlnrexgss(n=n,meanlog,sdlog,t0,tau,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
    is.in <- !is.na(sim.go$RT) # in case go failure
    plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  } else {TRIALS=NA;ts=0}

  # Simulate stop trials
  sim <- rlnrexgss(n=n,meanlog,sdlog,t0,tau,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)

  # Plot densities
  par(mfrow=c(1,2))
  dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
  dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
  x1c <- dns1$'2'$x; x2c <- dns2$'2'$x

  # Signal respond RT
  dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
  round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)

  if (do.trials) {
    tmp <- n1PDF.lnrexgss(sim$RT[!is.na(sim$RT)],meanlog[2:1],sdlog[2:1],t0,tau,
      ts=ts,TRIALS=TRIALS[!is.na(sim$RT)],SSD=SSD[!is.na(sim$RT)],Si=2,tf=tf,gf=gf)
    par(mfrow=c(1,2))
    # red=black?
    plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT")
    lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==.1],
     tmp[c(SSD==.1)[!is.na(sim$RT)]]),col="red")
    # red=black?
    plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT")
    lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==1],
      tmp[c(SSD==1)[!is.na(sim$RT)]]),col="red")
    print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
    tmp <- n1PDF.lnrss(rep(NA,n),meanlog,sdlog,t0,t0sg,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
    print(mean(tmp[SSD==.1]))
    print(mean(tmp[SSD==1]))
    plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
    lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
    lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  } else {
    # Save simulated densities
    r1 <- c(2,1)
    d.r1 <- n1PDF.lnrexgss(rt=c(x1c,x2c),meanlog=meanlog[r1],sdlog=sdlog[r1],t0,tau,
      SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
    # Plot simulated (black) and theoretical (red) densities
    par(mfrow=c(1,2))
    # red=black?
    plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
       ylim=c(0,max(dns1$'2'$y)))
    lines(x1c,d.r1[1:length(x1c)],col="red")
    # red=black?
    plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
       ylim=c(0,max(dns2$'2'$y)))
    lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")

    # p(Stop check)
    print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
    print(n1PDF.lnrexgss(NA,meanlog,sdlog,t0,tau,SSD=.1,Si=1,tf=tf,gf=gf))
    print(n1PDF.lnrexgss(NA,meanlog,sdlog,t0,tau,SSD=1,Si=1,tf=tf,gf=gf))
  }



  ########### THREE ACCUMULATOR CASE

  n=1e5
  meanlog=c(.5,-1,0); sdlog=c(.05,1,1); tau=.1; t0=.2
  SSD = rep(c(1,10)/10,each=n/2)

  do.trials=FALSE
  do.trials = TRUE # requires very differnet plotting check, can be SLOW!

  ### RUN ONE OF THE FOLLOWING FOUR LINES
  # Without trigger failure or go failure
  tf=0; gf=0
  # With trigger failure, no go failure
  tf=.1;gf=0
  # Without trigger failure, with go failure
  tf=0; gf=.1
  # With trigger failure and go failure
  tf=.1;gf=.1

  if (do.trials) {
    ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
    TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
    # Plot slowing in GO (usually nice and linear, up to smooting overfitting)
    sim.go <- rlnrexgss(n=n,meanlog,sdlog,t0,tau,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
    is.in <- !is.na(sim.go$RT) # in case go failure
    plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  } else {TRIALS=NA;ts=0}

  # Simulate stop trials
  sim <- rlnrexgss(n=n,meanlog,sdlog,t0,tau,SSD=SSD,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)

  par(mfrow=c(1,2))
  dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
  dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
  x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
  x1e <- dns1$'3'$x; x2e <- dns2$'3'$x

  # Signal respond RT
  dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
  round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)

  if (do.trials) {
    r1 <- c(2,1,3)
    is.in1 <- !is.na(sim$RT) & sim$R==2
    d.r1 <- n1PDF.lnrexgss(sim$RT[is.in1],meanlog[r1],ts=ts,TRIALS=TRIALS[is.in1],
              sdlog=sdlog[r1],t0=t0,tau=tau,SSD=SSD[is.in1],Si=2,tf=tf,gf=gf)
    r2 <- c(3,1,2)
    is.in2 <- !is.na(sim$RT) & sim$R==3
    d.r2 <- n1PDF.lnrexgss(sim$RT[is.in2],meanlog[r2],ts=ts,TRIALS=TRIALS[is.in2],
             sdlog=sdlog[r2],t0=t0,tau=tau,SSD=SSD[is.in2],Si=2,tf=tf,gf=gf)
    par(mfrow=c(1,3))
    # red=black?
    plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT",type="l")
    lines(x1e,dns1$'3'$y,lty=2)
    lines(smooth.spline(sim$RT[is.in1 & sim$SSD==.1],
                        d.r1[c(sim$SSD==.1)[is.in1]]),col="red")
    lines(smooth.spline(sim$RT[is.in2 & sim$SSD==.1],d.r2[c(sim$SSD==.1)[is.in2]]),
          lty=2,col="red")
    # red=black?
    plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT",type="l")
    lines(x2e,dns2$'3'$y,lty=2)
    lines(smooth.spline(sim$RT[is.in1 & sim$SSD==1],
                        d.r1[c(sim$SSD==1)[is.in1]]),col="red")
    lines(smooth.spline(sim$RT[is.in2 & sim$SSD==1],
                        d.r2[c(sim$SSD==1)[is.in2]]),col="red",lty=2)

    print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
    tmp <- n1PDF.lnrexgss(rep(NA,n),meanlog,sdlog,t0,tau,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
    print(mean(tmp[SSD==.1]))
    print(mean(tmp[SSD==1]))
    plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
    lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
    lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  } else {
    # Save simulated densities
    r1 <- c(2,1,3)
    d.r1 <- n1PDF.lnrexgss(rt=c(x1c,x2c),meanlog[r1],sdlog[r1],t0,tau,
      SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
    r2 <- c(3,1,2)
    d.r2 <- n1PDF.lnrexgss(rt=c(x1e,x2e),meanlog[r2],sdlog[r2],t0,tau,
      SSD=c(rep(.1,length(x1e)),rep(1,length(x2e))),Si=2,tf=tf,gf=gf)
    # Plot simulated (black) and theoretical (red) densities
    par(mfrow=c(1,2))
    # red=black?
    plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
       ylim=c(0,max(dns1$'2'$y)))
    lines(x1c,d.r1[1:length(x1c)],col="red")
    lines(x1e,dns1$'3'$y,lty=2)
    lines(x1e,d.r2[1:length(x1e)],col="red",lty=2)
    # red=black?
    plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
       ylim=c(0,max(dns2$'2'$y)))
    lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
    lines(x2e,dns2$'3'$y,lty=2)
    lines(x2e,d.r2[(length(x2e)+1):(2*length(x2e))],col="red",lty=2)

    # p(Stop check)
    print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
    print(n1PDF.lnrexgss(NA,meanlog,sdlog,t0,tau,SSD=.1,Si=1,tf=tf,gf=gf))
    print(n1PDF.lngextss(NA,meanlog,sdlog,t0,tau,SSD=1,Si=1,tf=tf,gf=gf))
  }
  }
  
