####################### Wald stop signal


######################## Wald ----

# n-choice uniformly varying start point (0-A) Wald 
#    race, with t0, v, A, b (boundary) parameterixaiton

### Single accumulator model

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with: 
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to 
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.

# Following functions use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence 
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A/2 

rWald <- function(n,B,v,A)
  # random function for single acumulator
{
  
  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss
    
    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }
    
    flag <- l>tiny
    x <- rep(NA,times=n)
    
    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2
    
    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]
    
    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)
    
    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }
  
  # Act as if negative v never terminates, cluge to do single accumulator
  # case by passing negative v
  if (length(v)!=n) v <- rep(v,length.out=n)
  if (length(B)!=n) B <- rep(B,length.out=n)
  if (length(A)!=n) A <- rep(A,length.out=n)
  
  # Kluge to return -Inf for negative rates, so can implment one accumulator case
  out <- numeric(n)
  ok <- v>0  
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}


dWald <- function(t,v,B,A)
  # density for single accumulator
{
  
  digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # pdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns digt.0 if a<1e-10 
    
    digt.0 <- function(t,k=1,l=1) {
      # pdf of inverse gaussian at t with no k variability
      # much faster than statmod's dinvgauss funciton
      
      lambda <- k^2
      l0 <- l==0
      e <- numeric(length(t))
      if ( any(!l0) ) {
        mu <- k[!l0]/l[!l0]
        e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
      }
      if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
      x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
      x[t<=0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- digt.0(t=t[atiny],k=k[atiny],l=l[atiny])
    
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
        term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
        
        term.2a <- log(.5)+log(l[notltiny])
        term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
        term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.2d <- term.2b+term.2c
        term.2 <- exp(term.2a)*term.2d
        
        term.3 <- term.1+term.2
        term.4 <- log(term.3)-log(2)-log(a[notltiny])
        x[notltiny] <- exp(term.4)
      }
      
      if ( any(ltiny) ) {  # rate zero
        log.t <- log(t[ltiny])
        term.1 <- -.5*(log(2)+log(pi)+log.t)
        term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
        term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
        term.4 <- (exp(-term.2)-exp(-term.3))
        term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
        x[ltiny] <- exp(term.5)
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- digt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out
}


pWald <- function(t,v,B,A)
  # cumulative density for single accumulator
{
  pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
    # cdf of inverse gaussian at t with k +/- a/2 uniform variability
    # returns pigt.0 if a<=0
    
    pigt.0 <- function(t,k=1,l=1) {
      # cdf of inverse gaussian at t with no k variability
      # much faster than statmod's pinvgauss funciton
      
      mu <- k/l
      lambda <- k^2
      
      e <- exp(log(2*lambda) - log(mu))
      add <- sqrt(lambda/t) * (1 + t/mu)
      sub <- sqrt(lambda/t) * (1 - t/mu)
      
      p.1 <- 1 - pnorm(add)
      p.2 <- 1 - pnorm(sub)
      x <- exp(e + log(p.1)) + p.2
      
      x[t<0] <- 0
      x
    }
    
    options(warn=-1)
    if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
    if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
    if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
    
    tpos <- t<=0
    
    atiny <- a<=tiny & !tpos
    a[atiny] <- 0
    
    ltiny <- (l<=tiny) & !atiny & !tpos
    notltiny <- (l>tiny) & !atiny & !tpos
    l[l<=tiny] <- 0
    
    x <- numeric(length(t))
    
    # No threshold variability
    if ( any(atiny) )
      x[atiny] <- pigt.0(t[atiny],k[atiny],l[atiny])
    
    # Threshold variability
    if ( any(!atiny) ) {
      
      if ( any(notltiny) ) { # rate non-zero
        
        log.t <- log(t[notltiny])
        sqr.t <- sqrt(t[notltiny])
        
        term.1a <- .5*log.t-.5*log(2*pi)
        term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
        term.1 <- exp(term.1a)*(term.1b-term.1c)
        
        term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) + 
                         log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) + 
                         log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
        term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
        
        term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
        term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
        term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
        term.4 <- term.4c*term.4a + term.4d*term.4b
        
        x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
      }
      
      if ( any(ltiny) ) {  # rate zero
        sqr.t <- sqrt(t[ltiny])
        log.t <- log(t[ltiny])
        term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
        term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
        term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
        
        term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
        term.6 <- 1 + exp(term.6b) - exp(term.6a)
        
        x[ltiny] <- term.5 + term.6
      }
      
    }
    
    x[x<0 | is.nan(x) ] <- 0
    x
  }
  
  out <- numeric(length(t))
  ok <- v>0
  out[ok] <- pigt(t[ok],k=B[ok]+A[ok]/2,l=v[ok],a=A[ok]/2)
  out[!ok] <- 0
  out
  
}



### Race model

rWaldRace <- function(n,v,B,A,t0,gf=0,return.ttf=FALSE) 
  # random function for Wald race.
{
  B[B<0] <- 0 # Protection for negatives 
  A[A<0] <- 0
  n_v  <- ifelse(is.null(dim(v)), length(v), dim(v)[1])
  ttf <- matrix(t0 + rWald(n*n_v,B=B,v=v,A=A),nrow=n_v)
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  out <- data.frame(RT = ttf[cbind(resp,1:n)], R = apply(ttf, 2, which.min))
  
  if (gf[1] > 0) {
    is.gf <- as.logical(rbinom(dim(out)[1],1,gf))
    out$RT[is.gf] <- NA
    out$R[is.gf] <- 1
  }
  
  out
}


n1Wald <- function(dt,v,B,A,t0=0,gf=0)
  # Generates defective PDF for responses on node=1, dt (decison time) is a vector of times
{
  B[B<0] <- 0 # Protection for negatives 
  A[A<0] <- 0
  n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  dt <- dt-t0
  
  is.go <- !is.na(dt[1,])
  n.go <- sum(is.go)
  
  if (!is.matrix(v)) v <- matrix(rep(v,n.go),nrow=n_acc)
  if (!is.matrix(B)) B <- matrix(rep(B,n.go),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,n.go),nrow=n_acc)
  
  # Winner
  dt[1,is.go] <- (1-gf[1])*dWald(dt[1,is.go],A=A[1,],v=v[1,],B=B[1,])
  if (n_acc > 1) for (i in 2:n_acc)
    dt[1,is.go] <- dt[1,is.go]*(1-pWald(dt[i,is.go],A=A[i,],v=v[i,],B=B[i,]))
  
  dt[1,!is.go] <- gf[1]
  
  dt[1,]
}


  # NB1: no st0
  
  # NB2: TRIALS effect on B (threshold), and values < 0 set to 0 
  
  rWaldss <- function (n, v, B, t0, t0sg, A=0, tf=0, gf=0, ts = 0,
                       SSD=Inf, TRIALS = NA, staircase=NA) 
    # Race among length(v) accumulators (if v a vector), 
    # or dim(v)[1] (if a matrix), first of which is a stop accumulator.
    # Acts the same as rWald except NA returned for RT when winner = 1. 
    # Optional SSD argument can be used to adjust start time for first
    # accumulator. SSD can be a scalar or vector length n; output has an SSD column 
    # For trials with winning first accumulator RT and R set to NA. 
    # tf = trigger failure probability, gf = go failure probability
    # If any !is.na in staircase runs a staircase
    # t0 is a scalar with the standard interpritaiton for GO accumulators. 
    # Definition of t0sg (a non-negative scalar) depends on whether response 
    # production time is assumed to be ballistic or not.  
  # If it is ballistic t0sg = stop encoding - go encoding time + t0 
  # If if it is NOT ballistic, t0sg = stop encoding 
  # IN BOTH CASES t0sg < 0 IS NOT ALLOWED
  # ts = slope of slowing (speeding if negative) over TRIALS, meanlog - ts*TRIALS
  # This has a linear effect on mean and sd
  
  {
    if ( t0sg < 0 ) stop("t0sg cannot be less than zero")  
    
    if ( length(SSD)==1 ) SSD <- rep(SSD,n)
    if ( any(is.na(SSD)) || length(SSD) != n )
      stop("SSD cannot have NAs and must be a scalar or same length as n")
    
    n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
    
    # t0sg -> sg for stop accumulator, relative to 0 for go accumulators
    t0S <- matrix(rep(c(t0sg-t0[1],rep(0,n_acc-1)),length.out=n*n_acc),nrow=n_acc)
    
    if (!is.matrix(v)) v <- matrix(rep(v,n),nrow=n_acc)
    if (!is.matrix(B)) B <- matrix(rep(B,n),nrow=n_acc)
    if (!is.matrix(A)) A <- matrix(rep(A,n),nrow=n_acc)
    
    if ( !any(is.na(TRIALS)) ) {
      if (length(TRIALS)!=n)
        stop("TRIALS must have length n")
      B[-1,] <- B[-1,] + rep(ts*TRIALS,each=n_acc-1)
    }
    
    if ( gf > 0 ) # Setup for GO failure
      is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
        is.gf <- logical(length(SSD))
    
    if ( all(!is.finite(SSD)) ) {              # ALL GO
      out <- rWaldRace(n,v=v[-1,,drop=FALSE],B=B[-1,,drop=FALSE],A=A[-1,,drop=FALSE],t0=t0)
      out$R <- out$R+1
    } else {                                   # SOME STOP
      if ( any(is.na(staircase)) ) {           # STOP fixed SSD
        # add SSD to stop accumulator
        t0S[1,] <- t0S[1,] + SSD 
        out <- rWaldRace(n,v=v,B=B,A=A,t0=t0S)
        if ( tf>0 ) {
          is.tf <- logical(length(SSD))
          is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE  
          if ( any(is.tf) ) { 
            out[is.tf,] <- rWaldRace(sum(is.tf),v=v[-1,is.tf,drop=FALSE],
                                     B=B[-1,is.tf,drop=FALSE],A=A[-1,is.tf,drop=FALSE],t0=0)
            out[is.tf,"R"] <- out[is.tf,"R"]+1
          }
        }
      } else {                                 # STOP, staircase
        if ( !is.numeric(staircase) | length(staircase)!=1 )
          stop("Staircase must be a numeric vector of length 1 specifying the step.")
        SSDi <- SSD[is.finite(SSD)][1] # begining SSD 
        dt <- rWaldRace(n,v=v,B=B,A=A,t0=t0S,return.ttf=TRUE)
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
      out$RT <- out$RT + t0 # Add t0 for go responses
    }
    
    # print(out)
    
    out[out$R==1,"RT"] <- NA
    if ( gf > 0 ) {
      out$RT[is.gf] <- NA
      out$R[is.gf] <- 1
    }
    if ( any(is.na(TRIALS)) ) cbind.data.frame(out,SSD=SSD) else
      cbind.data.frame(out,SSD=SSD,TRIALS=TRIALS)
  }
  
  
####################### Wald stop signal ----
  
  
  n1PDF.Waldss <- function(rt,v,B,t0,t0sg,A=0,tf=0,gf=0,ts=0,
                           SSD=Inf,TRIALS=NA,Si)
    # Same as n1Wald except SSD is either a scalar or vector of length(rt)
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
  #
  # Definition of t0sg (a non-negative scalar) depends on whether response 
  # production time is assumed to be ballistic or not.  
  # If it is ballistic t0sg = stop encoding - go encoding time + t0 
  # If if it is NOT ballistic, t0sg = stop encoding 
  # IN BOTH CASES t0sg < 0 IS NOT ALLOWED (zero likelihood returned)
  
  
  {
    
    stopfn <- function(t,vj,Bj,Aj,t0,t0sg,SSD,Si) 
    {
      # d = s-g = tD-t0(GO), then add SSDstart time to get 
      # start time go - start time stop, i.e., finishing time advantage for GO
      t0S <- rep(0,length(vj)) # Set relative to stop accumulator finish time
      t0S[-Si] <- t0S[-Si]+t0sg-t0+SSD   
      # subtracting min(t0S) keeps all times positive
      dt <- matrix(rep(t,each=length(vj)),nrow=length(vj))+t0S-min(t0S)
      i <- c(Si,c(1:length(vj))[-Si])
      n1Wald(dt[i,,drop=FALSE],v=vj[i],B=Bj[i],A=Aj[i])
    }
    
    # NOTE: t0 is not subtracted when making dt but passed to handle RT=NA case
    
    # Bad t0sg
    if (t0sg<0) return(rep(0,length(rt)))
    
    if ( length(SSD)==1 ) SSD <- rep(SSD,length(rt))
    if (length(SSD) != length(rt))
      stop("SSD must be a scalar or same length as rt")
    n_acc <- ifelse(is.null(dim(v)),length(v),dim(v)[1])
    
    rt <- matrix(rep(rt,each=n_acc),nrow=n_acc)
    is.stop <- is.na(rt[1,])  
    
    if (!is.matrix(v)) v <- matrix(rep(v,dim(rt)[2]),nrow=n_acc)
    if (!is.matrix(B)) B <- matrix(rep(B,dim(rt)[2]),nrow=n_acc)
    if (!is.matrix(A)) A <- matrix(rep(A,dim(rt)[2]),nrow=n_acc)
    
    if ( any(is.na(TRIALS)) | ts == 0 ) { 
      p <- SSD[is.stop]
      pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique SSD
    } else {
      B[-Si,] <- B[-Si,] + ts*TRIALS
      p <- apply(
        rbind(B[,is.stop,drop=FALSE],SSD[is.stop]),
        2,paste,collapse="")
      pj <- c(1:sum(is.stop))[!duplicated(p)] # index of unique p and SSD
    }
    
    if ( any(!is.stop) ) 
    {
      rt[Si,!is.stop] <- rt[Si,!is.stop] - t0sg - SSD[!is.stop]
      rt[-Si,!is.stop] <- rt[-Si,!is.stop]-t0
      if ( tf > 0 ) 
      {
        rt[1,!is.stop] <- (1-gf)*(
          tf*n1Wald(rt[-Si,!is.stop,drop=FALSE],v=v[-Si,!is.stop,drop=FALSE],
                    A=A[-Si,!is.stop,drop=FALSE],B=B[-Si,!is.stop,drop=FALSE]) +
            (1-tf)*n1Wald(rt[,!is.stop,drop=FALSE],v=v[,!is.stop,drop=FALSE],
                          B=B[,!is.stop,drop=FALSE],A=A[,!is.stop,drop=FALSE])
        )
      } else 
        rt[1,!is.stop] <- (1-gf)*n1Wald(dt=rt[,!is.stop,drop=FALSE],
                                        v=v[,!is.stop,drop=FALSE],B=B[,!is.stop,drop=FALSE],
                                        A=A[,!is.stop,drop=FALSE])
    }
    
    if ( any(is.stop) ) for (j in pj) {
      if ( !is.finite(SSD[j]) ) tmp <- 0 else
        tmp <- my.integrate(f=stopfn,lower=0,vj=v[,j],
                            Bj=B[,j],Aj=A[,j],t0=t0,t0sg=t0sg,SSD=SSD[j],Si=Si)
      rt[1,is.stop][p %in% p[j]] <- gf +(1-gf)*(1-tf)*tmp
    }
    rt[1,]
  }
  
  
  # # VERY EXTENSIVE TESTING WITH Two different SSDs
  # {
  # # ########### TWO ACCUMULATOR CASE
  # {
  # n=1e4
  # # n=10
  # v=c(.5,1); B=c(1,1); A=c(1,1)
  # SSD = rep(c(1,10)/10,each=n/2)
  # 
  # # Run one of the follwing two lines
  # do.trials=FALSE
  # do.trials = TRUE # requires very differnet plotting check, can be SLOW!
  # 
  # #### RUN ONE OF THE FOLLOWING THREE LINES, all assume .2s go Ter
  # t0=.2; t0sg= 0  # minimum possible value of t0sg, stop encoding .2 less than go
  # t0=.2; t0sg=.2  # equal stop and go enconding times
  # t0=.2; t0sg=.4  # stop .2 slower than go
  # 
  # ### RUN ONE OF THE FOLLOWING FOUR LINES
  # # Without trigger failure or go failure
  # tf=0; gf=0
  # # With trigger failure, no go failure
  # tf=.1;gf=0
  # # Without trigger failure, with go failure
  # tf=0; gf=.1
  # # With trigger failure and go failure
  # tf=.1;gf=.1
  # 
  # if (do.trials) {
  #   ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
  #   TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
  #   # Plot slowing in GO (usually nice and linear, up to smooting overfitting)
  #   sim.go <- rWaldss(n=n,v=v,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,SSD=SSD,TRIALS=TRIALS,ts=ts)
  #   is.in <- !is.na(sim.go$RT) # in case go failure
  #   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  # } else {TRIALS=NA;ts=0}
  # 
  # # Simulate stop trials
  # sim <- rWaldss(n=n,v=v,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,SSD=SSD,TRIALS=TRIALS,ts=ts)
  # 
  # # Plot densities
  # par(mfrow=c(1,2))
  # dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
  # dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
  # x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
  # 
  # # Signal respond RT
  # dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
  # round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
  # 
  # if (do.trials) {
  #   tmp <- n1PDF.Waldss(sim$RT[!is.na(sim$RT)],v=v[2:1],B=B[2:1],A=A[2:1],t0=t0,t0sg=t0sg,
  #     ts=ts,TRIALS=TRIALS[!is.na(sim$RT)],SSD=SSD[!is.na(sim$RT)],Si=2,tf=tf,gf=gf)
  #   par(mfrow=c(1,2))
  #   # red=black?
  #   plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT")
  #   lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==.1],
  #    tmp[c(SSD==.1)[!is.na(sim$RT)]]),col="red")
  #   # red=black?
  #   plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT")
  #   lines(smooth.spline(sim$RT[!is.na(sim$RT) & SSD==1],
  #     tmp[c(SSD==1)[!is.na(sim$RT)]]),col="red")
  #   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
  #   tmp <- n1PDF.Waldss(rep(NA,n),v=v,B=B,A=A,t0=t0,t0sg=t0sg,
  #     SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
  #   print(mean(tmp[SSD==.1]))
  #   print(mean(tmp[SSD==1]))
  #   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
  #   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
  #   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  # } else {
  #   # Save simulated densities
  #   r1 <- c(2,1)
  #   d.r1 <- n1PDF.Waldss(rt=c(x1c,x2c),v=v[r1],B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,
  #     SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
  #   # Plot simulated (black) and theoretical (red) densities
  #   par(mfrow=c(1,2))
  #   # red=black?
  #   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns1$'2'$y)))
  #   lines(x1c,d.r1[1:length(x1c)],col="red")
  #   # red=black?
  #   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns2$'2'$y)))
  #   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
  # 
  #   # p(Stop check)
  #   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
  #   print(n1PDF.Waldss(NA,v=v,A=A,B=B,t0=t0,t0sg=t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
  #   print(n1PDF.Waldss(NA,v=v,A=A,B=B,t0=t0,t0sg=t0sg,SSD=1,Si=1,tf=tf,gf=gf))
  # }
  # 
  # }
  # 
  # ########### THREE ACCUMULATOR CASE
  # {
  # n=1e4
  # v=c(.5,1,.5); B=c(1,1,1); A=c(1,1,1)
  # SSD = rep(c(1,10)/10,each=n/2)
  # 
  # do.trials=FALSE
  # do.trials = TRUE # requires very differnet plotting check, can be SLOW!
  # 
  # #### RUN ONE OF THE FOLLOWING THREE LINES, all assume .2s go Ter
  # t0=.2; t0sg= 0  # minimum possible value of t0sg, stop encoding .2 less than go
  # t0=.2; t0sg=.2 # equal stop and go enconding times
  # t0=.2; t0sg=.4 # stop .2 slower than go
  # 
  # ### RUN ONE OF THE FOLLOWING FOUR LINES
  # # Without trigger failure or go failure
  # tf=0; gf=0
  # # With trigger failure, no go failure
  # tf=.1;gf=0
  # # Without trigger failure, with go failure
  # tf=0; gf=.1
  # # With trigger failure and go failure
  # tf=.1;gf=.1
  # 
  # if (do.trials) {
  #   ts=.5; TRIALS=log10(seq(1,10,length.out=n)) # 1..10 natural so 0-1 on log
  #   TRIALS <- as.vector(t(matrix(TRIALS,nrow=2))) # interleave SSDs
  #   # Plot slowing in GO (usually nice and linear, up to smooting overfitting)
  #   sim.go <- rWaldss(n=n,v=v,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
  #   is.in <- !is.na(sim.go$RT) # in case go failure
  #   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  # } else {TRIALS=NA;ts=0}
  # 
  # # Simulate stop trials
  # sim <- rWaldss(n=n,v=v,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts,SSD=SSD)
  # 
  # 
  # par(mfrow=c(1,2))
  # dns1 <- plot.cell.density(sim[sim$SSD==.1,],xlim=c(0,7),save.density=TRUE,main="SSD=.1")
  # dns2 <- plot.cell.density(sim[sim$SSD!=.1,],xlim=c(0,7),save.density=TRUE,main="SSD=1")
  # x1c <- dns1$'2'$x; x2c <- dns2$'2'$x
  # x1e <- dns1$'3'$x; x2e <- dns2$'3'$x
  # 
  # # Signal respond RT
  # dat <- sim; dat <- dat[!is.na(dat$RT),]; dat$R <- factor(as.character(dat$R))
  # round(tapply(dat$RT,dat[,c("R","SSD")],mean),2)
  # 
  # if (do.trials) {
  #   r1 <- c(2,1,3)
  #   is.in1 <- !is.na(sim$RT) & sim$R==2
  #   d.r1 <- n1PDF.Waldss(sim$RT[is.in1],v=v[r1],ts=ts,TRIALS=TRIALS[is.in1],
  #           B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,SSD=SSD[is.in1],Si=2,tf=tf,gf=gf)
  #   r2 <- c(3,1,2)
  #   is.in2 <- !is.na(sim$RT) & sim$R==3
  #   d.r2 <- n1PDF.Waldss(sim$RT[is.in2],v=v[r2],ts=ts,TRIALS=TRIALS[is.in2],
  #                    B=B[r2],A=A[r2],t0,t0sg,SSD=SSD[is.in2],Si=2,tf=tf,gf=gf)
  #   par(mfrow=c(1,3))
  #   # red=black?
  #   plot(x1c,dns1$'2'$y,pch=".",main="SSD=.1",ylab="Density",xlab="RT",type="l")
  #   lines(x1e,dns1$'3'$y,lty=2)
  #   lines(smooth.spline(sim$RT[is.in1 & sim$SSD==.1],
  #                       d.r1[c(sim$SSD==.1)[is.in1]]),col="red")
  #   lines(smooth.spline(sim$RT[is.in2 & sim$SSD==.1],d.r2[c(sim$SSD==.1)[is.in2]]),
  #         lty=2,col="red")
  #   # red=black?
  #   plot(x2c,dns2$'2'$y,pch=".",main="SSD=1",ylab="Density",xlab="RT",type="l")
  #   lines(x2e,dns2$'3'$y,lty=2)
  #   lines(smooth.spline(sim$RT[is.in1 & sim$SSD==1],
  #                       d.r1[c(sim$SSD==1)[is.in1]]),col="red")
  #   lines(smooth.spline(sim$RT[is.in2 & sim$SSD==1],
  #                       d.r2[c(sim$SSD==1)[is.in2]]),col="red",lty=2)
  #   
  #   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
  #   tmp <- n1PDF.Waldss(rep(NA,n),v=v,A=A,B=B,t0=t0,t0sg=t0sg,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
  #   print(mean(tmp[SSD==.1]))
  #   print(mean(tmp[SSD==1]))
  #   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
  #   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
  #   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  # } else {
  #   # Save simulated densities
  #   r1 <- c(2,1,3)
  #   d.r1 <- n1PDF.Waldss(rt=c(x1c,x2c),v=v[r1],B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,
  #     SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
  #   r2 <- c(3,1,2)
  #   d.r2 <- n1PDF.Waldss(rt=c(x1e,x2e),v=v[r2],B=B[r2],A=A[r2],t0=t0,t0sg=t0sg,
  #     SSD=c(rep(.1,length(x1e)),rep(1,length(x2e))),Si=2,tf=tf,gf=gf)
  #   # Plot simulated (black) and theoretical (red) densities
  #   par(mfrow=c(1,2))
  #   # red=black?
  #   plot(x1c,dns1$'2'$y,type="l",main="SSD=.1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns1$'2'$y)))
  #   lines(x1c,d.r1[1:length(x1c)],col="red")
  #   lines(x1e,dns1$'3'$y,lty=2)
  #   lines(x1e,d.r2[1:length(x1e)],col="red",lty=2)
  #   # red=black?
  #   plot(x2c,dns2$'2'$y,type="l",main="SSD=1",ylab="Density",xlab="RT",
  #      ylim=c(0,max(dns2$'2'$y)))
  #   lines(x2c,d.r1[(length(x2c)+1):(2*length(x2c))],col="red")
  #   lines(x2e,dns2$'3'$y,lty=2)
  #   lines(x2e,d.r2[(length(x2e)+1):(2*length(x2e))],col="red",lty=2)
  # 
  #   # p(Stop check)
  #   print(tapply(is.na(sim$RT),sim$SSD,mean)) # empirical
  #   print(n1PDF.Waldss(NA,v=v,B=B,A=A,t0=t0,t0sg=t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
  #   print(n1PDF.Waldss(NA,v=v,B=B,A=A,t0=t0,t0sg=t0sg,SSD=1,Si=1,tf=tf,gf=gf))
  # }
  # }
