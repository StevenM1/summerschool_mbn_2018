####################### LBA stop signal

require(rtdists) 

##Simulate LBA trials
rlba.norm <- function (n,A,b,t0,mean_v,sd_v,st0=0,posdrift=TRUE,return.ttf=FALSE) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if (posdrift) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  make.r(drifts = drifts, b = b, A = A, n_v = n_v, t0 = t0, st0 = st0, n = n,
         return.ttf=return.ttf)
}


make.r <- function (drifts, b, A, n_v, t0, st0 = 0, n, return.ttf=FALSE) 
{
  drifts[drifts < 0] <- 0
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  ttf <- t0 + (b - starts)/drifts
  if (return.ttf) return(ttf)
  resp <- apply(ttf, 2, which.min)
  rt <- ttf[cbind(resp,1:n)]
  if (st0[1]>0) rt <- rt + runif(min = 0, max = st0[1], n = n)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad), "infinite RTs removed and less than", 
                  n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(rt = rt, response = resp)
}


dlba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like dlba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
  pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
    ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail), 
      ifelse(x < 0, 0, 1))
  
  dnormP <- function (x, mean = 0, sd = 1) 
    ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- rep(1, nn)
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmax(0, ((b[A_small]/t[A_small]^2) * 
                            dnorm1(b[A_small]/t[A_small], 
            mean_v[A_small], sd = sd_v[A_small]))/denom[A_small])
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A[!A_small])/zs
        out_o <- pmax(0, (mean_v[!A_small] * (pnorm1(chizu) - 
            pnorm1(chizumax)) + sd_v[!A_small] * (dnorm1(chizumax) - 
            dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
        out <- numeric(nn)
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        chizu <- chiminuszu/zs
        chizumax <- (chiminuszu - A)/zs
        return(pmax(0, (mean_v * (pnorm1(chizu) - pnorm1(chizumax)) + 
            sd_v * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
}


plba.norm.core <- function (t,A,b,mean_v,sd_v,posdrift=TRUE,robust=FALSE,nn) 
# like plba_norm_core but t0 dealt with outside (removed from t), and
# default added to nn
{
  
    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE) 
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail), 
        ifelse(x < 0, 0, 1))
  
    dnormP <- function (x, mean = 0, sd = 1) 
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
        pnorm1 <- pnormP
        dnorm1 <- dnormP
    }
    else {
        pnorm1 <- pnorm
        dnorm1 <- dnorm
    }
    if (posdrift) 
        denom <- pmax(pnorm1(mean_v/sd_v), 1e-10)
    else denom <- 1
    if (any(A < 1e-10)) {
        A_small <- A < 1e-10
        out_A <- pmin(1, pmax(0, (pnorm1(b[A_small]/t[A_small], 
            mean = mean_v[A_small], sd = sd_v[A_small], 
            lower.tail = FALSE))/denom[A_small]))
        zs <- t[!A_small] * sd_v[!A_small]
        zu <- t[!A_small] * mean_v[!A_small]
        chiminuszu <- b[!A_small] - zu
        xx <- chiminuszu - A[!A_small]
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        out_o <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]), 
            1)
        out <- numeric(length(mean_v))
        out[!A_small] <- out_o
        out[A_small] <- out_A
        return(out)
    }
    else {
        zs <- t * sd_v
        zu <- t * mean_v
        chiminuszu <- b - zu
        xx <- chiminuszu - A
        chizu <- chiminuszu/zs
        chizumax <- xx/zs
        tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
        tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
        return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
}


dlba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
#     A <- rep(A, length.out = nn)
#     b <- rep(b, length.out = nn)
#     mean_v <- rep(mean_v, length.out = nn)
#     sd_v <- rep(sd_v, length.out = nn)
#     if (any(b < A)) # b cannot be smaller than A! 
#         return(rep(0,nn))
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- dlba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], mean_v = mean_v[tpos], 
        sd_v = sd_v[tpos], posdrift = posdrift, robust = robust, nn = nn)
    out
}


plba.norm <- function (t, A, b, mean_v, sd_v, posdrift = TRUE, robust = FALSE) 
{
    nn <- length(t)
#     A <- rep(A, length.out = nn)
#     b <- rep(b, length.out = nn)
#     mean_v <- rep(mean_v, length.out = nn)
#     sd_v <- rep(sd_v, length.out = nn)
    tpos <- (t>0) & (b >= A)
    out <- numeric(length(t))
    out[tpos] <- plba.norm.core(t = t[tpos], A = A[tpos], b = b[tpos], mean_v = mean_v[tpos], 
        sd_v = sd_v[tpos], posdrift = posdrift, robust = robust, nn = nn)
    out
}


n1PDFfixedt0.norm=function(dt,A,b,mean_v,sd_v, 
                           posdrift=TRUE,robust = FALSE) 
# Generates defective PDF for responses on node= 1
# dt (decison time) is a matrix with length(mean_v) rows, one row for
# each accumulator to allow for different start times
{
  
  n_acc <- ifelse(is.null(dim(mean_v)),length(mean_v),dim(mean_v)[1])
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(mean_v)) mean_v <- matrix(rep(mean_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sd_v)) sd_v <- matrix(rep(sd_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)

  dt[1,] <- dlba.norm(dt[1,],A=A[1,],b=b[1,],mean_v=mean_v[1,],sd_v=sd_v[1,],
                      posdrift=posdrift,robust=robust)
  if (n_acc>1) for (i in 2:n_acc)
    dt[1,] <- dt[1,]*(1-plba.norm(dt[i,],
      A=A[i,],b=b[i,],mean_v=mean_v[i,],sd_v=sd_v[i,],
      posdrift=posdrift,robust=robust))
  dt[1,]
}

####################### LBA stop signal ----

  # NB1: no st0
  
  # NB2: TRIALS effect on B (threshold), and values < 0 set to 0 
  
  rLBAss <- function (n, v, sv, B, A, t0, t0sg, tf=0, gf=0, ts = 0,
                      SSD=Inf, TRIALS = NA, staircase=NA,posdrift=TRUE) 
    # Race among length(v) accumulators (if v a vector), 
    # or dim(v)[1] (if a matrix), first of which is a stop accumulator.
    # Acts the same as rlba.norm except NA returned for RT when winner = 1 
    # AND USES B PARAMETERIZATION
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
    
    if ( !is.matrix(v) ) v <- matrix(rep(v,n),nrow=n_acc)
    if ( !is.matrix(sv) ) sv <- matrix(rep(sv,n),nrow=n_acc)
    if ( !is.matrix(B) ) B <- matrix(rep(B,n),nrow=n_acc)
    if ( !is.matrix(A) ) A <- matrix(rep(A,n),nrow=n_acc)
    
    if ( !any(is.na(TRIALS)) ) {
      if (length(TRIALS)!=n)
        stop("TRIALS must have length n")
      B[-1,] <- B[-1,] + rep(ts*TRIALS,each=n_acc-1)
      B[B<0] <- 0
    }
    b <- B+A
    
    if ( gf > 0 ) # Setup for GO failure
      is.gf <- as.logical(rbinom(length(SSD),1,gf)) else
        is.gf <- logical(length(SSD))
    
    if ( all(!is.finite(SSD)) ) {              # ALL GO
      out <- rlba.norm(n,mean_v=v[-1,,drop=FALSE],sd_v=sv[-1,,drop=FALSE],
                       b=b[-1,,drop=FALSE],A=A[-1,,drop=FALSE],t0=t0,posdrift=posdrift)
      out$response <- out$response+1
    } else {                                   # SOME STOP
      if ( any(is.na(staircase)) ) {           # STOP fixed SSD
        # add SSD to stop accumulator
        t0S[1,] <- t0S[1,] + SSD 
        out <- rlba.norm(n,mean_v=v,sd_v=sv,b=b,A=A,t0=t0S,posdrift=posdrift)
        if ( tf>0 ) {
          is.tf <- logical(length(SSD))
          is.tf[is.finite(SSD)][as.logical(rbinom(sum(is.finite(SSD)),1,tf))] <- TRUE  
          if ( any(is.tf) ) { 
            out[is.tf,] <- rlba.norm(sum(is.tf),
                                     mean_v=v[-1,is.tf,drop=FALSE],sd_v=sv[-1,is.tf,drop=FALSE],
                                     b=b[-1,is.tf,drop=FALSE],A=A[-1,is.tf,drop=FALSE],
                                     t0=0,posdrift=posdrift)
            out[is.tf,"response"] <- out[is.tf,"response"]+1
          }
        }
      } else {                                 # STOP, staircase
        if ( !is.numeric(staircase) | length(staircase)!=1 )
          stop("Staircase must be a numeric vector of length 1 specifying the step.")
        SSDi <- SSD[is.finite(SSD)][1] # begining SSD 
        dt <- rlba.norm(n,mean_v=v,sd_v=sv,b=b,A=A,t0=t0S,posdrift=posdrift,return.ttf=TRUE)
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
        out <- data.frame(rt=dt[cbind(winner,1:n)],response=winner)
      }
      out$rt <- out$rt + t0 # Add t0 for go responses
    }
    names(out) <- c("RT","R")
    out[out$R==1,"RT"] <- NA
    if ( gf > 0 ) {
      out$RT[is.gf] <- NA
      out$R[is.gf] <- 1
    }
    if ( any(is.na(TRIALS)) ) cbind.data.frame(out,SSD=SSD) else
      cbind.data.frame(out,SSD=SSD,TRIALS=TRIALS)
  }
  
  
  
  
  n1PDF.LBAss <- function(rt,v,sv,B,A,t0,t0sg,tf=0,gf=0,ts=0,
                          SSD=Inf,TRIALS=NA,Si,posdrift=TRUE)
    # Same as n1PDFfixedt0.norm except SSD is either a scalar or vector of length(rt)
    # stop accumulator must have name "NR". SSD is subtracted stop accumulator time
    # and dt=NA done by integration.
    #
    # USES B PARAMETERIZATION
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
    
    stopfn <- function(t,vj,svj,bj,Aj,t0,t0sg,SSD,Si,posdrift=TRUE) 
    {
      # d = s-g = tD-t0(GO), then add SSDstart time to get 
      # start time go - start time stop, i.e., finishing time advantage for GO
      t0S <- rep(0,length(vj)) # Set relative to stop accumulator finish time
      t0S[-Si] <- t0S[-Si]+t0sg-t0+SSD   
      # subtracting min(t0S) keeps all times positive
      dt <- matrix(rep(t,each=length(vj)),nrow=length(vj))+t0S-min(t0S)
      i <- c(Si,c(1:length(vj))[-Si])
      n1PDFfixedt0.norm(dt[i,,drop=FALSE],mean_v=vj[i],sd_v=svj[i],b=bj[i],A=Aj[i],
                        posdrift=posdrift)
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
    if (!is.matrix(sv)) sv <- matrix(rep(sv,dim(rt)[2]),nrow=n_acc)
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
    b <- B+A
    
    if ( any(!is.stop) ) 
    {
      rt[Si,!is.stop] <- rt[Si,!is.stop] - t0sg - SSD[!is.stop] 
      rt[-Si,!is.stop] <- rt[-Si,!is.stop]-t0 
      if ( tf > 0 ) 
      {
        rt[1,!is.stop] <- (1-gf)*(
          tf*n1PDFfixedt0.norm(rt[-Si,!is.stop,drop=FALSE],
                               mean_v=v[-Si,!is.stop,drop=FALSE],sd_v=sv[-Si,!is.stop,drop=FALSE],
                               A=A[-Si,!is.stop,drop=FALSE],b=b[-Si,!is.stop,drop=FALSE],posdrift=posdrift) +
            (1-tf)*n1PDFfixedt0.norm(rt[,!is.stop,drop=FALSE],
                                     mean_v=v[,!is.stop,drop=FALSE],sd_v=sv[,!is.stop,drop=FALSE],
                                     b=b[,!is.stop,drop=FALSE],A=A[,!is.stop,drop=FALSE],posdrift=posdrift)
        )
      } else 
        rt[1,!is.stop] <- (1-gf)*n1PDFfixedt0.norm(dt=rt[,!is.stop,drop=FALSE],
                                                   mean_v=v[,!is.stop,drop=FALSE],sd_v=sv[,!is.stop,drop=FALSE],
                                                   b=b[,!is.stop,drop=FALSE],A=A[,!is.stop,drop=FALSE],posdrift=posdrift)
    }
    
    if ( any(is.stop) ) for (j in pj) {
      if ( !is.finite(SSD[j]) ) tmp <- 0 else
        tmp <- my.integrate(f=stopfn,lower=0,vj=v[,j],svj=sv[,j],
                            bj=b[,j],Aj=A[,j],t0=t0,t0sg=t0sg,SSD=SSD[j],Si=Si)
      rt[1,is.stop][p %in% p[j]] <- gf +(1-gf)*(1-tf)*tmp
    }
    rt[1,]
  }
  
  
  # VERY EXTENSIVE TESTING WITH Two different SSDs
  {
  # # ########### TWO ACCUMULATOR CASE
  # {
  # n=1e4
  # # n=10
  # v=c(.5,1); B=c(1,1); A=c(1,1); sv=c(1,1)
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
  #   sim.go <- rLBAss(n=n,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,SSD=SSD,TRIALS=TRIALS,ts=ts)
  #   is.in <- !is.na(sim.go$RT) # in case go failure
  #   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  # } else {TRIALS=NA;ts=0}
  # 
  # # Simulate stop trials
  # sim <- rLBAss(n=n,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,SSD=SSD,TRIALS=TRIALS,ts=ts)
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
  #   tmp <- n1PDF.LBAss(sim$RT[!is.na(sim$RT)],v=v[2:1],sv=sv[2:1],B=B[2:1],A=A[2:1],t0=t0,t0sg=t0sg,
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
  #   tmp <- n1PDF.LBAss(rep(NA,n),v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,
  #     SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
  #   print(mean(tmp[SSD==.1]))
  #   print(mean(tmp[SSD==1]))
  #   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
  #   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
  #   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  # } else {
  #   # Save simulated densities
  #   r1 <- c(2,1)
  #   d.r1 <- n1PDF.LBAss(rt=c(x1c,x2c),v=v[r1],sv=sv[r1],B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,
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
  #   print(n1PDF.LBAss(NA,v=v,sv=sv,A=A,B=B,t0=t0,t0sg=t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
  #   print(n1PDF.LBAss(NA,v=v,sv=sv,A=A,B=B,t0=t0,t0sg=t0sg,SSD=1,Si=1,tf=tf,gf=gf))
  # }
  # 
  # }
  # 
  # ########### THREE ACCUMULATOR CASE
  # {
  # n=1e4
  # v=c(.5,1,.5); B=c(1,1,1); A=c(1,1,1); sv=c(1,1,1)
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
  #   sim.go <- rLBAss(n=n,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts)
  #   is.in <- !is.na(sim.go$RT) # in case go failure
  #   plot(smooth.spline(TRIALS[is.in],sim.go$RT[is.in]),ylab="Smooth",xlab="TRIALS",type="l")
  # } else {TRIALS=NA;ts=0}
  # 
  # # Simulate stop trials
  # sim <- rLBAss(n=n,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,tf=tf,gf=gf,TRIALS=TRIALS,ts=ts,SSD=SSD)
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
  #   d.r1 <- n1PDF.LBAss(sim$RT[is.in1],v=v[r1],sv=sv[r1],ts=ts,TRIALS=TRIALS[is.in1],
  #           B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,SSD=SSD[is.in1],Si=2,tf=tf,gf=gf)
  #   r2 <- c(3,1,2)
  #   is.in2 <- !is.na(sim$RT) & sim$R==3
  #   d.r2 <- n1PDF.LBAss(sim$RT[is.in2],v=v[r2],sv=sv[r2],ts=ts,TRIALS=TRIALS[is.in2],
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
  #   tmp <- n1PDF.LBAss(rep(NA,n),v=v,sv=sv,A=A,B=B,t0=t0,t0sg=t0sg,SSD=SSD,Si=1,tf=tf,gf=gf,ts=ts,TRIALS=TRIALS)
  #   print(mean(tmp[SSD==.1]))
  #   print(mean(tmp[SSD==1]))
  #   plot(TRIALS,tmp,pch=".",xlab="TRIALS",ylab="p(NA)",ylim=c(0,1))
  #   lines(smooth.spline(TRIALS[SSD==.1],as.numeric(is.na(sim$RT)[SSD==.1])),col="red")
  #   lines(smooth.spline(TRIALS[SSD==1],as.numeric(is.na(sim$RT)[SSD==1])),col="red")
  # } else {
  #   # Save simulated densities
  #   r1 <- c(2,1,3)
  #   d.r1 <- n1PDF.LBAss(rt=c(x1c,x2c),v=v[r1],sv=sv[r1],B=B[r1],A=A[r1],t0=t0,t0sg=t0sg,
  #     SSD=c(rep(.1,length(x1c)),rep(1,length(x2c))),Si=2,tf=tf,gf=gf)
  #   r2 <- c(3,1,2)
  #   d.r2 <- n1PDF.LBAss(rt=c(x1e,x2e),v=v[r2],sv=sv[r2],B=B[r2],A=A[r2],t0=t0,t0sg=t0sg,
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
  #   print(n1PDF.LBAss(NA,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,SSD=.1,Si=1,tf=tf,gf=gf))
  #   print(n1PDF.LBAss(NA,v=v,sv=sv,B=B,A=A,t0=t0,t0sg=t0sg,SSD=1,Si=1,tf=tf,gf=gf))
  # }
  # 
  }
  
