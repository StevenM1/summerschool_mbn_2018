require(rtdists) 

#####################  LBA with Dont Know response

##################### Standard LBA funcitons
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


#####################  LBA with Dont Know response ----


rlba.DK <- function (n,A,b,d,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = nominal, 2 = DK
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  drifts[drifts < 0] <- 0
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  ttf <- (b - starts)/drifts
  resp <- apply(ttf, 2, which.min)
  looser <- resp; looser[resp==1] <- 2; looser[resp==2] <- 1
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  le <- rt*drifts[cbind(looser,1:n)]+starts[cbind(looser,1:n)]
  resp2 <- rep(1,length(resp))
  resp2[le>d[looser]] <- 2 # DK response
  if ( !(st0[1]>0) ) rt <- rt + t0 else
    rt <- rt + t0 + runif(min = 0, max = st0[1], n = n)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad), "infinite RTs removed and less than", 
                  n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(rt = rt, response = resp, response2 = resp2)
}


# tmp <- rlba.DK(n=10,A=1,b=2,d=1,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)
# 
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.DK <- function(dt,A,b,d,mean_v,sd_v, 
                            r2="",posdrift=TRUE,robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(mean_v)) mean_v <- matrix(rep(mean_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sd_v)) sd_v <- matrix(rep(sd_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d)) d <- matrix(rep(d,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.dk <- r2 =="DK"
  
  #   D1b <- dlba.norm(dt[1,!is.dk],A=A[1,!is.dk],b=b[1,!is.dk],mean_v=mean_v[1,!is.dk],
  #       sd_v=sd_v[1,!is.dk],posdrift=posdrift,robust=robust)
  #   S2d <- (1-plba.norm(dt[2,!is.dk],A=A[2,!is.dk],b=d[2,!is.dk],
  #       mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  # 
  #   if ( any(!is.dk) ) dt[1,!is.dk] <- D1b*S2d
  # 
  #   D1b <- dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],mean_v=mean_v[1,is.dk],
  #       sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)
  #   D2b <- dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],mean_v=mean_v[2,is.dk],
  #       sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)
  # 
  #   S1d <- (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
  #       mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))
  #   S1b <- (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
  #       mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))
  #   S2b <- (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
  #       mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))
  #   S2d <- (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
  #       mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))
  #   
  #   if ( any( is.dk) ) dt[1, is.dk] <- D1b*(S2b-S2d) + D2b*(S1b-S1d)
  
  
  if ( any(!is.dk) ) dt[1,!is.dk] <- 
    dlba.norm(dt[1,!is.dk],A=A[1,!is.dk],b=b[1,!is.dk],mean_v=mean_v[1,!is.dk],
              sd_v=sd_v[1,!is.dk],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,!is.dk],A=A[2,!is.dk],b=d[2,!is.dk],
                 mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  if ( any(is.dk) ) dt[1,is.dk] <- 
    dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
              mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
                  mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
                    mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))) +
    dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
              mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
                  mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
                    mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}


n1PDF.DK <- function (t,r2,A,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.DK(dt,r2=r2,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,mean_v,sd_v,scale,shape) 
      n1PDFfixedt0.DK(tau,r2,A,b,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],
                           A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
# # Check
# n=1e5
# t0=.2; A=c(.5,.5); d=c(1.5,1.5); b=c(2,2); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.DK(n=n,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# sim$response[sim$response2==2] <- 0
# names(sim) <- c("RT","R","R2")
# sim$R <- factor(sim$R,labels=c("DK","r1","r2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.DK,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="")$value+
# integrate(n1PDFfixedt0.DK,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="")$value+
# integrate(n1PDFfixedt0.DK,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="DK")$value
# 
# 
# # n1PDF check correct
# # red=black?
# d1 <- n1PDF.DK(dns$r1$x,r2="",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2 <- n1PDF.DK(dns$r2$x,r2="",A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# dDK <- n1PDF.DK(dns$DK$x,r2="DK",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# plot(dns$r1$x,dns$r1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1,d2,dDK))))
# lines(dns$r2$x,dns$r2$y,lty=3)
# lines(dns$DK$x,dns$DK$y,lty=1)
# lines(dns$r1$x,d1,col="red",lty=2)
# lines(dns$r2$x,d2,col="red",lty=3)
# lines(dns$DK$x,dDK,col="red",lty=1)
}

#####################  LBA with 2 level confidence response ----

rlba.2C <- function (n,A,b,d,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = low confidence, 2 = high confidence
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  drifts[drifts < 0] <- 0
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  ttf <- (b - starts)/drifts
  resp <- apply(ttf, 2, which.min)
  looser <- resp; looser[resp==1] <- 2; looser[resp==2] <- 1
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  le <- rt*drifts[cbind(looser,1:n)]+starts[cbind(looser,1:n)]
  resp2 <- rep(2,length(resp))  # High confidence response
  resp2[le>d[looser]] <- 1      # Low confidence response response
  if ( !(st0[1]>0) ) rt <- rt + t0 else
    rt <- rt + t0 + runif(min = 0, max = st0[1], n = n)
  bad <- !is.finite(rt)
  if (any(bad)) {
    warning(paste(sum(bad), "infinite RTs removed and less than", 
                  n, "rts returned"))
    resp <- resp[!bad]
    rt <- rt[!bad]
  }
  data.frame(rt = rt, response = resp, response2 = resp2)
}


# tmp <- rlba.2C(n=10,A=1,b=2,d=1,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)
# 
# dt = c(.5,1,1.5,.5,1,1.5); r2 = c(rep("1",3),rep("2",3))
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.2C <- function(dt,A,b,d,mean_v,sd_v, 
                            r2="2",posdrift=TRUE,robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(mean_v)) mean_v <- matrix(rep(mean_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sd_v)) sd_v <- matrix(rep(sd_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d)) d <- matrix(rep(d,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  
  if ( any(!is.lo) ) dt[1,!is.lo] <- 
    dlba.norm(dt[1,!is.lo],A=A[1,!is.lo],b=b[1,!is.lo],mean_v=mean_v[1,!is.lo],
              sd_v=sd_v[1,!is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,!is.lo],A=A[2,!is.lo],b=d[2,!is.lo],
                 mean_v=mean_v[2,!is.lo],sd_v=sd_v[2,!is.lo],posdrift=posdrift,robust=robust))
  if ( any(is.lo) ) dt[1,is.lo] <- 
    dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
              mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d[2,is.lo],
                    mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))
  dt[1,]
}


n1PDF.2C <- function (t,r2,A,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.2C(dt,r2=r2,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,mean_v,sd_v,scale,shape) 
      n1PDFfixedt0.2C(tau,r2,A,b,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],
                           A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
# # Check
# n=1e5
# t0=.2; A=c(.5,.5); d=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
# t0=.2; A=c(1,1); d=c(1.5,1.5); b=c(2.5,2.5); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.2C(n=n,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# names(sim) <- c("RT","R","R2")
# is.0 <- sim$R==1 & sim$R2==2
# is.3 <- sim$R==2 & sim$R2==2
# sim$R[is.0] <- 0
# sim$R[is.3] <- 3
# sim$R <- factor(sim$R,labels=c("hi1","lo1","lo2","hi2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="2")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value
# 
# # n1PDF check correct
# # red=black?
# d1lo <- n1PDF.2C(dns$lo1$x,r2="1",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2lo <- n1PDF.2C(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1hi <- n1PDF.2C(dns$hi1$x,r2="2",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2hi <- n1PDF.2C(dns$hi2$x,r2="2",A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1lo,d1hi,d2lo,d2hi))))
# lines(dns$hi1$x,d1hi,col="red",lty=2)
# lines(dns$lo1$x,dns$lo1$y,lty=3)
# lines(dns$lo1$x,d1lo,col="red",lty=3)
# lines(dns$lo2$x,dns$lo2$y,lty=1)
# lines(dns$lo2$x,d2lo,col="red",lty=1)
# lines(dns$hi2$x,dns$hi2$y,lty=4)
# lines(dns$hi2$x,d2hi,col="red",lty=4)
}