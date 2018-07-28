require(rtdists) 


##################### Standard LBA functions ----
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

##################### Lnorm Rate ----

dlba.lnorm <- function (t, A, b, meanlog_v, sdlog_v, robust = FALSE) 
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
  out[tpos] <- dlba.lnorm.core(t = t[tpos], A = A[tpos], b = b[tpos], meanlog_v = meanlog_v[tpos], 
                              sdlog_v = sdlog_v[tpos], robust = robust, nn = nn)
  out
}


plba.lnorm <- function (t, A, b, meanlog_v, sdlog_v, robust = FALSE) 
{
  nn <- length(t)
  #     A <- rep(A, length.out = nn)
  #     b <- rep(b, length.out = nn)
  #     mean_v <- rep(mean_v, length.out = nn)
  #     sd_v <- rep(sd_v, length.out = nn)
  tpos <- (t>0) & (b >= A)
  out <- numeric(length(t))
  out[tpos] <- plba.lnorm.core(t = t[tpos], A = A[tpos], b = b[tpos], meanlog_v = meanlog_v[tpos], 
                              sdlog_v = sdlog_v[tpos], robust = robust, nn = nn)
  out
}

dlba.lnorm.core <- function (t,A,b,meanlog_v,sdlog_v,robust=FALSE,nn) 
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
  min <- (b - A)/t
  max <- b/t
  zlognorm <- (exp(meanlog_v + (sdlog_v^2)/2) * (pnorm1((log(max) - 
                                                           meanlog_v - (sdlog_v^2))/sdlog_v) - pnorm1((log(min) - 
                                                                                                         meanlog_v - (sdlog_v^2))/sdlog_v)))/(pnorm1((log(max) - 
                                                                                                                                                        meanlog_v)/sdlog_v) - pnorm1((log(min) - meanlog_v)/sdlog_v))
  Gmax <- plnorm(max, meanlog = meanlog_v, sdlog = sdlog_v)
  Gmin <- plnorm(min, meanlog = meanlog_v, sdlog = sdlog_v)
  u <- (pnorm1((log(max) - meanlog_v - (sdlog_v)^2)/sdlog_v) - 
          pnorm1((log(min) - meanlog_v - (sdlog_v)^2)/sdlog_v))
  v <- (pnorm1((log(max) - meanlog_v)/sdlog_v) - pnorm1((log(min) - 
                                                           meanlog_v)/sdlog_v))
  udash <- (((-1/(sdlog_v * t)) * dnorm1((log(b/t) - meanlog_v - 
                                             (sdlog_v)^2)/sdlog_v)) - ((-1/(sdlog_v * t)) * dnorm1((log((b - 
                                                                                                            A)/t) - meanlog_v - (sdlog_v)^2)/sdlog_v)))
  vdash <- (((-1/(sdlog_v * t)) * dnorm1((log(b/t) - meanlog_v)/sdlog_v)) - 
              ((-1/(sdlog_v * t)) * dnorm1((log((b - A)/t) - meanlog_v)/sdlog_v)))
  const <- exp(meanlog_v + ((sdlog_v)^2)/2)
  diffzlognorm <- ((udash * v - vdash * u)/(v^2)) * const
  term1 <- (Gmax - Gmin) * (zlognorm + (t * diffzlognorm))
  term2 <- ((-b/(t^2)) * dlnorm(b/t, meanlog = meanlog_v, 
                                 sdlog = sdlog_v)) * ((zlognorm * t) - b)
  term3 <- (b - A - (zlognorm * t)) * ((-(b - A)/(t^2)) * 
                                          dlnorm((b - A)/t, meanlog = meanlog_v, sdlog = sdlog_v))
  out.value <- ((term1 + term2 + term3)/A)
  out.value[!is.finite(out.value)] <- 0
  return(pmax(0, out.value))

}


plba.lnorm.core <- function (t,A,b,meanlog_v,sdlog_v,robust=FALSE,nn) 
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
  min <- (b - A)/t
  max <- b/t
  zlognorm <- (exp(meanlog_v + (sdlog_v^2)/2) * (pnorm1((log(max) - 
                                                           meanlog_v - (sdlog_v^2))/sdlog_v) - pnorm1((log(min) - 
                                                                                                         meanlog_v - (sdlog_v^2))/sdlog_v)))/(pnorm1((log(max) - 
                                                                                                                                                        meanlog_v)/sdlog_v) - pnorm1((log(min) - meanlog_v)/sdlog_v))
  term1 <- ((t * zlognorm) - b)/A
  term2 <- (b - A - (t * zlognorm))/A
  pmax <- plnorm(max, meanlog = meanlog_v, sdlog = sdlog_v)
  pmin <- plnorm(min, meanlog = meanlog_v, sdlog = sdlog_v)
  out.value <- (1 + pmax * term1 + pmin * term2)
  out.value[t == Inf] <- 1
  out.value[!is.finite(out.value)] <- 0
  return(pmin(pmax(0, out.value), 1))

  
}

#####################  LBA with Dont Know response (no p_flc) ----

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
  is.dk <- r2 =="2"
  
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
    tmpf <- function(tau,ti0,r2,A,b,d,mean_v,sd_v) 
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
# Check
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

#####################  LBA with Dont Know response Threshold var. (no p_flc) ----

rlba.DK_A <- function (n,A,A1,b,d,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = nominal, 2 = DK
{
  
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  drifts[drifts < 0] <- 0
  ends<-matrix(runif(min = b-A, max = b, n = n * n_v), nrow = n_v)
  dkthresh<-matrix(runif(min = d-A1, max = d, n = n * n_v), nrow = n_v)
  ttf <- ends/drifts
  ttd<- dkthresh/drifts
  resp <- apply(ttf, 2, which.min)
  loser<- 3-resp
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  resp2 <- rep(1,n)
  
  resp2[ttd[cbind(loser,1:n)]< ttf[cbind(resp,1:n)]] <- 2 #DK response
  
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

n1PDFfixedt0.DK_A <- function(dt,A,A1,b,d,mean_v,sd_v, 
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
  if (!is.matrix(A1)) A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  
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
    (1-plba.norm(dt[2,!is.dk],A=A1[2,!is.dk],b=d[2,!is.dk],
                 mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  if ( any(is.dk) ) dt[1,is.dk] <- 
    dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
              mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
                  mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.dk],A=A1[2,is.dk],b=d[2,is.dk],
                    mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))) +
    dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
              mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
                  mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.dk],A=A1[1,is.dk],b=d[1,is.dk],
                    mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}


n1PDF.DK_A <- function (t,r2,A,A1,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.DK_A(dt,r2=r2,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                             posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,A1,b,d,mean_v,sd_v) 
      n1PDFfixedt0.DK_A(tau,r2,A,b,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],
                           A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
#Check
# n=1e5
# t0=.2; A=c(.5,.5);A1=c(.1,.1); d=c(1.5,1.5); b=c(2,2); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.DK_A(n=n,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# sim$response[sim$response2==2] <- 0
# names(sim) <- c("RT","R","R2")
# sim$R <- factor(sim$R,labels=c("DK","r1","r2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.DK_A,lower=0,upper=Inf,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="")$value+
#   integrate(n1PDFfixedt0.DK_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="")$value+
#   integrate(n1PDFfixedt0.DK_A,lower=0,upper=Inf,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="DK")$value
# 
# 
# # n1PDF check correct
# # red=black?
# d1 <- n1PDF.DK_A(dns$r1$x,r2="",A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2 <- n1PDF.DK_A(dns$r2$x,r2="",A=A[2:1],A1=A1[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# dDK <- n1PDF.DK_A(dns$DK$x,r2="DK",A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# plot(dns$r1$x,dns$r1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1,d2,dDK))))
# lines(dns$r2$x,dns$r2$y,lty=3)
# lines(dns$DK$x,dns$DK$y,lty=1)
# lines(dns$r1$x,d1,col="red",lty=2)
# lines(dns$r2$x,d2,col="red",lty=3)
# lines(dns$DK$x,dDK,col="red",lty=1)
}

#####################  LBA with Dont Know TC response (no p_flc) ----

rlba.DKTC <- function (n,A,b,d,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
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
  ttd<- (d-starts)/drifts
  minb <- apply(ttf, 2, which.min)
  maxd<- apply(ttd,2,which.max)
  rt <- ttf[cbind(minb,1:n)] # actually decision time, add t0 later
  resp2 <- rep(1,n)
  resp<-minb
  
  is.dk<- ttd[cbind(maxd,1:n)]<rt
  rt[is.dk]<-ttd[cbind(maxd,1:n)][is.dk]
  resp[is.dk]<-maxd[is.dk]
  # looser evidence = decison time * looser drift
  resp2[is.dk] <- 2
  
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

n1PDFfixedt0.DKTC <- function(dt,A,b,d,mean_v,sd_v, 
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
    dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
              mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
                  mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
                    mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))) +
    dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
              mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
                  mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
                    mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}


n1PDF.DKTC <- function (t,r2,A,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.DKTC(dt,r2=r2,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                             posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,mean_v,sd_v) 
      n1PDFfixedt0.DKTC(tau,r2,A,b,mean_v,sd_v)/st0
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
# ## Check
# n=1e5
# t0=.2; A=c(.5,.5); d=c(1.5,1.5); b=c(2,2); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.DKTC(n=n,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# sim$response[sim$response2==2] <- 0
# names(sim) <- c("RT","R","R2")
# sim$R <- factor(sim$R,labels=c("DK","r1","r2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.DKTC,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="")$value+
#   integrate(n1PDFfixedt0.DKTC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="")$value+
#   integrate(n1PDFfixedt0.DKTC,lower=0,upper=Inf,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="DK")$value
# 
# 
# # n1PDF check correct
# # red=black?
# d1 <- n1PDF.DKTC(dns$r1$x,r2="",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2 <- n1PDF.DKTC(dns$r2$x,r2="",A=A[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# dDK <- n1PDF.DKTC(dns$DK$x,r2="DK",A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# plot(dns$r1$x,dns$r1$y,type="l",ylab="density",x)
# lines(dns$r1$x,d1,col="red",lty=2)
# lines(dns$r2$x,d2,col="red",lty=3)lab="RT",lty=2,ylim=c(0,max(c(d1,d2,dDK))))
# lines(dns$r2$x,dns$r2$y,lty=3)
# lines(dns$DK$x,dns$DK$y,lty=1
# lines(dns$DK$x,dDK,col="red",lty=1)
}

#####################  LBA with Dont Know response TC Threshold var. (no p_flc) ----

rlba.DKTC_A <- function (n,A,A1,b,d,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = nominal, 2 = DK
{
  
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  drifts[drifts < 0] <- 0
  ends<-matrix(runif(min = b-A, max = b, n = n * n_v), nrow = n_v)
  dkthresh<-matrix(runif(min = d-A1, max = d, n = n * n_v), nrow = n_v)
  ttf <- ends/drifts
  ttd<- dkthresh/drifts
  resp <- apply(ttf, 2, which.min)
  maxd<-apply(ttd, 2, which.max)
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  resp2 <- rep(1,n)
  
  is.dk<-ttd[cbind(maxd,1:n)]< ttf[cbind(resp,1:n)]
  
  resp2[is.dk] <- 2 #DK response
  resp[is.dk]<-maxd[is.dk]
  rt[is.dk]<-ttd[cbind(maxd,1:n)][is.dk]
  
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

n1PDFfixedt0.DKTC_A <- function(dt,A,A1,b,d,mean_v,sd_v, 
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
  if (!is.matrix(A1)) A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  
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
    (1-plba.norm(dt[2,!is.dk],A=A1[2,!is.dk],b=d[2,!is.dk],
                 mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  if ( any(is.dk) ) dt[1,is.dk] <- 
    dlba.norm(dt[1,is.dk],A=A1[1,is.dk],b=d[1,is.dk],
              mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
                  mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.dk],A=A1[2,is.dk],b=d[2,is.dk],
                    mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))) +
    dlba.norm(dt[2,is.dk],A=A1[2,is.dk],b=d[2,is.dk],
              mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
                  mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.dk],A=A1[1,is.dk],b=d[1,is.dk],
                    mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}


n1PDF.DKTC_A <- function (t,r2,A,A1,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.DKTC_A(dt,r2=r2,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                               posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,A1,b,d,mean_v,sd_v) 
      n1PDFfixedt0.DKTC_A(tau,r2,A,b,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],
                           A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
#Check
# n=1e5
# t0=.2; A=c(.5,.5);A1=c(.1,.1); d=c(1.4,1.4); b=c(2,2); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.DKTC_A(n=n,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# sim$response[sim$response2==2] <- 0
# names(sim) <- c("RT","R","R2")
# sim$R <- factor(sim$R,labels=c("DK","r1","r2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.DKTC_A,lower=0,upper=Inf,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="")$value+
#   integrate(n1PDFfixedt0.DKTC_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="")$value+
#   integrate(n1PDFfixedt0.DKTC_A,lower=0,upper=Inf,A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,r2="DK")$value
# 
# 
# # n1PDF check correct
# # red=black?
# d1 <- n1PDF.DKTC_A(dns$r1$x,r2="",A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2 <- n1PDF.DKTC_A(dns$r2$x,r2="",A=A[2:1],A1=A1[2:1],b=b[2:1],d=d[2:1],mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# dDK <- n1PDF.DKTC_A(dns$DK$x,r2="DK",A=A,A1=A1,b=b,d=d,mean_v=mean_v,sd_v=sd_v,t0=t0)
# plot(dns$r1$x,dns$r1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1,d2,dDK))))
# lines(dns$r2$x,dns$r2$y,lty=3)
# lines(dns$DK$x,dns$DK$y,lty=1)
# lines(dns$r1$x,d1,col="red",lty=2)
# lines(dns$r2$x,d2,col="red",lty=3)
# lines(dns$DK$x,dDK,col="red",lty=1)
}

#####################  LBA with 2 level confidence response ----

rlba.2C <- function (n,A,b,d,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
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
  ttd<- (d-starts)/drifts  #low confidence  threshold
  
  nflc<-rbinom(1,n,p_flc)
  if(nflc>0){
    
    resp<-rep(1,times=n)
    resp2<-rep(1,times=n)
    rt<- rep(0,times=n)
    looser<-rep(2,times=n)
    
    resp[1:nflc]<-apply(ttd[,1:nflc], 2, which.min)
    resp2[1:nflc]<-1  #fast low confidence
    rt[1:nflc]<-ttd[cbind(resp[1:nflc],1:nflc)]
  
  resp[-(1:nflc)] <- apply(ttf[,-(1:nflc)], 2, which.min)
  looser[-(1:nflc)] <- apply(ttf[,-(1:nflc)], 2, which.max)
  #looser[-(1:nflc)] <- -1*looser[-(1:nflc)]+3  #faster?
  
  rt[-(1:nflc)] <- ttf[cbind(resp[-(1:nflc)],(nflc+1):n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  resp2[-(1:nflc)] <- rep(2,length(n-nflc))  # High confidence response
  resp2[-(1:nflc)][ttd[cbind(looser[-(1:nflc)],(nflc+1):n)]<rt[-(1:nflc)]] <- 1      # Low confidence response response
  
  } else
  {
    resp <- apply(ttf, 2, which.min)
    resp2<-rep(2,times=n)
    rt<-ttf[cbind(resp,1:n)]
    
    
    looser <- -1*resp+3
    # looser evidence = decison time * looser drift
    resp2 <- rep(2,times=n)  # High confidence response
    resp2[ttd[cbind(looser,1:n)]<rt] <- 1      # Low confidence response response
    
    
    
  }
  
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

n1PDFfixedt0.2C <- function(dt,A,b,d,p_flc, mean_v,sd_v, 
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
    (1-p_flc)*dlba.norm(dt[1,!is.lo],A=A[1,!is.lo],b=b[1,!is.lo],mean_v=mean_v[1,!is.lo],
              sd_v=sd_v[1,!is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,!is.lo],A=A[2,!is.lo],b=d[2,!is.lo],
                 mean_v=mean_v[2,!is.lo],sd_v=sd_v[2,!is.lo],posdrift=posdrift,robust=robust))
  
  if ( any(is.lo) ) dt[1,is.lo] <-
    ((p_flc)*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=d[1,is.lo],
                     mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))+
    ((1-p_flc)*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
              mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d[2,is.lo],
               mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)-
      plba.norm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))

   
  dt[1,]
  
  ##Note for Andrew, aren't we doing the (1-plba)-(1-plba) 
  #is plba.norm the cdf of an accumulator finishing before time t,
  #or after a time t?
  #Unesasarily? 
}


n1PDF.2C <- function (t,r2,A,b,d,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is the low crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.2C(dt,r2=r2,A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,p_flc,mean_v,sd_v) 
      n1PDFfixedt0.2C(tau,r2,A,b,d,p_flc,mean_v,sd_v)/st0
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
# Check
# n=1e5
# t0=.2; A=c(.5,.5); d=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
# t0=.2; A=c(1,1); d=c(1.5,1.5); b=c(2.5,2.5);p_flc=.9; mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.2C(n=n,A=A,b=b,d=d,p_flc=p_flc, mean_v=mean_v,sd_v=sd_v,t0=t0)
# names(sim) <- c("RT","R","R2")
# is.0 <- sim$R==1 & sim$R2==2
# is.3 <- sim$R==2 & sim$R2==2
# sim$R[is.0] <- 0
# sim$R[is.3] <- 3
# sim$R <- factor(sim$R,labels=c("hi1","lo1","lo2","hi2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A,b=b,d=d,p_flc=p_flc, mean_v=mean_v,sd_v=sd_v,r2="2")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
# integrate(n1PDFfixedt0.2C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value
# 
# 
# # n1PDF check correct
# # red=black?
# d1lo <- n1PDF.2C(dns$lo1$x,r2="1",A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2lo <- n1PDF.2C(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1hi <- n1PDF.2C(dns$hi1$x,r2="2",A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2hi <- n1PDF.2C(dns$hi2$x,r2="2",A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1lo,d1hi,d2lo,d2hi))))
# lines(dns$hi1$x,d1hi,col="red",lty=2)
# lines(dns$lo1$x,dns$lo1$y,lty=3)
# lines(dns$lo1$x,d1lo,col="red",lty=3)
# lines(dns$lo2$x,dns$lo2$y,lty=1)
# lines(dns$lo2$x,d2lo,col="red",lty=1)
# lines(dns$hi2$x,dns$hi2$y,lty=4)
# lines(dns$hi2$x,d2hi,col="red",lty=4)
}

#####################  LBA with 3 level confidence response ----

rlba.3C <- function (n,A,b,d1,d2,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = low confidence, 2 = medium confidence, 3 = high confidence
{
#   if (any(b < A)) 
#     stop("b cannot be smaller than A!")
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
  resp2 <- rep(3,n)                # High confidence response
  resp2[le>d1[looser] & le<=d2[looser]] <- 2  # Medium confidence response response
  resp2[le>d2[looser]] <- 1                   # Low confidence response response
  
  nflc<-rbinom(1,n,p_flc)
  
  if(nflc>0)
  {
    ttd1 <- (d1 - starts)/drifts
    # ADDED DROP
    respf<- apply(ttd1[,1:nflc,drop=FALSE], 2, which.min)

    rt[1:nflc]<-ttd1[cbind(respf,1:nflc)]
    resp[1:nflc]<-respf
    resp2[1:nflc]<-1
  }
  
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


# tmp <- rlba.3C(n=10,A=1,b=3,d1=1,d2=2,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)

# dt = c(.5,1,1.5,.5,1,1.5); r2 = c(rep("1",3),rep("2",3))
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.3C <- function(dt,A,b,d1,d2,p_flc,mean_v,sd_v, 
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
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    dlba.norm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],mean_v=mean_v[1,is.high],
              sd_v=sd_v[1,is.high],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.high],A=A[2,is.high],b=d1[2,is.high],
                 mean_v=mean_v[2,is.high],sd_v=sd_v[2,is.high],posdrift=posdrift,robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.med) ) dt[1,is.med] <- 
    dlba.norm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
              mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.med],A=A[2,is.med],b=d1[2,is.med],
                 mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust)-
       plba.norm(dt[2,is.med],A=A[2,is.med],b=d2[2,is.med],
                  mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust))*
    (1-p_flc[1])
       
  if ( any(is.lo) ) dt[1,is.lo] <- 
    ((1-p_flc[1])*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
              mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d2[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)-
      plba.norm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))+
    (p_flc[1]*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=d1[1,is.lo],
                    mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))
       
  dt[1,]
}


n1PDF.3C <- function (t,r2,A,b,d1,d2,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.3C(dt,r2=r2,A=A,b=b,d1=d1,d2=d2,p_flc,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d1,d2,p_flc,mean_v,sd_v) 
      n1PDFfixedt0.3C(tau,r2,A,b,d1,d2,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    # ADDED p_flc
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,b=b,d1=d1,d2=d2,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

 {
# #Check
# n=1e6
# t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
# t0=.2; A=c(1,1); d1=c(1.25,1.25); d2=c(1.75,1.75);p_flc=.2; b=c(2.5,2.5); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.3C(n=n,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# names(sim) <- c("RT","R","R2")
# is.1 <- sim$R==1 & sim$R2==3
# is.2 <- sim$R==1 & sim$R2==2
# is.3 <- sim$R==1 & sim$R2==1
# is.4 <- sim$R==2 & sim$R2==1
# is.5 <- sim$R==2 & sim$R2==2
# is.6 <- sim$R==2 & sim$R2==3
# sim$R[is.1] <- 1
# sim$R[is.2] <- 2
# sim$R[is.3] <- 3
# sim$R[is.4] <- 4
# sim$R[is.5] <- 5
# sim$R[is.6] <- 6
# sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="3")$value+
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="2")$value+
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value+
# integrate(n1PDFfixedt0.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="3")$value
# 
# # n1PDF check correct
# # red=black?
# d1lo <- n1PDF.3C(dns$lo1$x,r2="1",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2lo <- n1PDF.3C(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1med <- n1PDF.3C(dns$med1$x,r2="2",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2med <- n1PDF.3C(dns$med2$x,r2="2",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1hi <- n1PDF.3C(dns$hi1$x,r2="3",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2hi <- n1PDF.3C(dns$hi2$x,r2="3",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# 
# plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#   ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
# lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
# lines(dns$med1$x,dns$med1$y,lty=2)
# lines(dns$med1$x,d1med,col="red",lty=2)
# 
# lines(dns$lo1$x,dns$lo1$y,lty=3)
# lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
# lines(dns$lo2$x,dns$lo2$y,lty=4)
# lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
# lines(dns$med2$x,dns$med2$y,lty=5)
# lines(dns$med2$x,d2med,col="red",lty=5)
# 
# lines(dns$hi2$x,dns$hi2$y,lty=6)
# lines(dns$hi2$x,d2hi,col="red",lty=6)
}

#####################  LBA with 2 confidence levels and threshold counting ----

rlba.2TC <- function (n,A,b,d,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
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
  
  #RT at minimum of first cross of b, or 2nd cross of d
  ttb <- (b - starts)/drifts  #time to b
  ttd<- (d-starts)/drifts   # time to d
  
  resphi <- apply(ttb, 2, which.min)
  resplo<- apply(ttd,2,which.min)
  
  #looserhi <- resphi;looserhi[resphi==1] <- 2; looserhi[resphi==2] <- 1
  #looserhi redundant
  
  looserlo <- resplo;looserlo[resplo==1] <- 2; looserlo[resplo==2] <- 1
  rthi <- ttb[cbind(resphi,1:n)] # actaully decision time, add t0 later
  rtlo  <- ttd[cbind(looserlo,1:n)] # actaully decision time, add t0 lat
  #rt lo finishes when the losing accumulator reaches d.
  
  ## Calculate finishing time for each possible count.  
  
  rt<- apply(rbind(rthi,rtlo), 2, min)
  resp<- resphi
  resp[rtlo<rthi]<-resplo[rtlo<rthi]
  
  resp2 <- rep(2,length(resp))  # High confidence response
  resp2[rt==rtlo]<-1  #lo confidence response
  
  #Possibly not the most efficient way to get here.
  #Is there an if (both(rtd)< min (rtb)) style logic that could be used? 
  
  nflc<-rbinom(1,n,p_flc)
  if(nflc>0)
  {
    
    rt[1:nflc]<-ttd[cbind(resplo[1:nflc],1:nflc)]
    resp[1:nflc]<-resplo[1:nflc]
    resp2[1:nflc]<-1
  }
  
  
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


n1PDFfixedt0.2TC <- function(dt,A,b,d,p_flc,mean_v,sd_v, 
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
    (1-p_flc[1])*dlba.norm(dt[1,!is.lo],A=A[1,!is.lo],b=b[1,!is.lo],mean_v=mean_v[1,!is.lo],
              sd_v=sd_v[1,!is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,!is.lo],A=A[2,!is.lo],b=d[2,!is.lo],
                 mean_v=mean_v[2,!is.lo],sd_v=sd_v[2,!is.lo],posdrift=posdrift,robust=robust))
  if ( any(is.lo) ) dt[1,is.lo] <- 
    (1-p_flc[1])*dlba.norm(dt[2,is.lo],A=A[2,is.lo],b=d[2,is.lo],
              mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[1,is.lo],A=A[1,is.lo],b=d[1,is.lo],
                 mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)-
    plba.norm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
                  mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust))+
    
    (p_flc[1])*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=d[1,is.lo],
                          mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d[2,is.lo],
               mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust))
       
  dt[1,]     
}


n1PDF.2TC <- function (t,r2,A,b,d,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is the low crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.2TC(dt,r2=r2,A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,p_flc=p_flc,mean_v,sd_v) 
      n1PDFfixedt0.2TC(tau,r2,A,b,d,mean_v,sd_v)/st0
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
  # Check
  # n=1e5
  # t0=.2; A=c(.5,.5); d=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
  # t0=.2; A=c(1,1); d=c(1.5,1.5); b=c(2.5,2.5);p_flc=.1; mean_v=c(2,1); sd_v=c(1,1)
  # sim <- rlba.2TC(n=n,A=A,b=b,d=d,p_flc=p_flc, mean_v=mean_v,sd_v=sd_v,t0=t0)
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
  # integrate(n1PDFfixedt0.2TC,lower=0,upper=Inf,A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="2")$value+
  # integrate(n1PDFfixedt0.2TC,lower=0,upper=Inf,A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
  # integrate(n1PDFfixedt0.2TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
  # integrate(n1PDFfixedt0.2TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value
  # 
  # # n1PDF check correct
  # # red=black?
  # d1lo <- n1PDF.2TC(dns$lo1$x,r2="1",A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
  # d2lo <- n1PDF.2TC(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
  # d1hi <- n1PDF.2TC(dns$hi1$x,r2="2",A=A,b=b,d=d,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
  # d2hi <- n1PDF.2TC(dns$hi2$x,r2="2",A=A[2:1],b=b[2:1],d=d[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
  # plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=2,ylim=c(0,max(c(d1lo,d1hi,d2lo,d2hi))))
  # lines(dns$hi1$x,d1hi,col="red",lty=2)
  # lines(dns$lo1$x,dns$lo1$y,lty=3)
  # lines(dns$lo1$x,d1lo,col="red",lty=3)
  # lines(dns$lo2$x,dns$lo2$y,lty=1)
  # lines(dns$lo2$x,d2lo,col="red",lty=1)
  # lines(dns$hi2$x,dns$hi2$y,lty=4)
  # lines(dns$hi2$x,d2hi,col="red",lty=4)
}

#####################  LBA with 3 confidence levels, threshold counting. ----

rlba.3TC <- function (n,A,b,d1,d2,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                          lower = 0), nrow = n_v) else 
    drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                      nrow = n_v)
drifts[drifts < 0] <- 0
starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)

ttb <- (b - starts)/drifts  #time to b
ttd1<- (d1-starts)/drifts   # time to d1
ttd2<- (d2-starts)/drifts   # time to d2
resphi <- apply(ttb, 2, which.min)
respme <- apply(ttd2,2,which.min) 
# resplo<- apply(ttd2,2,which.min)  
looserlo<- apply(ttd2,2,which.max)


rt<-ttb[cbind(resphi,1:n)]
resp<-resphi

resp2 <- rep(3,n)  # High confidence response
#Set everything to High confidence and work down.


#med conf logical
tmp<-(rt>ttd1[cbind(looserlo,1:n)])&
  (ttd1[cbind(looserlo,1:n)]>ttd2[cbind(respme,1:n)])

#med conf responses and RTs
resp2[tmp]<-2   
resp[tmp]<-respme[tmp]   
rt[tmp]<-ttd1[cbind(looserlo,1:n)][tmp]


#low confidence logical
tmp<-(rt>ttd1[cbind(looserlo,1:n)])&
  (ttd1[cbind(looserlo,1:n)]<ttd2[cbind(respme,1:n)])

#low confidence responses and RTs
resp2[tmp]<-1    
resp[tmp]<-respme[tmp]
rt[tmp]<-ttd2[cbind(respme,1:n)][tmp]


nflc<-rbinom(1,n,p_flc)
if(nflc>0)
{
  resplo<- apply(ttd2[,1:nflc,drop=FALSE],2,which.min)  
  rt[1:nflc]<-ttd1[cbind(resplo,1:nflc)]
  resp[1:nflc]<-resplo
  resp2[1:nflc]<-1
  
}

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

n1PDFfixedt0.3TC <- function(dt,A,b,d1,d2,p_flc,mean_v,sd_v, 
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
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    (1-p_flc[1])*dlba.norm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],mean_v=mean_v[1,is.high],
              sd_v=sd_v[1,is.high],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.high],A=A[2,is.high],b=d1[2,is.high],
                 mean_v=mean_v[2,is.high],sd_v=sd_v[2,is.high],posdrift=posdrift,robust=robust))
  if ( any(is.med) ) dt[1,is.med] <- 
    (1-p_flc[1])*dlba.norm(dt[2,is.med],A=A[2,is.med],b=d1[2,is.med],
              mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
                  mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.med],A=A[1,is.med],b=d2[1,is.med],
                    mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust)))
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    (1-p_flc[1])*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=d2[1,is.lo],
              mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
               mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)-
      plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d2[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust))+
    (p_flc[1])*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=d1[1,is.lo],
                       mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust))
  
  dt[1,]
}


n1PDF.3TC <- function (t,r2,A,b,d1,d2,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.3TC(dt,r2=r2,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d1,d2,p_flc,mean_v,sd_v,scale,shape) 
      n1PDFfixedt0.3TC(tau,r2,A,b,d1,d2,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,b=b,d1=d1,d2=d2,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}


{
# Check
#   n=1e5
# 
#   t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
#   t0=.2; A=c(1,1); d1=c(1.5,1.5); d2=c(1.75,1.75);p_flc=.2; b=c(2,2); mean_v=c(2,1); sd_v=c(1,1)
#   sim <- rlba.3TC(n=n,A=A,b=b,d1=d1,d2=d2,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   names(sim) <- c("RT","R","R2")
#   is.1 <- sim$R==1 & sim$R2==3
#   is.2 <- sim$R==1 & sim$R2==2
#   is.3 <- sim$R==1 & sim$R2==1
#   is.4 <- sim$R==2 & sim$R2==1
#   is.5 <- sim$R==2 & sim$R2==2
#   is.6 <- sim$R==2 & sim$R2==3
#   sim$R[is.1] <- 1
#   sim$R[is.2] <- 2
#   sim$R[is.3] <- 3
#   sim$R[is.4] <- 4
#   sim$R[is.5] <- 5
#   sim$R[is.6] <- 6
#   sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#   dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
#   tapply(sim$RT,sim$R,mean)
# 
#   # Sum to 1?
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="3")$value+
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="2")$value+
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value+
#   integrate(n1PDFfixedt0.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="3")$value
# 
#   # n1PDF check correct
#   # red=black?
#   d1lo <- n1PDF.3TC(dns$lo1$x,r2="1",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2lo <- n1PDF.3TC(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
#   d1med <- n1PDF.3TC(dns$med1$x,r2="2",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2med <- n1PDF.3TC(dns$med2$x,r2="2",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
#   d1hi <- n1PDF.3TC(dns$hi1$x,r2="3",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2hi <- n1PDF.3TC(dns$hi2$x,r2="3",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# 
#   plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#     ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#   lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#   lines(dns$med1$x,dns$med1$y,lty=2)
#   lines(dns$med1$x,d1med,col="red",lty=2)
# 
#   lines(dns$lo1$x,dns$lo1$y,lty=3)
#   lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#   lines(dns$lo2$x,dns$lo2$y,lty=4)
#   lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#   lines(dns$med2$x,dns$med2$y,lty=5)
#   lines(dns$med2$x,d2med,col="red",lty=5)
# 
#   lines(dns$hi2$x,dns$hi2$y,lty=6)
#   lines(dns$hi2$x,d2hi,col="red",lty=6)

}

#####################  UL 3TC ----
rlba.lnorm.3TC <- function (n,A,b,d1,d2,p_flc,t0,meanlog_v,sdlog_v,st0=0) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(meanlog_v)), length(meanlog_v), dim(meanlog_v)[1])
  drifts <- exp(matrix(rnorm(n = n * n_v, mean = meanlog_v, sd = sdlog_v), 
                        nrow = n_v))

    starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  
  ttb <- (b - starts)/drifts  #time to b
  ttd1<- (d1-starts)/drifts   # time to d1
  ttd2<- (d2-starts)/drifts   # time to d2
  resphi <- apply(ttb, 2, which.min)
  respme <- apply(ttd2,2,which.min) 
  # resplo<- apply(ttd2,2,which.min)  
  looserlo<- apply(ttd2,2,which.max)
  
  
  rt<-ttb[cbind(resphi,1:n)]
  resp<-resphi
  
  resp2 <- rep(3,n)  # High confidence response
  #Set everything to High confidence and work down.
  
  
  #med conf logical
  tmp<-(rt>ttd1[cbind(looserlo,1:n)])&
    (ttd1[cbind(looserlo,1:n)]>ttd2[cbind(respme,1:n)])
  
  #med conf responses and RTs
  resp2[tmp]<-2   
  resp[tmp]<-respme[tmp]   
  rt[tmp]<-ttd1[cbind(looserlo,1:n)][tmp]
  
  
  #low confidence logical
  tmp<-(rt>ttd1[cbind(looserlo,1:n)])&
    (ttd1[cbind(looserlo,1:n)]<ttd2[cbind(respme,1:n)])
  
  #low confidence responses and RTs
  resp2[tmp]<-1    
  resp[tmp]<-respme[tmp]
  rt[tmp]<-ttd2[cbind(respme,1:n)][tmp]
  
  
  nflc<-rbinom(1,n,p_flc)
  if(nflc>0)
  {
    resplo<- apply(ttd2[,1:nflc,drop=FALSE],2,which.min)  
    rt[1:nflc]<-ttd1[cbind(resplo,1:nflc)]
    resp[1:nflc]<-resplo
    resp2[1:nflc]<-1
    
  }
  
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


n1PDFfixedt0.lnorm.3TC <- function(dt,A,b,d1,d2,p_flc,meanlog_v,sdlog_v, 
                             r2="3",robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(meanlog_v)) meanlog_v <- matrix(rep(meanlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog_v)) sdlog_v <- matrix(rep(sdlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    dlba.lnorm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],meanlog_v=meanlog_v[1,is.high],
                           sdlog_v=sdlog_v[1,is.high],robust=robust)*
    (1-plba.lnorm(dt[2,is.high],A=A[2,is.high],b=d1[2,is.high],
                 meanlog_v=meanlog_v[2,is.high],sdlog_v=sdlog_v[2,is.high],robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.med) ) dt[1,is.med] <- 
    dlba.lnorm(dt[2,is.med],A=A[2,is.med],b=d1[2,is.med],
                           meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust)*
    ((1-plba.lnorm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
                  meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust))-
       (1-plba.lnorm(dt[1,is.med],A=A[1,is.med],b=d2[1,is.med],
                    meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust)))*
    (1-p_flc[1])
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    (dlba.lnorm(dt[1,is.lo],A=A[1,is.lo],b=d2[1,is.lo],
                           meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
    (plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
               meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)-
       plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=d2[2,is.lo],
                 meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust))*(1-p_flc[1]))+
    (p_flc[1])*dlba.lnorm(dt[1,is.lo],A=A[1,is.lo],b=d1[1,is.lo],
                         meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
    (1-plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
                 meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust))
  
  dt[1,]
}


n1PDF.lnorm.3TC <- function (t,r2,A,b,d1,d2,p_flc,meanlog_v,sdlog_v,t0,st0=0,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.lnorm.3TC(dt,r2=r2,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,
      meanlog_v=meanlog_v,sdlog_v=sdlog_v, robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d1,d2,p_flc,meanlog_v,sdlog_v) 
      n1PDFfixedt0.lnorm.3TC(tau,r2,A,b,d1,d2,meanlog_v,sdlog_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,b=b,d1=d1,d2=d2,meanlog_v=meanlog_v,sdlog_v=sdlog_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
#  #Check
#     n=1e5
# 
#     t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     t0=.2; A=c(1,1); d1=c(1.5,1.5); d2=c(1.75,1.75);p_flc=.2; b=c(2,2); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     sim <- rlba.lnorm.3TC(n=n,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     names(sim) <- c("RT","R","R2")
#     is.1 <- sim$R==1 & sim$R2==3
#     is.2 <- sim$R==1 & sim$R2==2
#     is.3 <- sim$R==1 & sim$R2==1
#     is.4 <- sim$R==2 & sim$R2==1
#     is.5 <- sim$R==2 & sim$R2==2
#     is.6 <- sim$R==2 & sim$R2==3
#     sim$R[is.1] <- 1
#     sim$R[is.2] <- 2
#     sim$R[is.3] <- 3
#     sim$R[is.4] <- 4
#     sim$R[is.5] <- 5
#     sim$R[is.6] <- 6
#     sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#     dns <- plot.cell.density(sim,xlim=c(0,2),save.density=TRUE)
#     tapply(sim$RT,sim$R,mean)
# 
#     # Sum to 1?
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="3")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="3")$value
# 
#     # n1PDF check correct
#     # red=black?
#     d1lo <- n1PDF.lnorm.3TC(dns$lo1$x,r2="1",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2lo <- n1PDF.lnorm.3TC(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1med <- n1PDF.lnorm.3TC(dns$med1$x,r2="2",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2med <- n1PDF.lnorm.3TC(dns$med2$x,r2="2",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1hi <- n1PDF.lnorm.3TC(dns$hi1$x,r2="3",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2hi <- n1PDF.lnorm.3TC(dns$hi2$x,r2="3",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
# 
#     plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#       ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#     lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#     lines(dns$med1$x,dns$med1$y,lty=2)
#     lines(dns$med1$x,d1med,col="red",lty=2)
# 
#     lines(dns$lo1$x,dns$lo1$y,lty=3)
#     lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#     lines(dns$lo2$x,dns$lo2$y,lty=4)
#     lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#     lines(dns$med2$x,dns$med2$y,lty=5)
#     lines(dns$med2$x,d2med,col="red",lty=5)
# 
#     lines(dns$hi2$x,dns$hi2$y,lty=6)
#     lines(dns$hi2$x,d2hi,col="red",lty=6)
  
}

#####################  UL 3 ----

rlba.lnorm.3C <- function (n,A,b,d1,d2,p_flc,t0,meanlog_v,sdlog_v,st0=0) 
  # response2: 1 = low confidence, 2 = medium confidence, 3 = high confidence
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(meanlog_v)), length(meanlog_v), dim(meanlog_v)[1])
  drifts <- exp(matrix(rnorm(n = n * n_v, mean = meanlog_v, sd = sdlog_v), 
                       nrow = n_v))
  
  starts <- matrix(runif(min = 0, max = A, n = n * n_v), nrow = n_v)
  
  ttf <- (b - starts)/drifts
  resp <- apply(ttf, 2, which.min)
  looser <- resp; looser[resp==1] <- 2; looser[resp==2] <- 1
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  le <- rt*drifts[cbind(looser,1:n)]+starts[cbind(looser,1:n)]
  resp2 <- rep(3,n)                # High confidence response
  resp2[le>d1[looser] & le<=d2[looser]] <- 2  # Medium confidence response response
  resp2[le>d2[looser]] <- 1                   # Low confidence response response
  
  nflc<-rbinom(1,n,p_flc)
  
  if(nflc>0)
  {
    ttd1 <- (d1 - starts)/drifts
    respf<- apply(ttd1[,1:nflc,drop=FALSE], 2, which.min)
    rt[1:nflc]<-ttd1[cbind(respf,1:nflc)]
    resp[1:nflc]<-respf
    resp2[1:nflc]<-1
  }
  
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


# tmp <- rlba.3C(n=10,A=1,b=3,d1=1,d2=2,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)

# dt = c(.5,1,1.5,.5,1,1.5); r2 = c(rep("1",3),rep("2",3))
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.lnorm.3C <- function(dt,A,b,d1,d2,p_flc,meanlog_v,sdlog_v, 
                                  r2="3",robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(meanlog_v)) meanlog_v <- matrix(rep(meanlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog_v)) sdlog_v <- matrix(rep(sdlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    dlba.lnorm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],meanlog_v=meanlog_v[1,is.high],
               sdlog_v=sdlog_v[1,is.high])*
    (1-plba.lnorm(dt[2,is.high],A=A[2,is.high],b=d1[2,is.high],
                  meanlog_v=meanlog_v[2,is.high],sdlog_v=sdlog_v[2,is.high]))*
    (1-p_flc[1])
  
  if ( any(is.med) ) dt[1,is.med] <- 
    dlba.lnorm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
               meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust)*
    (plba.lnorm(dt[2,is.med],A=A[2,is.med],b=d1[2,is.med],
                meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust)-
       plba.lnorm(dt[2,is.med],A=A[2,is.med],b=d2[2,is.med],
                  meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    ((1-p_flc[1])*dlba.lnorm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
                             meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
       (plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=d2[2,is.lo],
                   meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)-
          plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                     meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)))+
    (p_flc[1]*dlba.lnorm(dt[1,is.lo],A=A[1,is.lo],b=d1[1,is.lo],
                         meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
       (1-plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=d1[2,is.lo],
                     meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)))
  
  dt[1,]
}


n1PDF.lnorm.3C <- function (t,r2,A,b,d1,d2,p_flc,meanlog_v,sdlog_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.lnorm.3C(dt,r2=r2,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d1,d2,p_flc,meanlog_v,sdlog_v) 
      n1PDFfixedt0.3C(tau,r2,A,b,d1,d2,meanlog_v,sdlog_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,b=b,d1=d1,d2=d2,meanlog_v=meanlog_v,sdlog_v=sdlog_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}


{
# Check
#     n=1e5
# 
#     t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     t0=.2; A=c(1,1); d1=c(1.5,1.5); d2=c(1.75,1.75);p_flc=.2; b=c(2,2); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     sim <- rlba.lnorm.3C(n=n,A=A,b=b,d1=d1,d2=d2,p_flc=.2,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     names(sim) <- c("RT","R","R2")
#     is.1 <- sim$R==1 & sim$R2==3
#     is.2 <- sim$R==1 & sim$R2==2
#     is.3 <- sim$R==1 & sim$R2==1
#     is.4 <- sim$R==2 & sim$R2==1
#     is.5 <- sim$R==2 & sim$R2==2
#     is.6 <- sim$R==2 & sim$R2==3
#     sim$R[is.1] <- 1
#     sim$R[is.2] <- 2
#     sim$R[is.3] <- 3
#     sim$R[is.4] <- 4
#     sim$R[is.5] <- 5
#     sim$R[is.6] <- 6
#     sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#     dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
#     tapply(sim$RT,sim$R,mean)
# 
#     # Sum to 1?
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="3")$value+
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3C,lower=0,upper=Inf,A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="3")$value
# 
#     # n1PDF check correct
#     # red=black?
#     d1lo <- n1PDF.lnorm.3C(dns$lo1$x,r2="1",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2lo <- n1PDF.lnorm.3C(dns$lo2$x,r2="1",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1med <- n1PDF.lnorm.3C(dns$med1$x,r2="2",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2med <- n1PDF.lnorm.3C(dns$med2$x,r2="2",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1hi <- n1PDF.lnorm.3C(dns$hi1$x,r2="3",A=A,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2hi <- n1PDF.lnorm.3C(dns$hi2$x,r2="3",A=A[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
# 
#     plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#       ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#     lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#     lines(dns$med1$x,dns$med1$y,lty=2)
#     lines(dns$med1$x,d1med,col="red",lty=2)
# 
#     lines(dns$lo1$x,dns$lo1$y,lty=3)
#     lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#     lines(dns$lo2$x,dns$lo2$y,lty=4)
#     lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#     lines(dns$med2$x,dns$med2$y,lty=5)
#     lines(dns$med2$x,d2med,col="red",lty=5)
# 
#     lines(dns$hi2$x,dns$hi2$y,lty=6)
#     lines(dns$hi2$x,d2hi,col="red",lty=6)
#   
}

#####################  LBA 3 A  ---- 

rlba.3C_A <- function (n,A,A1,A2,b,d1,d2,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
  # response2: 1 = low confidence, 2 = medium confidence, 3 = high confidence
{

  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) 
    drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                            lower = 0), nrow = n_v)
  else drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                        nrow = n_v)
  drifts[drifts < 0] <- 0
  
  ends <- matrix(runif(min = b-A, max = b, n = n * n_v), nrow = n_v)
  d2t<- matrix(runif(min = d2-A2, max = d2, n = n * n_v), nrow = n_v)
  d1t<- matrix(runif(min = d1-A1, max = d1, n = n * n_v), nrow = n_v)

  ttf <- (ends)/drifts
  resp <- apply(ttf, 2, which.min)
  looser <- 3-resp
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  le <- rt*drifts[cbind(looser,1:n)]
  resp2 <- rep(3,n)                # High confidence response
  resp2[le>d1t[cbind(looser,1:n)] & le<=d2t[cbind(looser,1:n)]] <- 2  # Medium confidence response response
  resp2[le>d2t[cbind(looser,1:n)]] <- 1                   # Low confidence response response
  
  nflc <- rbinom(1,n,p_flc)
  
  if(nflc>0)
  {
    ttd1 <- (d1t)/drifts
    # ADDED DROP
    respf<- apply(ttd1[,1:nflc,drop=FALSE], 2, which.min)
    
    rt[1:nflc]<-ttd1[cbind(respf,1:nflc)]
    resp[1:nflc]<-respf
    resp2[1:nflc]<-1
  }
  
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


# tmp <- rlba.3C(n=10,A=1,b=3,d1=1,d2=2,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)

# dt = c(.5,1,1.5,.5,1,1.5); r2 = c(rep("1",3),rep("2",3))
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.3C_A <- function(dt,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v, 
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
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A1)) A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A2)) A2 <- matrix(rep(A2,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    dlba.norm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],mean_v=mean_v[1,is.high],
              sd_v=sd_v[1,is.high],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.high],A=A1[2,is.high],b=d1[2,is.high],
                 mean_v=mean_v[2,is.high],sd_v=sd_v[2,is.high],posdrift=posdrift,robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.med) ) dt[1,is.med] <- 
    dlba.norm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
              mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.med],A=A1[2,is.med],b=d1[2,is.med],
               mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust)-
       plba.norm(dt[2,is.med],A=A2[2,is.med],b=d2[2,is.med],
                 mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    ((1-p_flc[1])*dlba.norm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
                            mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
       (plba.norm(dt[2,is.lo],A=A2[2,is.lo],b=d2[2,is.lo],
                  mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)-
          plba.norm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                    mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))+
    (p_flc[1]*dlba.norm(dt[1,is.lo],A=A1[1,is.lo],b=d1[1,is.lo],
                        mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
       (1-plba.norm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
                    mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}


n1PDF.3C_A <- function (t,r2,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.3C_A(dt,r2=r2,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,r2,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v) 
      n1PDFfixedt0.3C_A(tau,r2,A,A1,A2,b,d1,d2,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    # ADDED p_flc
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
#  #Check
# n=1e5
# t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); mean_v=c(2,1); sd_v=c(1,1)
# t0=.2; A=c(1,1);A1=c(.1,.1);A2=c(.1,.1); d1=c(1.25,1.25); d2=c(1.75,1.75);p_flc=.1; b=c(2.5,2.5); mean_v=c(2,1); sd_v=c(1,1)
# sim <- rlba.3C_A(n=n,A=A,A1=A1,A2=A2, b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# names(sim) <- c("RT","R","R2")
# is.1 <- sim$R==1 & sim$R2==3
# is.2 <- sim$R==1 & sim$R2==2
# is.3 <- sim$R==1 & sim$R2==1
# is.4 <- sim$R==2 & sim$R2==1
# is.5 <- sim$R==2 & sim$R2==2
# is.6 <- sim$R==2 & sim$R2==3
# sim$R[is.1] <- 1
# sim$R[is.2] <- 2
# sim$R[is.3] <- 3
# sim$R[is.4] <- 4
# sim$R[is.5] <- 5
# sim$R[is.6] <- 6
# sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
# dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
# tapply(sim$RT,sim$R,mean)
# 
# # Sum to 1?
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="3")$value+
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="2")$value+
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value+
# integrate(n1PDFfixedt0.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="3")$value
# 
# # n1PDF check correct
# # red=black?
# d1lo <- n1PDF.3C_A(dns$lo1$x,r2="1",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2lo <- n1PDF.3C_A(dns$lo2$x,r2="1",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1med <- n1PDF.3C_A(dns$med1$x,r2="2",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2med <- n1PDF.3C_A(dns$med2$x,r2="2",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# d1hi <- n1PDF.3C_A(dns$hi1$x,r2="3",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
# d2hi <- n1PDF.3C_A(dns$hi2$x,r2="3",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# 
# plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#   ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
# lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
# lines(dns$med1$x,dns$med1$y,lty=2)
# lines(dns$med1$x,d1med,col="red",lty=2)
# 
# lines(dns$lo1$x,dns$lo1$y,lty=3)
# lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
# lines(dns$lo2$x,dns$lo2$y,lty=4)
# lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
# lines(dns$med2$x,dns$med2$y,lty=5)
# lines(dns$med2$x,d2med,col="red",lty=5)
# 
# lines(dns$hi2$x,dns$hi2$y,lty=6)
# lines(dns$hi2$x,d2hi,col="red",lty=6)
}

#####################  LBA 3TC A. ----

rlba.3TC_A <- function (n,A,A1,A2,b,d1,d2,p_flc,t0,mean_v,sd_v,st0=0,posdrift=TRUE) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  if ( posdrift ) drifts <- matrix(rtnorm(n = n * n_v, mean = mean_v, sd = sd_v, 
                                   lower = 0), nrow = n_v) else { 
    drifts <- matrix(rnorm(n = n * n_v, mean = mean_v, sd = sd_v), 
                          nrow = n_v)
    drifts[drifts < 0] <- 0
  } 
  
  ttb  <- matrix(runif(min = b-A,   max = b,  n = n * n_v), nrow = n_v)/drifts  #time to b
  ttd2 <- matrix(runif(min = d2-A2, max = d2, n = n * n_v), nrow = n_v)/drifts  # time to d2
  ttd1 <- matrix(runif(min = d1-A1, max = d1, n = n * n_v), nrow = n_v)/drifts  # time to d1
 
  resphi   <- apply(ttb, 2, which.min)                            # min at hi
  respme   <- apply(ttd2,2, which.min)                            # min at med
  # resplo <- apply(ttd2,2,which.min)                             
  looserlo <- 3-respme  # equivilent to apply(ttd2,2,which.max)   # max at med
  
  #Set everything to High confidence and work down.
  resp2 <- rep(3,n)  # High confidence response
  rt <- ttb[cbind(resphi,1:n)]
  resp <- resphi
  
  looser_gt_lo <- (rt > ttd1[cbind(looserlo,1:n)])
  loose_lo_after_win_med <- (ttd1[cbind(looserlo,1:n)] > ttd2[cbind(respme,1:n)])

  # vmed conf logical
  is.med <- looser_gt_lo & loose_lo_after_win_med  
  #med conf responses and RTs
  resp2[is.med] <- 2   
  rt[is.med] <- ttd1[cbind(looserlo,1:n)][is.med]
  resp[is.med] <- respme[is.med]   
  
  #low confidence logical
  is.lo <- looser_gt_lo & !loose_lo_after_win_med
  #low confidence responses and RTs
  resp2[is.lo] <- 1    
  rt[is.lo] <- ttd2[cbind(respme,1:n)][is.lo]
  resp[is.lo] <- respme[is.lo]
  
  nflc <- rbinom(1,n,p_flc)
  if ( nflc>0 )
  {
    resplo <- apply(ttd2[,1:nflc,drop=FALSE],2,which.min)  
    rt[1:nflc] <- ttd1[cbind(resplo,1:nflc)]
    resp[1:nflc] <- resplo
    resp2[1:nflc] <- 1
  }
  
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


n1PDFfixedt0.3TC_A <- function(dt,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v, 
                               r2,posdrift=TRUE,robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt)))   dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(mean_v)) mean_v <- matrix(rep(mean_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sd_v))   sd_v <- matrix(rep(sd_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b))      b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d1))     d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2))     d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A))      A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A1))     A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A2))     A2 <- matrix(rep(A2,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  if ( any(is.high) ) dt[1,is.high] <- 
    (1-p_flc[1])*dlba.norm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],mean_v=mean_v[1,is.high],
                           sd_v=sd_v[1,is.high],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.high],A=A1[2,is.high],b=d1[2,is.high],
                 mean_v=mean_v[2,is.high],sd_v=sd_v[2,is.high],posdrift=posdrift,robust=robust))
  
  if ( any(is.med) ) dt[1,is.med] <- 
    (1-p_flc[1])*dlba.norm(dt[2,is.med],A=A1[2,is.med],b=d1[2,is.med],
                           mean_v=mean_v[2,is.med],sd_v=sd_v[2,is.med],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[1,is.med],A=A2[1,is.med],b=d2[1,is.med],
               mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust) -
       plba.norm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
                  mean_v=mean_v[1,is.med],sd_v=sd_v[1,is.med],posdrift=posdrift,robust=robust))
    
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    (1-p_flc[1])*dlba.norm(dt[1,is.lo],A=A2[1,is.lo],b=d2[1,is.lo],
                           mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (plba.norm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
               mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust) -
       plba.norm(dt[2,is.lo],A=A2[2,is.lo],b=d2[2,is.lo],
                 mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust)) +
    (p_flc[1])*dlba.norm(dt[1,is.lo],A=A1[1,is.lo],b=d1[1,is.lo],
                         mean_v=mean_v[1,is.lo],sd_v=sd_v[1,is.lo],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
                 mean_v=mean_v[2,is.lo],sd_v=sd_v[2,is.lo],posdrift=posdrift,robust=robust))
  
  dt[1,]
}


n1PDF.3TC_A <- function (t,r2,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.3TC_A(dt,r2=r2,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v, 
                            posdrift=posdrift,robust = robust))
  else {
#     tmpf <- function(tau,ti0,r2,A,A1,A2,b,d1,d2,p_flc,mean_v,sd_v) 
#       n1PDFfixedt0.3TC(tau,ti0=dt[i],r2,A,A1,A2,b,d1,d2,mean_v,sd_v)/st0
#     outs <- numeric(length(dt))
#     for (i in c(1:length(outs))) {
#       tmp <- try(integrate(f=tmpf,
#                            lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
#                            A=A,b=b,d1=d1,d2=d2,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
#       if (is.numeric(tmp)) 
#         outs[i] <- tmp else
#           outs[i] <- 0
#     }
#     return(outs)
    #fix it
  }
}


{
#  #Check
#   n=1e5
# 
#   t0=.2; A=c(.5,.5);A1=c(.5,.5);A2=c(.3,.3); d1=c(1.5,1.5); d2=c(1.75,1.75);p_flc=.5; b=c(2.5,2.5); mean_v=c(2,1); sd_v=c(1,1)
#   sim <- rlba.3TC_A(n=n,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   names(sim) <- c("RT","R","R2")
#   is.1 <- sim$R==1 & sim$R2==3
#   is.2 <- sim$R==1 & sim$R2==2
#   is.3 <- sim$R==1 & sim$R2==1
#   is.4 <- sim$R==2 & sim$R2==1
#   is.5 <- sim$R==2 & sim$R2==2
#   is.6 <- sim$R==2 & sim$R2==3
#   sim$R[is.1] <- 1
#   sim$R[is.2] <- 2
#   sim$R[is.3] <- 3
#   sim$R[is.4] <- 4
#   sim$R[is.5] <- 5
#   sim$R[is.6] <- 6
#   sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#   dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
#   tapply(sim$RT,sim$R,mean)
# 
#   # Sum to 1?
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="3")$value+
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="2")$value+ #
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,r2="1")$value+
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1,A2=A2,b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="1")$value+
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1,A2=A2,b=b[2:1],d1=d1[2:1],d2=d2[2:1],
#     p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="2")$value+ #
#   integrate(n1PDFfixedt0.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1,A2=A2,b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],r2="3")$value
# 
#   # n1PDF check correct
#   # red=black?
#   d1lo <- n1PDF.3TC_A(dns$lo1$x,r2="1",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2lo <- n1PDF.3TC_A(dns$lo2$x,r2="1",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
#   d1med <- n1PDF.3TC_A(dns$med1$x,r2="2",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2med <- n1PDF.3TC_A(dns$med2$x,r2="2",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
#   d1hi <- n1PDF.3TC_A(dns$hi1$x,r2="3",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,mean_v=mean_v,sd_v=sd_v,t0=t0)
#   d2hi <- n1PDF.3TC_A(dns$hi2$x,r2="3",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,mean_v=mean_v[2:1],sd_v=sd_v[2:1],t0=t0)
# 
#   plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#     ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#   lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#   lines(dns$med1$x,dns$med1$y,lty=2)
#   lines(dns$med1$x,d1med,col="red",lty=2)
# 
#   lines(dns$lo1$x,dns$lo1$y,lty=3)
#   lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#   lines(dns$lo2$x,dns$lo2$y,lty=4)
#   lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#   lines(dns$med2$x,dns$med2$y,lty=5)
#   lines(dns$med2$x,d2med,col="red",lty=5)
# 
#   lines(dns$hi2$x,dns$hi2$y,lty=6)
#   lines(dns$hi2$x,d2hi,col="red",lty=6)

}



#####################  UL 3 A To be created ----
rlba.lnorm.3C_A <- function (n,A,A1,A2,b,d1,d2,p_flc,t0,meanlog_v,sdlog_v,st0=0) 
  # response2: 1 = low confidence, 2 = medium confidence, 3 = high confidence
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(meanlog_v)), length(meanlog_v), dim(meanlog_v)[1])
  drifts <- exp(matrix(rnorm(n = n * n_v, mean = meanlog_v, sd = sdlog_v), 
                       nrow = n_v))
  
  ends <- matrix(runif(min = b-A, max = b, n = n * n_v), nrow = n_v)
  d2t<- matrix(runif(min = d2-A2, max = d2, n = n * n_v), nrow = n_v)
  d1t<- matrix(runif(min = d1-A1, max = d1, n = n * n_v), nrow = n_v)
  
  ttf <- (ends)/drifts
  resp <- apply(ttf, 2, which.min)
  looser <- 3-resp
  rt <- ttf[cbind(resp,1:n)] # actaully decision time, add t0 later
  
  # looser evidence = decison time * looser drift
  le <- rt*drifts[cbind(looser,1:n)]
  resp2 <- rep(3,n)                # High confidence response
  resp2[le>d1t[cbind(looser,1:n)] & le<=d2t[cbind(looser,1:n)]] <- 2  # Medium confidence response response
  resp2[le>d2t[cbind(looser,1:n)]] <- 1                   # Low confidence response response
  
  nflc <- rbinom(1,n,p_flc)
  
  if(nflc>0)
  {
    ttd1 <- (d1t)/drifts
    # ADDED DROP
    respf<- apply(ttd1[,1:nflc,drop=FALSE], 2, which.min)
    
    rt[1:nflc]<-ttd1[cbind(respf,1:nflc)]
    resp[1:nflc]<-respf
    resp2[1:nflc]<-1
  }
  
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


# tmp <- rlba.3C(n=10,A=1,b=3,d1=1,d2=2,t0=.1,mean_v=c(1,0),sd_v=1,st0=0,posdrift=TRUE)

# dt = c(.5,1,1.5,.5,1,1.5); r2 = c(rep("1",3),rep("2",3))
# d=A=c(1,1); b=d+1; mean_v=c(1,0); sd_v=c(1,1)

n1PDFfixedt0.lnorm.3C_A <- function(dt,A,A1,A2,b,d1,d2,p_flc,meanlog_v,sdlog_v, 
                                  r2="3",robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(meanlog_v)) meanlog_v <- matrix(rep(meanlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog_v)) sdlog_v <- matrix(rep(sdlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A1)) A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  
  if (!is.matrix(A2)) A2 <- matrix(rep(A2,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    dlba.lnorm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],meanlog_v=meanlog_v[1,is.high],
               sdlog_v=sdlog_v[1,is.high])*
    (1-plba.lnorm(dt[2,is.high],A=A1[2,is.high],b=d1[2,is.high],
                  meanlog_v=meanlog_v[2,is.high],sdlog_v=sdlog_v[2,is.high]))*
    (1-p_flc[1])
  
  if ( any(is.med) ) dt[1,is.med] <- 
    dlba.lnorm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
               meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust)*
    (plba.lnorm(dt[2,is.med],A=A1[2,is.med],b=d1[2,is.med],
                meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust)-
       plba.lnorm(dt[2,is.med],A=A2[2,is.med],b=d2[2,is.med],
                  meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust))*
    (1-p_flc[1])
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    ((1-p_flc[1])*dlba.lnorm(dt[1,is.lo],A=A[1,is.lo],b=b[1,is.lo],
                             meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
       (plba.lnorm(dt[2,is.lo],A=A2[2,is.lo],b=d2[2,is.lo],
                   meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)-
          plba.lnorm(dt[2,is.lo],A=A[2,is.lo],b=b[2,is.lo],
                     meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)))+
    (p_flc[1]*dlba.lnorm(dt[1,is.lo],A=A1[1,is.lo],b=d1[1,is.lo],
                         meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
       (1-plba.lnorm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
                     meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)))
  
  dt[1,]
}


n1PDF.lnorm.3C_A <- function (t,r2,A,A1,A2,b,d1,d2,p_flc,meanlog_v,sdlog_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.lnorm.3C_A(dt,r2=r2,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v))
  else {
    tmpf <- function(tau,ti0,r2,A,A1=A1,A2=A2,b,d1,d2,p_flc,meanlog_v,sdlog_v) 
      n1PDFfixedt0.3C_A(tau,r2,A,A1,A2,b,d1,d2,meanlog_v,sdlog_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,meanlog_v=meanlog_v,sdlog_v=sdlog_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}


{
#     #Check
#     n=1e5
# 
#     t0=.2; A=c(.5,.5); d1=c(1.25,.75); d2=c(1.5,1); b=c(2,1.5); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     t0=.2; A=c(.1,.1);A1=c(.1,.1);A2=c(.1,.1); d1=c(1.5,1.5); d2=c(1.75,1.75);p_flc=0; b=c(2,2); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     sim <- rlba.lnorm.3C_A(n=n,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     names(sim) <- c("RT","R","R2")
#     is.1 <- sim$R==1 & sim$R2==3
#     is.2 <- sim$R==1 & sim$R2==2
#     is.3 <- sim$R==1 & sim$R2==1
#     is.4 <- sim$R==2 & sim$R2==1
#     is.5 <- sim$R==2 & sim$R2==2
#     is.6 <- sim$R==2 & sim$R2==3
#     sim$R[is.1] <- 1
#     sim$R[is.2] <- 2
#     sim$R[is.3] <- 3
#     sim$R[is.4] <- 4
#     sim$R[is.5] <- 5
#     sim$R[is.6] <- 6
#     sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#     dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
#     tapply(sim$RT,sim$R,mean)
# 
#     # Sum to 1?
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="3")$value+
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3C_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="3")$value
# 
#     # n1PDF check correct
#     # red=black?
#     d1lo <- n1PDF.lnorm.3C_A(dns$lo1$x,r2="1",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2lo <- n1PDF.lnorm.3C_A(dns$lo2$x,r2="1",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1med <- n1PDF.lnorm.3C_A(dns$med1$x,r2="2",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2med <- n1PDF.lnorm.3C_A(dns$med2$x,r2="2",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1hi <- n1PDF.lnorm.3C_A(dns$hi1$x,r2="3",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2hi <- n1PDF.lnorm.3C_A(dns$hi2$x,r2="3",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
# 
#     plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#       ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#     lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#     lines(dns$med1$x,dns$med1$y,lty=2)
#     lines(dns$med1$x,d1med,col="red",lty=2)
# 
#     lines(dns$lo1$x,dns$lo1$y,lty=3)
#     lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#     lines(dns$lo2$x,dns$lo2$y,lty=4)
#     lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#     lines(dns$med2$x,dns$med2$y,lty=5)
#     lines(dns$med2$x,d2med,col="red",lty=5)
# 
#     lines(dns$hi2$x,dns$hi2$y,lty=6)
#     lines(dns$hi2$x,d2hi,col="red",lty=6)
  
}
#####################  UL 3TC A  . ----

rlba.lnorm.3TC_A <- function (n,A,A1,A2,b,d1,d2,p_flc,t0,meanlog_v,sdlog_v,st0=0) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(meanlog_v)), length(meanlog_v), dim(meanlog_v)[1])
  drifts <- exp(matrix(rnorm(n = n * n_v, mean = meanlog_v, sd = sdlog_v), 
                       nrow = n_v))
  
  ttb  <- matrix(runif(min = b-A,   max = b,  n = n * n_v), nrow = n_v)/drifts  #time to b
  ttd2 <- matrix(runif(min = d2-A2, max = d2, n = n * n_v), nrow = n_v)/drifts  # time to d2
  ttd1 <- matrix(runif(min = d1-A1, max = d1, n = n * n_v), nrow = n_v)/drifts  # time to d1
  
  resphi   <- apply(ttb, 2, which.min)                            # min at hi
  respme   <- apply(ttd2,2, which.min)                            # min at med
  # resplo <- apply(ttd2,2,which.min)                             
  looserlo <- 3-respme  # equivilent to apply(ttd2,2,which.max)   # max at med
  
  #Set everything to High confidence and work down.
  resp2 <- rep(3,n)  # High confidence response
  rt <- ttb[cbind(resphi,1:n)]
  resp <- resphi
  
  looser_gt_lo <- (rt > ttd1[cbind(looserlo,1:n)])
  loose_lo_after_win_med <- (ttd1[cbind(looserlo,1:n)] > ttd2[cbind(respme,1:n)])
  
  # vmed conf logical
  is.med <- looser_gt_lo & loose_lo_after_win_med  
  #med conf responses and RTs
  resp2[is.med] <- 2   
  rt[is.med] <- ttd1[cbind(looserlo,1:n)][is.med]
  resp[is.med] <- respme[is.med]   
  
  #low confidence logical
  is.lo <- looser_gt_lo & !loose_lo_after_win_med
  #low confidence responses and RTs
  resp2[is.lo] <- 1    
  rt[is.lo] <- ttd2[cbind(respme,1:n)][is.lo]
  resp[is.lo] <- respme[is.lo]
  
  
  
  nflc<-rbinom(1,n,p_flc)
  if(nflc>0)
  {
    resplo<- apply(ttd2[,1:nflc,drop=FALSE],2,which.min)  
    rt[1:nflc]<-ttd1[cbind(resplo,1:nflc)]
    resp[1:nflc]<-resplo
    resp2[1:nflc]<-1
    
  }
  
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


n1PDFfixedt0.lnorm.3TC_A <- function(dt,A,A1,A2,b,d1,d2,p_flc,meanlog_v,sdlog_v, 
                                     r2="3",robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(meanlog_v)) meanlog_v <- matrix(rep(meanlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sdlog_v)) sdlog_v <- matrix(rep(sdlog_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d1)) d1 <- matrix(rep(d1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d2)) d2 <- matrix(rep(d2,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A1)) A1 <- matrix(rep(A1,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A2)) A2 <- matrix(rep(A2,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.lo <- r2 =="1"
  is.med <- r2 =="2"
  is.high <- r2 =="3"
  
  
  if ( any(is.high) ) dt[1,is.high] <- 
    (1-p_flc[1])*dlba.lnorm(dt[1,is.high],A=A[1,is.high],b=b[1,is.high],meanlog_v=meanlog_v[1,is.high],
                            sdlog_v=sdlog_v[1,is.high],robust=robust)*
    (1-plba.lnorm(dt[2,is.high],A=A1[2,is.high],b=d1[2,is.high],
                  meanlog_v=meanlog_v[2,is.high],sdlog_v=sdlog_v[2,is.high],robust=robust))
  
  if ( any(is.med) ) dt[1,is.med] <- 
    (1-p_flc[1])*dlba.lnorm(dt[2,is.med],A=A1[2,is.med],b=d1[2,is.med],
                            meanlog_v=meanlog_v[2,is.med],sdlog_v=sdlog_v[2,is.med],robust=robust)*
    (plba.lnorm(dt[1,is.med],A=A2[1,is.med],b=d2[1,is.med],
                meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust)-
       plba.lnorm(dt[1,is.med],A=A[1,is.med],b=b[1,is.med],
                  meanlog_v=meanlog_v[1,is.med],sdlog_v=sdlog_v[1,is.med],robust=robust))
  
  if ( any(is.lo) ) dt[1,is.lo] <- 
    (1-p_flc[1])*dlba.lnorm(dt[1,is.lo],A=A2[1,is.lo],b=d2[1,is.lo],
                            meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
    (plba.lnorm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
                meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust)-
       plba.lnorm(dt[2,is.lo],A=A2[2,is.lo],b=d2[2,is.lo],
                  meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust))+
    (p_flc[1])*dlba.lnorm(dt[1,is.lo],A=A1[1,is.lo],b=d1[1,is.lo],
                          meanlog_v=meanlog_v[1,is.lo],sdlog_v=sdlog_v[1,is.lo],robust=robust)*
    (1-plba.lnorm(dt[2,is.lo],A=A1[2,is.lo],b=d1[2,is.lo],
                  meanlog_v=meanlog_v[2,is.lo],sdlog_v=sdlog_v[2,is.lo],robust=robust))
  
  dt[1,]
}


n1PDF.lnorm.3TC_A <- function (t,r2,A,A1,A2,b,d1,d2,p_flc,meanlog_v,sdlog_v,t0,st0=0,robust = FALSE)
  # r2 is second respond, d1 is low crierion, d2 is medium
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.lnorm.3TC_A(dt,r2=r2,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,
                                    meanlog_v=meanlog_v,sdlog_v=sdlog_v, robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,A1,A2,b,d1,d2,p_flc,meanlog_v,sdlog_v) 
      n1PDFfixedt0.lnorm.3TC_A(tau,r2,A,A1,A2,b,d1,d2,meanlog_v,sdlog_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],p_flc=p_flc,
                           A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,meanlog_v=meanlog_v,sdlog_v=sdlog_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

{
#     #Check
#     n=1e5
# 
#     t0=.2; A=c(.2,.2);A1=c(.25,.25);A2=c(.25,.25); d1=c(1,1); d2=c(1.5,1.5);p_flc=0; b=c(2,2); meanlog_v=c(2,1); sdlog_v=c(1,1)
#     sim <- rlba.lnorm.3TC_A(n=n,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     names(sim) <- c("RT","R","R2")
#     is.1 <- sim$R==1 & sim$R2==3
#     is.2 <- sim$R==1 & sim$R2==2
#     is.3 <- sim$R==1 & sim$R2==1
#     is.4 <- sim$R==2 & sim$R2==1
#     is.5 <- sim$R==2 & sim$R2==2
#     is.6 <- sim$R==2 & sim$R2==3
#     sim$R[is.1] <- 1
#     sim$R[is.2] <- 2
#     sim$R[is.3] <- 3
#     sim$R[is.4] <- 4
#     sim$R[is.5] <- 5
#     sim$R[is.6] <- 6
#     sim$R <- factor(sim$R,labels=c("hi1","med1","lo1","lo2","med2","hi2"))
#     dns <- plot.cell.density(sim,xlim=c(0,5),save.density=TRUE)
#     tapply(sim$RT,sim$R,mean)
# 
#     # Sum to 1?
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="3")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="1")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="2")$value+
#     integrate(n1PDFfixedt0.lnorm.3TC_A,lower=0,upper=Inf,A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],r2="3")$value
# 
#     # n1PDF check correct
#     # red=black?
#     d1lo <- n1PDF.lnorm.3TC_A(dns$lo1$x,r2="1",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2lo <- n1PDF.lnorm.3TC_A(dns$lo2$x,r2="1",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1med <- n1PDF.lnorm.3TC_A(dns$med1$x,r2="2",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2med <- n1PDF.lnorm.3TC_A(dns$med2$x,r2="2",A=A[2:1],A1=A1[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
#     d1hi <- n1PDF.lnorm.3TC_A(dns$hi1$x,r2="3",A=A,A1=A1,A2=A2,b=b,d1=d1,d2=d2,p_flc=p_flc,meanlog_v=meanlog_v,sdlog_v=sdlog_v,t0=t0)
#     d2hi <- n1PDF.lnorm.3TC_A(dns$hi2$x,r2="3",A=A[2:1],A1=A[2:1],A2=A2[2:1],b=b[2:1],d1=d1[2:1],d2=d2[2:1],p_flc=p_flc,meanlog_v=meanlog_v[2:1],sdlog_v=sdlog_v[2:1],t0=t0)
# 
#     plot(dns$hi1$x,dns$hi1$y,type="l",ylab="density",xlab="RT",lty=1,
#       ylim=c(0,max(c(d1lo,d1med,d1hi,d2lo,d2med,d2hi))))
#     lines(dns$hi1$x,d1hi,col="red",lty=1)
# 
#     lines(dns$med1$x,dns$med1$y,lty=2)
#     lines(dns$med1$x,d1med,col="red",lty=2)
# 
#     lines(dns$lo1$x,dns$lo1$y,lty=3)
#     lines(dns$lo1$x,d1lo,col="red",lty=3)
# 
#     lines(dns$lo2$x,dns$lo2$y,lty=4)
#     lines(dns$lo2$x,d2lo,col="red",lty=4)
# 
#     lines(dns$med2$x,dns$med2$y,lty=5)
#     lines(dns$med2$x,d2med,col="red",lty=5)
# 
#     lines(dns$hi2$x,dns$hi2$y,lty=6)
#     lines(dns$hi2$x,d2hi,col="red",lty=6)
  
}


