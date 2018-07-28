### SEARCH Models (using LBA)

### Self terminiating and exhuastive of a set of objects either all lures or
### one target.

library(rtdists)

rST <- function(n, A, b, vT, vL, svT, svL, tp, nStim, t0)
  # RANDOM FUNCTION FOR SELF-TERMINATING SEARCH
  # n = number of trials to simulate
  # PARAMETERS (2-vectors, yes,no accumulator)
  # A, b:   start point/threshold for "yes" and "no" accumulators 
  # vT, vL: mean_v target/lure for "yes" and "no" accumulators 
  # svT, svL: md_v target/lure for "yes" and "no" accumulators
  # t0 = non-decision time (SCALAR!).
  #   DEFINE DISPLAY  
  # tp: target position, if is.na(tp) all lure stimuli
  # nStim: number of stimuli

{
  
  # ASSUME THIS IS OK AND COMMENT OUT IN FINAL VERSION
  # Make sure A, b, svT, and svL all have length 2 and t0 length 1, A and b ordered
  A <- rep(A,length.out=2)
  b <- rep(b,length.out=2)
  svT <- rep(svT,length.out=2)
  svL <- rep(svL,length.out=2)
  t0 <- t0[1]
  if (any(b < A)) stop("b cannot be smaller than A!")
  
  # Arrays of decision times for yes and no accumulators for each stimulus.
  yesAccumulators <- array(0, dim=c(n, nStim))
  noAccumulators <- array(0, dim=c(n, nStim))
  
  if ( !is.na(tp) ) { # target present
    yesAccumulators[,tp] <- rlba_norm(n=n, A=A[1], b=b[1], t0=0, 
                                      mean_v=vT[1], sd_v=svT[1], st0=0)$rt

    noAccumulators[,tp] <- rlba_norm(n=n, A=A[2], b=b[2], t0=0, 
                                      mean_v=vT[2], sd_v=svT[2], st0=0)$rt
    lure_cols <- c(1:nStim)[-tp]
    nLure <- nStim-length(tp)
  } else {
    lure_cols <- 1:nStim
    nLure <- nStim
  }
  yesAccumulators[,lure_cols] <- rlba_norm(n=n*nLure, A=A[1], b=b[1], t0=0, 
                                           mean_v=vL[1], sd_v=svL[1], st0=0)$rt

  noAccumulators[,lure_cols] <- rlba_norm(n=n*nLure, A=A[2], b=b[2], t0=0, 
                                          mean_v=vL[2], sd_v=svL[2], st0=0)$rt

  out <- data.frame("rt"=rep(0, n), "response"=rep(1, n))
  # A yes response is produced when the fastest "yes" accumulator finishes first
  yesmin <- apply(yesAccumulators, 1, min) 
  # A no response is produced when the slowest "no" accumulator finishes before
  # the fastest yes
  nomax <- apply(noAccumulators, 1, max)
  
  yes <- yesmin < nomax
  out$rt[yes] <- yesmin[yes]+t0
  out$rt[!yes] <- nomax[!yes]+t0
  out$response <- ifelse(yes, 1, 2)
  
  out
}

# # More likley to say yes, need v for lure twice as strong as for target
# # to get about equal proportions of correct responses here (~85%)
# A=0.3; b=1; vT=c(2, 1); vL=c(1, 4); svT=1; svL=1; nStim=3
# n=1e5; t0=0
# simST <- rST(n, A, b, vT, vL, svT, svL, tp=1, nStim, t0)
# table(simST$response)/dim(simST)[1]
# simST <- rST(n, A, b, vT, vL, svT, svL, tp=NA, nStim, t0)
# table(simST$response)/dim(simST)[1]



dST <- function(t, response, A, b, vT, vL, svT, svL, tp, nStim)
  #   DENSITY FOR SELF-TERMINATING SEARCH
  # t = decision time vector (rt-t0)
  # response: 1 = yes (target present), 2 = no (target absent)
  #   PARAMETERS (2-vectors, 1st value matches response, e.g., response=2 => no,yes)
  # A, b:   start point/threshold for "yes" and "no" accumulators 
  # vT, vL: mean_v target/lure for "yes" and "no" accumulators 
  # svT, svL: md_v target/lure for "yes" and "no" accumulators
  #   DEFINE DISPLAY  
  # tp: target position, if is.na(tp) all lure stimuli
  # nStim: number of stimuli
{

  # ASSUME THIS IS OK AND COMMENT OUT IN FINAL VERSION
  # Make sure A, b, svT, and svL all have length 2 and A and b ordered
  A <- rep(A,length.out=2)
  b <- rep(b,length.out=2)
  svT <- rep(svT,length.out=2)
  svL <- rep(svL,length.out=2)
  if (any(b < A)) stop("b cannot be smaller than A!")

  # CDF matrix for Y(es) and N(o) accumulators for each t and stimulus.
  cdfsY <- array(0, dim=c(length(t), nStim)) 
  cdfsN <- array(0, dim=c(length(t), nStim)) 
  # PDF matrix of response (i.e., Y if response=1, N if response = 2)
  pdfs <- array(0, dim=c(length(t), nStim))

  # Functions for calculating the CDF and PDF values depending on 
  #   response (1 or 2, picks element of parameter vector),
  #   stimulus (target boolean, picks vT/svT or vL/svL) and 
  #   type: Y or N accumulator (N flips parameter vector order).
  
  CDF <- function(t, A, b, vT, vL, svT, svL, target, index)
  {
    if ( target )
        plba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vT[index], sd_v=svT[index], t0=0) else
        plba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vL[index], sd_v=svL[index], t0=0)
  }
  

  PDF <- function(t, A, b, vT, vL, svT, svL, target, index){
    if ( target )
        dlba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vT[index], sd_v=svT[index], t0=0) else
        dlba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vL[index], sd_v=svL[index], t0=0)
  }
  

  # Fill cdf and pdf matrices
  if (!is.na(tp)) { # Only calcualte target stimulus if needed
      cdfsY[,tp] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(1:2)[response])
      cdfsN[,tp] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(2:1)[response])
      if (response==1) # Need pdfsY else pdfsN
        pdfs[,tp] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(1:2)[response]) else
        pdfs[,tp] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(2:1)[response])
      lure_cols <- c(1:nStim)[-tp]
  } else lure_cols <- c(1:nStim)
  # lure stimuli, calcualte only once and fill redundant columns
  cdfsY[,lure_cols] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(1:2)[response])
  cdfsN[,lure_cols] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(2:1)[response])
  if (response==1) # Need pdfsY else pdfsN
    pdfs[,lure_cols] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(1:2)[response]) else
    pdfs[,lure_cols] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(2:1)[response])
  
  
  term2 <- 0
  if ( response==1 ) { # L_YES
    term1 <- 1-apply(cdfsN,1,prod)
    for (i in 1:nStim) {
      term3 <- apply(1-cdfsY[,-i,drop=FALSE],1,prod)
      term2 <- term2+(pdfs[,i,drop=FALSE]*term3)
    }
  } else {           # L_NO
    term1 <- apply(1-cdfsY,1,prod) 
    for (i in 1:nStim) {
      term3 <- apply(cdfsN[,-i,drop=FALSE],1,prod)
      term2 <- term2+(pdfs[,i,drop=FALSE]*term3)
    }
  }  
  as.vector(term1*term2)
}

# # Checks:
# 
# A=0.3; b=1; vT=c(2, 1); vL=c(1, 4); svT=1; svL=1; nStim=3
# 
# # 1. Integrate to 1?
# tp=NA
# ps <- c(
#   integrate(dST, lower=0, upper=Inf, A=A, b=b, vT=vT, vL=vL, 
#     svT=svT, svL=svL, tp=tp,nStim=nStim,response=1)$value,
#   integrate(dST, lower=0, upper=Inf, A=A, b=b, vT=vT[2:1], vL=vL[2:1],
#     svT=svT, svL=svL, tp=tp, nStim=nStim,response=2)$value
# )
# c(ps,sum(ps))
# 
# tp=1
# ps <- c(
#   integrate(dST, lower=0, upper=Inf, A=A, b=b, vT=vT, vL=vL, 
#     svT=svT, svL=svL, tp=tp,nStim=nStim,response=1)$value,
#   integrate(dST, lower=0, upper=Inf, A=A, b=b, vT=vT[2:1], vL=vL[2:1],
#     svT=svT, svL=svL, tp=tp, nStim=nStim,response=2)$value
# )
# c(ps,sum(ps))
# 
# # Plot check
# par(mfrow=c(1,2))
# 
# tp=NA
# simST <- rST(n=1e5, A, b, vT, vL, svT, svL, tp, nStim, t0=0)
# ps <- table(simST$response)/dim(simST)[1]
# dns1 <- density(simST$rt[simST$response==1]); dns1$y <- dns1$y*ps[1]
# dns2 <- density(simST$rt[simST$response==2]); dns2$y <- dns2$y*ps[2]
# plot(dns1,main="ST(non-target)",ylim=c(0,max(c(dns1$y,dns2$y))))
# lines(x=dns1$x, y=dST(t=dns1$x, response=1, A, b, vT, vL, svT, svL, tp, nStim), col="red")
# lines(dns2,lty=2)
# lines(x=dns2$x, y=dST(t=dns2$x, response=2, A, b, vT[2:1], vL[2:1], svT, svL, tp, nStim), 
#   col="red",lty=2)
# legend("topright",c("Yes","No"),lty=1:2,bty="n")
# 
# tp=1
# simST <- rST(n=1e5, A, b, vT, vL, svT, svL, tp, nStim, t0=0)
# ps <- table(simST$response)/dim(simST)[1]
# dns1 <- density(simST$rt[simST$response==1]); dns1$y <- dns1$y*ps[1]
# dns2 <- density(simST$rt[simST$response==2]); dns2$y <- dns2$y*ps[2]
# plot(dns1,main="ST(target)",ylim=c(0,max(c(dns1$y,dns2$y))))
# lines(x=dns1$x, y=dST(t=dns1$x, response=1, A, b, vT, vL, svT, svL, tp, nStim), col="red")
# lines(dns2,lty=2)
# lines(x=dns2$x, y=dST(t=dns2$x, response=2, A, b, vT[2:1], vL[2:1], svT, svL, tp, nStim), 
#   col="red",lty=2)
# legend("topright",c("Yes","No"),lty=1:2,bty="n")
# 
# # TIMING
# 
# # 2e6 random numbers
# system.time({
#   rST(n=1e6, A, b, vT, vL, svT, svL, tp=NA, nStim, t0=0)
#   rST(n=1e6, A, b, vT, vL, svT, svL, tp=1, nStim, t0=0)
# })
# #    user  system elapsed 
# #  83.505   1.338  84.892
#  
# # Negligible speedup over Peters routine by a factor of
# 99.717/84.892
# # [1] 1.174634
# 
# # 400,000 density evaluations
# n <- 1e5
# system.time({
#   dST(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=1, A, b, vT, vL, svT, svL, tp=NA, nStim)
#   dST(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=2, A, b, vT, vL, svT, svL, tp=NA, nStim)
#   dST(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=1, A, b, vT, vL, svT, svL, tp=1, nStim)
#   dST(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=2, A, b, vT, vL, svT, svL, tp=1, nStim)
# })
# #    user  system elapsed 
# #   3.832   0.029   3.862
# 
# # Speedup over Peters routine by a factor of
# 100*16.087/3.862
# # [1] 416.5458


rEX <- function(n, A, b, vT, vL, svT, svL, tp, nStim, t0)
  # RANDOM FUNCTION FOR EXHAUSTIVE SEARCH
  # n = number of trials to simulate
  # PARAMETERS (2-vectors, yes,no accumulator)
  # A, b:   start point/threshold for "yes" and "no" accumulators 
  # vT, vL: mean_v target/lure for "yes" and "no" accumulators 
  # svT, svL: md_v target/lure for "yes" and "no" accumulators
  # t0 = non-decision time (SCALAR!).
  #   DEFINE DISPLAY  
  # tp: target position, if is.na(tp) all lure stimuli
  # nStim: number of stimuli
{
  
  # ASSUME THIS IS OK AND COMMENT OUT IN FINAL VERSION
  # Make sure A, b, svT, and svL all have length 2 and t0 length 1, A and b ordered
  A <- rep(A,length.out=2)
  b <- rep(b,length.out=2)
  svT <- rep(svT,length.out=2)
  svL <- rep(svL,length.out=2)
  t0 <- t0[1]
  if (any(b < A)) stop("b cannot be smaller than A!")
  
  # Arrays of decision times for yes and no accumulators for each stimulus.
  yesAccumulators <- array(0, dim=c(n, nStim))
  noAccumulators <- array(0, dim=c(n, nStim))
  
  if (!is.na(tp)) { # target present
    yesAccumulators[,tp] <- rlba_norm(n=n, A=A[1], b=b[1], t0=0, 
                                      mean_v=vT[1], sd_v=svT[1], st0=0)$rt

    noAccumulators[,tp] <- rlba_norm(n=n, A=A[2], b=b[2], t0=0, 
                                      mean_v=vT[2], sd_v=svT[2], st0=0)$rt
    lure_cols <- c(1:nStim)[-tp]
    nLure <- nStim-length(tp)
  } else {
    lure_cols <- 1:nStim
    nLure <- nStim
  }
  yesAccumulators[,lure_cols] <- rlba_norm(n=n*nLure, A=A[1], b=b[1], t0=0, 
                                           mean_v=vL[1], sd_v=svL[1], st0=0)$rt

  noAccumulators[,lure_cols] <- rlba_norm(n=n*nLure, A=A[2], b=b[2], t0=0, 
                                          mean_v=vL[2], sd_v=svL[2], st0=0)$rt
  
  data.frame(
    rt=apply(pmin(yesAccumulators, noAccumulators), 1, max)+t0,
    # RT is always the slowest winner  
    response=as.numeric(apply(noAccumulators < yesAccumulators, 1, all))+1
    # A yes response is produced when at least one "yes" accumulator finished 
    # before its paired "no" accumulator
    # A no response is produced when each "no" accumulators won the race with 
    # its paired "yes" accumulator (used above to get a boolean)
  )

}

# # More likley to say yes, need v for lure stronger than for target
# # to get about equal proportions of correct responses here (~77%)
# A=0.3; b=1; vT=c(2, 1); vL=c(1, 3.1); svT=1; svL=1; nStim=3
# n=1e5; t0=0
# simST <- rEX(n, A, b, vT, vL, svT, svL, tp=1, nStim, t0)
# table(simST$response)/dim(simST)[1]
# simST <- rEX(n, A, b, vT, vL, svT, svL, tp=NA, nStim, t0)
# table(simST$response)/dim(simST)[1]


dEX <- function(t, response, A, b, vT, vL, svT, svL, tp, nStim)
  #   DENSITY FOR SELF-TERMINATING SEARCH
  # t = decision time vector (rt-t0)
  # response: 1 = yes (target present), 2 = no (target absent)
  #   PARAMETERS (2-vectors, 1st value matches response, e.g., response=2 => no,yes)
  # A, b:   start point/threshold for "yes" and "no" accumulators 
  # vT, vL: mean_v target/lure for "yes" and "no" accumulators 
  # svT, svL: md_v target/lure for "yes" and "no" accumulators
  #   DEFINE DISPLAY  
  # tp: target position, if is.na(tp) all lure stimuli
  # nStim: number of stimuli
{
  
  # ASSUME THIS IS OK AND COMMENT OUT IN FINAL VERSION
  # Make sure A, b, svT, and svL all have length 2 and A and b ordered
  A <- rep(A,length.out=2)
  b <- rep(b,length.out=2)
  svT <- rep(svT,length.out=2)
  svL <- rep(svL,length.out=2)
  if (any(b < A)) stop("b cannot be smaller than A!")
  
  # CDF matrix for Y(es) and N(o) accumulators for each t and stimulus.
  cdfsY <- array(0, dim=c(length(t), nStim)) 
  cdfsN <- array(0, dim=c(length(t), nStim)) 
  # PDF matrix of no response 
  pdfsN <- array(0, dim=c(length(t), nStim))
  # Only need pdfsY for target display but still need pdfsN
  if (response==1) pdfsY <- array(0, dim=c(length(t), nStim))

  # Functions for calculating the CDF and PDF values depending on 
  #   response (1 or 2, picks element of parameter vector),
  #   stimulus (target boolean, picks vT/svT or vL/svL) and 
  #   type: Y or N accumulator (N flips parameter vector order).
  
  CDF <- function(t, A, b, vT, vL, svT, svL, target, index)
  {
    if ( target )
        plba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vT[index], sd_v=svT[index], t0=0) else
        plba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vL[index], sd_v=svL[index], t0=0)
  }
  

  PDF <- function(t, A, b, vT, vL, svT, svL, target, index){
    if ( target )
        dlba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vT[index], sd_v=svT[index], t0=0) else
        dlba_norm(rt=t, A=A[index], b=b[index], 
                  mean_v=vL[index], sd_v=svL[index], t0=0)
  }
  

  # Fill cdf and pdf matrices
  if ( !is.na(tp) ) { # Only calcualte target stimulus if needed
      cdfsY[,tp] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(1:2)[response])
      cdfsN[,tp] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(2:1)[response])
      if (response==1) # Dont need for response==2
        pdfsY[,tp] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                          target=TRUE, index=c(1:2)[response]) else
      pdfsN[,tp] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                        target=TRUE, index=c(2:1)[response])
      lure_cols <- c(1:nStim)[-tp]
  } else lure_cols <- c(1:nStim)
  # lure stimuli, calcualte only once and fill redundant columns
  cdfsY[,lure_cols] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(1:2)[response])
  cdfsN[,lure_cols] <- CDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(2:1)[response])
  if (response==1) # Dont need for response==2 
    pdfsY[,lure_cols] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                             target=FALSE, index=c(1:2)[response]) else
  pdfsN[,lure_cols] <- PDF(t, A=A, b=b, vT=vT, vL=vL, svT=svT, svL=svL,
                           target=FALSE, index=c(2:1)[response])
  
  term1 <- 0
  if ( response==1 ) { # L_YES
    term2 <- 0
    for (i in 1:nStim){
      term1b <- apply(1-(1-cdfsN[,-i,drop=FALSE])*(1-cdfsY[,-i,drop=FALSE]),1,prod)
      term1 <- term1+(pdfsY[,i]*(1-cdfsN[,i])*term1b)
      term2c <- apply(cdfsN[,-i,drop=FALSE]*(1-cdfsY[,-i,drop=FALSE]),1,prod)
      term2 <- term2+(pdfsN[,i]*(1-cdfsY[,i])*(term1b-term2c))
    }
    term1+term2
  } else {           # L_NO
    for (i in 1:nStim){
      term2 <- apply(cdfsN[,-i]*(1-cdfsY[,-i]),1,prod)
      term1 <- term1+(pdfsN[,i]*(1-cdfsY[,i])*term2)
    }
    term1
  }  
}


# Checks:

A=0.3; b=1; vT=c(2, 1); vL=c(1, 3.1); svT=1; svL=1; nStim=3

# 1. Integrate to 1?
tp=NA
ps <- c(
  integrate(dEX, lower=0, upper=Inf, A=A, b=b, vT=vT, vL=vL, 
    svT=svT, svL=svL, tp=tp,nStim=nStim,response=1)$value,
  integrate(dEX, lower=0, upper=Inf, A=A, b=b, vT=vT[2:1], vL=vL[2:1],
    svT=svT, svL=svL, tp=tp, nStim=nStim,response=2)$value
)
c(ps,sum(ps))

tp=1
ps <- c(
  integrate(dEX, lower=0, upper=Inf, A=A, b=b, vT=vT, vL=vL, 
    svT=svT, svL=svL, tp=tp,nStim=nStim,response=1)$value,
  integrate(dEX, lower=0, upper=Inf, A=A, b=b, vT=vT[2:1], vL=vL[2:1],
    svT=svT, svL=svL, tp=tp, nStim=nStim,response=2)$value
)
c(ps,sum(ps))

# Plot check
par(mfrow=c(1,2))
to=2 # usually to=NULL

tp=NA
simEX <- rEX(n=1e5, A, b, vT, vL, svT, svL, tp, nStim, t0=0)
ps <- table(simEX$response)/dim(simEX)[1]
dns1 <- density(simEX$rt[simST$response==1],to=to); dns1$y <- dns1$y*ps[1]
dns2 <- density(simEX$rt[simST$response==2],to=to); dns2$y <- dns2$y*ps[2]
plot(dns1,main="EX(non-target)",ylim=c(0,max(c(dns1$y,dns2$y))))
lines(x=dns1$x, y=dEX(t=dns1$x, response=1, A, b, vT, vL, svT, svL, tp, nStim), col="red")
lines(dns2,lty=2)
lines(x=dns2$x, y=dEX(t=dns2$x, response=2, A, b, vT[2:1], vL[2:1], svT, svL, tp, nStim), 
  col="red",lty=2)
legend("topright",c("Yes","No"),lty=1:2,bty="n")

tp=1
simEX <- rEX(n=1e5, A, b, vT, vL, svT, svL, tp, nStim, t0=0)
ps <- table(simEX$response)/dim(simST)[1]
dns1 <- density(simEX$rt[simST$response==1],to=to); dns1$y <- dns1$y*ps[1]
dns2 <- density(simEX$rt[simST$response==2],to=to); dns2$y <- dns2$y*ps[2]
plot(dns1,main="EX(target)",ylim=c(0,max(c(dns1$y,dns2$y))))
lines(x=dns1$x, y=dEX(t=dns1$x, response=1, A, b, vT, vL, svT, svL, tp, nStim), col="red")
lines(dns2,lty=2)
lines(x=dns2$x, y=dEX(t=dns2$x, response=2, A, b, vT[2:1], vL[2:1], svT, svL, tp, nStim), 
  col="red",lty=2)
legend("topright",c("Yes","No"),lty=1:2,bty="n")

# TIMING

# 2e6 random numbers
system.time({
  rEX(n=1e6, A, b, vT, vL, svT, svL, tp=NA, nStim, t0=0)
  rEX(n=1e6, A, b, vT, vL, svT, svL, tp=1, nStim, t0=0)
})
#    user  system elapsed 
#  85.068   1.136  86.219 
 
# Negligible speedup over Peters routine by a factor of
86.074/86.219
# [1] 0.9983182

# 400,000 density evaluations
n <- 1e5
system.time({
  dEX(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=1, A, b, vT, vL, svT, svL, tp=NA, nStim)
  dEX(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=2, A, b, vT, vL, svT, svL, tp=NA, nStim)
  dEX(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=1, A, b, vT, vL, svT, svL, tp=1, nStim)
  dEX(t=seq(min(dns1$x),max(dns1$x),length.out=n), response=2, A, b, vT, vL, svT, svL, tp=1, nStim)
})

# Speedup over Peters routine by a factor of
100*16.087/3.862
# [1] 416.5458
