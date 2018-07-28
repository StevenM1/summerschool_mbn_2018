
require(rtdists) 

#Simulate correlated rate LBA
rlba.cnorm <- function (n,A,b,t0,mean_v,sd_v,st0=0,corr_v=0,return.ttf=FALSE) 
{
  if (any(b < A)) 
    stop("b cannot be smaller than A!")
  n_v  <- ifelse(is.null(dim(mean_v)), length(mean_v), dim(mean_v)[1])
  D=diag(sd_v)
  Cor=matrix(nrow=n_v,ncol=n_v)
  diag(Cor)<-1
  Cor[lower.tri(Cor)]<-corr_v[1]
  Cor<-t(Cor)
  Cor[lower.tri(Cor)]<-corr_v[1]
  Sigma=matrix(nrow=n,ncol=n)
  Sigma=D%*%Cor%*%D
  drifts <- t(rtmvnorm(n = n, mean = mean_v, sigma = Sigma,lower=c(0,rep(-Inf,(n_v-1)))))
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


n1PDFcorrelated <- function (rt, A, b, t0, mean_v,sd_v,corr, st0 = 0, distribution = c("norm","gamma", "frechet", "lnorm"), 
                             args.dist = list(), silent = FALSE) 
{
  dots <- list(...)
  if (is.null(names(dots))) 
    stop("... arguments need to be named.")
  if (any(names(dots) == "")) 
    stop("all ... arguments need to be named.")
  n_v <- max(vapply(dots, length, 0))
  if (!silent) 
    message(paste("Results based on", n_v, "accumulators/drift rates."))
  if (n_v < 2) 
    stop("There need to be at least two accumulators/drift rates.")
  #distribution <- match.arg(distribution)
  nn <- length(rt)
  A <- check_n1_arguments(A, nn = nn, n_v = n_v)
  b <- check_n1_arguments(b, nn = nn, n_v = n_v)
  t0 <- check_n1_arguments(t0, nn = nn, n_v = n_v)
  st0 <- rep(unname(st0), length.out = nn)
  
  pdf <- dlba_norm_core
  cdf <- plba_norm_core
  if (any(!(c("mean_v", "sd_v") %in% names(dots)))) stop("mean_v and sd_v need to be passed for distribution = \"norm\"")
  dots$mean_v <- check_n1_arguments(dots$mean_v, nn = nn, 
                                    n_v = n_v, dots = TRUE)
  dots$sd_v <- check_n1_arguments(dots$sd_v, nn = nn, n_v = n_v, 
                                  dots = TRUE)
  dots <- dots[c("mean_v", "sd_v")]
  
  for (i in seq_len(length(dots))) {
    if (length(dots[[i]]) < n_v) 
      dots[[i]] <- rep(dots[[i]], length.out = n_v)
  }
  do.call(n1PDF_core, args = c(rt = list(rt), A = list(A), 
                               b = list(b), t0 = list(t0), st0 = list(st0), dots, pdf = pdf, 
                               cdf = cdf, args.dist = list(args.dist)))
}


###


