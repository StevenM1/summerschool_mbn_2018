
# Here is the update package.  Specifically, in 'rcddm1.cpp'
# 
#   - Change unit circle to scaled circle. That is, in rcddm1.cpp, the kappa (k)
#   parameter in von Mise distribution is scaled by a step size (ie p in the
#   function argument); hence changing code from
#   'double k = std::sqrt(pVec[2]*pVec[2] + pVec[1]*pVec[1]) / pVec[4];' to
#   'double k = p * std::sqrt(pVec[2]*pVec[2] + pVec[1]*pVec[1]) / pVec[4];'.
# 
# Note there is no calculation of k (ie von Mise's kappa) in rcddm2_internal.cpp, 
# because as Peter's email explains 'the angles give both the magnitude and angle information ...'.


  
## dcddm example
require(CircularDDM)

# Make some RTs and Response (radians, angle = radians*180/pi)
x <- cbind(
RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
)
# a=threshold, x and y drift rates, non decision time, within trial variability
pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)

# Likeliood  (one column matrix)
dcddm(x, pVec)

## rcddm example
pVec <- c(a=2, vx=1.5, vy=1.25, t0=.25, s=1)
# 2 column matrix RT, angle in radians.
den  <- rcddm(1e3, pVec)
par(mfrow=c(1,2))
hist(den[,1], breaks = "fd", xlab="Response Time",  main="Density")
hist(den[,2], breaks = "fd", xlab="Response Angle", main="Density")
