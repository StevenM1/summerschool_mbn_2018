rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("MTR-LBA","lba_DK.R")

# Note that responses2 must have a second value called "DK". The first value can
# be anything including empty as here.
model <- model.dmc(p.map=list(A="1",B="SA",D="SA",mean_v=c("D","M"),sd_v="1",t0="SA",st0="1"),
  constants=c(st0=0,sd_v=1),match.map=list(M=list(left="LEFT",right="RIGHT")),
  factors=list(S=c("left","right"),D=c("easy","hard"),SA=c("accuracy","speed")),
  responses=c("LEFT","RIGHT"),responses2=c("","DK"),type="normDK")
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "A"                 "B.accuracy"        "B.speed"          
#  [4] "D.accuracy"        "D.speed"           "mean_v.easy.true" 
#  [7] "mean_v.hard.true"  "mean_v.easy.false" "mean_v.hard.false"
# [10] "t0.accuracy"       "t0.speed"         
# 
# Constants are (see attr(,"constants") ):
#  st0 sd_v 
#    0    1 
# 
# Model type = normDK 

p.vector  <- c(A=.5,D.accuracy=.5,D.speed=.25,B.accuracy=.5,B.speed=.25,
  mean_v.hard.true=1,mean_v.hard.false=.5,mean_v.easy.true=1.5,mean_v.easy.false=.75,
  t0.speed=.2,t0.accuracy=.3)
data <- simulate.dmc(p.vector,model,n=1e4)
data.model <- data.model.dmc(data,model)

# Likelihood works
head(likelihood.dmc(p.vector,data.model))


# R2 colmn indicates "DK" for dont know responses but is otherwise blank 
head(data)

# some checks on DK behaviour 
tmp <- table(data$R,data$S,data$R2)
# DK, more likley for errors 
round(100*tmp[,,"DK"]/apply(tmp,1:2,sum))

# To plot collapse DK over respone
data$R <- as.character(data$R); data$R[data$R2=="DK"] <- "DK"
data$R <- factor(data$R,levels=c("LEFT","DK","RIGHT"))

is.in <- data$D=="easy" & data$SA=="accuracy" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")

### Do sampling.

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","beta","beta"),
  p1=c(A=.3,B.accuracy=.3,B.speed=.3,D.accuracy=.3,D.speed=.3,
    mean_v.hard.true=1,mean_v.hard.false=0,mean_v.easy.true=1,mean_v.easy.false=0,
    t0.accuracy=1,t0.speed=1),                           
  p2=c(1,1,1,1,1,3,3,3,3,1,1),
  lower=c(0,0,0,0,0,NA,NA,NA,NA,.1,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,1,1)
)
par(mfcol=c(2,6)); for (i in names(p.prior.real)) plot.prior(i,p.prior.real)

samples <- samples.dmc(nmc=100,p.prior,data=data.model)
samples <- run.unstuck.dmc(samples,p.migrate=.05,verbose=TRUE,cores=11)
samples1 <- run.converge.dmc(samples,minN=500,nmc=50,verbose=TRUE,max.try=20,cores=11)

# Chains look OK
plot.dmc(samples1,layout=c(2,6))

# Good recovery
check.recovery.dmc(samples1,p.vector)
#                   A D.accuracy D.speed B.accuracy B.speed mean_v.hard.true mean_v.hard.false
# True           0.50       0.50    0.25       0.50    0.25             1.00              0.50
# 2.5% Estimate  0.45       0.47    0.22       0.50    0.25             0.98              0.49
# 50% Estimate   0.52       0.52    0.26       0.51    0.26             1.03              0.55
# 97.5% Estimate 0.57       0.58    0.31       0.52    0.27             1.08              0.60
# Median-True    0.02       0.02    0.01       0.01    0.01             0.03              0.05
#                mean_v.easy.true mean_v.easy.false t0.speed t0.accuracy
# True                       1.50              0.75     0.20        0.30
# 2.5% Estimate              1.50              0.75     0.18        0.28
# 50% Estimate               1.54              0.80     0.20        0.29
# 97.5% Estimate             1.57              0.85     0.21        0.31
# Median-True                0.04              0.05     0.00       -0.01

# Updating weakest on A and D. but OK.
plot.dmc(samples1,p.prior=p.prior,layout=c(2,6),show.obs=FALSE)

# Suffers from some very strong correlations (particularly D.speed and 
# D.accruacy, also these with t0 and A, and among the v parameters)
pairs.dmc(samples1)


save_data(data,data.model,samples,samples1,file="dmc_5_DK.RData")


