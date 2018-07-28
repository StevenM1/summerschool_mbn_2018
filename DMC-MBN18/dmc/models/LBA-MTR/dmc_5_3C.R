rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA-MTR","lba_3C.R")

# Note that responses2 must have values 1 (low confidence) and 2 (high confidence).
model <- model.dmc(p.map=list(A="1",B="SA",D1="SA",D2="SA",p_flc="SA", mean_v=c("D","M"),sd_v="1",t0="SA",st0="1"),
  constants=c(st0=0,sd_v=1),match.map=list(M=list(left="LEFT",right="RIGHT")),
  factors=list(S=c("left","right"),D=c("easy","hard"),SA=c("accuracy","speed")),
  responses=c("LEFT","RIGHT"),responses2=c("1","2"),type="norm3C")
# Parameter vector names are: ( see attr(,"p.vector") )
#  [1] "A"                 "B.accuracy"       
#  [3] "B.speed"           "D.accuracy"       
#  [5] "D.speed"           "mean_v.easy.true" 
#  [7] "mean_v.hard.true"  "mean_v.easy.false"
#  [9] "mean_v.hard.false" "t0.accuracy"      
# [11] "t0.speed"         
# 
# Constants are (see attr(,"constants") ):
#  st0 sd_v 
#    0    1 
# 
# Model type = norm2C 


p.vector  <- c(A=.5,D1.accuracy=.5,D1.speed=.25,D2.accuracy=.5,D2.speed=.25,
               p_flc.speed=.1,p_flc.accuracy=0, B.accuracy=.5,B.speed=.25,
  mean_v.hard.true=1,mean_v.hard.false=.5,mean_v.easy.true=1.5,mean_v.easy.false=.75,
  t0.speed=.2,t0.accuracy=.3)
data <- simulate.dmc(p.vector,model,n=1e4)
data.model <- data.model.dmc(data,model)

# Likelihood works
head(likelihood.dmc(p.vector,data.model))
length(likelihood.dmc(p.vector,data.model))

# Check profiles
par(mfrow=c(2,3))
profile.dmc("A",               .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("B.accuracy",               .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("p_flc.accuracy",               0,  .5,p.vector,data.model,ylim=NA)
profile.dmc("p_flc.speed",               0,  .5,p.vector,data.model,ylim=NA)

profile.dmc("mean_v.hard.true",  .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.easy.false", .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("D1.speed",    .1,  1,p.vector,data.model,ylim=NA)
profile.dmc("D2.speed",    .1,  1,p.vector,data.model,ylim=NA)
profile.dmc("D2.accuracy",    .1,  1,p.vector,data.model,ylim=NA)
profile.dmc("t0.accuracy",             .01,.4,p.vector,data.model,ylim=NA)

#Profiles look good.

# R2 colmn indicates "1" for low and "2" for high confidence 
head(data)

# some checks on C behaviour 
tmp <- table(data$R,data$S,data$R2)
# low confidnece is more likley for errors 
round(100*tmp[,,"1"]/apply(tmp,1:2,sum))

# To plot map hi lo lo hi
is.1 <- data$R=="LEFT" & data$R2==3
is.2 <- data$R=="LEFT" & data$R2==2
is.3 <- data$R=="LEFT" & data$R2==1
is.4 <- data$R=="RIGHT" & data$R2==1
is.5 <- data$R=="RIGHT" & data$R2==2
is.6 <- data$R=="RIGHT" & data$R2==3



data$R <- as.numeric(data$R)
data$R[is.1] <- 1
data$R[is.2] <- 2
data$R[is.3] <- 3
data$R[is.4] <- 4
data$R[is.5] <- 5
data$R[is.6] <- 6
data$R <- factor(data$R,labels=c("LEFThi","LEFTmed", "LEFTlo","RIGHTlo","RIGHTmed", "RIGHThi"))

tapply(data$RT,data$R,quantile,probs=c(.1,.3,.5,.7,.9))

# $LEFThi
# 10%       30%       50%       70%       90% 
# 0.5737420 0.7725961 0.9855159 1.2704909 1.9865876 
# 
# $LEFTmed
# 10%       30%       50%       70%       90% 
# 0.6388880 0.8725386 1.1012628 1.3937252 2.1186934 
# 
# $LEFTlo
# 10%       30%       50%       70%       90% 
# 0.6904364 0.9638853 1.2211972 1.5426796 2.3837170 
# 
# $RIGHTlo
# 10%       30%       50%       70%       90% 
# 0.6919002 0.9580817 1.2127155 1.5382075 2.3453741 
# 
# $RIGHTmed
# 10%       30%       50%       70%       90% 
# 0.6405807 0.8761171 1.0973017 1.3843727 2.1084720 
# 
# $RIGHThi
# 10%       30%       50%       70%       90% 
# 0.5692835 0.7690748 0.9826878 1.2582866 1.9637624 


dat<-matrix(unlist(tapply(data$RT,data$R,quantile,probs=c(.1,.3,.5,.7,.9))),ncol=6)
     


rownames(dat)<-c(.1,.3,.5,.7,.9)
colnames(dat)<-c("LEFThi","LEFTmed", "LEFTlo","RIGHTlo","RIGHTmed", "RIGHThi")
boxplot(dat)     ##Kind of shows what we want.

par(mfrow=c(1,1))
plot(dat[1,],type="b",ylim=c(min(dat),max(dat)),xlab="Response",ylab="RT",
main="RT quantiles by Response  ",sub=".10,.30,.50,.70,.90 quantiles")
for (i in 2:5)
{lines(dat[i,],col=1)
  points(dat[i,],col=1)}

##Think I need to use ggplot to label the responses.



is.in <- data$D=="easy" & data$SA=="accuracy" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")

### Do sampling.

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","beta","beta"),
  p1=c(A=.3,B.accuracy=.3,B.speed=.3,D1.accuracy=.3,D1.speed=.3,D2.accuracy=.3,D2.speed=.3,
    mean_v.hard.true=1,mean_v.hard.false=0,mean_v.easy.true=1,mean_v.easy.false=0,
    t0.accuracy=1,t0.speed=1),                           
  p2=c(1,1,1,1,1,1,1,3,3,3,3,1,1),
  lower=c(0,0,0,0,0,0,0,NA,NA,NA,NA,.1,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1,1)
)
par(mfcol=c(2,7)); for (i in names(p.prior)) plot.prior(i,p.prior)

samples <- samples.dmc(nmc=300,p.prior,data=data.model)
samples <- run.unstuck.dmc(samples,p.migrate=.05,verbose=FALSE,cores=13)
samples1 <- run.converge.dmc(samples,minN=500,nmc=50,verbose=TRUE,
  max.try=20,cores=39)
# [1] "Final multivariate psrf = 1.09585102305701"
# Effective sample size
#                 A        B.accuracy           B.speed        D.accuracy           D.speed 
#               612               604               614               603               618 
#  mean_v.easy.true  mean_v.hard.true mean_v.easy.false mean_v.hard.false       t0.accuracy 
#               553               573               562               601               612 
#          t0.speed 
#               630 
              
# Chains look OK
plot.dmc(samples1,layout=c(2,6))

# Good recovery
check.recovery.dmc(samples1,p.vector)
#                   A D.accuracy D.speed B.accuracy B.speed mean_v.hard.true
# True           0.50       0.50    0.25       0.50    0.25             1.00
# 2.5% Estimate  0.43       0.46    0.22       0.49    0.25             0.98
# 50% Estimate   0.50       0.51    0.26       0.50    0.25             1.03
# 97.5% Estimate 0.55       0.57    0.32       0.52    0.26             1.08
# Median-True    0.00       0.01    0.01       0.00    0.00             0.03
#                mean_v.hard.false mean_v.easy.true mean_v.easy.false t0.speed t0.accuracy
# True                        0.50             1.50              0.75     0.20        0.30
# 2.5% Estimate               0.47             1.47              0.72     0.18        0.28
# 50% Estimate                0.53             1.51              0.77     0.20        0.30
# 97.5% Estimate              0.59             1.55              0.83     0.21        0.31
# Median-True                 0.03             0.01              0.02     0.00        0.00

# All nicely updated 
plot.dmc(samples1,p.prior=p.prior,layout=c(2,6),show.obs=FALSE)

# 
pairs.dmc(samples1)


save_data(data,data.model,samples,samples1,file="dmc_5_2C.RData")


