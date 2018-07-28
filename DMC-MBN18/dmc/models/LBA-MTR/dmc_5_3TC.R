rm(list=ls()) 

# Current working directory must be set to the top-level folder  
# containing the dmc and tutorial subfolders 
source ("dmc/dmc.R")
load_model ("LBA-MTR","lba_3TC.R")

# Note that responses2 must have values 1 (low confidence) and 2 (high confidence).
model <- model.dmc(p.map=list(A="1",B="SA",D1="SA",D2="SA",p_flc="SA",mean_v=c("D","M"),sd_v="1",t0="SA",st0="1"),
                   constants=c(st0=0,sd_v=1),match.map=list(M=list(left="LEFT",right="RIGHT")),
                   factors=list(S=c("left","right"),D=c("easy","hard"),SA=c("accuracy","speed")),
                   responses=c("LEFT","RIGHT"),responses2=c("1","2"),type="norm3TC")
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
               B.accuracy=.5,B.speed=.25,p_flc.speed=.1,p_flc.accuracy=0,
               mean_v.hard.true=1,mean_v.hard.false=.5,mean_v.easy.true=1.5,mean_v.easy.false=.75,
               t0.speed=.2,t0.accuracy=.3)
data <- simulate.dmc(p.vector,model,n=1e3)
data.model <- data.model.dmc(data,model)

# Likelihood works
head(likelihood.dmc(p.vector,data.model))


# Check profiles
par(mfrow=c(2,3))
profile.dmc("A",               .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("B.accuracy",               .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.hard.true",  .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("mean_v.easy.false", .1,  2,p.vector,data.model,ylim=NA)
profile.dmc("p_flc.speed", 0,  .5,p.vector,data.model,ylim=NA)
profile.dmc("p_flc.accuracy", 0,  .5,p.vector,data.model,ylim=NA)
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
# , ,  = 1
# 
# 
# left right
# LEFT   7206  5737
# RIGHT  5718  7214
# 
# , ,  = 2
# 
# 
# left right
# LEFT   4408  2940
# RIGHT  2872  4334
# 
# , ,  = 3
# 
# 
# left right
# LEFT  13590  6214
# RIGHT  6206 13561

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
# 0.5728004 0.7651145 0.9823042 1.2573125 1.9659448 
# 
# $LEFTmed
# 10%       30%       50%       70%       90% 
# 0.5569913 0.7576259 0.9480723 1.1971229 1.7979855 
# 
# $LEFTlo
# 10%       30%       50%       70%       90% 
# 0.5547452 0.7501902 0.9349991 1.1586043 1.7350080 
# 
# $RIGHTlo
# 10%       30%       50%       70%       90% 
# 0.5584179 0.7583267 0.9379532 1.1600406 1.7448519 
# 
# $RIGHTmed
# 10%       30%       50%       70%       90% 
# 0.5575534 0.7554769 0.9473176 1.1912423 1.8342338 
# 
# $RIGHThi
# 10%       30%       50%       70%       90% 
# 0.5726832 0.7733311 0.9903970 1.2711019 1.9541322 


dat<-matrix(unlist(tapply(data$RT,data$R,quantile,probs=c(.1,.3,.5,.7,.9))),ncol=6)

rownames(dat)<-c(.1,.3,.5,.7,.9)
colnames(dat)<-c("LEFThi","LEFTmed", "LEFTlo","RIGHTlo","RIGHTmed", "RIGHThi")

par(mfrow=c(1,1))
plot(dat[1,],type="b",ylim=c(min(dat,.5),max(dat)),xaxt="n",xlab="Response",ylab="RT",
     main="RT quantiles by Response  ",sub=".10,.30,.50,.70,.90 quantiles")
for (i in 2:5)
{lines(dat[i,],col=1)
  points(dat[i,],col=1)}
axis(1,c(1,2,3,4,5,6))

##So threshold count can give you faster low confidence responses.
##

dat<-matrix(unlist(tapply(data[data$D=="easy",]$RT,data[data$D=="easy",]$R,quantile,probs=c(.1,.3,.5,.7,.9))),ncol=6)

rownames(dat)<-c(.1,.3,.5,.7,.9)
colnames(dat)<-c("LEFThi","LEFTmed", "LEFTlo","RIGHTlo","RIGHTmed", "RIGHThi")

par(mfrow=c(1,1))
plot(dat[1,],type="b",ylim=c(min(dat,.5),max(dat)),xlab="Response",ylab="RT",
     main="RT quantiles by Response  ",sub=".10,.30,.50,.70,.90 quantiles")
for (i in 2:5)
{lines(dat[i,],col=1)
  points(dat[i,],col=1)}





par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left",],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right",],xlim=c(0,4),main="RIGHT=correct")


is.in <- data$D=="easy" & data$SA=="accuracy" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")

{
is.in <- data$D=="easy" & data$SA=="speed" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")


is.in <- data$D=="hard" & data$SA=="accuracy" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")

is.in <- data$D=="hard" & data$SA=="speed" # Change to look at other three cells.
par(mfrow=c(1,2))
plot.cell.density(data.cell=data[data$S=="left" & is.in,],xlim=c(0,4),main="LEFT=correct")
plot.cell.density(data.cell=data[data$S=="right" & is.in,],xlim=c(0,4),main="RIGHT=correct")

}
### Do sampling.

p.prior <- prior.p.dmc(
  dists = c("tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","tnorm","beta","beta"),
  p1=c(A=.3,B.accuracy=.3,B.speed=.3,D1.accuracy=.3,D1.speed=.3,D2.accuracy=.3,D2.speed=.3,
       p_flc.speed=0,p_flc.accuracy=0,mean_v.hard.true=1,mean_v.hard.false=0,mean_v.easy.true=1,mean_v.easy.false=0,
       t0.accuracy=1,t0.speed=1),                           
  p2=c(1,1,1,1,1,1,1,.5,.5,3,3,3,3,1,1),
  lower=c(0,0,0,0,0,0,0,0,0,NA,NA,NA,NA,.1,.1),upper=c(NA,NA,NA,NA,NA,NA,NA,1,1,NA,NA,NA,NA,1,1)
)
par(mfcol=c(2,5)); for (i in names(p.prior)) plot.prior(i,p.prior)

samples <- samples.dmc(nmc=300,p.prior,data=data.model)
samples <- run.unstuck.dmc(samples,p.migrate=.05,verbose=FALSE,cores=1)
samples1 <- run.converge.dmc(samples,minN=500,nmc=50,verbose=TRUE,
                             max.try=20,cores=1)


plot.dmc(samples1,layout=c(2,6))

# Good recovery
check.recovery.dmc(samples1,p.vector)

plot.dmc(samples1,p.prior=p.prior,layout=c(2,6),show.obs=FALSE)

# 
pairs.dmc(samples1)


save_data(data,data.model,samples,samples1,file="dmc_5_3TC.RData")



##Let's look at a simpler model.

model <- model.dmc(p.map=list(A="1",B="1",D1="1",D2="1",mean_v="M",sd_v="M",t0="1",st0="1"),
                   constants=c(st0=0,sd_v.false=1),match.map=list(M=list(left="LEFT",right="RIGHT")),
                   factors=list(S=c("left","right")),
                   responses=c("LEFT","RIGHT"),responses2=c("1","2"),type="norm3TC")



p.vector  <- c(A=.3,D1=.6,D2=.3,
               B=.1,mean_v.true=1,mean_v.false=.7,
               sd_v.true=.6,t0=.2)
data <- simulate.dmc(p.vector,model,n=1e5)
data.model <- data.model.dmc(data,model)


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


dat<-matrix(unlist(tapply(data$RT,data$R,quantile,probs=c(.1,.3,.5,.7,.9))),ncol=6)

rownames(dat)<-c(.1,.3,.5,.7,.9)
colnames(dat)<-c("LEFThi","LEFTmed", "LEFTlo","RIGHTlo","RIGHTmed", "RIGHThi")

par(mfrow=c(1,1))
plot(dat[1,],type="b",ylim=c(min(dat,.5),max(dat)),xlab="Response",ylab="RT",
     main="RT quantiles by Response  ",sub=".10,.30,.50,.70,.90 quantiles")
for (i in 2:5)
{lines(dat[i,],col=1)
  points(dat[i,],col=1)}






