##################  DMC Stop Signal Model Lesson
#
rm(list=ls())
source ("dmc/dmc.R")
#
# This lesson illustrates the use of 
# 1) waldSSes:  NON-BALLISTIC ASSUMPTION (can stop response production)
#               es is encoding time for stop (must be positive). Note that es
#               is just t0sg interprited and bounded in a particualar way.
# 2) waldSSesN: The same as (1) but with a variable number of accumulators
# 3) waldSSsg:  BALLISTIC ASSUMPTION (cant stop response production), where 
#               sg = stop encoding time - go encoding time. Note if sg < -t0 
#               likelihood is set to zero as this has no process interpritation). 
#               This minimum occurs when response production time and encoding 
#               time for stop = 0 so sg = 0 - t0 
# 4) waldSSsgN: Ths same as (3) but with a variable number of accumulators

# Use a pure Wald model (no start point variability, so A=0, AS=0).


# is.tf <- TRUE
# is.gf <- TRUE
# is.ts <- FALSE
# use.staircase <- TRUE

########################################################## 1)

load_model ("WALD-SS","waldSSes.R") 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),              # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             vS="1",BS="1",AS="1",             # STOP accumulator parameters
             t0="1",es="1",                    # GO non-decision, stop encoding time
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,AS=0,ts=0)) # No start point noise or TRIAL effects
# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "v.true"  "v.false" "B"       "vS"      "BS"      "t0"     
# [7] "es"      "tf"      "gf"     
# 
# Constants are (see attr(,"constants") ):
#  A AS ts 
#  0  0  0 

# Here we set t0 = .2, s-g =-.1 and hence t0sg = .1 
p.vector  <- c(v.true=2.5,v.false=.5,vS=2,B=1,BS=0.5,
               t0=.2,es=.1,tf=.1,gf=.1)

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
n=6
data <- data.model.dmc(simulate.dmc(p.vector,model,
                       n=n,SSD=c(Inf,Inf,.25,.25)),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 
# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("vS",1.5,2.5,p.vector,data)
profile.dmc("BS",.25,.75,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("es",.05,.15,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)

########################################################## 2)

load_model ("WALD-SS","waldSSesN.R") 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2","s3"),SS=c("GO","SS"), # Go stimulus, go and stop trials
               NC=c("a2","a3")),               # 2 choice vs. 3 choice case
  responses <- c("NR","r1","r2","r3"),         # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             vS="1",BS="1",AS="1",             # STOP accumulator parameters
             t0="1",es="1",                    # GO non-decision, stop encoding time
             tf="1",gf="1",                    # trigger failure, go failure
             N="NC",                           # Number of accumulators
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")), # NR mapping ignored
  constants = c(A=0,AS=0,ts=0,N.a2=3,N.a3=4)) # No start point noise or TRIAL effects

# Here we set t0 = .2, s-g =-.1 and hence t0sg = .1 
p.vector  <- c(v.true=3,v.false=.5,vS=2,B=1,BS=0.5,
               t0=.2,es=.1,tf=.1,gf=.1)

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
# Full factorial 12 cells, but remove s3 in a2 for both GO and SS
SSD <- c(rep(Inf,3),rep(.25,3),rep(Inf,3),rep(.25,3)) # Per cell, two ignored
# SSD <- c(rep(Inf,2*6),rep(.25,2*6),rep(Inf,3*6),rep(.25,3*6)) # Per value

# Check likelihood
n <- c(6,6,0,6,6,0,6,6,6,6,6,6)  
data <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=SSD),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,0,2.5e3,2.5e3,0,7.5e3,7.5e3,7.5e3,2.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=SSD),model) 

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2.5,3.5,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("vS",1.5,2.5,p.vector,data)
profile.dmc("BS",.25,.75,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("es",.05,.15,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)

########################################################## 3)

load_model ("WALD-SS","waldSSsg.R") 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2"),SS=c("GO","SS")),# Go stimulus, go and stop trials
  responses = c("NR","r1","r2"),              # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             vS="1",BS="1",AS="1",             # STOP accumulator parameters
             t0="1",sg="1",                    # GO non-decision, stop-go encoding time
             tf="1",gf="1",                    # trigger failure, go failure
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s1="NR")), # NR mapping ignored
  constants = c(A=0,AS=0,ts=0)) # No start point noise or TRIAL effects

# Here we set t0 = .2, s-g =-.1 and hence t0sg = .1 
p.vector  <- c(v.true=2.5,v.false=.5,vS=2,B=1,BS=0.5,
               t0=.2,sg=-.1,tf=.1,gf=.1)

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
n=6
data <- data.model.dmc(simulate.dmc(p.vector,model,
                       n=n,SSD=c(Inf,Inf,.25,.25)),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=c(Inf,Inf,.25,.25)),model) 
# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2,3,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("vS",1.5,2.5,p.vector,data)
profile.dmc("BS",.25,.75,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("sg",-.15,-.05,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)

########################################################## 4)

load_model ("WALD-SS","waldSSsgN.R") 

model <- model.dmc(type="waldss",              # Wald stop signal model
  factors=list(S=c("s1","s2","s3"),SS=c("GO","SS"), # Go stimulus, go and stop trials
               NC=c("a2","a3")),               # 2 choice vs. 3 choice case
  responses <- c("NR","r1","r2","r3"),         # NR=non-response & 2 choices
  p.map=list(v="M",B="1",A="1",                # GO accumulator's parameters 
             vS="1",BS="1",AS="1",             # STOP accumulator parameters
             t0="1",sg="1",                    # GO non-decision, stop-go encoding time
             tf="1",gf="1",                    # trigger failure, go failure
             N="NC",                           # Number of accumulators
             ts="1"),                          # TRIAL covariate slope
  match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")), # NR mapping ignored
  constants = c(A=0,AS=0,ts=0,N.a2=3,N.a3=4)) # No start point noise or TRIAL effects

# Here we set t0 = .2, s-g =-.1 and hence t0sg = .1 
p.vector  <- c(v.true=3,v.false=.5,vS=2,B=1,BS=0.5,
               t0=.2,sg=-.1,tf=.1,gf=.1)

# check.p.vector(p.vector,model)
# print.cell.p(p.vector,model)
# Full factorial 12 cells, but remove s3 in a2 for both GO and SS
SSD <- c(rep(Inf,3),rep(.25,3),rep(Inf,3),rep(.25,3)) # Per cell, two ignored
# SSD <- c(rep(Inf,2*6),rep(.25,2*6),rep(Inf,3*6),rep(.25,3*6)) # Per value

# Check likelihood
n <- c(6,6,0,6,6,0,6,6,6,6,6,6)  
data <- data.model.dmc(simulate.dmc(p.vector,model,n=n,SSD=SSD),model)
likelihood.dmc(p.vector,data)

n <- c(7.5e3,7.5e3,0,2.5e3,2.5e3,0,7.5e3,7.5e3,7.5e3,2.5e3,2.5e3,2.5e3)/10
data <- data.model.dmc(simulate.dmc(p.vector,model,staircase=.05,
    n=n,SSD=SSD),model) 

# SSDs
sort(tapply(as.character(data$SSD),data[,c("SS")],unique)$SS) 
# Probability of non-reponse
tapply(is.na(data$RT),data[,c("SS")],mean)
# Broken down by SSD
tapply(is.na(data$RT),data[,c("SS","SSD")],mean)
# Accuracy
tapply(as.numeric(data$S)==(as.numeric(data$R)-1),data$SS,mean,na.rm=TRUE)["GO"]

# Check the profiles
par(mfrow=c(2,5))
profile.dmc("v.true",2.5,3.5,p.vector,data)
profile.dmc("v.false",.25,0.75,p.vector,data)
profile.dmc("B",.5,1.5,p.vector,data)
profile.dmc("vS",1.5,2.5,p.vector,data)
profile.dmc("BS",.25,.75,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("sg",-.15,-.05,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)
