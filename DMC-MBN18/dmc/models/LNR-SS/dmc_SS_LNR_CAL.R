### An LNR Stop-signal example with TRIGGER and GO FAILURE
# Context Independent Parameterization (i.e., seperate meanlogS and sdlogS for stop).
# On some proprotion of trials, tf, the stop signal is ignored.
#
# NB1: Unlike lnr st0 is not available
#
# NB2: given g = go encoding time, s = stop encoding time and 
#      r = go response production time, can only identify two functions of these:
# t0 = g+r (identified by go RTs)
# t0sg = t0+s-g (identified by the probability of stopping, adding t0 keeps
#        always positive, zero when r=0,s=0,g=t0).
#
# This lesson illustrates the use of 
# 1) lnrSSes:  NON-BALLISTIC ASSUMPTION (can stop response production)
#              es is encoding time for stop (must be positive). Note that es
#              is just t0sg interprited and bounded in a particualar way.
# 2) lnrSSesN: The same as (1) but with a variable number of accumulators
# 3) lnrSSsg:  BALLISTIC ASSUMPTION (cant stop response production), where 
#              sg = stop encoding time - go encoding time. Note if sg < -t0 
#              likelihood is set to zero as this has no process interpritation). 
#              This minimum occurs when response production time and encoding 
#              time for stop = 0 so sg = 0 - t0 
# 4) lnrSSsgN: Ths same as (3) but with a variable number of accumulators

rm(list=ls())
source ("dmc/dmc.R")

# is.tf <- TRUE
# is.gf <- TRUE
# is.ts <- FALSE
# use.staircase <- TRUE


load_model ("LNR-SS","lnrSSsgCAL.R") 

model <- model.dmc(type="lnrss",               # LNR stop signal model
  factors=list(S=c("s1","s2","s3","s4"),       # Go stimulus, two SG stimuli
               SS=c("GO","SS","SG","GG")),     # Go and stop, simple, choice 
  responses <- c("NR","r1","r2"),              # NR=non-response & 2 choices
  p.map=list(meanlog="M",sdlog="M",            # GO accumulator's parameters 
             meanlogS="1",sdlogS="1",          # STOP accumulator parameters
             meanlogGG="M",sdlogGG="M",        # Choice block parameters
             meanlogSG="1",sdlogSG="1",        # Simple block parameters
             t0SG,t0GG,                        # simple and GO/GG non-decision
             sg="1",                           # Stop encoding time
             tfSS,                             # Trigger failure for stop signal, 
             gf="1",                           # Go failure
             N1="SS",                          # Start number of accumulators
             N2="SS",                          # End number of accumulators
             ts="1"),                          # TRIAL covariate slope
  cvs=c("isSS","isSG","isGG"),
  match.map=list(M=list(s1="r1",s2="r2",s3="r3",s1="NR")), # NR mapping ignored
  constants = c(ts=0,sg=0,                      # No TRIAL, sg computed
                N1.GO=1,N2.GO=3,N1.SS=1,N2.SS=3,
                N1.SG=1,N2.SG=1,N1.GG=2,N2.GG=3))

p.vector  <- c(meanlog.true=-.5,meanlog.false=1,meanlogS=-1,
               sdlog.true=0.5,sdlog.false=1,sdlogS=0.75,
               t0=.2,sg=-.1,gf=.1,tf=.1) 

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
profile.dmc("meanlog.true",-1,1,p.vector,data)
profile.dmc("meanlog.false",0,2,p.vector,data)
profile.dmc("meanlogS",-2,0,p.vector,data)
profile.dmc("sdlog.true",.4,.6,p.vector,data)
profile.dmc("sdlog.false",.75,1.25,p.vector,data)
profile.dmc("sdlogS",.5,1,p.vector,data)
profile.dmc("t0",.1,.3,p.vector,data)
profile.dmc("sg",-.15,0,p.vector,data)
profile.dmc("gf",.05,.15,p.vector,data)
profile.dmc("tf",.05,.15,p.vector,data)


