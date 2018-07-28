# Stop-signal context independent parameterization n-choice Wald, 
# External parameters types: v, B, A, t0, vS, BS, AS, t0sg, tf, gf, ts
# Internal parameters types: v, B, A, t0, vS, BS, AS, t0sg, tf, gf, ts
#
# NB1: st0 is not available
# NB2: ts is slope of TRIALS covariate on B (NB: Threshold = B+A)
# 
# Given g = go encoding time, s = stop encoding time and 
#       r = go response production time can only identify two functions of these:
# t0 = g+r (identified by go RTs)
# Assuming a BALLISTIC (cant be stopped) response production stage
# t0sg = t0+s-g (identified by the probability of stopping, adding t0 keeps
#        always positive, zero when r=0,s=0,g=t0).
# Assuming NON-BALLISTIC (can be stopped) response production stage
# t0sg = s-g (identified by the probability of stopping, must be positive).
#
# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all appliciatons

my.integrate <- function(...,big=10)
# Avoids but in integrate upper=Inf that uses only 1  subdivision
# Use of  big=10 is arbitary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (class(out)=="try-error") 0 else 
  {
    if (out$subdivisions==1) 
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (class(out)=="try-error") 0 else
      {
         if (out$subdivisions==1) 0 else out$value   
      }
    } else out$value
  }
}

# source("rtdists_extras.R")


transform.dmc <- function(par.df) 
# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used is present.
{
  # Context independence: seperate go and stop accumulator parameterization.
  par.df["NR",c("v","B","A")] <- par.df["NR",c("vS","BS","AS")]
  par.df[,c("v","B","A","vS","BS","AS","t0","t0sg","tf","gf","ts")]
}

random.dmc<- function(n,p.df,model,SSD=Inf,staircase=NULL,TRIALS=NULL)
{
  rWaldss(n,v=p.df$v,B=p.df$B,A=p.df$A,t0=p.df$t0[1],t0sg=p.df$t0sg[1],
  tf=p.df$tf[1],gf=p.df$gf[1],ts=p.df$ts[1],TRIALS=TRIALS,
  SSD=SSD,staircase=staircase)
}

likelihood.dmc <- function(p.vector,data,min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
{
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    if ( is.null(data$TRIALS[attr(data,"cell.index")[[i]]]) )
        TRIALS <- NA else TRIALS <- data$TRIALS[attr(data,"cell.index")[[i]]]
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      n1PDF.Waldss(rt=data$RT[attr(data,"cell.index")[[i]]],
          v=p.df$v,
          B=p.df$B,
          A=p.df$A,
          t0=p.df$t0[1], 
          t0sg=p.df$t0sg[1],
          tf=p.df$tf[1],
          gf=p.df$gf[1], 
          ts=p.df$ts[1],
          # Stop-signal delays
          SSD=data$SSD[attr(data,"cell.index")[[i]]],
          # TRIAL regression
          TRIALS=TRIALS,# In case no TRIALS
          # Index of stop signal accumulator
          Si=c(1:dim(p.df)[1])[row.names(p.df)=="NR"]
      )
 }
 pmax(likelihood,min.like)
}


