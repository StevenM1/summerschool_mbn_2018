# Stop-signal context independent parameterization n-choice LNR, 
# External parameters types: meanlog, sdlog, t0, meanlogS, sdlogS, es, tf, gf, ts
# Internal parameters types: meanlog, sdlog, t0, t0sg, tf, gf, ts

# VARIANT of lnrSS.R:  NON-BALLISTIC ASSUMPTION (can stop response production)
#  es is encoding time for stop (must be positive) 
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
  par.df["NR",c("meanlog","sdlog")] <- par.df["NR",c("meanlogS","sdlogS")]
  par.df$t0sg <- par.df$es
  par.df[,c("meanlog","sdlog","meanlogS","sdlogS","t0","t0sg","tf","gf","ts")]
}

random.dmc<- function(n,p.df,model,SSD=Inf,staircase=NULL,TRIALS=NULL)
{
  rlnrss(n,meanlog=p.df$meanlog,sdlog=p.df$sdlog,
         t0=p.df$t0[1],t0sg=p.df$t0sg[1],tf=p.df$tf[1],gf=p.df$gf[1],
         ts=p.df$ts[1],TRIALS=TRIALS,SSD=SSD,staircase=staircase)
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
      n1PDF.lnrss(rt=data$RT[attr(data,"cell.index")[[i]]],
          meanlog=p.df$meanlog,
          sdlog=p.df$sdlog,
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


