# Stop-signal context independent parameterization n-choice LNR, N accumulators
# External parameters types: meanlog/S/SG/GG, sdlog/S/SG/GG, t0, sg, tf, gf, ts, N1, N2
# Internal parameters types: meanlog, sdlog, t0, t0sg, tf, gf, ts

# VARIANT of lnrSS.R:  BALLISTIC ASSUMPTION (cant stop response production)
#  This version allows sg to be direclty set (but note if sg < -t0 likelihood
#  is set to zero as this has no process interpritation). This minimum occurs
#  when response production time and encoding time for stop = 0 so sg = 0 - t0

# CALIBRATION paradigm SS factor has 4 levels: 
#  GO, SS: as in ordinary stop-signal paradigm
#  SG: simple RT to only stop signal
#  GG: choice RT with no SS trials
# 
# PARAMTER CONSTRAINTS
# 
# tf=0 for SG and GG
# For all but SG meanlog/sdlog.true/false estimated as usual. For SG no errors
# so fix meanlog.false to large value/sdlog.false to small value
# ter~SS, all but ter.SS sampled (ter.SS set to arbitary constant) 
# sg = ter.SG - ter.GG for SS
# 
# Hence must set constants tf.SG=tf.GG=0 (and not sampled), as is conventional
# with LNR2 models meanlog.false=1e6,sdlog.false=1e-10, and 
# ter.SS = arbitray, and sg=aribtary and replace with difference of ter samples. 
#
# Also becasue must assume that meanlog and sdlog vary with SS, must then
# equate them over GO and SS
# 
# Assume 4 levels of stimulus, s1, s2, two are choice stimuli present in all but
# SG, and s3, s4, specific to SG and mapped to r1 and r2 respectively, really
# the same stop signal stimulus in both cases, but corresponding to blocks where
# it is mapped to r1 and r2. 
# 
# As usuall map s1=r1, s2=r2, s1=NR (the usual dummy mapping for stop-signal
# trials) and also s3=r1 and s4=r2 for SG condition.
# 
# In simulation/data not s1 or s2 in SG
# 
# Define two "N" parameters N1 = start accumulator and N2 = end accumulator, 
# used to pick out approriate race architecture for each condition:
# GO, SS : use N1=1, N2=3 (i.e., stop and two go accumualators)
# GG: N1=2, N2=3 (i.e., two go accumulators)
# SG: N1=1, N2=1 (i.e., only stop accumulator)


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
  par.df$t0sg <- par.df$sg + par.df$t0
  par.df[,c("meanlog","sdlog","meanlogS","sdlogS","meanlogSG","sdlogSG",
            "meanlogGG","sdlogGG","t0","t0sg","tf","gf","ts","N")]
}

random.dmc<- function(n,p.df,model,SSD=Inf,staircase=NULL,TRIALS=NULL)
{
  rlnrss(n,meanlog=p.df$meanlog[p.df$N1[1]:p.df$N1[2]],
           sdlog=p.df$sdlog[p.df$N1[1]:p.df$N1[2]],
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
          meanlog=p.df$meanlog[p.df$N1[1]:p.df$N1[2]],
          sdlog=p.df$sdlog[p.df$N1[1]:p.df$N1[2]],
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


