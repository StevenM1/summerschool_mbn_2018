# Template setup for 2-choice LBA 2 Confidence choices, identical to DK,
# except different pdf.
# Two thresholds, b (upper) and d (lower). 
#   External parameters types: A, B, D, t0, mean_v, sd_v, st0 = 0 (optional)
#   Internal parameters types: A, b, d, t0, mean_v, sd_v, st0 = 0 (optional)

# User edited functions for the DMC (Dynamic Models of Choice)
#    Source in all applications

# source("rtdists_extras.R")


# This function transfroms parameters to a form suitbale for the model 
#   being used. Called inside of get.par.mat. 
# "par.df" is a data frame of parameters types , some of which may need to be 
#   transformed, or new columns created, so that the full set of internal 
#   parameter types, specified in "type.par.names", required by the type of 
#   evidence accumulation model being used ("norm" etc.) is present.
transform.dmc <- function(par.df) 
{
  # User supplied tranforms go here
  par.df$d <- par.df$D+par.df$A
  par.df$b <- par.df$B+par.df$d
  
#   # COMMENT OUT this check for speed after debugging
#   if ( !all(type.par.names %in% names(par.df)) )
#     stop("Trasform has not created parameter(s) required by the model.")
  
  par.df[,c("A","b","d","p_flc","t0","mean_v","sd_v","st0")]
}

random.dmc<- function(n,p.df,model)
{
  rlba.2C(n,A=p.df$A,d=p.df$d,b=p.df$b,p_flc=p.df$p_flc, t0=p.df$t0[1], 
          mean_v=p.df$mean_v,sd_v=p.df$sd_v,st0=p.df$st0[1],
          posdrift = attr(model,"posdrift"))
}


likelihood.dmc <- function(p.vector,data,ok.types=c("norm2C"),min.like=1e-10) 
# Returns vector of likelihoods for each RT in data (in same order)
# !!! TO DO: types other than norm
{

#   # COMMENT OUT this check for speed after debugging
#   if ( !all(attr(model,"type") %in% ok.types) )
#     stop("Distribution function type not supported by likelihood.dmc")
  
  likelihood <- numeric(dim(data)[1])
  for ( i in row.names(attr(data,"model")) ) if ( !attr(data,"cell.empty")[i] )
  {
    p.df <- p.df.dmc(p.vector,i,attributes(data)$model,n1order=TRUE)
    likelihood[ attr(data,"cell.index")[[i]] ] <-
      switch(attr(attributes(data)$model,"type"),
        norm2C=n1PDF.2C(t=data$RT[attr(data,"cell.index")[[i]]],
          r2=data$R2[attr(data,"cell.index")[[i]]],
          A=p.df$A,
          b=p.df$b,
          d=p.df$d,
          p_flc=p.df$p_flc,
          t0=p.df$t0[1], 
          mean_v=p.df$mean_v,
          sd_v=p.df$sd_v,
          st0=p.df$st0[1],
          posdrift = attr(attr(data,"model"),"posdrift"))
      )
 }
 pmax(likelihood,min.like)
}


