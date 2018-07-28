# Template setup for 2-choice LBA 3 Confidence choices,
# except different pdf.
# Two thresholds, b (upper) and d (lower). 
#   External parameters types: A, B, D1, D2, t0, meanlog_v, sdlog_v, st0 = 0 (optional)
#   Internal parameters types: A, b, d1, d2, t0, meanlog_v, sdlog_v, st0 = 0 (optional)

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
  par.df$d1 <- par.df$D1+par.df$A
  par.df$d2 <- par.df$D2+par.df$d1
  par.df$b <- par.df$B+par.df$d2
  
  
  par.df[,c("A","b","d1","d2","p_flc", "t0","meanlog_v","sdlog_v","st0")]
}

random.dmc<- function(n,p.df,model)
{
  rlba.lnorm.3TC(n,A=p.df$A,d1=p.df$d1,d2=p.df$d2,b=p.df$b,t0=p.df$t0[1], 
           p_flc=p.df$p_flc,meanlog_v=p.df$meanlog_v,sdlog_v=p.df$sdlog_v,st0=p.df$st0[1])
}


likelihood.dmc <- function(p.vector,data,ok.types=c("lnorm3TC"),min.like=1e-10) 
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
             lnorm3TC=n1PDF.lnorm.3TC(t=data$RT[attr(data,"cell.index")[[i]]],
                               r2=data$R2[attr(data,"cell.index")[[i]]],
                               A=p.df$A,
                               b=p.df$b,
                               d1=p.df$d1,
                               d2=p.df$d2,
                               p_flc=p.df$p_flc,
                               t0=p.df$t0[1], 
                               meanlog_v=p.df$meanlog_v,
                               sdlog_v=p.df$sdlog_v,
                               st0=p.df$st0[1])
      )
  }
  pmax(likelihood,min.like)
}


