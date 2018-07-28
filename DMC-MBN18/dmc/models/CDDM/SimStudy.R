# Hey Yishin & Andrew,
# 
# So following up on what we talked about today, I wrote down roughly the things to consider in 
# the tutorial. The main things to vary across parameter recovery studies are:
# - Whether we use the exact likelihood or PDA in the sampler
# - The size of the simulated sample data (for a reasonably-sized study vs. to see 
#   asymptotic recovery)
# - Starting point bias (can be done with PDA)
# - Different sampling distributions or starting point variability that are 
#   specified as an input vector to the random CDDM function (again with PDA)
# 
# I recall we had said there would be approximately 10 parameter recovery studies 
# to run by varying these four things.  I'm not sure if this is what you had in 
# mind, but a list of 10 could be:
# - Likelihood method / reasonable sample size
# - PDA / reasonable sample size
# - Likelihood method / asymptotic sample size
# - PDA / asymptotic sample size
# - PDA / reasonable sample + start point bias
# - PDA / asymptotic sample + start point bias
# - PDA / reasonable + input sampling distribution 
#         (e.g., bimodal, two Von Meiss opposite directions)
# - PDA / asymptotic + input sampling distribution
# - PDA / reasonable + input starting point variability (e.g., circular uniform 
#   around a biased mean)
# - PDA / asymptotic + input starting point variability
# 
# So I think the things that we will have to do next is get these included in 
# the models that we use:
# - Threshold variability
# - Pre-drawn random numbers for sampling distribution
# - Pre-drawn random numbers for starting point
# 
# reasonable parameter values to use in the recovery study.
# 
