#############################################################################
# This file is part of the RSV modelling project.
# 
# => MISCELLANEOUS HELP FUNCTION FOR COST-EFFECTIVENESS ANALYSES
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################


# calculate the net benefit, given a price per DALY 'p'
get_net_benefit <- function(daly_averted,incr_cost_disc,wtp_level)
{
  NB <-   wtp_level*daly_averted - incr_cost_disc
  return(NB)
}


# NB <- NB_tmp; input_parameter_values <- country_param[,i]
evppi_gam <- function(NB,model_parameter_values){
  
  # note: current should have NB == 0
  D_opt  <- which(!colSums(NB) == 0)
  D <- ncol(NB)
  N <- nrow(NB)
  
  g.hat_new <- matrix(0,nrow=N,ncol=D)
  for(d in D_opt)
  {
    model <- gam(NB[,d] ~ model_parameter_values)
    g.hat_new[,d] <- model$fitted
  }   
  
  perfect.info  <- mean(apply(g.hat_new,1,max))
  baseline      <- max(colSums(g.hat_new)/N)
  
  partial.evpi  <- round(perfect.info - baseline, digits=4) ## estimate EVPPI 
  
  return(partial.evpi)
}

