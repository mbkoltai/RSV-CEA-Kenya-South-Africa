#############################################################################
# This file is part of the RSV modelling project.
# 
# => GET COUNTRY SPECIFIC OUTPATIENT AND INPATIENT RSV COSTS
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################
#rm(list=ls())
# f_country_iso <- 'KEN' ; f_num_sim <- num_sim; filename_costs <- './input/cost_data_outpatient.csv'
# main function to sample cost data for RSV
get_cost_data <- function(f_country_iso,f_num_sim,filename_costs)
{
  # load
  df_cost <- read.table(filename_costs,sep=',')
  
  # select country
  df_cost <- df_cost[rownames(df_cost)==f_country_iso,]
  
  # select random samples
   if(f_num_sim <= length(df_cost)){
     # note: this fix is to obtain an identical random sample as in the initial analysis, 
     # but without the meta-analysis
    df_cost <- df_cost[1:f_num_sim]  
  } else{ 
    # note: this is the normal procedure, which is used is the requested number of samples
    # is larger than the available samples.
    df_cost <- df_cost[sample(length(df_cost),f_num_sim,replace = T)]
  }
  
  # change format
  df_cost <- as.numeric(df_cost)
  
  # return
  return(df_cost)
}

