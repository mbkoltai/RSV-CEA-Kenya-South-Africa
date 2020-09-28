#############################################################################
# This file is part of the RSV modelling project.
# 
# => WRITE TABLE WITH SUMMARY STATISTICS
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

################################################################
# Function header
################################################################

#run_tag <- 'RSV_gavi72'
write_global_summary_tables <- function(sim_output_filename){
  
  cli_print('Create global summary table for',sim_output_filename)
  
  ##################################################################################
  ## LOAD SIMULATION DATA
  ##################################################################################
  
  # load data
  load(paste0(sim_output_filename,'.RData'))
  
  # add discounded total cost
  sim_output$total_cost_0to1y_disc              <- sim_output$total_medical_cost_0to1y_disc
  sim_output$total_cost_0to1y_disc_intervention <- sim_output$total_medical_cost_0to1y_disc - sim_output$total_medical_cost_0to1y_averted + sim_output$intervention_cost_0to1y
  sim_output$total_cost_1to5y_disc              <- sim_output$total_medical_cost_1to5y_disc
  sim_output$total_cost_1to5y_disc_intervention <- sim_output$total_medical_cost_1to5y_disc - sim_output$total_medical_cost_1to5y_averted + sim_output$intervention_cost_1to5y

  # hack for incremental cost => add negative incremental AVERTED and 0 for "current"
  col_cost      <- grepl('incremental_cost',names(sim_output))
  colnames_cost_averted <- paste0(names(sim_output)[col_cost],'_averted')
  sim_output[,colnames_cost_averted] <- -sim_output[,col_cost]
  sim_output[,col_cost] <- 0
  
  # hack for intervention cost => add negative intervention AVERTED and 0 for "current"
  col_cost      <- grepl('intervention_cost',names(sim_output))
  colnames_cost_averted <- paste0(names(sim_output)[col_cost],'_averted')
  sim_output[,colnames_cost_averted] <- -sim_output[,col_cost]
  sim_output[,col_cost] <- 0
  
  col_burden_averted      <- grepl('averted',names(sim_output))
  colnames_burden_averted <- names(sim_output)[col_burden_averted]
  colnames_burden         <- gsub('_averted','',colnames_burden_averted)
  colnames_burden_intervention <- paste0(colnames_burden,'_intervention')
  names(sim_output)
  
  # calculate burden with intervention
  sim_output_intervention <- sim_output[,colnames_burden] - sim_output[,colnames_burden_averted]
  sim_output[,colnames_burden_intervention] <- sim_output_intervention
  
  # realisation id
  num_sim <- unique(sim_output$num_sim)
  sim_output$sim_id <- 1:num_sim
 
  # convert country_iso into levels => semi-numeric
  sim_output$country_iso <- factor(sim_output$country_iso)
  any(table(sim_output$sim_id)!= 9*72) # check the size of all the tables
  
  if(any(is.na(sim_output$efficacy_maternal))){
    sim_output$efficacy_maternal <- 0
  }
 
  # aggregate: sum over all countries (by sim_id)
  sim_data_global <- aggregate(. ~ config_tag + scenario + intervention + sim_id + outputFileDir, data=sim_output, sum, na.rm=TRUE, na.action = NULL)
  
  # aggregate: summary statistics
  sim_data_mean      <- aggregate(. ~ config_tag + scenario + intervention + outputFileDir, data=sim_data_global, mean, na.rm = TRUE, na.action = NULL) 
  sim_data_CI_LR     <- aggregate(. ~ config_tag + scenario + intervention + outputFileDir, data=sim_data_global, quantile,0.025, na.rm = TRUE, na.action = NULL)
  sim_data_CI_HR     <- aggregate(. ~ config_tag + scenario + intervention + outputFileDir, data=sim_data_global, quantile,0.975, na.rm = TRUE, na.action = NULL)
  
  # round => mean
  sim_data_mean  <- round_sim_data_scale(sim_data_mean,digits_x = 0, scale = 1)
  
  # write table with mean and CI
  sim_data_all <- generate_mean_CI_matrix(sim_data_mean,
                                          sim_data_CI_LR,
                                          sim_data_CI_HR)
  
  # write table with factor 1000 for mean and CI
  sim_data_all_k <- generate_mean_CI_matrix(sim_data_mean,
                                            sim_data_CI_LR,
                                            sim_data_CI_HR,
                                            f_digits = 0,
                                            f_scale = 1e3)
  
  burden_base_flag        <- 1
  global_data_burden      <- sim_data_mean[burden_base_flag,colnames_burden]
  global_data_burden_ci   <- sim_data_all[burden_base_flag,colnames_burden]
  global_data_burden_ci_k <- sim_data_all_k[burden_base_flag,colnames_burden]
  
  intervention_scenario_opt             <- unique(sim_data_mean[,c('intervention','scenario')])
  row.names(intervention_scenario_opt)  <- apply(intervention_scenario_opt,1,paste,collapse='_')
  
  global_data_all   <- NULL 
  global_data_all_k <- NULL
  
  for(i in 1:nrow(intervention_scenario_opt)){
    intervention_scenario_flag <- sim_data_mean$intervention==intervention_scenario_opt$intervention[i] & 
                                  sim_data_mean$scenario==intervention_scenario_opt$scenario[i]
    global_data_all  <- rbind(global_data_all,
                               sim_data_all[intervention_scenario_flag,colnames_burden_intervention])
    global_data_all_k <- rbind(global_data_all_k,
                             sim_data_all_k[intervention_scenario_flag,colnames_burden_intervention])
    }
  row.names(global_data_all)    <- apply(intervention_scenario_opt,1,paste,collapse='_')
  row.names(global_data_all_k)  <- apply(intervention_scenario_opt,1,paste,collapse='_')

 names(global_data_all) <- colnames_burden
 global_data_all <- rbind(no_intervention = global_data_burden_ci,
                            global_data_all)
 
 names(global_data_all_k) <- colnames_burden
 global_data_all_k <- rbind(no_intervention = global_data_burden_ci_k,
                          global_data_all_k)
 
 # add column with parameter names
 table_all   <- data.frame(output=names(global_data_all),t(global_data_all),stringsAsFactors = F)
 table_all_k <- data.frame(output=names(global_data_all_k),t(global_data_all_k),stringsAsFactors = F)
 
 # write summary table with all results
 write.table(table_all,file.path(unique(sim_output$outputFileDir),paste0(run_tag,'_summary_global_statistics.csv')),sep=',',row.names=F)

 #______________________________________________________________________
 ### manuscript results table: burden in 0-1 year and 1-5 years of age
 row_names_results_table_0to1y=c('rsv_cases_0to1y',
                           'hosp_cases_0to1y',
                           'rsv_deaths_0to1y',
                           
                           'total_YLD_0to1y_disc',
                           'total_YLL_0to1y_disc',
                           'total_DALY_0to1y_disc',
                           
                           'intervention_cost_0to1y_disc',
                           'cost_rsv_outpatient_0to1y_disc',
                           'cost_rsv_hosp_0to1y_disc',
                           'total_cost_0to1y_disc',
                           
                           'rsv_cases_0to1y_averted',
                           'hosp_cases_0to1y_averted',
                           'rsv_deaths_0to1y_averted',
                           
                           'total_DALY_0to1y_disc_averted',
                           'incremental_cost_0to1y_disc')
 
 row_names_intervention_0to1y<-c('rsv_cases_0to1y_intervention',
                                  'hosp_cases_0to1y_intervention',
                                  'rsv_deaths_0to1y_intervention',
    
                                  'total_YLD_0to1y_disc_intervention',#added by Xiao
                                  'total_YLL_0to1y_disc_intervention',
                                  'total_DALY_0to1y_disc_intervention',
    
                                  'intervention_cost_0to1y_disc_intervention',
                                  'cost_rsv_outpatient_0to1y_disc_intervention',
                                  'cost_rsv_hosp_0to1y_disc_intervention',
                                  'total_cost_0to1y_disc_intervention',
    
                                  'rsv_cases_0to1y_averted',
                                  'hosp_cases_0to1y_averted',
                                  'rsv_deaths_0to1y_averted',
    
                                  'total_DALY_0to1y_disc_averted',
                                  'incremental_cost_0to1y_disc_intervention')
 
 row_names_results_table_general <- gsub('_0to1y','',row_names_results_table_0to1y)
 row_names_results_table_1to5y   <- gsub('0to1y','1to5y',row_names_results_table_0to1y)
 row_names_intervention_1to5y    <- gsub('0to1y','1to5y',row_names_intervention_0to1y)
 
 Table2 <- data.frame(RSV_associated_burden = row_names_results_table_general,
                      NoIntervention_0to1y = unlist(sim_data_all_k[sim_data_all_k$config_tag == 'mAb_basecase',row_names_results_table_0to1y]),
                      NoIntervention_1to5y = unlist(sim_data_all_k[sim_data_all_k$config_tag == 'mAb_basecase',row_names_results_table_1to5y]),
                      
                      mAb_0to1y = unlist(sim_data_all_k[sim_data_all_k$config_tag == 'mAb_basecase',row_names_intervention_0to1y]),
                      maternal_0to1y = unlist(sim_data_all_k[sim_data_all_k$config_tag == 'maternal_basecase',row_names_intervention_0to1y]))
 
 # fix 'averted' for NoIntervention
 Table2[grepl('averted',Table2$RSV_associated_burden),grepl('NoIntervention',names(Table2))] <- '0  [0 - 0]'
 
 # write table2
 write.table(Table2,file.path(unique(sim_output$outputFileDir),paste0(run_tag,'_Table2.csv')),sep=',',row.names=F)
 }


# round simulation data with given number of digits and rescale (optional) 
round_sim_data_scale <- function(sim_data_x,digits_x=0,scale = 1){
  
  flag_numeric              <- unlist(lapply(sim_data_x[1,],is.numeric))
  sim_data_x[,flag_numeric] <- round(sim_data_x[,flag_numeric],digits=2)
  
  flag_a <- !is.na(sim_data_x[1,]) & flag_numeric & sim_data_x[1,] > 1
  sim_data_x[,flag_a] <- round(sim_data_x[,flag_a]/scale,digits=digits_x)
  
  return(sim_data_x)
}

# combine mean with confidence intervals (CI)
generate_mean_CI_matrix <- function(f_mean,f_CI_LR,f_CI_HR,f_digits=0,f_scale = 1){
  
  # round
  f_mean  <- round_sim_data_scale(f_mean,f_digits,f_scale)
  f_CI_LR <- round_sim_data_scale(f_CI_LR,f_digits,f_scale)
  f_CI_HR <- round_sim_data_scale(f_CI_HR,f_digits,f_scale)
  
  # get output columns with numeric content
  flag_numeric     <- unlist(lapply(f_mean[1,],is.numeric))
  colnames_output  <- colnames(f_mean)[flag_numeric]
  
  # start from the mean...
  f_data_all <- f_mean  
  
  # help function to format the numbers (number of decimals, but no leading spaces)
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k, big.mark=","))

  # for each output variable => add median and CI
  col_out <- colnames_output[1]
  for(col_out in colnames_output){
    f_data_all[,col_out] <- paste0(specify_decimal(f_mean[,col_out],f_digits),
                                         '  [',specify_decimal(f_CI_LR[,col_out],f_digits),
                                         ' - ',specify_decimal(f_CI_HR[,col_out],f_digits),']')
  }
  
  # return result
  return(f_data_all)
}




