#############################################################################
# This file is part of the RSV modelling project.
# 
# => DEFAULT CONFIGURATION FOR THE RSV ANALYSIS
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

# Main function to retrieve model settings
get_rsv_ce_config <- function(configList){
  # Initialize config variable
  config <- list(efficacy=NA)
  
  # Add directory names for input and output
  config$inputFileDir              <- "./input/"
  config$outputFileDir             <- configList$outputFileDir

  # time horizon and param
  config$nYearsOfAges   <- 5
  config$monthsInYear   <- 12
  
  # add default config parameters
  config$num_sim                     <- configList$num_sim
  config$rng_seed                    <- NA
  config$efficacy_maternal           <- NA
  config$efficacy_infant             <- NA
  config$dur_protection_maternal     <- NA
  config$dur_protection_infant       <- NA
  config$year                        <- configList$year
  config$country_iso                 <- NA 
  config$config_tag                  <- NA
  config$scenario                    <- NA
  config$intervention                <- NA
  config$scenario_id                 <- NA
  config$coverage_maternal           <- NA
  config$coverage_infant             <- NA
  config$mortality_adjustment_fctr       <- 1  # sensitivity analysis
  config$price_mAb_sens_factor           <- 1  # sensitivity analysis
  config$protection_infant_sens_factor   <- NA # sensitivity analysis
  config$protection_maternal_sens_factor <- NA # sensitivity analysis
  
  # set country name
  config$country_iso <- configList$country_iso
  
  # DALY parameters
  config$duration_illness_mean         <- 11.2/365 # duration of illness is 11.2 days  
  config$duration_illness_stdev        <- (12.3-10.1)/3.92/365 # ref: Briggs Book
  config$severe_LRTI_DALY_mean         <- 0.21
  config$severe_LRTI_DALY_stdev        <- (0.298-0.139)/3.92
  config$moderate_LRTI_DALY_mean       <- 0.053  
  config$moderate_LRTI_DALY_stdev      <- (0.074-0.032)/3.92
  
  config$severe_rsv_DALYloss      <- config$severe_LRTI_DALY   * config$duration_illness    
  config$non_severe_rsv_DALYloss  <- config$moderate_LRTI_DALY * config$duration_illness 
  config$hosp_CFR_DALYloss        <- NA  # simulated ages and period can vary, so computed in get_burden function  
  config$hosp_CFR_DALYloss_disc   <- NA  # discount rate can vary, so computed in get_burden function  
  
  # BIRTHS
  all_country_data <- read.table('./input/country_details_gavi72.csv',sep=',',header=T,stringsAsFactors=F)
  if (!any(all_country_data$country_iso3 %in% configList$country_iso)){
    print(paste0('loading data for ',configList$country_iso, ' (not in original gavi 72)'))
    all_country_data <- data.frame(read_csv('./input/country_details_gavi72_expanded.csv'))
  }
  
  # check if required country (and year) is present in the data base  
  flag  <- all_country_data$country_iso3 == configList$country_iso & all_country_data$year == configList$year 
  if(!any(flag))
  { flag_country            <- all_country_data$country_iso3 == configList$country_iso
    country_years           <- all_country_data$year[flag_country]
    country_year_selection  <- country_years[which.min(abs(country_years - as.numeric(configList$year)))]
    flag_year               <- all_country_data$year == country_year_selection
    flag                    <- flag_country & flag_year
    
    if(!exists('warning_population')){
      cli_print('WARNING: NO TARGET POPULATION AND BIRTH RATE DATA AVAILABLE FOR ',configList$country_iso,configList$year)
      cli_print('=======> USE DEMOGRAPHIC DATA FROM ',country_year_selection)
      warning_population <<- TRUE
    }
    
  }
  config$target_population  <- all_country_data$target_population[flag]
  config$pre_term_rate      <- all_country_data$incomplete_maternal_transfer_rate[flag]   
  config$stillbirth_rate    <- all_country_data$stillbirth_rate[flag]
  config$income_region      <- as.factor(all_country_data$income_region[flag])
  
  # BURDEN OF DISEASE
  # Get country specific incidence
  df_country <- get_incidence(config$country_iso,config$outputFileDir)

  # Set config values
  config$rsv_rate            <- df_country[[1]]
  config$hosp_prob           <- df_country[[2]]
  config$hosp_CFR            <- df_country[[3]]
  
  # sample community CFR ratios
  comm_CFR_base            <- sample((150:290)/100 , configList$num_sim,replace = T)  # uniform distribution
  comm_CFR_flu_factor      <- sample((90:100)/100 , configList$num_sim,replace = T)   # uniform distribution
  comm_CFR_adj             <- comm_CFR_base * comm_CFR_flu_factor
  
  # replicate these values to get the same value for each age
  config$comm_CFR_adj      <- matrix(rep(comm_CFR_adj,each=nrow(config$hosp_CFR)),nrow=nrow(config$hosp_CFR),byrow=F)
  
  # Set discounting values
  config$disc_rate_cost      <- 0.03
  config$disc_rate_effect    <- 0.03
  
  ########################
  # COST parameters      #
  ########################
  
  filename_cost_outpatient='./input/cost_data_outpatient.csv'; filename_cost_inpatient='./input/cost_data_inpatient.csv'
  config$sample_outpatient_cost <- get_cost_data(configList$country_iso,config$num_sim, filename_cost_outpatient)
  config$sample_inpatient_cost  <- get_cost_data(configList$country_iso,config$num_sim, filename_cost_inpatient)
  
  config$sample_admin_cost      <- 0
  
  config$price_dose_maternal    <- 3     
  config$price_dose_mAb         <- 6    
  
  # dose wastage
  config$wastage_maternal    <- 0.05
  config$wastage_mAb         <- 0.05
  
  
  #############################
  # EFFICACY parameters       #
  #############################
  
  # default: point estimate for primary endpoint and no additional efficacy for hospitalization or mortality
  config$efficacy_maternal_primary    <- rep(configList$efficacy_maternal,configList$num_sim)
  config$efficacy_maternal_hospital   <- rep(configList$efficacy_maternal,configList$num_sim)
  config$efficacy_maternal_cfr        <- rep(configList$efficacy_maternal,configList$num_sim)
  
  # option: differential efficacy for hospitalization and mortality
  config$efficacy_maternal_mean       <- NA
  config$efficacy_maternal_stdev      <- NA
  config$efficacy_maternal_hosp_mean  <- NA
  config$efficacy_maternal_hosp_stdev <- NA
  config$efficacy_maternal_cfr_mean   <- NA
  config$efficacy_maternal_cfr_stdev  <- NA
  config$efficacy_maternal_lognormal  <- FALSE  # boolean to sample from (log)Normal distrubtion
  
  # update on 20190916:
  # Sample from (log)normal distribution in 'get_burden' script, 
  # to NOT interfere with random sampling of other (burden) estimates
  # do not (occasionally) sample efficacy values here here!!
  
  # monoclonal antibodies
  config$efficacy_infant_primary    <- rep(configList$efficacy_infant,configList$num_sim)
  config$efficacy_infant_hospital   <- rep(configList$efficacy_infant,configList$num_sim)
  config$efficacy_infant_cfr        <- rep(configList$efficacy_infant,configList$num_sim)
  
  
  #############################
  ## TRANSFER CONFIG DETAILS ##
  #############################
  # Check if the provided configuration, has only valid parameters
  # valid if all column names are present in the default configuration
  if(any(!names(configList) %in% names(config))){
    stop(paste(c('!! STOP EXECUTION => UNUSED CONFIG PARAMETER(S): ',
                paste(names(configList)[!names(configList) %in% names(config)],collapse=', '))))
  }
  
  # transfer the function parameter to 'config'
  config[names(configList)] <- configList
  
  #############################
  ## CONSISTENCY CHECKS      ##
  #############################
  
  flag <- grepl('efficacy',names(config)) |
          grepl('prob',names(config)) |
          grepl('rate',names(config)) |
          grepl('coverage',names(config))
  
  if(any(unlist(config[flag])>1,na.rm=T)){
    cli_print('INCONSISTENT CONFIGURATION: EFFICACY, RATE, PROBABILITY OR COVERAGE > 1 WITH', 
              configList$config_tag,WARNING=T,FORCED=T)
  }
  
  # return config list
  return(config)
}
#########################################################
# return the filename of the burden files
get_burden_filename <- function()
{ # Country adaptation input => take data from Shi 2017 BoD paper, 
  filename_burden_country <- "input/RSV_burden_Shi_2017.csv"
  # return filenames
  return(filename_burden_country)
}

# sample from (log)Normal distribution
sample_normal_dist <- function(f_num_samples,f_mean,f_sd,f_bool_lognormal=FALSE){
  # sample from (log)Normal distrubtion
  if(any(f_bool_lognormal)){
    sample_all <- 1-exp(rnorm(f_num_samples,f_mean,f_sd))
  } else {
    sample_all <- rnorm(f_num_samples,f_mean,f_sd)
  }
  # truncate negative values
  sample_all[sample_all<0] <- 0 
  #return samples
  return(sample_all)
}