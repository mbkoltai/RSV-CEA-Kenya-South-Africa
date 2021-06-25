###### get_burden_flexible ------------------------------------
get_burden_flexible <- function(configList,incid_matrix,hosp_prob_matrix,dose_price) {
  # set RNG seed
  set.seed(configList$rng_seed)
  ###############################
  # Configuration               #
  ###############################
  select_1to5y <- function(x){ x_1to5y=x*0; x_1to5y[13:60,]=x[13:60,]; return(x_1to5y) }
  select_0to1y <- function(x){ x_0to1y<-x*0; x_0to1y[1:12,] <- x[1:12,];     return(x_0to1y)  }

  # get model parameters (based on country, year, num_sim, etc... as specified in configList)
  config <- get_rsv_ce_config(configList)
  # this loads costs/efficacy but also
  
  # load own data if it was supplied to the functions
  flag_own_data=""
  if (!is.null(dim(incid_matrix)) & !is.null(dim(hosp_prob_matrix))){
    print('check if first matrix is RSV disease episode/person year, second is rate of hospitalisation/person year!')
    config$rsv_rate=data.frame(incid_matrix); config$hosp_prob=data.frame(hosp_prob_matrix); flag_own_data="yes"}
  
  ###############################
  # Pre-processing              #
  ###############################
  
  # get 'infant' target population
  config$target_population_lifebirth <- config$target_population*(1-config$stillbirth_rate)
  
  # use the nYearsOfAges to construct age-dependent vectors
  config$nMonthsOfAges <- config$nYearsOfAges*config$monthsInYear
  
  ####################################
  # Set DALY CFR parameters          #
  ####################################
  
  # births => cohort size
  # get the life table for the given country and period
  life_table <- get_life_table(config$country_iso,config$year,config$outputFileDir,config$disc_rate_effect)
  # select the target ages	
  life_table <- life_table[1:config$nMonthsOfAges,]
  # adjust for target population PER MONTH = (% ppl left)*(livebirth)/12
  life_table$pop <- life_table$lx_rate * config$target_population_lifebirth / config$monthsInYear #fix 2018-01-11
  config$hosp_CFR_DALYloss      <- life_table$life_expectancy
  config$hosp_CFR_DALYloss_disc <- life_table$life_expectancy_disc
  
  ##################################
  # Parameter uncertainty: BURDEN  #
  ##################################
  # sample num_sim indices
  # config$rsv_rate is the dataframe of rsv cases/person adjusted to the given country
  # iterations to select
  # selecting num_sim samples from 
  sample_rsv_rate_ind      <- sample(ncol(config$rsv_rate),config$num_sim,replace=T)
  sample_hosp_prob_ind     <- sample(ncol(config$hosp_prob),config$num_sim,replace = T)
  sample_hosp_CFR_ind      <- sample(ncol(config$hosp_CFR),config$num_sim,replace = T)
  sample_comm_CFR_adj      <- sample(ncol(config$comm_CFR_adj),config$num_sim,replace = T)
  
  # select age-specific rsv rates from the provided options
  # this takes n samples from the iterations
  config$rsv_rate          <- config$rsv_rate[,sample_rsv_rate_ind] 
  config$hosp_prob         <- config$hosp_prob[,sample_hosp_prob_ind]
  # monot declining from 2.5 to 0.75
  config$hosp_CFR          <- config$hosp_CFR[,sample_hosp_CFR_ind]
  # community cfr adjustment factor=2.121466
  config$comm_CFR_adj      <- config$comm_CFR_adj[,sample_comm_CFR_adj] 

  # aggregate CFR 
  config$hosp_comm_CFR     <- config$hosp_CFR * config$comm_CFR_adj 

  if(config$num_sim == 1)  {
    config$rsv_rate      <- as.matrix(config$rsv_rate,ncol=1)
    config$hosp_prob     <- as.matrix(config$hosp_prob,ncol=1)
    config$hosp_comm_CFR <- as.matrix(config$hosp_comm_CFR,ncol=1)   }
  
  # SAMPLE DALY VALUES: gamma distribution. mean=alpha/beta
  sample_rgamma <- function(mean,stdev,num_sim){
    alpha <- (mean^2) / (stdev^2) # shape parameter
    beta  <- (stdev^2) / (mean) # rate parameter (beta=1/theta)
    sample <- rgamma(num_sim,shape=alpha,rate = 1/beta)
  }
  
  # duration illness, LRTI DALY sampled from a gamma distribution
  config$duration_illness <- sample_rgamma(config$duration_illness_mean,config$duration_illness_stdev,config$num_sim)
  config$severe_LRTI_DALY <- sample_rgamma(config$severe_LRTI_DALY_mean,config$severe_LRTI_DALY_stdev,config$num_sim)
  config$moderate_LRTI_DALY <- sample_rgamma(config$moderate_LRTI_DALY_mean,config$moderate_LRTI_DALY_stdev,config$num_sim)
  # DALY loss
  config$severe_rsv_DALYloss      <- config$severe_LRTI_DALY*config$duration_illness
  config$non_severe_rsv_DALYloss  <- config$moderate_LRTI_DALY*config$duration_illness
  
  ###############################
  # Infant protection (this is for antibodies)
  ###############################
  # convert duration of protection from years to months
  dur_prot_infant <- round(config$monthsInYear*config$dur_protection_infant)
  
  # Calculate the effective protection
  # init matrices
  iAgeEffectiveProtection_primary=matrix(0,ncol=config$num_sim,nrow=config$nMonthsOfAges)
  iAgeEffectiveProtection_hospital=iAgeEffectiveProtection_primary; iAgeEffectiveProtection_cfr=iAgeEffectiveProtection_primary
  
  # Fill matrices
  # i_eff <- 1
  for(i_eff in 1:config$num_sim){
    iAgeEffectiveProtection_primary[1:dur_prot_infant,i_eff]  <- config$efficacy_infant_primary[i_eff]*config$coverage_infant
    iAgeEffectiveProtection_hospital[1:dur_prot_infant,i_eff] <- config$efficacy_infant_hospital[i_eff]*config$coverage_infant
    iAgeEffectiveProtection_cfr[1:dur_prot_infant,i_eff]      <- config$efficacy_infant_cfr[i_eff]*config$coverage_infant
    }
  
###############################
# Maternal protection         #
###############################
# convert duration of protection from years to months
  dur_prot_maternal <- round(config$monthsInYear*config$dur_protection_maternal)
  
# Calculate the effective protection # init matrices
mAgeEffectiveProtection_primary <- matrix(0,ncol=config$num_sim,nrow=config$nMonthsOfAges)  #dummy
mAgeEffectiveProtection_hospital <- mAgeEffectiveProtection_primary
mAgeEffectiveProtection_cfr <- mAgeEffectiveProtection_primary # dummy
  
# if no maternal efficacy present yet, sample from (log)normal distribution
  if(any(is.na(config$efficacy_maternal_primary))){ 
    # dependent samples for all-cause, hospital and cfr efficacy ==>> use identical random-stream
    set.seed(configList$rng_seed)
    config$efficacy_maternal_primary <- sample_normal_dist(config$num_sim,configList$efficacy_maternal_mean,
                                                           configList$efficacy_maternal_stdev,
                                                           configList$efficacy_maternal_lognormal)
    set.seed(configList$rng_seed)
    config$efficacy_maternal_hospital <- sample_normal_dist(config$num_sim,configList$efficacy_maternal_hosp_mean,
                                                            configList$efficacy_maternal_hosp_stdev,
                                                            configList$efficacy_maternal_lognormal)
    set.seed(configList$rng_seed)
    config$efficacy_maternal_cfr <- sample_normal_dist(config$num_sim,configList$efficacy_maternal_cfr_mean,
                                                       configList$efficacy_maternal_cfr_stdev,
                                                       configList$efficacy_maternal_lognormal)
  }

  # Fill matrices # i_eff <- 1
  for(i_eff in 1:config$num_sim){
    mAgeEffectiveProtection_primary[1:dur_prot_maternal,i_eff] <- config$efficacy_maternal_primary[i_eff]*config$coverage_maternal
    mAgeEffectiveProtection_hospital[1:dur_prot_maternal,i_eff] <- config$efficacy_maternal_hospital[i_eff]*config$coverage_maternal
    mAgeEffectiveProtection_cfr[1:dur_prot_maternal,i_eff] <- config$efficacy_maternal_cfr[i_eff]*config$coverage_maternal }
  
  # Decrease the efficacy factor for the pre-term infants
  mAgeEffectiveProtection_primary  <- mAgeEffectiveProtection_primary*(1-config$pre_term_rate)
  mAgeEffectiveProtection_hospital <- mAgeEffectiveProtection_hospital*(1-config$pre_term_rate)
  mAgeEffectiveProtection_cfr      <- mAgeEffectiveProtection_cfr*(1-config$pre_term_rate)
  
  ###############################
  # Total protection         #
  ###############################
  # Combine the Maternal and Infant protection
  totalAgeEffectiveProtection_primary  <-  mAgeEffectiveProtection_primary  + iAgeEffectiveProtection_primary 
  totalAgeEffectiveProtection_hospital <-  mAgeEffectiveProtection_hospital + iAgeEffectiveProtection_hospital
  totalAgeEffectiveProtection_cfr      <-  mAgeEffectiveProtection_cfr      + iAgeEffectiveProtection_cfr
  config$outpatient_cost
  
  # Differential efficacy problem if number averted hospital cases > total number averted cases
  # Check if the non_hosp protection is negative ===>>> a negative number of averted cases ??!! ===>>> set NA
  if(any(totalAgeEffectiveProtection_primary < config$hosp_prob * totalAgeEffectiveProtection_hospital)){
    flag_sim <- colSums(totalAgeEffectiveProtection_primary < config$hosp_prob * totalAgeEffectiveProtection_hospital)
    table(flag_sim)
    dim(totalAgeEffectiveProtection_primary)
    totalAgeEffectiveProtection_primary[,flag_sim]   <- NA
    totalAgeEffectiveProtection_hospital[,flag_sim]  <- NA
    totalAgeEffectiveProtection_cfr[,flag_sim]       <- NA
    warning_message <- paste0(c("SCENARIO '", config$scenario,"' ENCOUNTERS A DIFFERENTIAL EFFICACY ISSUE:",
              "\n=================>> averted hospital cases > total averted cases",
              "\n=================>> disabled", sum(flag_sim), "simulations [",round(sum(flag_sim)/length(flag_sim)*100),"%]"))
    # print warning
    cli_print(warning_message,WARNING=T,FORCED =T)
    # store in 'config'
    config$total_efficacy_num_na <- sum(flag_sim)
  } else{
    config$total_efficacy_num_na <- 0
  }
  ##################################
  # Parameter uncertainty: COSTS   #
  ##################################
  print("calculating costs")
  # reshape the cost parameter into a [age, sim] structure
  config$outpatient_cost  <- matrix(rep(config$sample_outpatient_cost,each=config$nMonthsOfAges),ncol=config$num_sim)
  # ggplot(data.frame(value=config$sample_outpatient_cost)) + geom_density(aes(x=value)) + theme_bw() + xlab("outpatient cost")
  config$hosp_cost <- matrix(rep(config$sample_inpatient_cost,each=config$nMonthsOfAges),ncol=config$num_sim)
  # ggplot(data.frame(value=config$hosp_cost[1,])) + geom_density(aes(x=value)) + theme_bw() + xlab("inpatient cost")
  config$admin_cost_maternal <- matrix(0,nrow=config$nMonthsOfAges,ncol=config$num_sim)
  # input price data
  if (is.numeric(dose_price) & length(dose_price)>1) {config$price_dose_maternal<-dose_price["mat_vacc"];config$price_dose_mAb<-dose_price["mAb"]}
  config$admin_cost_maternal[1,]  <- config$sample_admin_cost + (config$price_dose_maternal*(1+config$wastage_maternal))
  
  config$admin_cost_mAb <- matrix(0,nrow=config$nMonthsOfAges,ncol=config$num_sim)
  config$admin_cost_mAb[1,] <- config$sample_admin_cost + (config$price_dose_mAb*(1+config$wastage_mAb)*config$price_mAb_sens_factor)
  ##############################
  # Discounting                #
  ##############################
  # setup discounting-over-time vector
  years_after_vaccin     <- rep(0:(config$nYearsOfAges-1),each=config$monthsInYear)
  disc_time_effect       <- 1/((1+config$disc_rate_effect)^years_after_vaccin)
  disc_time_cost         <- 1/((1+config$disc_rate_cost)^years_after_vaccin)
  ###############################
  # Burden of Disease           #
  ###############################	
  # Cases
  rsv_cases <- life_table$pop*config$rsv_rate 
  non_hosp_cases       <- rsv_cases* (1-config$hosp_prob)
  hosp_cases           <- rsv_cases* config$hosp_prob
  rsv_deaths           <- hosp_cases* config$hosp_comm_CFR
  # discounted
  rsv_cases_disc            <- colSums(rsv_cases       * disc_time_effect)
  non_hosp_cases_disc       <- colSums(non_hosp_cases  * disc_time_effect)
  hosp_cases_disc           <- colSums(hosp_cases      * disc_time_effect)
  rsv_deaths_disc           <- colSums(rsv_deaths      * disc_time_effect)
  
  # DALY
  # YLD = years lived with disability
  non_hosp_YLD         <- non_hosp_cases    * config$non_severe_rsv_DALYloss
  hosp_YLD             <- hosp_cases        * config$severe_rsv_DALYloss
  non_hosp_YLD_disc    <- non_hosp_YLD      * disc_time_effect
  hosp_YLD_disc        <- hosp_YLD          * disc_time_effect
  total_YLD            <- non_hosp_YLD + hosp_YLD
  total_YLD_disc       <- (non_hosp_YLD_disc + hosp_YLD_disc)
  
  # Life years lost => YLL 
  total_YLL      <- rsv_deaths  *  config$hosp_CFR_DALYloss
  total_YLL_disc <- (rsv_deaths  *  config$hosp_CFR_DALYloss_disc  * disc_time_effect)
  
  # DALY loss  
  total_DALY      <- total_YLD + total_YLL
  total_DALY_disc <- total_YLD_disc + total_YLL_disc	
  
  ###############################
  # Averted Burden              #
  ###############################	
  print("calculating averted burden")
  # Averted Cases
  rsv_cases_averted            <- rsv_cases             * totalAgeEffectiveProtection_primary
  hosp_cases_averted           <- hosp_cases            * totalAgeEffectiveProtection_hospital
  rsv_deaths_averted           <- rsv_deaths            * totalAgeEffectiveProtection_cfr
  non_hosp_cases_averted       <- rsv_cases_averted - hosp_cases_averted
  # Averted DALY loss => YLD averted
  # YLD = years lived with disability
  non_hosp_YLD_averted         <- non_hosp_cases_averted    * config$non_severe_rsv_DALYloss
  hosp_YLD_averted             <- hosp_cases_averted        * config$severe_rsv_DALYloss
  
  # Averted Life years lost => YLL averted
  # YLL = years of life lost due to premature mortality
  total_YLL_averted      <- rsv_deaths_averted       *  config$hosp_CFR_DALYloss 
  
  # discounted
  rsv_cases_disc_averted            <- colSums(rsv_cases_averted        * disc_time_effect)
  non_hosp_cases_disc_averted       <- colSums(non_hosp_cases_averted   * disc_time_effect)
  hosp_cases_disc_averted           <- colSums(hosp_cases_averted       * disc_time_effect)
  rsv_deaths_disc_averted           <- colSums(rsv_deaths_averted       * disc_time_effect)
  
  # YLD (note: this cannot be aggregated yet because it is used to calculate the "Total Averted Burden")
  non_hosp_YLD_disc_averted         <- (non_hosp_YLD_averted     * disc_time_effect)   # note
  hosp_YLD_disc_averted             <- (hosp_YLD_averted         * disc_time_effect)
  
  # Prevented Life years lost using the dicounted YLL averted (note: idem as YLD)
  total_YLL_disc_averted       <- (rsv_deaths_averted *  config$hosp_CFR_DALYloss_disc  * disc_time_effect)
  
  ###############################
  # Total Averted Burden        #
  ###############################	
  
  # YLD averted
  total_YLD_averted      <- non_hosp_YLD_averted + hosp_YLD_averted
  total_YLD_disc_averted <- (non_hosp_YLD_disc_averted + hosp_YLD_disc_averted)
  
  # YLL averted
  total_YLL_averted      <- total_YLL_averted
  total_YLL_disc_averted <- total_YLL_disc_averted
  
  # DALY averted  
  total_DALY_averted      <- total_YLD_averted + total_YLL_averted	
  total_DALY_disc_averted <- total_YLD_disc_averted + total_YLL_disc_averted
  
  ###############################
  # Costs                       #
  ###############################
  # medical costs: baseline
  print("calculating total costs")
  cost_rsv_outpatient <- non_hosp_cases * config$outpatient_cost 
  cost_rsv_hosp       <- hosp_cases* config$hosp_cost   
  # total medical cost
  total_medical_cost <- cost_rsv_outpatient + cost_rsv_hosp
  # medical costs: averted
  cost_rsv_outpatient_averted <- non_hosp_cases_averted * config$outpatient_cost 
  cost_rsv_hosp_averted       <- hosp_cases_averted     * config$hosp_cost   
  total_medical_cost_averted  <- cost_rsv_outpatient_averted + cost_rsv_hosp_averted
  # program costs maternal vaccination => target pop = live babies + stillbirths
  pop_maternal  <- config$target_population
  cost_maternal <- pop_maternal * config$coverage_maternal * config$admin_cost_maternal
  # program costs mAb => target pop = live babies
  cost_mAb      <- config$target_population_lifebirth * config$coverage_infant   * config$admin_cost_mAb
  # total intervention cost
  intervention_cost <- (cost_maternal + cost_mAb)
  # total intervention cost = averted costs over time & program costs for first month
  incremental_cost     <- intervention_cost - total_medical_cost_averted
  # discounting
  cost_rsv_outpatient_disc         <- (cost_rsv_outpatient         * disc_time_cost)
  cost_rsv_hosp_disc               <- (cost_rsv_hosp               * disc_time_cost)
  cost_rsv_outpatient_disc_averted <- (cost_rsv_outpatient_averted * disc_time_cost)
  cost_rsv_hosp_disc_averted       <- (cost_rsv_hosp_averted       * disc_time_cost)   
  total_medical_cost_disc_averted  <- (total_medical_cost_averted  * disc_time_cost)
  cost_maternal_disc               <- cost_maternal                * disc_time_cost
  cost_mAb_disc                    <- cost_mAb                     * disc_time_cost
  total_medical_cost_disc          <- (total_medical_cost  * disc_time_cost)
  intervention_cost_disc           <- (intervention_cost   * disc_time_cost) # this should be the same as the undiscounted value, because the intervention is given at the beginning
  incremental_cost_disc            <- (incremental_cost    * disc_time_cost)
  #################################
  ## OUTPUT                      ##
  #################################
  ## FUNCTION TO AGGREGATE BURDEN OUTPUT
  print("aggregating output")
  # note: fix to get age-specific results, withough duplicating code
  aggregate_output <- function(age_function,age_tag){
    output_all <- data.frame(
      rsv_cases      = colSums(age_function(rsv_cases)),
      non_hosp_cases = colSums(age_function(non_hosp_cases)),
      hosp_cases     = colSums(age_function(hosp_cases)),
      rsv_deaths     = colSums(age_function(rsv_deaths)),
      rsv_cases_disc,
      non_hosp_cases_disc,
      hosp_cases_disc,
      rsv_deaths_disc,
      # years of life with disab
      non_hosp_YLD            = colSums(age_function(non_hosp_YLD)),
      hosp_YLD                = colSums(age_function(hosp_YLD)),
      non_hosp_YLD_disc       = colSums(age_function(non_hosp_YLD_disc)),
      hosp_YLD_disc           = colSums(age_function(hosp_YLD_disc)),
      # cases averted
      rsv_cases_averted       = colSums(age_function(rsv_cases_averted)),
      non_hosp_cases_averted  = colSums(age_function(non_hosp_cases_averted)),
      hosp_cases_averted      = colSums(age_function(hosp_cases_averted)),
      rsv_deaths_averted      = colSums(age_function(rsv_deaths_averted)),
      
      rsv_cases_disc_averted,
      non_hosp_cases_disc_averted,
      hosp_cases_disc_averted,
      rsv_deaths_disc_averted,
      
      non_hosp_YLD_averted        = colSums(age_function(non_hosp_YLD_averted)),
      hosp_YLD_averted            = colSums(age_function(hosp_YLD_averted)),
      non_hosp_YLD_disc_averted   = colSums(age_function(non_hosp_YLD_disc_averted)),
      hosp_YLD_disc_averted       = colSums(age_function(hosp_YLD_disc_averted)),
      
      cost_rsv_outpatient         = colSums(age_function(cost_rsv_outpatient)),
      cost_rsv_hosp               = colSums(age_function(cost_rsv_hosp)),
      cost_rsv_outpatient_averted = colSums(age_function(cost_rsv_outpatient_averted)),
      cost_rsv_hosp_averted       = colSums(age_function(cost_rsv_hosp_averted)),
      total_medical_cost          = colSums(age_function(total_medical_cost)),
      total_medical_cost_averted  = colSums(age_function(total_medical_cost_averted)),
      intervention_cost           = colSums(age_function(intervention_cost)),
      incremental_cost            = colSums(age_function(incremental_cost)),
      
      cost_rsv_outpatient_disc         = colSums(age_function(cost_rsv_outpatient_disc)),
      cost_rsv_hosp_disc               = colSums(age_function(cost_rsv_hosp_disc)),
      total_medical_cost_disc          = colSums(age_function(total_medical_cost_disc)),
      cost_rsv_outpatient_disc_averted = colSums(age_function(cost_rsv_outpatient_disc_averted)),
      cost_rsv_hosp_disc_averted       = colSums(age_function(cost_rsv_hosp_disc_averted)),
      total_medical_cost_disc_averted  = colSums(age_function(total_medical_cost_disc_averted)),
      intervention_cost_disc           = colSums(age_function(intervention_cost_disc)),
      incremental_cost_disc            = colSums(age_function(incremental_cost_disc)),
      
      total_YLD                        = colSums(age_function(total_YLD)),
      total_YLL                        = colSums(age_function(total_YLL)),
      total_DALY                       = colSums(age_function(total_DALY)),
      
      total_YLD_disc                   = colSums(age_function(total_YLD_disc)),
      total_YLL_disc                   = colSums(age_function(total_YLL_disc)),
      total_DALY_disc                  = colSums(age_function(total_DALY_disc)),
      
      total_YLD_averted                = colSums(age_function(total_YLD_averted)),
      total_YLL_averted                = colSums(age_function(total_YLL_averted)),
      total_DALY_averted               = colSums(age_function(total_DALY_averted)),
      
      total_YLD_disc_averted           = colSums(age_function(total_YLD_disc_averted)),
      total_YLL_disc_averted           = colSums(age_function(total_YLL_disc_averted)),
      total_DALY_disc_averted          = colSums(age_function(total_DALY_disc_averted))
    )
    message("adjusting column names")
    # adjust column names
    if(nchar(age_tag)>0){
      # names(output_all) <- paste(names(output_all),age_tag,sep='_')
      # get base names
      b_names <- gsub('_disc','',names(output_all))
      b_names <- unique(gsub('_averted','',b_names))
      # fix: remove "non_hosp" names
      b_names <- b_names[!grepl('non_hosp',b_names)]
      # add age category
      b_names_age <- paste(b_names,age_tag,sep='_')
      # adapt column names
      for(i_name in 1:length(b_names)){ names(output_all) <- gsub(b_names[i_name],b_names_age[i_name],names(output_all)) }
    }
    # return
    return(output_all)
  }
  
  output_all <- data.frame(aggregate_output(select_0to5y,''),
                           aggregate_output(select_0to1y,'0to1y'),
                           aggregate_output(select_1to5y,'1to5y'))
  
  # add population
  output_all$population <- sum(life_table$pop)
  # reformat: round
  output_all[] <- round(output_all,digits=0)
  # add incidence
  output_all$rsv_incidence  <- round(colSums(rsv_cases) / output_all$population * 1000,digits=2)
  # add input parameters
  output_all$outpatient_cost    <- config$outpatient_cost[1,]
  output_all$hosp_cost          <- config$hosp_cost[1,]
  output_all$admin_cost         <- config$sample_admin_cost
  output_all$rsv_rate           <- colSums(config$rsv_rate)
  output_all$hosp_prob          <- colSums(config$hosp_prob)
  output_all$hosp_CFR           <- colSums(config$hosp_CFR)
  output_all$comm_CFR_adj       <- config$comm_CFR_adj[1,]
  output_all$hosp_comm_CFR      <- colSums(config$hosp_comm_CFR)
  output_all$duration_illness   <- config$duration_illness
  output_all$severe_LRTI_DALY   <- config$severe_LRTI_DALY
  output_all$moderate_LRTI_DALY <- config$moderate_LRTI_DALY
  
  output_all$FVP_maternal          <- config$target_population * config$coverage_maternal
  output_all$FVP_mAb               <- config$target_population_lifebirth * config$coverage_infant
  output_all$population_lifebirth  <- round(config$target_population_lifebirth)
  output_all$population_target     <- config$target_population
  
  output_all$price_dose_maternal   <- config$price_dose_maternal
  output_all$price_dose_mAb        <- config$price_dose_mAb * config$price_mAb_sens_factor
  
  output_all$efficacy_maternal_primary   <- config$efficacy_maternal_primary
  output_all$efficacy_maternal_hospital  <- config$efficacy_maternal_hospital
  output_all$efficacy_maternal_cfr       <- config$efficacy_maternal_cfr
  output_all$efficacy_total_num_na       <- config$total_efficacy_num_na
  
  output_all$income_region               <- config$income_region
  output_all$stochastic_proc_id  <- 1:config$num_sim
  
  # return results
  return(output_all)
}
