#############################################################################
# This file is part of the RSV modelling project.
# 
# => FIT INCIDENCE AND MORTALITY
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

#f_country_iso <- 'GBR' ; f_outputFileDir <- './output'; 
get_incidence <- function(f_country_iso,f_outputFileDir){

  # filename
  filenames    <- file.path(get_temp_output_folder(f_outputFileDir),paste0('df_',f_country_iso,'_all','.RData'))

    if(!file.exists(filenames)){
    
      cli_print('GET COUNTRY SPECIFIC BURDEN FOR:',filenames)
      
      filename_country_BoD <- 'input/RSV_burden_Shi_2017.csv'
      
      #################################
      ##  UNCERTAINTY SPLINES        ##
      #################################
      
      # original data format: long table
      spline_datafiles <- data.frame(incidence ='./input/incidence_lmic_ts_n5000.csv',   # LMIC only,  5000x
                                     hosp_prob ='./input/hosp_prob_ts_n5000.csv',        # Nokes et al. 5000x, stdev = 2x standard error
                                     hosp_cfr  ='./input/cfr_lmic_ts_n5000.csv',         # LMIC only,  5000x
                                     stringsAsFactors = F) # hospital CRF
      
      #######################
      ##  CHECK FILES      ##
      #######################
      
      if(any(!file.exists(unlist(spline_datafiles)))){
        cli_print("Spline file(s) missing:", spline_datafiles[!file.exists(unlist(spline_datafiles))],WARNING=T)
        stop(paste("Spline file(s) missing:", spline_datafiles[!file.exists(unlist(spline_datafiles))]))
      }
      
      #######################
      ##  RSV CASES        ##
      #######################
      
      # load RSV incidence by monthly age
      # => from uncertainty analysis Marina Antillon
      incidence_mat <- convert_pred_into_model_input(spline_datafiles$incidence)

      #######################
      ##  Hospitalizations ##
      #######################
      
      # load rsv hospital probability by monthly age
      # => from uncertainty analysis Marina Antillon
      hosp_prob_mat <- convert_pred_into_model_input(spline_datafiles$hosp_prob)
      
      ###########################
      ##  Deaths
      ###########################
      
      # load rsv hospital rate by monthly age
      # => from uncertainty analysis Marina Antillon
      cfr_mat <- convert_pred_into_model_input(spline_datafiles$hosp_cfr)
      
      
      #####################
      # Country specific
      #####################
      
      # Country adaptations input
      # => take the data from Shi 2017 BoD paper, 
      burden_country <- read.csv(filename_country_BoD,stringsAsFactors = FALSE)
      
      # select details for the given country
      flag <- burden_country$country_iso == f_country_iso
      
      if(!any(flag)){
        cli_print('!! Reference data missing in ', filename_country_BoD, 'for: ',f_country_iso,'  !!')
        cli_print('!! USE KENYA !!')
        flag <- burden_country$country_iso == 'KEN'
      }
      
      # get reference incidence 
      burden_country_reference <- burden_country$incidence_RSV_associated_ALRI[flag]
      
      # Life table
      life_table <- get_life_table(f_country_iso,2015,f_outputFileDir,0)
      under5_pop <- life_table$lx[1:60]
      
      # adapt for the population (per month) and reshape parameters to handle uncertainty matrices
      RSV_incidence <- incidence_mat
      population    <- under5_pop
      
      # cases = incidence per 1000 * population by age
      country_rsv_cases      <- RSV_incidence / 1000 * population
      # rsv rate = cases / total cases (per column)
      country_rsv_rate       <- country_rsv_cases / rep(colSums(country_rsv_cases),each=dim(country_rsv_cases)[1])
      # country rate = rsv rate * country_reference per 1000 * country population
      country_rsv_pred_cases <- country_rsv_rate * burden_country_reference / 1000 * sum(population)
      # rsv rate = rsv cases / population
      country_rsv_rate       <- data.frame(country_rsv_pred_cases  / population)
      
      
      ###########################
      ## CHECK
      ###########################
  
      pdf(sub('_gavi.RData','.pdf',filenames))
        # baseline RSV incidence
        plot(0:59,incidence_mat[,1],type='l',
             ylim=range(incidence_mat),
             ylab='RSV incidence (baseline)',
             xlab='age',col=0)
        for(i in 1:dim(incidence_mat)[2]){
          lines(0:59,incidence_mat[,i],col=alpha(1,0.1))
        }
        
        # Hospitalisation probability
        plot(0:59,hosp_prob_mat[,1],type='l',
             ylim=range(hosp_prob_mat),
             ylab='RSV hospitalisation probability',
             xlab='age', col=0)
        for(i in 1:dim(hosp_prob_mat)[2]){
          lines(0:59,hosp_prob_mat[,i],col=alpha(1,0.1))
        }
  
        # HOSP CFR
        plot(0:59,cfr_mat[,1],type='l',
             ylim=range(cfr_mat),
             ylab='CFR',
             xlab='age', col=0)
        for(i in 1:dim(cfr_mat)[2]){
          lines(0:59,cfr_mat[,i],col=alpha(1,0.1))
        }

        par(mfrow=c(1,3))
        boxplot(colSums(incidence_mat),ylab='sum(RSV incidence per age)')
        boxplot(colSums(hosp_prob_mat),ylab='sum(RSV hopitalisaion probability per age)')
        boxplot(colSums(cfr_mat),ylab='sum(RSV CFR per age)')
        
      dev.off()
  
      ###########################
      ## SAVE
      ###########################
      
      # incidence rates and probabilties
      df_country  <- list(rsv_rate  = country_rsv_rate,
                          hosp_prob = hosp_prob_mat,
                          hCFR_prob = cfr_mat)
   
      save(df_country,file=filenames)
      
      cli_print('CALCULATE BURDEN COMPLETE:',f_country_iso)
  } 

  # load data and return
  load(filenames)
  return(df_country)

}

####################################################################
# HELP FUNCTION TO CONVERT SPLINE PREDICTIONS INTO MODEL OUTPUT    #
####################################################################

convert_pred_into_model_input <- function(prediction_file){
  
  # load data
  raw_data <- read.table(prediction_file,sep=',',header=T)
  dim(raw_data)
  num_ages_raw <- length(unique(raw_data$mos))
  num_sim_raw  <- length(unique(raw_data$iter)) 
  
  # convert into matrix [sim,age]
  clean_data <- matrix(raw_data$pred,ncol=num_ages_raw,byrow=T)
  
  # select first 60 age groups
  clean_data <- clean_data[,1:60]
  
  # remove names and transpose
  names(clean_data) <- NULL
  clean_data <- data.frame(t(clean_data))

  # return results  
  return(clean_data)
  
}


