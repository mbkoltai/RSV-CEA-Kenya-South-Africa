#############################################################################
# This file is part of the RSV modelling project.
# 
# => VISUAL PRESENTATION OF THE CEA RESULTS BY COUNTRY
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################


plot_CEA_country_results <- function(sim_output_filename)
{
  
  ###############################
  # LOAD DATA                   #
  ###############################
  
  load(paste0(sim_output_filename,'.RData'))
  dim(sim_output)
  
  # set pdf dimensions
  pdf_w_h <- c(4,4)
  
  # get plot output directory
  plot_output_folder <- get_plot_output_folder(sim_output,'CEA')
  
  # get all scenario-country-year combinations
  sim_output$scen_country_year <- paste(sim_output$scenario,
                                        sim_output$country_iso,
                                        sim_output$year,
                                        sep='_')
  # get unique combinations
  scen_country_year_opt <- unique(sim_output$scen_country_year)
  
  # safety check
  check_parallel_workers()
  
  # run all scenario-country_year combinations (in parallel)
  i_scen_country_year <- scen_country_year_opt[1]
  foreach(i_scen_country_year = scen_country_year_opt,
          .export=c('plot_country_results','evppi_gam','get_intervention_legend',
                    'get_net_benefit','get_all_net_benefit_results','get_intervention_legend',
                    'get_temp_output_folder','plot_evppi','cli_print'),
          .packages=all_packages,.combine=rbind) %dopar% 
    {
      plot_country_results(sim_output,i_scen_country_year,plot_output_folder,pdf_w_h)
    } -> dummy
  
  # plot a summary with all EVPPI results
  plot_all_evppi_results(sim_output_filename)
}


# Function to plot CEA results per country and year
# - CE plane
# - ICER => temporarly switched off
# - CEAC/CEAF
# - EVPPI
if(0==1) {
  scen_country_year <- 'basecase_KEN_2020'
}
plot_country_results <- function(sim_output,scen_country_year,plot_output_folder,pdf_w_h){
  
  country_data <- sim_output[sim_output$scen_country_year == scen_country_year,]
  cli_print('Create country CEA plots for:',scen_country_year)
  
  # set scenario colors and legend
  country_data$intervention_factor <- as.factor(country_data$intervention)
  scenario_code_legend             <- get_intervention_legend(country_data$intervention_factor)
  
  # open PDF stream
  pdf(paste0(plot_output_folder,'CEA_',scen_country_year,'.pdf'), width = 3*pdf_w_h[1], height = 1*pdf_w_h[2])
  par(mfrow=c(1,3),mar=c(5.1,4.1,1.2,1.2))
  
  # CE-PLANE
  x_lim <- quantile(country_data$total_DALY_disc_averted,c(0,0.99),na.rm=T)
  plot(country_data$total_DALY_disc_averted,country_data$incremental_cost_disc,
       col=alpha(as.numeric(country_data$intervention_factor)+1,0.2),
       xlab='DALY averted, discounted',
       ylab='Incremental cost, discounted (USD)',
       xlim=x_lim,
       pch=20)
  abline(h=0,lty=3)
  legend('bottomright',scenario_code_legend$name,
         col=scenario_code_legend$color,
         pch=rep(20,nrow(scenario_code_legend)), 
         bg='white',
         lty=rep(NA,nrow(scenario_code_legend)),
         lwd=rep(1,nrow(scenario_code_legend)))
  legend('topright',scen_country_year,bty='n')
  
  # willingness-to-pay values
  num_wtp <- 70
  wtp_max <- 30000 #4000
  opt_wtp <- seq(0,wtp_max,length.out=num_wtp)
  
  # get CEAC_CEAF
  CEA_out                    <- get_all_net_benefit_results(country_data,opt_wtp)
  net_benefit_all            <- CEA_out$net_benefit_all
  prob_high_net_benefit      <- CEA_out$prob_high_net_benefit        # CEAC
  prob_high_mean_net_benefit <- CEA_out$prob_high_mean_net_benefit   # CEAF
  net_benefit_legend         <- CEA_out$net_benefit_legend
  
  # CEAC - CEAF
  plot(range(opt_wtp),c(0,1.12),col=0,
       xlab = 'Willingness to pay for a DALY averted (USD)',
       ylab = 'Probability highest net benefit')
  abline(v=seq(0,max(opt_wtp),1000),lty=3,col=alpha(1,0.5))
  # CEAC
  for(i in 1:dim(prob_high_net_benefit)[1])
  {
    lines(opt_wtp,prob_high_net_benefit[i,],col=alpha(net_benefit_legend$color[i],0.8),lwd=4)
  }
  # CEAF
  for(i in 1:dim(prob_high_mean_net_benefit)[1])
  {
    lines( opt_wtp,prob_high_mean_net_benefit[i,],col=net_benefit_legend$color[i],lwd=4) # add line
    points(opt_wtp,prob_high_mean_net_benefit[i,],col=1,lwd=2,pch=1)
  }
  legend_ncol <- ceiling((length(net_benefit_legend$color)+1)/2)
  legend('topleft',
         c(net_benefit_legend$name_legend,'CEAF'),
         pch=c(rep(NA,length(net_benefit_legend$color)),1),
         lty=c(rep(1,length(net_benefit_legend$color)),0),
         lwd=2,
         col=c(net_benefit_legend$color,1),
         bg='white',
         ncol=legend_ncol,
         cex=0.9)
  
  legend('topright',scen_country_year,bty='n')
  
  ## EVPPI
  ## Based on code from Joke for Typhoid
  param_opt <- c('outpatient_cost','hosp_cost','rsv_rate','hosp_prob','hosp_CFR','comm_CFR_adj')
  num_param <- length(param_opt)
  country_param <- unique(country_data[,names(country_data) %in% param_opt])
  
  dim(country_param)
  evppi <- matrix(NA,num_wtp,num_param)  
  j <- 1; i <- 1
  for(j in 1:num_wtp){
    NB_tmp <- t(net_benefit_all[j,,])
    for(i in 1:num_param){
      evppi[j,i] <- evppi_gam(NB_tmp,country_param[,i])
    }
  }  
  colnames(evppi) <- param_opt
  
  # plot EVPPI
  plot_evppi(evppi,opt_wtp)
  
  # save EVPPI data
  output_dir       <- unique(country_data$outputFileDir)
  evppi_filename   <- file.path(get_temp_output_folder(output_dir,'plot_EVPPI'),paste0('CEA_',scen_country_year,'.RData'))
  save(evppi,opt_wtp,scen_country_year,file=evppi_filename)
  
  # close pdf stream
  dev.off()
  
}

get_all_net_benefit_results <- function(country_data,opt_wtp){
  
  # get standardized scenario names
  country_data$intervention_factor <- as.factor(country_data$intervention)
  
  # get the number of wtp values
  num_wtp <- length(opt_wtp)
  
  # get dimensions for the "highest net benefit probability" parameter
  num_sim  <- unique(country_data$num_sim)
  num_scen <- nlevels(country_data$intervention_factor) + 1
  
  # initiate the output parameter
  prob_high_net_benefit      <- matrix(NA,num_scen,num_wtp)  # CEAC
  prob_high_mean_net_benefit <- prob_high_net_benefit        # CEAF
  net_benefit_all            <- array(NA,dim=c(num_wtp,num_scen,num_sim))
  net_benefit_fctr           <- factor(c(as.character(country_data$intervention_factor),rep("comparator",num_sim)),
                                       levels=c(levels(country_data$intervention_factor),"comparator"))
  net_benefit_legend         <- get_intervention_legend(net_benefit_fctr)
  
  # get the probability of the highest NB per stochastic simulation
  for(i_wtp in 1:num_wtp){
    
    # calculate the net benefit
    net_benefit <- get_net_benefit(country_data$total_DALY_disc_averted,
                                   country_data$incremental_cost_disc,
                                   opt_wtp[i_wtp])
    
    # reshape
    net_benefit        <- data.frame(matrix(net_benefit,nrow=num_sim))
    names(net_benefit) <- unique(country_data$intervention_factor)
    
    # add "comparator" == no intervention
    net_benefit$comparator <- 0
    
    # reorder columns
    net_benefit <- net_benefit[,net_benefit_legend$name]
    
    # get highest net benefit per simulation
    high_net_benefit <- (apply(X = net_benefit, MARGIN = 1,FUN=function(X){X==max(X,na.rm=T)}))
    # aggregate to get a probability
    prob <- rowSums(high_net_benefit,na.rm=T) / num_sim
    
    # store probability of highest net benefit => CEAC
    prob_high_net_benefit[,i_wtp] <- prob
    
    # if not the highest "mean net benefit", set NA => CEAF
    mean_nb   <- colMeans(net_benefit,na.rm = T)
    prob_mean <- prob
    prob_mean[mean_nb != max(mean_nb,na.rm = T)] <- NA
    
    # store
    prob_high_mean_net_benefit[,i_wtp] <- prob_mean
    
    # store net_benefit
    dim(t(net_benefit))
    dim(net_benefit_all)
    net_benefit_all[i_wtp,,] <- t(net_benefit)
  }
  
  # return results
  return(list(net_benefit_all            = net_benefit_all,
              prob_high_net_benefit      = prob_high_net_benefit,
              prob_high_mean_net_benefit = prob_high_mean_net_benefit,
              net_benefit_legend         = net_benefit_legend))
}

#factor_codes <- net_benefit_fctr
get_intervention_legend <- function(factor_codes){
  
  factor_legend  <- data.frame(name=levels(factor_codes),
                               color=seq(1,nlevels(factor_codes))+1,
                               stringsAsFactors = F)
  
  # manual adjustments (to make plots uniform)
  factor_legend$color[factor_legend$name == 'maternal']    <- 3
  factor_legend$color[factor_legend$name == 'comparator']  <- 4
  
  factor_legend$name_legend <- factor_legend$name
  factor_legend$name_legend[factor_legend$name == 'comparator']  <- 'current'
  
  return(factor_legend)
}


plot_evppi <- function(evppi,opt_wtp,plot_main=''){
  
  # modify evppi colnames for figure legend 
  evpi_legend <- colnames(evppi)
  evpi_legend <- gsub('_',' ',evpi_legend)
  evpi_legend <- gsub('hosp','hosp.',evpi_legend)
  evpi_legend <- gsub('prob','probability',evpi_legend)
  evpi_legend <- gsub('rsv','RSV',evpi_legend)
  evpi_legend <- gsub('comm','comm.',evpi_legend)
  evpi_legend
  
  # rescale or transform the EVPPI
  evppi_plot   <- evppi
  ylim_plot <- range(c(0,evppi*1.15))
  
  plot(range(opt_wtp),range(evppi_plot),col=0,type='l',lwd=2,ylim=ylim_plot,
       xlab='Willingness to pay for a DALY averted (USD)',ylab=paste0('EVPPI (USD)'),
       main=plot_main)
  
  for(i in 1:ncol(evppi_plot)){
    lines(opt_wtp,evppi_plot[,i],col=i,lwd=2)
  }
  
  legend('top',evpi_legend,
         col=(1:ncol(evppi_plot)),
         lty=1,
         lwd=2,
         pch=NA,
         ncol=3,
         cex=0.9)
}

plot_all_evppi_results <- function(sim_output_filename){
  
  cli_print("AGGREGATE ALL EVPPI RESULTS")
  
  # load data
  load(paste0(sim_output_filename,'.RData'))
  dim(sim_output)
  
  # get number of realisations
  num_sim <- unique(sim_output$num_sim)
  
  # get EVPPI files
  output_dir           <- unique(sim_output$outputFileDir)
  temp_evppi_directory <- get_temp_output_folder(output_dir,'plot_EVPPI')
  temp_evppi_files     <- dir(temp_evppi_directory,full.names = T)
  
  # get plot output directory 
  plot_output_folder <- get_plot_output_folder(sim_output,'CEA_EVPPI')
  
  # get scenarios
  opt_scenario <- unique(sim_output$scenario)
  
  # create summary per senario
  sel_scen <- opt_scenario[1]
  for(sel_scen in opt_scenario){
    # open pdf stream
    pdf(file.path(plot_output_folder,paste0('EVPPI_all_',sel_scen,'_n',num_sim,'.pdf')),width=12,height=16)
    par(mfrow=c(4,3))
    
    # select evppi files from scenario 'sel_scen'
    temp_evppi_files_scen <- temp_evppi_files[grepl(sel_scen,temp_evppi_files)]
    
    # add plots
    i_evppi_file <-temp_evppi_files_scen[1]
    for(i_evppi_file in temp_evppi_files_scen){
      load(i_evppi_file)
      plot_evppi(evppi,opt_wtp,scen_country_year)
    }
    
    # close pdf stream
    dev.off()
    
  } # end for-loopt i_scen
  
}

