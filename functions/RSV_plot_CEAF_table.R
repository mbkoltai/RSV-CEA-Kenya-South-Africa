#############################################################################
# This file is part of the RSV modelling project.
# 
# => GET AGGREGATED CEAF STATISTICS OVER MULTIPLE COUNTRIES
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

## GENERAL FUNCTION TO DISPLAY ALL RESULTS
plot_CEAF_table <- function(sim_output_filename){
  
  # LOAD DATA
  load(paste0(sim_output_filename,'.RData'))
  sim_output_all <- sim_output
  
  # get plot output directory
  plot_output_folder <- get_plot_output_folder(sim_output,'CEA_table')
  
  if(length(unique(sim_output$country_iso))==1){
    cli_print("CEAF TABLE REQUIRES AT LEAST TWO COUNTRIES => NO TABLE FIGURES")
  } else{
  
  # safety check
  check_parallel_workers()
  
  # SELECT SCENARIO AND GENERATE RESULTS
  scenario_opt <- unique(sim_output_all$scenario)
  
  # note: parallel loop caused issues with very large samples...
  i_scenario <- 1
  foreach(i_scenario = 1:length(scenario_opt),
          .combine='rbind') %do%{

            # select simulation output
            sim_output <- sim_output_all[sim_output_all$scenario==scenario_opt[i_scenario],]
            
            # create CEAF table
            plot_CEAF_table_scenario(sim_output,plot_output_folder)
    
    } -> dummy # end foreach
  
  } # end if-clause to check the number of countries
  
} # end function

## FUNCTION TO GENERATE RESULTS FOR ONE SCENARIO
plot_CEAF_table_scenario <- function(sim_output, plot_output_folder)
{
  
  # create scenario tag          
  scenario_tag <- paste(unique(sim_output$scenario))
  
  # check the given number of cournties
  if(length(unique(sim_output$country_iso))==1){
    #cli_print("CEAF TABLE NOT POSSIBLE FOR:", scenario_tag)
    return(-1)
  }
  
  # print scenario tag if we continue        
  cli_print('CEAF TABLE START: ',scenario_tag)
  
  # COUNTRY DETAILS
  opt_countries <- unique(sim_output$country_iso)
  num_countries <- length(opt_countries)
  
  ## SELECT ONE YEAR
  if(length(unique(sim_output$year))>1){
    sim_output <- sim_output[sim_output$year == sim_output$year[1],]
    cli_print('SELECT ONE YEAR FOR THE CEAF TABLE')  
  }
  
  # with differential efficacy => use primary efficacy as 'base'
  if(any(is.na(sim_output$efficacy_maternal))){
    sim_output$efficacy_maternal <- sim_output$efficacy_maternal_primary
  }
  
  # aggregate the burden per country (select the first intervention type)
  country_summary <- aggregate(. ~ scenario_id + config_tag + scenario + intervention + outputFileDir + country_iso ,  
                               data=sim_output[sim_output$intervention == sim_output$intervention[1],],mean, na.omit = NULL)

  # calculate RSV incidence (absolute)
  country_summary$rsv_inc    <- country_summary$rsv_cases/country_summary$population
  
  # set income region back to character value
  country_summary$income_region <- levels(sim_output$income_region)[country_summary$income_region]
  
  ##############################
  ## Willingness-to-pay values #
  ##############################
  
  # settings
  wtp_stepsize       <- 500
  wtp_max            <- 30000 
  plot_x_unit        <- 'USD'
  plot_grid_unit     <- 5000

  # get local copy of WTP summary statistics
  opt_wtp            <- seq(0,wtp_max,wtp_stepsize)
  num_wtp            <- length(opt_wtp)
  
  # summary matrix 
  wtp_country_data        <- data.frame(matrix(NA,num_wtp,num_countries))
  names(wtp_country_data) <- opt_countries
  wtp_country_data_prob   <- wtp_country_data
   
  ## PER COUNTRY
  i_country <- 1
  for(i_country in 1:num_countries)
  {
    # select country data
    sel_country  <- opt_countries[i_country]
    country_data <- sim_output[sim_output$country_iso == sel_country,]
 
    # get all net benefit results
    out_tmp                    <- get_all_net_benefit_results(country_data,opt_wtp)
    prob_high_mean_net_benefit <- out_tmp$prob_high_mean_net_benefit  # CEAF
    net_benefit_legend         <- out_tmp$net_benefit_legend          # legend
    
    # select preferred strategy
    flag                               <- !is.na(prob_high_mean_net_benefit)
    intervention_id                    <- apply(flag,2,which)
    tmp                                <- unlist(lapply(intervention_id,length))
    wtp_country_data[tmp==1,i_country] <- intervention_id
    
    # store the probability of the preferred strategy
    intervention_prob                 <- prob_high_mean_net_benefit[flag]
    wtp_country_data_prob[,i_country] <- intervention_prob
  }
  
  # add CEA results
  bool_comparator  <- colSums(wtp_country_data ==  which(net_benefit_legend$name=='comparator'))
  bool_comparator[bool_comparator==0] <- NA
  country_summary$wtp_maternal        <- opt_wtp[bool_comparator]

  # sort on RSV incidence
  country_order <- data.frame(index            = order(country_summary$income_region,country_summary$rsv_inc,country_summary$wtp_maternal),
                              incidence        = round(country_summary$rsv_cases/country_summary$population*1000,digits=1),
                              country_iso      = country_summary$country_iso,
                              income_region    = country_summary$income_region,
                              stringsAsFactors = F)

  # plot title
  flag_maternal <- sim_output$intervention == "maternal"
  maternal_eff_tag <- paste0('efficacy ',round(mean(sim_output$efficacy_maternal_primary[flag_maternal])*100),'%')
  
  # add quantiles if sampled from distribution
  if(any(sim_output$efficacy_maternal_primary[flag_maternal]>0)) {
    maternal_eff_tag <- paste0(maternal_eff_tag,' [',paste0(round(quantile(sim_output$efficacy_maternal_primary[flag_maternal],c(0.025,0.975))*100),collapse=';'),']')
  }
  
  # add hospital efficacy if specified
  if(any(sim_output$efficacy_maternal_hospital>0)){
    maternal_eff_tag <- paste0('overall ',maternal_eff_tag,', hosp. efficacy ',
                               round(mean(sim_output$efficacy_maternal_hospital[flag_maternal])*100), '% [',
                               paste0(round(quantile(sim_output$efficacy_maternal_hospital[flag_maternal],c(0.025,0.975))*100),collapse=';'),']')
  }
  
  plot_title <- paste0('maternal: ', round(max(sim_output$dur_protection_maternal)*12,digits=2),'m protection, ',
                       maternal_eff_tag,', ',
                       max(sim_output$price_dose_maternal),'USD',
                       '\n',
                       'mAb: ', round(max(sim_output$dur_protection_infant)*12,digits=2),'m protection, ',
                       'efficacy ',max(sim_output$efficacy_infant)*100,'%, ',
                       max(sim_output$price_dose_mAb),'USD')
  
  # plot
  plot_filename <- paste0('preferred_strategy_country_table_',unique(sim_output$scenario),'.pdf')
  pdf(file=file.path(plot_output_folder,plot_filename),width=9,height=12)
  par(mar=c(5,7,5,2.5))
  plot_prefered_strategy_table(wtp_country_data,
                               wtp_country_data_prob,
                               net_benefit_legend,
                               opt_wtp,
                               country_order,
                               plot_x_unit,
                               plot_grid_unit,
                               plot_title)
  
  dev.off()
  

 }

# f_wtp_country_data <- wtp_country_data
# f_wtp_country_data_prob <- wtp_country_data_prob
# f_net_benefit_legend <- net_benefit_legend
# f_opt_wtp <- opt_wtp
# f_country_order <- country_order
# f_plot_x_unit <- plot_x_unit
# f_plot_grid_unit <- plot_grid_unit
# f_plot_title <- 'debug'
plot_prefered_strategy_table <- function(f_wtp_country_data,
                                         f_wtp_country_data_prob,
                                         f_net_benefit_legend,
                                         f_opt_wtp,
                                         f_country_order,
                                         f_plot_x_unit,
                                         f_plot_grid_unit,
                                         f_plot_title){
 
  # add additional info to country name
  names(f_wtp_country_data) <- paste0(names(f_wtp_country_data),' [',format(f_country_order$incidence,digits=3),']')
  
  # reorder
  f_wtp_country_data      <- f_wtp_country_data[,f_country_order$index]
  f_wtp_country_data_prob <- f_wtp_country_data_prob[,f_country_order$index]
  f_country_income_region <- f_country_order$income_region[f_country_order$index]

  # convert list into matrix
  num_countries        <- ncol(f_wtp_country_data)
  num_wtp              <- nrow(f_wtp_country_data)
  f_wtp_country_matrix <- matrix(unlist(f_wtp_country_data),ncol=num_countries,nrow=num_wtp)
  f_wtp_country_prob_matrix <- matrix(unlist(f_wtp_country_data_prob),ncol=num_countries,nrow=num_wtp)
  
  country_y_value <- 1:num_countries
  country_y_value <- country_y_value + (f_country_income_region =='LMIC')*2
  
  # store a copy of the intervention colors and legend
  image_colors <- f_net_benefit_legend$color
  image_legend <- f_net_benefit_legend$name_legend
  
  # plot configuration
  plot_lwd        <- cut(f_wtp_country_prob_matrix,seq(0,1,0.25))
  plot_lwd_value  <- as.numeric(plot_lwd)/2.5 
  plot_lwd_pch    <- seq(nlevels(plot_lwd))/2.5 
  plot_lwd_levels <- levels(plot_lwd)
  plot_lwd_levels[c(1,4)] <- c(expression(""<="0.25"),expression("">"0.75"))
  
  plot(x=rep(f_opt_wtp,num_countries),
       y=rep(country_y_value,each=length(f_opt_wtp)),
       col=image_colors[c(f_wtp_country_matrix)],
       pch=15,
       cex=plot_lwd_value,
       yaxt='n',
       ylim=c(2,max(country_y_value)+1),
       ylab='',
       xlab=paste0('Willingness to pay for a DALY averted (',f_plot_x_unit,')'),
       main=f_plot_title)
  
  # add axis with country ISO codes
  axis(2,at=country_y_value,labels=names(f_wtp_country_data),las=2,cex.axis=0.8)
  mtext('Country [RSV incidence per 1000 person years]',side=2,padj=-8)
  
  # add tag for LIC and LMIC
  mean_y_values <- c(mean(country_y_value[f_country_income_region=='LIC']),
                     mean(country_y_value[f_country_income_region=='LMIC']))
  axis(4,mean_y_values,
       c('Low-Income Countries (LIC)','Lower-Middle-Income Countries (LMIC)'),
       tick=F,line=-0.5)
  
  # add grid lines
  abline(v=seq(0,max(f_opt_wtp),f_plot_grid_unit),lty=3,col=9,lwd=2)

  # enable to plot outside the figure box
  par(xpd=TRUE)
  
  # add legend with interventions (colors)
  legend('topleft',
         rev(image_legend),
         ncol=nrow(f_net_benefit_legend),
         x.intersp = 0.5,
         pch=15,
         col=rev(image_colors),
         title='Preferred strategy',
         bg='white',
         inset=c(-0.02,-0.02))
  
  # add legend with probability coding (size)
  legend('topright', 
         plot_lwd_levels, 
         col=1, 
         pch=15, 
         pt.cex=plot_lwd_pch, 
         ncol=4, 
         x.intersp = 0.5,
         title='Degree of certainty',
         bg='white',
         inset=c(-0.02,-0.02)) 
  
  # disable to plot outside the figure box
  par(xpd=FALSE)
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

