
for (k_par in 1:nrow(par_table)) {
  n_cntr_output <- par_table$n_cntr_output[k_par]; n_interv <- par_table$n_interv[k_par]
  # intervention config table
  sel_interv <- sim_config_matrix[which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])[n_interv],]
  if (cntrs_cea[n_cntr_output]=="ZAF"){
    sel_interv$country_iso=cntrs_cea[n_cntr_output]}
  # modify duration
  if (n_interv==1) { sel_interv$dur_protection_maternal=3/12 
  } else { sel_interv$dur_protection_infant<-5/12 }
  # half-life
  half_life <- ifelse(n_interv==1,36.5,59.3)
  # lower coverage
  if (lower_cov) {if (n_interv==1) { sel_interv$coverage_maternal <- lower_cov_val } else {
    sel_interv$coverage_infant=lower_cov_val} }
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # calculate
  if (grepl("calc",read_calc_flag)) {
    # loop through PRICE levels
    for (k_price in 1:length(pricelist[[n_interv]])) {
      doseprice=c("mat_vacc"=ifelse(n_interv==1,pricelist$mat_vacc[k_price],pricelist$mat_vacc[1]),
                  "mAb"=ifelse(n_interv==2,pricelist$mAb[k_price],pricelist$mAb[1]))
      # SARI, ARI distinguished, hosp/nonhosp distinguished
      if (cntrs_cea[n_cntr_output]=="ZAF") { cost_input <- list_SA_costs } else { cost_input<-list_KEN_costs }
      # input death data if available
      if (kenya_deaths_input) { print("using propr death data for Kenya")
        kenya_nonhosp_hosp_incid_ari_sari$deaths$hosp=kenya_deaths_incid$hosp
        kenya_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=kenya_deaths_incid$nonhosp } else {
          kenya_nonhosp_hosp_incid_ari_sari$deaths=NULL}
      if (SA_deaths_input) { 
        print("using proprietary death data for SA")
        sa_nonhosp_hosp_incid_ari_sari$deaths$hosp=SA_deaths_incid$hosp
        sa_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=SA_deaths_incid$nonhosp 
        } else {
          sa_nonhosp_hosp_incid_ari_sari$deaths=NULL}
      # RUN CALCULATIONS
      sim_output_user_input_ari_sari <- get_burden_flexible_ari_sari(sel_interv,
                                         list(kenya_nonhosp_hosp_incid_ari_sari,sa_nonhosp_hosp_incid_ari_sari)[[n_cntr_output]],
                                         efficacy_figures,effic_prob=T,effic_distr=effic_dist_fit,list_effic_fit=list_effic_betafit,
                                         exp_wane=exp_wane_val,list_exp_waning_param,dose_price=doseprice,cost_data=cost_input)
      # calculation with projections from mcmarcel (community-based)
      sim_output=sim_output_user_input_ari_sari # get_burden_flexible(sel_interv,NA,NA,exp_wane=exp_wane_val,doseprice)
      sim_output[1:nrow(sim_output_user_input_ari_sari),1:ncol(sim_output_user_input_ari_sari)]=NA
      # sim_output=data.frame()
      ### processing samples -> mean, median, CIs
      # ICER with discounted DALYs
      with_discounted_DALYs <- fcn_process_burden_output(
        user_output=sim_output_user_input_ari_sari,default_output=sim_output, # 
        sel_cntr=sel_interv$country_iso,cols_burden_sel="",
        plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",
                      own="new data (hospital-based + HUS)"),
        icercolname="incremental_cost/DALY_averted",
        icercols=c("incremental_cost","total_DALY_disc_averted")) %>%
        filter(variable %in% "incremental_cost/DALY_averted") %>% 
        mutate(variable="incremental_cost/DALY_disc_averted")
      # every other outcome (appending `with_discounted_DALYs` to it)
      burden_mcmarcel_owndata <- bind_rows(fcn_process_burden_output(
        user_output=sim_output_user_input_ari_sari,default_output=sim_output,
        sel_cntr=sel_interv$country_iso,cols_burden_sel="",
        plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",
                      own="new data (hospital-based + HUS)"),
        icercolname="incremental_cost/DALY_averted",
        icercols=c("incremental_cost","total_DALY_averted")),
        with_discounted_DALYs) %>%
        mutate(name_root=gsub("_averted","",variable)) %>% group_by(source,iter,name_root) %>% 
        mutate(value_norm=ifelse(grepl("avert",variable)&!grepl("incremental",variable),
                                 value[grepl("avert",variable)]/value[!grepl("avert",variable)],NA)) %>% ungroup() %>% select(!name_root)
      
      # calculate mean, median, CI95
      x=burden_mcmarcel_owndata %>% group_by(source,variable) %>% 
        summarise(mean=mean(value,na.rm=T),median=median(value,na.rm=T),
                  CI50_low=quantile(value,probs=ci50_range,na.rm=T)[1],
                  CI50_high=quantile(value,probs=ci50_range,na.rm=T)[2],
                  CI95_low=quantile(value,probs=ci95_range,na.rm=T)[1],
                  CI95_high=quantile(value,probs=ci95_range,na.rm=T)[2],
                  norm_mean=mean(value_norm,na.rm=T),norm_median=median(value_norm,na.rm=T),
                  norm_CI50_low=quantile(value_norm,probs=ci50_range,na.rm=T)[1],
                  norm_CI50_high=quantile(value_norm,probs=ci50_range,na.rm=T)[2],
                  norm_CI95_low=quantile(value_norm,probs=ci95_range,na.rm=T)[1],
                  norm_CI95_high=quantile(value_norm,probs=ci95_range,na.rm=T)[2]) %>%
        mutate(price=doseprice[n_interv],source_num=c(0.2,2)[as.numeric(source)],
               country_iso=sel_interv$country_iso,intervention=sel_interv$intervention)
      # collect all outputs
      if (k_price==1) {cea_summary=x} else {cea_summary=bind_rows(cea_summary,x)}  } # END price scan
    # save CEA summary table
    folder_path=paste0("output/cea_plots/",subfolder_name)
    if (!dir.exists(folder_path)) {dir.create(folder_path)}
    write_csv(cea_summary,paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,
                                 "_cea_summary_mean_CI_",gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv") ) } 
  # LOAD EXISTING DATA IF already available
  else { 
    cea_summary <- read_csv(paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,
                                   "_cea_summary_mean_CI_",gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv")) 
  }
} # loop country
# stop cluster