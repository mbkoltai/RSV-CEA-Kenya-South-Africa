source("functions/get_burden_flexible_ari_sari.R")
# loop thru: cntrs * interventions * dose prices
# n_cntr_output=1:length(cntrs_cea); n_interv=1:2
par_table=expand_grid(n_cntr_output=1:length(cntrs_cea),n_interv=1:2); read_calc_flag=c("calc","read")[2]
subfolder_name="new_price_efficacy/"
cl=parallel::makeCluster(8); registerDoParallel(cl)
foreach (k_par=1:nrow(par_table),.packages=c("dplyr","ggplot2","tidyr","readr")) %dopar% {
    n_cntr_output=par_table$n_cntr_output[k_par]; n_interv=par_table$n_interv[k_par]
    # intervention config table
    sel_interv=sim_config_matrix[which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])[n_interv],]
    if (cntrs_cea[n_cntr_output]=="ZAF"){sel_interv$country_iso=cntrs_cea[n_cntr_output]}
    ######
    # calculate
    if (grepl("calc",read_calc_flag)) {
      for (k_price in 1:length(pricelist[[n_interv]])) {
        doseprice=c("mat_vacc"=ifelse(n_interv==1,pricelist$mat_vacc[k_price],pricelist$mat_vacc[1]),
                    "mAb"=ifelse(n_interv==2,pricelist$mAb[k_price],pricelist$mAb[1]))
        # calculation with data from mcmarcel (community-based)
        sim_output=get_burden_flexible(sel_interv,NA,NA,doseprice)
        # SARI, ARI distinguished, hosp/nonhosp distinguished
        if (cntrs_cea[n_cntr_output]=="ZAF") {
          cost_input <- list("inpatient"=s_afr_inpatient_cost %>% filter(name %in% "total"),"outpatient"=s_afr_outpatient_cost)} else {
            cost_input<-NA }
    sim_output_user_input_ari_sari=get_burden_flexible_ari_sari(sel_interv,
          list(kenya_nonhosp_hosp_incid_ari_sari,sa_nonhosp_hosp_incid_ari_sari)[[n_cntr_output]],
          efficacy_figures,doseprice,cost_data=cost_input)
  ### Plot results -------------------------
  # are nonhosp SARIs accounted as SARIs? (or as nonhosp = mild cases)
  burden_mcmarcel_owndata=fcn_process_burden_output(user_output=sim_output_user_input_ari_sari,default_output=sim_output,
    sel_cntr=sel_interv$country_iso,cols_burden_sel="",
    plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",own="new data (hospital-based + HUS)"),
    icercolname="incremental_cost/DALY_averted",icercols=c("incremental_cost_disc","total_DALY_disc_averted"))
    # calculate mean, median, CI95
    x=burden_mcmarcel_owndata %>% group_by(source,variable) %>% summarise(mean=mean(value,na.rm=T),median=median(value,na.rm=T),
          CI50_low=quantile(value,probs=c(25,75)/1e2,na.rm=T)[1],CI50_high=quantile(value,probs=c(25,75)/1e2,na.rm=T)[2],
          CI95_low=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[1],CI95_high=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[2]) %>%
          mutate(price=doseprice[n_interv],source_num=c(0.2,2)[as.numeric(source)],
                 country_iso=sel_interv$country_iso,intervention=sel_interv$intervention)
        if (k_price==1) {cea_summary=x} else {cea_summary=bind_rows(cea_summary,x)}  } # END price scan
      # save CEA summary table
      write_csv(cea_summary,paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                                   gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv") )
    } else { cea_summary<-read_csv(paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                                          gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv")) }
    # plot
    for (k_col in 1:4){
      df_plot=subset(cea_summary, variable %in% list(selvars,all_cols,burden_cols,cost_cols)[[k_col]])
      if (k_col==3) {df_plot=subset(df_plot,price==min(price) ); x_dodge_val=0.1 } else {x_dodge_val=0.35}
      # PLOT
      p <- ggplot(df_plot,aes(x=factor(price))) + 
        # geom_linerange(aes(ymin=CI95_low,ymax=CI95_high,group=interaction(source,price),color=factor(source)),
        # alpha=0.3,position=position_dodge(width=x_dodge_val),size=3) + 
        geom_linerange(aes(ymin=CI50_low,ymax=CI50_high,group=interaction(source,price),color=factor(source)),alpha=0.5,
          position=position_dodge(width=x_dodge_val),size=ifelse(k_col==2,3,4)) +
        geom_point(aes(y=mean,group=interaction(source,price)),pch="-",size=ifelse(k_col==2,12,14),position=position_dodge(width=x_dodge_val)) + 
        facet_wrap(~variable,scales="free") +
        geom_rect(data=subset(df_plot,grepl("cost",variable)),fill=NA,colour="blue",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
 # geom_rect(data=subset(df_plot,!grepl("cost|averted",variable)),fill=NA,colour="red",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
        geom_rect(data=subset(df_plot,grepl("averted",variable)),fill=NA,colour="green",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
        theme_bw() + standard_theme + labs(color="data source") + scale_x_discrete(expand=c(0,ifelse(k_col!=3,0.4,0.1))) + 
 theme(legend.position="bottom",axis.text.x=element_text(vjust=0.5,size=12,angle=0),axis.text.y=element_text(size=ifelse(k_col==2,8,11)),
  strip.text=element_text(size=ifelse(k_col==2,8,12)),legend.text=element_text(size=11)) + xlab("dose price") + ylab("mean (CI50, CI95)") +
        geom_text(aes(x=factor(price),y=mean,group=interaction(source,price),
            label=ifelse(mean>0,paste0(round(mean/(10^floor(log10(mean))),1),"e",floor(log10(mean))),
          paste0(round(mean/(10^floor(log10(abs(mean)))),1),"e",floor(log10(abs(mean))))) ),size=ifelse(k_col==2,3,4.5),#check_overlap=T,
            position=position_dodge(width=ifelse(k_col==3,0.13,ifelse(k_col!=2,0.9,1.05))),angle=ifelse(k_col!=2,90,90),show.legend=F) +
        ggtitle(paste0(cntrs_cea[n_cntr_output]," cost effectiveness for: ",gsub("maternal","maternal vaccination",sel_interv$intervention)))
      if (min(df_plot$mean,na.rm=T)>0) {p=p+scale_y_log10(expand=expansion(0.19,0));p} else {p=p+scale_y_continuous(expand=expansion(0.3,0));p}
      if (length(unique(df_plot$price))==1) {p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + xlab("")+
        ggtitle(paste0(cntrs_cea[n_cntr_output]," disease burden")); p}
      # filename
      cea_summ_plot_filename=paste("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                                   gsub("maternal","Mat_Vacc",sel_interv$intervention),".png",sep="") # calc_tag
      if (identical(unique(df_plot$variable), sort(selvars))) {cea_summ_plot_filename=gsub(".png","_SELVARS.png",cea_summ_plot_filename)}
      if (identical(unique(df_plot$variable), sort(burden_cols))){
        cea_summ_plot_filename=paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_burden.png")}
      if (identical(unique(df_plot$variable), sort(cost_cols))){cea_summ_plot_filename=gsub(".png","_COSTS.png",cea_summ_plot_filename)}
      # SAVE
      ggsave(cea_summ_plot_filename,width=36,height=18,units="cm")
    }
} # loop country

# stop cluster
stopCluster(cl)
