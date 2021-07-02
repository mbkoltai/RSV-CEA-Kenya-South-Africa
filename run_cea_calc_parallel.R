source("functions/get_burden_flexible_ari_sari.R"); source("functions/get_burden_flexible.R")
# loop thru: cntrs * interventions * dose prices
# n_cntr_output=1:length(cntrs_cea); n_interv=1:2
par_table=expand_grid(n_cntr_output=1:length(cntrs_cea),n_interv=1:2); read_calc_flag=c("calc","read")[1]
subfolder_name="new_price_efficacy_kenyadeaths_CIs_expwaning/" # "new_price_efficacy_CIs/"; 
kenya_deaths_input=TRUE; exp_wane_val=TRUE
cl=parallel::makeCluster(8); registerDoParallel(cl)
foreach (k_par=1:nrow(par_table),.packages=c("dplyr","ggplot2","tidyr","readr","rriskDistributions")) %dopar% {
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
        sim_output=get_burden_flexible(sel_interv,NA,NA,exp_wane=exp_wane_val,doseprice)
        # SARI, ARI distinguished, hosp/nonhosp distinguished
        if (cntrs_cea[n_cntr_output]=="ZAF") {
          cost_input <- list("inpatient"=s_afr_inpatient_cost %>% filter(name %in% "total"),"outpatient"=s_afr_outpatient_cost)} else {
            cost_input<-NA }
        # input death data if available
        if (kenya_deaths_input) {print("using propr death data"); kenya_nonhosp_hosp_incid_ari_sari$deaths$hosp=kenya_deaths_incid$hosp
          kenya_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=kenya_deaths_incid$nonhosp } else {kenya_nonhosp_hosp_incid_ari_sari$deaths=NULL}
    sim_output_user_input_ari_sari=get_burden_flexible_ari_sari(sel_interv,
          list(kenya_nonhosp_hosp_incid_ari_sari,sa_nonhosp_hosp_incid_ari_sari)[[n_cntr_output]],
          efficacy_figures,effic_prob=T,exp_wane=exp_wane_val,doseprice,cost_data=cost_input)
  ### Plot results -------------------------
  # are nonhosp SARIs accounted as SARIs? (or as nonhosp = mild cases)
  burden_mcmarcel_owndata=fcn_process_burden_output(user_output=sim_output_user_input_ari_sari,default_output=sim_output,
    sel_cntr=sel_interv$country_iso,cols_burden_sel="",
    plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",own="new data (hospital-based + HUS)"),
    icercolname="incremental_cost/DALY_averted",icercols=c("incremental_cost","total_DALY_averted"))
    # calculate mean, median, CI95
    x=burden_mcmarcel_owndata %>% group_by(source,variable) %>% summarise(mean=mean(value,na.rm=T),median=median(value,na.rm=T),
          CI50_low=quantile(value,probs=c(25,75)/1e2,na.rm=T)[1],CI50_high=quantile(value,probs=c(25,75)/1e2,na.rm=T)[2],
          CI95_low=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[1],CI95_high=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[2]) %>%
          mutate(price=doseprice[n_interv],source_num=c(0.2,2)[as.numeric(source)],
                 country_iso=sel_interv$country_iso,intervention=sel_interv$intervention)
        if (k_price==1) {cea_summary=x} else {cea_summary=bind_rows(cea_summary,x)}  } # END price scan
      # save CEA summary table
      folder_path=paste0("output/cea_plots/",subfolder_name); if (!dir.exists(folder_path)) {dir.create(folder_path)}
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
      plot_folder_path=paste0("output/cea_plots/",subfolder_name,"comparisons/old/")
      if (!dir.exists(plot_folder_path)) {dir.create(paste0("output/cea_plots/",subfolder_name,"comparisons/")); dir.create(plot_folder_path)}
      cea_summ_plot_filename=paste(plot_folder_path,sel_interv$country_iso,"_cea_summary_mean_CI_",
                                   gsub("maternal","Mat_Vacc",sel_interv$intervention),".png",sep="") # calc_tag
      if (identical(unique(df_plot$variable), sort(selvars))) {cea_summ_plot_filename=gsub(".png","_SELVARS.png",cea_summ_plot_filename)}
      if (identical(unique(df_plot$variable), sort(burden_cols))){
        cea_summ_plot_filename=paste0(plot_folder_path,sel_interv$country_iso,"_cea_summary_mean_CI_burden.png")}
      if (identical(unique(df_plot$variable), sort(cost_cols))){cea_summ_plot_filename=gsub(".png","_COSTS.png",cea_summ_plot_filename)}
      # SAVE
      ggsave(cea_summ_plot_filename,width=36,height=18,units="cm")
    }
} # loop country

# stop cluster
stopCluster(cl)

### ### ###
# save in one file
filenames=list.files(paste0("output/cea_plots/",subfolder_name),pattern=".csv"); filenames=filenames[!filenames %in% "cea_summary_all.csv"]
for (k_filename in 1:4) {
  x=read_csv(paste0("output/cea_plots/",subfolder_name,filenames[k_filename]))
  if (k_filename==1){ cea_summary_all=x } else {cea_summary_all=bind_rows(cea_summary_all,x)} }
write_csv(cea_summary_all,paste0("output/cea_plots/",subfolder_name,"cea_summary_all.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Make summary plots (After running calcs in parallel by "run_cea_calc_parallel.R")
# composition of total burden
# total_DALY <- total_YLD + total_YLL
# total_YLL <- rsv_deaths*config$hosp_CFR_DALYloss
# total_YLD <- hosp_med_att_YLD + non_hosp_YLD
# hosp_med_att_YLD <- hosp_SARI * config$severe_rsv_DALYloss + med_att_ARI*config$non_severe_rsv_DALYloss
# hosp_YLD <- hosp_SARI* config$severe_rsv_DALYloss
# non_hosp_YLD <- non_hosp_SARI*config$severe_rsv_DALYloss + non_med_att_ARI*config$non_severe_rsv_DALYloss
for (k_plot in 1:3) {
  sel_vars <- list(c("total_YLD","total_YLL","hosp_YLD","hosp_med_att_YLD","non_hosp_YLD","total_DALY"), # ,"ARI_YLD","SARI_YLD"
                   c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),
                   c("admin_cost","cost_rsv_hosp","hosp_cost","cost_rsv_outpatient","outpatient_cost","total_medical_cost"))[[k_plot]]
  df_plot <- cea_summary_all %>% filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) & # intervention=="mAb" & 
                                          ((price==3&intervention=="maternal")|(price==6&intervention=="mAb")) & grepl("new",source)) %>%
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
           burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
           source=ifelse(grepl("projection",as.character(source)),"projection (from [Shi 2017])",as.character(source)),
           vartype=gsub("_averted","",variable)) %>% group_by(source,vartype,price,country_iso,intervention) %>% 
    mutate(norm_mean=mean/mean[!grepl("averted",variable)], norm_median=median/median[!grepl("averted",variable)],
           norm_CI50_low=CI50_low/mean[!grepl("averted",variable)],norm_CI50_high=CI50_high/mean[!grepl("averted",variable)],
         norm_CI95_low=CI95_low/mean[!grepl("averted",variable)],norm_CI95_high=CI95_high/mean[!grepl("averted",variable)]) %>% ungroup() %>% 
    mutate(vec=as.character((10^floor(log10(median/norm_median)))*round(median/norm_median/(10^floor(log10(median/norm_median))),3)))
  df_plot$orig_burden_round=gsub("^\\.","",
                                 sapply(df_plot$vec, function(vec) paste0(substring(vec,first=c(1,seq(nchar(vec)-floor(nchar(vec)/3)*3,
                          nchar(vec)-1,by=3)+1),last=c(seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by =3),nchar(vec))),collapse=".")))
  df_plot <- df_plot %>% mutate(orig_burden_round=ifelse(as.numeric(vec)<1e4,gsub("\\.","",orig_burden_round),orig_burden_round))
  if (any(df_plot$vartype %in% "total_DALY")) {
    df_plot$vartype=factor(df_plot$vartype,levels=
        unique(df_plot$vartype)[c(which(!grepl("total_DALY",unique(df_plot$vartype))),which(grepl("total_DALY",unique(df_plot$vartype))))])}
  # plot for 2 cntrs, 2 intervents
  dodge_val=1; round_val=2; caption_txt=paste0("Numbers above medians are pre-intervention",ifelse(any(grepl("YLL",sel_vars))," DALYs",
        ifelse(any(grepl("death",sel_vars))," case/death numbers"," costs (USD)") ),
        ifelse(any(grepl("YLL",sel_vars)),". YLD=years lived with disability. YLL=years of life lost",""))
  ylab_txt=paste0("% reduction in ",ifelse(any(grepl("YLL",sel_vars)),"DALYs",
                                           ifelse(any(grepl("death",sel_vars)),"cases/deaths","cost"))," (mean, CI50)")
  # plot with cntr on x-axis, MV/mAb as colors
  ggplot(df_plot %>% filter(grepl("averted",burden_interv)) ) +
    geom_hpline(aes(x=country_iso,y=norm_median*1e2,group=intervention,color=intervention), # ,linetype=source
                position=position_dodge(width=dodge_val),width=0.42,size=1) + # scale_linetype_manual(values=c("solid","longdash"))+
    geom_linerange(aes(x=country_iso,ymin=norm_CI50_low*1e2,ymax=norm_CI50_high*1e2,group=,color=intervention),
                   alpha=0.35,position=position_dodge(width=dodge_val),size=30,show.legend=F) +
    facet_wrap(~vartype) + scale_color_manual(values=c("red","blue")) +
    geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
    geom_text(aes(x=country_iso,y=norm_CI50_high*1e2+2,group=intervention,label=ifelse(intervention!="MV",orig_burden_round,"")),
              position=position_dodge(width=dodge_val),size=5) + labs(color="",linetype="",caption=caption_txt) +
    scale_x_discrete(expand=expansion(0.02,0)) + scale_y_continuous(breaks=(0:10)*10) + 
    theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=15),
          strip.text=element_text(size=14),legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=18))
  # save
  ggsave(paste0("output/cea_plots/",subfolder_name,ifelse(any(grepl("YLL",sel_vars)),"DALY",
      ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),width=36,height=18,units="cm")
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables  # "total_medical_cost_averted"
df_plot <- cea_summary_all %>% filter(variable %in% c("intervention_cost","incremental_cost", "incremental_cost/DALY_averted") &
                                        grepl("new",source)) %>% mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
                vec=as.character(abs(round((10^floor(log10(abs(median))))*round(median/(10^floor(log10(abs(median)))),3)))),
                price_interv=factor(paste0(price,"$ (",intervention,")"),levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")"))),
                variable=factor(variable,levels=c("intervention_cost","incremental_cost", "incremental_cost/DALY_averted")))
df_plot$orig_burden_round=df_plot$vec; df_plot$orig_burden_round[as.numeric(df_plot$vec)>1e4]=sapply(df_plot$vec[as.numeric(df_plot$vec)>1e4], 
      function(vec) paste0(substring(vec, first=c(1,seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by=3)+1),
            last=c(seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by=3),nchar(vec))),collapse="."))
df_plot <- df_plot %>% dplyr::select(!vec) %>% mutate(orig_burden_round=gsub("^\\.","",orig_burden_round)) %>% 
  mutate(orig_burden_round=ifelse(median<0,paste0("-",orig_burden_round),orig_burden_round))
# plot
ylab_txt="cost in USD (mean, CI50)"
ggplot(df_plot) + geom_hpline(aes(x=country_iso,y=median,group=intervention,color=price_interv), # 
                    position=position_dodge(width=dodge_val),width=0.43,size=1) + #scale_linetype_manual(values=c("solid","longdash"))+
  geom_linerange(aes(x=country_iso,ymin=CI50_low,ymax=CI50_high,group=intervention,color=price_interv),
                 alpha=0.35,position=position_dodge(width=dodge_val),size=28,show.legend=F) +
  facet_wrap(~variable,scales = "free") +
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=ifelse(abs(CI50_high)>1e6,CI50_high+2e6,CI50_high+80),group=intervention,label=orig_burden_round),
            position=position_dodge(width=dodge_val)) + labs(color="",linetype="",caption="Numbers show median values.") + 
  scale_x_discrete(expand=expansion(0.02,0)) + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=15),
            strip.text=element_text(size=14),legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=18)) + 
  guides(color=guide_legend(ncol=2))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"incremental_costs_KEN_ZAF.png"),width=36,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plots COMPARING CEA results with projected data vs new data
# subfolder_name<-"new_price_efficacy_kenyadeaths_CIs/"
for (k_plot in 1:3) {
  sel_vars <- list(c("total_YLD","total_YLL","hosp_YLD","hosp_med_att_YLD","non_hosp_YLD","total_DALY"), # ,"ARI_YLD","SARI_YLD"
                   c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),
                   c("admin_cost","cost_rsv_hosp","hosp_cost","cost_rsv_outpatient","outpatient_cost","total_medical_cost"))[[k_plot]]
  df_plot <- cea_summary_all %>% filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) &
                                          ((price==3&intervention=="maternal")|(price==6&intervention=="mAb"))) %>% #  &grepl("new",source)
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
           burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
           source=ifelse(grepl("projection",as.character(source)),"projection (from [Shi 2017])",as.character(source)),
           vartype=gsub("_averted","",variable)) %>% group_by(source,vartype,price,country_iso,intervention) %>% 
    mutate(norm_mean=mean/mean[!grepl("averted",variable)], norm_median=median/median[!grepl("averted",variable)],
           norm_CI50_low=CI50_low/mean[!grepl("averted",variable)],norm_CI50_high=CI50_high/mean[!grepl("averted",variable)],
           norm_CI95_low=CI95_low/mean[!grepl("averted",variable)],norm_CI95_high=CI95_high/mean[!grepl("averted",variable)]) %>% ungroup() %>%
    mutate(vec=(10^floor(log10(median/norm_median)))*round(median/norm_median/(10^floor(log10(median/norm_median))),3)) %>%
    mutate(orig_burden_round=ifelse(vec>1e6,paste0(round(vec/1e6,2),"e6"),ifelse(vec>1e3,paste0(round(vec/1e3,2),"e3"),vec)))
  if (any(df_plot$vartype %in% "total_DALY")) {
    df_plot$vartype=factor(df_plot$vartype,levels=
                             unique(df_plot$vartype)[c(which(!grepl("total_DALY",unique(df_plot$vartype))),which(grepl("total_DALY",unique(df_plot$vartype))))])}
  #
  dodge_val=1; round_val=2; caption_txt=paste0("Numbers above medians are pre-intervention",ifelse(any(grepl("YLL",sel_vars))," DALYs",
                                                                                                   ifelse(any(grepl("death",sel_vars))," case/death numbers"," costs (USD)") ),
                                               ifelse(any(grepl("YLL",sel_vars)),". YLD=years lived with disability. YLL=years of life lost",""))
  ylab_txt=paste0("% reduction in ",ifelse(any(grepl("YLL",sel_vars)),"DALYs",
                                           ifelse(any(grepl("death",sel_vars)),"cases/deaths","cost"))," (mean, CI50)")
  # plot with cntr on x-axis, MV/mAb as colors
  ggplot(df_plot %>% filter(grepl("averted",burden_interv)) ) +
    geom_hpline(aes(x=country_iso,y=norm_median*1e2,group=interaction(source,intervention),color=intervention,linetype=source), # 
                position=position_dodge(width=dodge_val),width=0.22,size=1) + # scale_linetype_manual(values=c("solid","longdash"))+
    geom_linerange(aes(x=country_iso,ymin=norm_CI50_low*1e2,ymax=norm_CI50_high*1e2,group=interaction(source,intervention),
                       color=intervention),alpha=0.35,position=position_dodge(width=dodge_val),size=17,show.legend=F) +
    facet_wrap(~vartype) + scale_y_continuous(breaks=(0:10)*10) + scale_color_manual(values=c("red","blue")) +
    geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
    geom_text(aes(x=country_iso,y=norm_CI50_high*1e2+2,group=interaction(source,intervention),
                  label=ifelse(intervention!="MV",orig_burden_round,"")),position=position_dodge(width=dodge_val)) +
    labs(color="",linetype="",caption=caption_txt)+scale_x_discrete(expand=expansion(0.02,0)) +
    theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=15),strip.text=element_text(size=14),
          legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=18)) + 
    guides(linetype=guide_legend(nrow=2),color=guide_legend(nrow=2))
  # save
  ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/",ifelse(any(grepl("YLL",sel_vars)),"DALY",
                                                                         ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),width=36,height=18,units="cm") 
} # end of for-loop for diff vars
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables  # "total_medical_cost_averted"
df_plot <- cea_summary_all %>% filter(variable %in% c("total_DALY_averted","total_medical_cost_averted","incremental_cost")) %>% #
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention), # grepl("new",source)
         vec=sign(median)*abs(round((10^floor(log10(abs(median))))*round(median/(10^floor(log10(abs(median)))),3))),
         price_interv=factor(paste0(price,"$ (",intervention,")"),levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")"))),
         variable=factor(variable,levels=c("total_DALY_averted","total_medical_cost_averted","incremental_cost"))) %>%
  mutate(orig_burden_round=ifelse(abs(vec)>1e6,paste0(round(vec/1e6,2),"e6"),ifelse(abs(vec)>1e3,paste0(round(vec/1e3,2),"e3"),vec))) %>%
  filter(!((grepl("total_medical",variable)|grepl("total_DALY",variable)) & price>6)) # ,"incremental_cost/DALY_averted"
# plot
ylab_txt="cost in USD (mean, CI50)"
ggplot(df_plot) + geom_hpline(aes(x=country_iso,y=median,group=interaction(source,intervention),color=price_interv,linetype=source),
                              position=position_dodge(width=dodge_val),width=0.22,size=1) + # scale_linetype_manual(values=c("solid","longdash")) +
  geom_linerange(aes(x=country_iso,ymin=CI50_low,ymax=CI50_high,group=interaction(source,intervention),color=price_interv),
                 alpha=0.25,position=position_dodge(width=dodge_val),size=20,show.legend=F) +
  facet_wrap(~variable,scales="free",nrow=2) + 
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=c(1.5),linetype="dashed",size=0.5) + geom_vline(xintercept=c(0.75,1,1.25,1.75,2,2.25),linetype="dashed",size=0.1) +
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=ifelse(abs(median)>1e6,median+3e6,ifelse(median>1e5,median+1e6,median+1e4)),
                group=interaction(source,intervention),label=orig_burden_round),position=position_dodge(width=dodge_val)) + 
  labs(color="",linetype="",caption="Numbers show median values.") + scale_x_discrete(expand=expansion(0.02,0)) + 
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=13),strip.text=element_text(size=14),
        legend.position="top",legend.text=element_text(size=14)) + guides(color=guide_legend(ncol=2),linetype=guide_legend(nrow=2))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/incremental_costs_KEN_ZAF.png"),width=30,height=22,units="cm")

### ### ### ### ### ### ### ### ### ###
# ICER plot
ylab_txt="cost in USD (mean, CI50)"
ggplot(cea_summary_all %>% filter(variable %in% "incremental_cost/DALY_averted") %>%
         mutate(intervention=ifelse(intervention=="maternal","MV",intervention), # grepl("new",source)
                vec=sign(median)*abs(round((10^floor(log10(abs(median))))*round(median/(10^floor(log10(abs(median)))),3))),
                price_interv=factor(paste0(price,"$ (",intervention,")"),levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")")))) %>%
         mutate(orig_burden_round=ifelse(abs(vec)>1e6,paste0(round(vec/1e6,2),"e6"),ifelse(abs(vec)>1e3,paste0(round(vec/1e3,2),"e3"),vec)))) +
  geom_hpline(aes(x=country_iso,y=median,group=interaction(price,source,intervention),color=price_interv,linetype=source),
              position=position_dodge(width=dodge_val),width=0.22/3,size=1) + # scale_linetype_manual(values=c("solid","longdash")) +
  geom_linerange(aes(x=country_iso,ymin=CI50_low,ymax=CI50_high,group=interaction(price,source,intervention),color=price_interv),
                 alpha=0.25,position=position_dodge(width=dodge_val),size=45/3,show.legend=F) +
  facet_wrap(~variable,scales="free",nrow=2) + 
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=c(1.5),linetype="dashed",size=0.5) + geom_vline(xintercept=c(0.75,1,1.25,1.75,2,2.25),linetype="dashed",size=0.1) +
  geom_hline(yintercept=0,linetype="dashed",size=0.25) +
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + labs(color="",linetype="",caption="Numbers show median values.") + 
  geom_text(aes(x=country_iso,y=CI50_high+1e3,group=interaction(price,source,intervention),label=orig_burden_round),
            position=position_dodge(width=dodge_val)) + 
  scale_x_discrete(expand=expansion(0.05,0)) + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=13),
                                                     strip.text=element_text(size=14),legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=15)) +
  guides(color=guide_legend(ncol=2),linetype=guide_legend(nrow=2))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/ICER_KEN_ZAF.png"),width=36,height=22,units="cm")