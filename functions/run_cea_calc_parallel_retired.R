# plot
#    for (k_col in 1:4){
#      df_plot=subset(cea_summary, variable %in% list(selvars,all_cols,burden_cols,cost_cols)[[k_col]])
#      if (k_col==3) {df_plot=subset(df_plot,price==min(price) ); x_dodge_val=0.1 } else {x_dodge_val=0.35}
#      # PLOT
#      p <- ggplot(df_plot,aes(x=factor(price))) + 
#        # geom_linerange(aes(ymin=CI95_low,ymax=CI95_high,group=interaction(source,price),color=factor(source)),
#        # alpha=0.3,position=position_dodge(width=x_dodge_val),size=3) + 
#        geom_linerange(aes(ymin=CI95_low,ymax=CI95_high,group=interaction(source,price),color=factor(source)),alpha=0.5,
#          position=position_dodge(width=x_dodge_val),size=ifelse(k_col==2,3,4)) +
#    geom_point(aes(y=mean,group=interaction(source,price)),pch="-",size=ifelse(k_col==2,12,14),position=position_dodge(width=x_dodge_val)) +
#        facet_wrap(~variable,scales="free") +
#        geom_rect(data=subset(df_plot,grepl("cost",variable)),fill=NA,colour="blue",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
# # geom_rect(data=subset(df_plot,!grepl("cost|averted",variable)),fill=NA,colour="red",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
#        geom_rect(data=subset(df_plot,grepl("averted",variable)),fill=NA,colour="green",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
#        theme_bw() + standard_theme + labs(color="data source") + scale_x_discrete(expand=c(0,ifelse(k_col!=3,0.4,0.1))) + 
# theme(legend.position="bottom",axis.text.x=element_text(vjust=0.5,size=12,angle=0),axis.text.y=element_text(size=ifelse(k_col==2,8,11)),
#  strip.text=element_text(size=ifelse(k_col==2,8,12)),legend.text=element_text(size=11)) + xlab("dose price") + ylab("mean (CI50, CI95)") +
#        geom_text(aes(x=factor(price),y=mean,group=interaction(source,price),
#            label=ifelse(mean>0,paste0(round(mean/(10^floor(log10(mean))),1),"e",floor(log10(mean))),
#          paste0(round(mean/(10^floor(log10(abs(mean)))),1),"e",floor(log10(abs(mean))))) ),size=ifelse(k_col==2,3,4.5),#check_overlap=T,
#            position=position_dodge(width=ifelse(k_col==3,0.13,ifelse(k_col!=2,0.9,1.05))),angle=ifelse(k_col!=2,90,90),show.legend=F) +
#        ggtitle(paste0(cntrs_cea[n_cntr_output]," cost effectiveness for: ",gsub("maternal","maternal vaccination",sel_interv$intervention)))
#   if (min(df_plot$mean,na.rm=T)>0) {p=p+scale_y_log10(expand=expansion(0.19,0));p} else {p=p+scale_y_continuous(expand=expansion(0.3,0));p}
#   if (length(unique(df_plot$price))==1) {p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + xlab("")+
#        ggtitle(paste0(cntrs_cea[n_cntr_output]," disease burden")); p}
#      plot_folder_path=paste0("output/cea_plots/",subfolder_name,"comparisons/old/")
#    if (!dir.exists(plot_folder_path)) {dir.create(paste0("output/cea_plots/",subfolder_name,"comparisons/")); dir.create(plot_folder_path)}
#      cea_summ_plot_filename=paste(plot_folder_path,sel_interv$country_iso,"_cea_summary_mean_CI_",
#                                   gsub("maternal","Mat_Vacc",sel_interv$intervention),".png",sep="") # calc_tag
#      if (identical(unique(df_plot$variable), sort(selvars))) {cea_summ_plot_filename=gsub(".png","_SELVARS.png",cea_summ_plot_filename)}
#      if (identical(unique(df_plot$variable), sort(burden_cols))){
#        cea_summ_plot_filename=paste0(plot_folder_path,sel_interv$country_iso,"_cea_summary_mean_CI_burden.png")}
#      if (identical(unique(df_plot$variable), sort(cost_cols))){cea_summ_plot_filename=gsub(".png","_COSTS.png",cea_summ_plot_filename)}
#      # SAVE
#      ggsave(cea_summ_plot_filename,width=36,height=18,units="cm")     
#      }