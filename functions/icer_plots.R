levels_except_icer = gsub("disc","(disc)",
                          gsub("_"," ",c("intervention cost","incremental cost",
                          "total medical cost averted","total_DALY_disc_averted","total_DALY_averted")))

for (k_daly in c("incremental_cost/DALY_averted","incremental_cost/DALY_disc_averted")) {
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c("intervention_cost","incremental_cost","total_medical_cost_averted",
                           "total_DALY_disc_averted","total_DALY_averted",k_daly) & grepl("new",source)) %>%
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
           vec=as.character(abs(round((10^floor(log10(abs(median))))*
                                        round(median/(10^floor(log10(abs(median)))),3)))),
           price_interv=factor(paste0(price,"$ (",intervention,")"),
                               levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")"))),
           variable=factor(gsub("disc","(disc)",gsub("_"," ",variable)),
                           levels=c(levels_except_icer,gsub("disc","(disc)",gsub("_"," ",k_daly))) )) %>%
    select(!(contains("norm")|plot_variable)) %>%
    mutate(cnt_int=paste0(country_iso,"\n$",price,"\n",intervention), 
           mean=round(ifelse(!grepl("DALY",variable),mean/1e6,mean),2),
           median=round(ifelse(!grepl("DALY",variable),median/1e6,median),2),
           CI50_low=round(ifelse(!grepl("DALY",variable),CI50_low/1e6,CI50_low),2),
           CI50_high=round(ifelse(!grepl("DALY",variable),CI50_high/1e6,CI50_high),2),
           CI95_low=round(ifelse(!grepl("DALY",variable),CI95_low/1e6,CI95_low),2),
           CI95_high=round(ifelse(!grepl("DALY",variable),CI95_high/1e6,CI95_high),2),
           variable=ifelse(!grepl("DALY",variable),
                           paste0(variable," (million US$)"),paste0(variable, " (US$)"))) %>%
    mutate(cnt_int=factor(cnt_int,levels=unique(cnt_int[order(country_iso,intervention,price)]),ordered=T),
    variable=factor(variable,levels=unique(variable))) 
  
  # save table
  if (k_daly %in% "incremental_cost/DALY_averted"){ 
    df_interv_incremcosts_icer <- df_plot 
    } else {
    df_interv_incremcosts_icer <- bind_rows(df_plot,df_interv_incremcosts_icer) %>% 
      mutate(across(where(is.numeric),round,3)) }
}
 
### ### ### ### ### ### ### ### ###


# scan with different prices
list_scaled <- list()
for (k_scale in 1:length(price_scaling_vect)) {
  price_scaling <- price_scaling_vect[k_scale] 
  # replace values
  xx <- df_interv_incremcosts_icer[!duplicated(df_interv_incremcosts_icer),] %>% 
    select(!c(source_num,source)) %>%
    pivot_longer(c(mean,median,CI50_low,CI50_high,CI95_low,CI95_high)) %>%
    group_by(country_iso,intervention) %>% 
    mutate(value=ifelse(grepl("intervention cost",variable),price_scaling*value,value),
           price=price_scaling*price, price_interv=paste0(price,"$(",intervention,")"),
           cnt_int=paste0(country_iso,"\n$",price,"\n",intervention) ) %>% 
    # paste0(country_iso ," $",price," ",intervention)
    group_by(country_iso,intervention,name) %>%
    mutate(value=ifelse(grepl("incremental cost",variable) & grepl("million",variable),
      value[grepl("intervention cost",variable)]-
        value[grepl("total medical cost averted",variable)],value)) %>%
    mutate(row=row_number()) %>% pivot_wider(names_from=c(name),values_from=value) %>% 
    select(!row) %>% group_by(country_iso,intervention) %>% # ,variable
    # incremental costs: ci50/95 high and low needs to be swapped
    mutate(CI50_low_store=CI50_low,CI95_low_store=CI95_low) %>% 
    mutate(CI50_low=ifelse(grepl("incremental cost \\(million",variable),CI50_high,CI50_low),
           CI50_high=ifelse(grepl("incremental cost \\(million",variable),CI50_low_store,CI50_high),
           CI95_low=ifelse(grepl("incremental cost \\(million",variable),CI95_high,CI95_low),
           CI95_high=ifelse(grepl("incremental cost \\(million",variable),CI95_low_store,CI95_high)) %>%
    select(!c(CI50_low_store,CI95_low_store)) %>%
    mutate(mean=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                       1e6*mean[grepl("incremental cost \\(million",variable)]/
                         mean[grepl("total DALY \\(disc",variable)],mean),
           median=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                         1e6*median[grepl("incremental cost \\(million",variable)]/
                           median[grepl("total DALY \\(disc",variable)],median), 
           CI50_low=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                           1e6*CI50_low[grepl("incremental cost \\(million",variable)]/
                             CI50_high[grepl("total DALY \\(disc",variable)],CI50_low), 
           CI50_high=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                            1e6*CI50_high[grepl("incremental cost \\(million",variable)]/
                              CI50_low[grepl("total DALY \\(disc",variable)],CI50_high),
           CI95_low=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                           1e6*CI95_low[grepl("incremental cost \\(million",variable)]/
                             CI95_high[grepl("total DALY \\(disc",variable)],CI95_low),
           CI95_high=ifelse(grepl("incremental cost/DALY \\(disc",variable),
                            1e6*CI95_high[grepl("incremental cost \\(million",variable)]/
                              CI95_low[grepl("total DALY \\(disc",variable)],CI95_high)) %>%
    # incremental cost/DALY averted
    mutate(mean=ifelse(grepl("incremental cost/DALY averted",variable),
                       1e6*mean[grepl("incremental cost \\(million",variable)]/
                         mean[grepl("total DALY averted",variable)],mean),
           median=ifelse(grepl("incremental cost/DALY averted",variable),
                         1e6*median[grepl("incremental cost \\(million",variable)]/
                           median[grepl("total DALY averted",variable)],median), 
           CI50_low=ifelse(grepl("incremental cost/DALY averted",variable),
                           1e6*CI50_low[grepl("incremental cost \\(million",variable)]/
                             CI50_high[grepl("total DALY averted",variable)],CI50_low), 
           CI50_high=ifelse(grepl("incremental cost/DALY averted",variable),
                            1e6*CI50_high[grepl("incremental cost",variable) & grepl("million",variable)]/
                              CI50_low[grepl("total DALY averted",variable)],CI50_high),
           CI95_low=ifelse(grepl("incremental cost/DALY averted",variable),
                           1e6*CI95_low[grepl("incremental cost",variable) & grepl("million",variable)]/
                             CI95_high[grepl("total DALY averted",variable)],CI95_low),
           CI95_high=ifelse(grepl("incremental cost/DALY averted",variable),
                            1e6*CI95_high[grepl("incremental cost",variable) & grepl("million",variable)]/
                              CI95_low[grepl("total DALY averted",variable)],CI95_high)) %>%
# mutate(cnt_int=factor(cnt_int,levels=unique(cnt_int[order(country_iso,intervention,price)]),ordered=T) ) %>%
        arrange(country_iso,intervention,price) # %>%
# mutate(cnt_int=factor(cnt_int,ordered=T) ) # levels=unique(cnt_int[order(country_iso,intervention,price)]),
    
  list_scaled[[k_scale]] <- xx
  # "6$ (mAb)" "ZAF $6 mAb" # list_scaled[[k_scale]]
}

# plot
# ylab_txt <- "cost in USD (median, CI95)"
# for (k_plot in 1:2) {
#   if (k_plot==2) { 
#     df_plot <- df_plot %>% filter(grepl("averted",variable) & grepl("incremental",variable)) }
#   df_n_price = df_plot %>% group_by(intervention) %>% summarise(price=unique(price)) %>% 
#     group_by(intervention) %>% summarise(n_price=n())
#   df_n_price$n_price
#   
#   # create plot
#     p <- ggplot(df_plot) + geom_boxplot(aes(x=cnt_int,middle=median, color=price_interv,
#                                 lower=CI50_low,upper=CI50_high,ymin=CI95_low,ymax=CI95_high),
#                    position=position_dodge(width=dodge_val),stat="identity",width=0.85) + # ,size=1.1
#     facet_wrap(~variable,scales="free_y",nrow=3) +
#     scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(df_n_price$n_price[1]),
#                                 colorRampPalette(colors=c("blue","blueviolet"))(df_n_price$n_price[2]))) +
#     geom_vline(xintercept=c(3.5,6.5,9.5),linetype="dashed",size=0.3) + 
#     geom_vline(xintercept=c(6.5),size=1/2) + geom_hline(yintercept=0,linetype="dashed",size=1/2) +
#     xlab("") + ylab(ylab_txt) + labs(color="",linetype="") +guides(color=guide_legend(ncol=2)) + 
#     scale_x_discrete(expand=expansion(0.05,0)) + 
#     theme_bw() + standard_theme + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=18),
#                                         axis.text.y=element_text(size=17),strip.text=element_text(size=14),
#                                         legend.text=element_text(size=20),legend.position="top",
#                                         axis.title.y=element_text(size=20),strip.text.x=element_text(size=18)) 
#   if (k_plot==2) { p <- p + scale_y_continuous(breaks=(-3:20)*2*1e3) } # 2.5e3
#   p
# save
# full_filename <- paste0("output/cea_plots/",subfolder_name,ifelse(k_plot==2,"icer","interv_incremcosts_icer"),
#                         ifelse(grepl("disc",k_daly),"_disc",""),"_KEN_ZAF_3rows",
#                         ifelse(any(is.na(df_plot$CI95_high)),"_ci95_removed",""),".png")

# ggsave(full_filename,width=25,height=40/ifelse(k_plot==2,2,1),units="cm")
# }