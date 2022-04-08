price_interv_vals <- gsub("maternal","MV",unique(paste0(unique(cea_summary_all$price),
                                                        "$ (",unique(cea_summary_all$intervention),")")))
# process dataframe for plotting
df_total_DALY_medcost_averted_KEN_ZAF <- cea_summary_all %>% select(!(contains("norm")|plot_variable)) %>%
  filter(variable %in% c("total_DALY_averted","total_DALY_disc_averted",
                         "total_medical_cost_averted","incremental_cost")) %>% #
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention), # grepl("new",source)
         vec=sign(median)*abs(round((10^floor(log10(abs(median))))*round(
           median/(10^floor(log10(abs(median)))),3))),
         price_interv=paste0(price,"$ (",intervention,")"), # factor(,levels=price_interv_vals), 
         variable=factor(gsub("_"," ",variable),
                         levels=gsub("_"," ",c("total_DALY_averted","total_DALY_disc_averted",
                                               "total_medical_cost_averted","incremental_cost")) )) %>%
  mutate(orig_burden_round=ifelse(abs(vec)>1e6,paste0(round(vec/1e6,1),"e6"),
                                  ifelse(abs(vec)>1e3,paste0(round(vec/1e3,1),"e3"),vec)),
         price_interv=factor(price_interv,levels=unique(price_interv))) %>%
  filter(!((grepl("total_medical",variable)|grepl("total_DALY",variable)) & price>6)) %>%
  pivot_longer(!c(price,source_num,country_iso,intervention,source,variable,
                  country_plot,vec,price_interv,orig_burden_round)) %>%
  group_by(variable) %>% mutate(money_unit=ifelse(max(value)>1e6,"million USD","USD")) %>% ungroup() %>%
  mutate(value=ifelse(grepl("mill",money_unit),value/1e6,value/1e3),
         variable=ifelse(grepl("mill",money_unit),paste0(as.character(variable)," (million US$)"),
                         paste0(as.character(variable)," (thousand US$)"))) %>% 
  mutate(variable=factor(variable,levels=unique(variable)[c(2,3,4,1)])) %>%
  pivot_wider(names_from=name,values_from=value) %>%
  mutate(orig_burden_round=ifelse(grepl("mill",money_unit),paste0(round(median,1),"m"),paste0(round(median,1),"k")),
         source=gsub(", community-based","",source)) %>%
  mutate(across(where(is.numeric),round,3))
ylab_txt="cost in USD (median, CI50, CI95)"
# plot total DALY AVERTED, total medical cost averted
p_compar_plot_cost_red <- ggplot(df_total_DALY_medcost_averted_KEN_ZAF %>% filter(!grepl("incremental",variable))) + 
  geom_hpline(aes(x=country_iso,y=median,group=interaction(source,intervention),
                  linetype=source),position=position_dodge(width=dodge_val),width=0.23,size=2) + # color=price_interv,
  geom_linerange(aes(x=country_iso,ymin=CI95_low,ymax=CI95_high,group=interaction(source,intervention),
                     color=intervention),alpha=0.25,position=position_dodge(width=dodge_val),size=42) +
  facet_wrap(~variable,scales="free_y",ncol=1) + scale_color_manual(values=c("red","blue")) +
  geom_vline(xintercept=1.5,size=0.5) + geom_vline(xintercept=c(1,2),linetype="dashed",size=0.75) +
  geom_vline(xintercept=c(0.75,1.25,1.75,2.25),linetype="dashed",size=0.25) + 
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=CI95_high+2,
                group=interaction(source,intervention),label=orig_burden_round),
            position=position_dodge(width=dodge_val),size=8) +
  scale_x_discrete(expand=expansion(0.02,0)) + labs(color="",linetype="") +#,caption="Numbers show median values."
  guides(color=guide_legend(ncol=2,override.aes=list(size=5)),
         linetype=guide_legend(nrow=2,override.aes=list(size=1.25))) +
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=18),axis.text.y=element_text(size=18),
        strip.text=element_text(size=22),legend.position="top",legend.text=element_text(size=21),
        plot.caption=element_text(10),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24))
# save
if (SAVE_FLAG) {
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/total_DALY_medcost_averted_KEN_ZAF.png"),
       width=30,height=30,units="cm")
  }