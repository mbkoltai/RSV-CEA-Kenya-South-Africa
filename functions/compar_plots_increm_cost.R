ggplot(df_total_DALY_medcost_averted_KEN_ZAF %>% filter(grepl("incremental",variable))) + 
  geom_hpline(aes(x=country_iso,y=median,group=interaction(source,intervention),
                  linetype=source),position=position_dodge(width=dodge_val),width=0.23,size=2) + 
  geom_linerange(aes(x=country_iso,ymin=CI95_low,ymax=CI95_high,group=interaction(source,intervention),
                     color=price_interv),alpha=0.4,position=position_dodge(width=dodge_val),size=44) +
  facet_wrap(~variable,scales="free_y",ncol=1) + # ,nrow=2
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),
                              colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=1.5,size=0.5) + geom_vline(xintercept=c(1,2),linetype="dashed",size=1) +
  geom_vline(xintercept=c(0.75,1.25,1.75,2.25),linetype="dashed",size=0.25) + 
  geom_hline(yintercept=0,linetype="dashed",size=0.5) +
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=ifelse(median>=0,CI95_high+2,CI95_low-2),
                group=interaction(source,intervention),label=orig_burden_round),
            position=position_dodge(width=dodge_val),size=7) +
  scale_x_discrete(expand=expansion(0.02,0)) + labs(color="",linetype="") +
  guides(color=guide_legend(ncol=2,override.aes=list(size=4)),
         linetype=guide_legend(nrow=2,override.aes=list(size=1.4))) + #
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=20),axis.text.y=element_text(size=20),
        strip.text=element_text(size=22),legend.position="top",legend.text=element_text(size=17),
        axis.title.x=element_text(21),axis.title.y=element_text(size=18),
        plot.caption=element_text(size=10))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/incremental cost_KEN_ZAF.png"),
       width=35,height=25,units="cm")