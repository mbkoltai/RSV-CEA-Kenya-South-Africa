ylab_txt<-"cost in USD (mean, CI95)"
price_interv_vals <- gsub("maternal","MV",
                          unique(paste0(cea_summary_all$price,"$ (",cea_summary_all$intervention,")")))
df_plot_icer_comp <- cea_summary_all %>% select(!(contains("norm")|plot_variable)) %>% 
  filter(grepl("incremental_cost/DALY",variable)) %>%
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
         vec=sign(median)*abs(round((10^floor(log10(abs(median))))*
                                      round(median/(10^floor(log10(abs(median)))),3))),
         price_interv=factor(paste0(price,"$ (",intervention,")"),levels=price_interv_vals),
         variable=gsub("disc","(discounted)",gsub("_"," ",variable)) ) %>%
  mutate(orig_burden_round=paste0(round(median/1e3,2),"K"),source=gsub(", community-based","",source)) %>%
  mutate(across(where(is.numeric),round,1))

n_price= length(unique(df_plot_icer_comp$price))/2
# PLOT
p_icer_comp <- ggplot(df_plot_icer_comp) + 
  geom_hpline(aes(x=country_iso,y=median/1e3,group=interaction(price,source),linetype=source),
              position=position_dodge(width=dodge_val),width=ifelse(exists("width_val"),width_val,0.46),size=1) +
  geom_linerange(aes(x=country_iso,ymin=CI50_low/1e3,ymax=CI50_high/1e3,
                     group=interaction(price,source),color=price_interv),
                 alpha=0.5,position=position_dodge(width=dodge_val),size=ifelse(exists("linerange_val"),linerange_val,42)) +
  facet_grid(intervention~variable,scales="free") + 
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(n_price),
                              colorRampPalette(colors=c("blue","blueviolet"))(n_price))) +
  geom_vline(xintercept=c(1.5),size=0.5) + 
  # geom_vline(xintercept=c(1,2),size=0.75,linetype="dashed") +
  # geom_vline(xintercept=setdiff(seq(2/3,4,by=1/6),c(1,1.5,2)),linetype="dashed",size=0.2) +
  geom_hline(yintercept=0,size=0.25,linetype="dashed") +
  xlab("") + ylab("incremental cost (thousand US$) per DALY averted") + 
  labs(color="",linetype="",caption="*cost-saving") + # 
  geom_text(aes(x=country_iso,y=(CI50_high+2e3)/1e3,group=interaction(price,source),
                label=ifelse(median>=0,orig_burden_round,"*")),
            size=ifelse(exists("geom_text_font_size"),geom_text_font_size,6),position=position_dodge(width=dodge_val)) + 
  scale_x_discrete(expand=expansion(0.2,0)) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=20),
        axis.text.y=element_text(size=18),axis.title.y=element_text(size=19),
        strip.text=element_text(size=21),plot.caption=element_text(size=15),
        legend.position="top",legend.text=element_text(size=16)) + guides(color=guide_legend(ncol=2),linetype=guide_legend(nrow=2))

if (CI95_FLAG) {
  p_icer_comp <- p_icer_comp + 
    geom_linerange(aes(x=country_iso,ymin=CI95_low/1e3,ymax=CI95_high/1e3,group=interaction(price,source),color=price_interv),
                 alpha=0.2,position=position_dodge(width=dodge_val),size=ifelse(exists("linerange_val"),linerange_val,42),
                 show.legend = F)
}