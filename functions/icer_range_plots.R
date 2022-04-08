# plot range of ICERs as a function of price
# plot theme
plot_theme <- theme(axis.text.x=element_text(angle=90,vjust=1/2,size=17),axis.text.y=element_text(size=17),
                    strip.text=element_text(size=14),strip.text.x=element_text(size=18),
                    legend.text=element_text(size=20),legend.position="top",
                    axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) 

p_icer_plotlist <- list(); k_plot=0

for (k_icer in 1:2) {
  sel_icer_var_pattern <- c("incremental cost/DALY \\(disc","incremental cost/DALY averted")[k_icer]
  df_ylim <- ICER_sensit_price %>% filter(grepl(sel_icer_var_pattern,variable)) %>% 
    filter(!(intervention %in% "MV" & price>price_limit_MV)) %>%
    group_by(intervention,country_plot) %>%
    summarise(CI95_low_sep_cntr=min(CI95_low),CI95_high_sep_cntr=max(CI95_high),
              CI50_low_sep_cntr=min(CI50_low),CI50_high_sep_cntr=max(CI50_high)) %>% 
    group_by(intervention) %>% mutate(CI95_low=min(CI95_low_sep_cntr),CI95_high=max(CI95_high_sep_cntr),
                                      CI50_low=min(CI50_low_sep_cntr),CI50_high=max(CI50_high_sep_cntr))
  # horizontal lines
  hline_vals <- ICER_sensit_price %>% group_by(intervention,country_plot) %>% 
    summarise(zero_limit=0,midlim=3e3,highlim=5e3) %>% pivot_longer(!c(intervention,country_plot)) %>% 
    filter(!(intervention %in% "MV" & value>3e3)) 
  
  # plot
  for (k_ylim in 1:2) {
    ylim_var_name=list(c("CI50_low","CI50_high"),c("CI95_low","CI95_high"))[[k_ylim]]
    p <- ICER_sensit_price %>% filter(grepl(sel_icer_var_pattern,variable)) %>%
      filter(!(intervention %in% "MV" & price>price_limit_MV)) %>%
      ggplot() + geom_line(aes(x=price,y=median)) + 
      geom_ribbon(aes(x=price,ymin=CI50_low,ymax=CI50_high,fill=intervention),alpha=1/2) + 
      geom_point(data=df_ylim,aes(x=10,y=get(ylim_var_name[1])),color="white") + 
      geom_point(data=df_ylim,aes(x=10,y=get(ylim_var_name[2])),color="white") +
      facet_wrap(intervention~country_plot,scales="free",nrow=2) + 
      geom_hline(data=hline_vals,aes(yintercept=value),size=1/3,linetype="dashed") + # ,5e3
      scale_x_continuous(breaks=(0:20)*5,expand=expansion(0.03,0)) + scale_y_continuous(expand=expansion(0.01,0)) +
      xlab("dose price (2019 USD)") + 
      ylab(paste0("ICER (incremental cost per ",ifelse(grepl("disc",sel_icer_var_pattern),"disc. ",""),"DALY averted)")) +
      theme_bw() + standard_theme + plot_theme + theme(legend.title=element_blank())
    if (any(grepl("CI50",ylim_var_name))) { 
      p } else {
        p <- p + geom_ribbon(aes(x=price,ymin=CI95_low,ymax=CI95_high,fill=intervention),alpha=1/5)
      }
    # SAVE
    k_plot=k_plot+1
    p_icer_plotlist[[k_plot]] <- p
    if (SAVE_FLAG){
    ggsave(paste0("output/cea_plots/",subfolder_name,"ICER_",
                   ifelse(grepl("disc",sel_icer_var_pattern),"disc_DALY_",""),
                   "price_scan",ifelse(any(grepl("CI50",ylim_var_name)),"_CI50",""),".png"),width=40,height=30,units="cm")}
  }
}
