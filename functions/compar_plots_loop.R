list_compar_plots <- list()
for (k_plot in 1:3) {
  sel_vars <- list(c("total_YLD","total_YLL","hosp_YLD","non_hosp_YLD","total_DALY"),
                   c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),
                   c("admin_cost","cost_rsv_hosp","cost_rsv_outpatient", # ,"hosp_cost"
                     "outpatient_cost","total_medical_cost"))[[k_plot]]
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted"))) %>% 
    #  & ((price==3&intervention=="maternal")|(price==6&intervention=="mAb"))
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
           burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
           source=ifelse(grepl("projection",as.character(source)),
                         "projection (from [Shi 2017])",as.character(source)),
           vartype=gsub("hosp","hospitalised",
                        gsub("rsv","RSV",gsub("_"," ",gsub("_averted","",variable)))) ) %>% 
    mutate(vec=(10^floor(log10(median/norm_median)))*round(
      median/norm_median/(10^floor(log10(median/norm_median))),3)) %>%
    mutate(orig_burden_round=ifelse(vec>1e6,paste0(round(vec/1e6,1),"e6"),
                                    ifelse(vec>1e3,paste0(round(vec/1e3,1),"K"),vec))) %>% 
    filter(grepl("averted",burden_interv))
  if (any(grepl("cost",sel_vars))) {
    df_plot <- df_plot %>% mutate(vartype=gsub("outpatient","outpatient care",gsub("cost RSV","cost of RSV",
                                                    gsub("hospitalised","hospitalisation",vartype))))
  }
  # labels
  if (any(df_plot$vartype %in% "total_DALY")) {
    df_plot$vartype=factor(df_plot$vartype,levels=
                             unique(df_plot$vartype)[c(which(!grepl("total_DALY",unique(df_plot$vartype))),
                                                       which(grepl("total_DALY",unique(df_plot$vartype))))])}
  #
  dodge_val=1; round_val=2; 
  caption_txt=paste0("Numbers above bars are pre-intervention values",
                     ifelse(any(grepl("YLL",sel_vars))," DALYs",
                            ifelse(any(grepl("death",sel_vars))," case/death numbers"," costs (USD)") ),
                     ifelse(any(grepl("YLL",sel_vars)),". YLD=years lived with disability. YLL=years of life lost",""))
  ylab_txt=paste0("% reduction in ",ifelse(any(grepl("YLL",sel_vars)),"DALYs",
                                           ifelse(any(grepl("death",sel_vars)),"cases/deaths","cost"))," (mean, CI95)")
  ### ### ### ### ###
  
  if (k_plot==1) { df_comparison_reductions_KEN_ZAF <- df_plot } else {
    df_comparison_reductions_KEN_ZAF <- bind_rows(df_plot,df_comparison_reductions_KEN_ZAF) }
  ### ### ### ### ###
  # plot with cntr on x-axis, MV/mAb as colors
  list_compar_plots[[k_plot]] <- ggplot(df_plot) + 
    geom_hpline(aes(x=country_iso,y=norm_median*1e2,group=interaction(source,intervention),
                color=intervention,linetype=source),position=position_dodge(width=dodge_val),
                width=ifelse(exists("width_val"),width_val,0.22),size=1.2) + 
    geom_linerange(aes(x=country_iso,ymin=norm_CI95_low*1e2,ymax=norm_CI95_high*1e2,
                group=interaction(source,intervention),color=intervention),
      alpha=0.35,position=position_dodge(width=dodge_val),size=ifelse(exists("linerange_val"),linerange_val,17),show.legend=F) +
    facet_wrap(~vartype) + scale_y_continuous(breaks=(0:10)*10) + scale_color_manual(values=c("red","blue")) +
    geom_vline(xintercept=1.5,size=0.3) + geom_vline(xintercept=c(1,2),linetype="dashed",size=0.3) +
    geom_text(aes(x=country_iso,y=norm_CI95_high*1e2+2,group=interaction(source,intervention),
                  label=orig_burden_round ), # ifelse(intervention!="MV",orig_burden_round,"")
              position=position_dodge(width=dodge_val),size=ifelse(exists("geom_text_font_size"),geom_text_font_size,4)) +
    labs(color="",linetype="",caption=caption_txt) + scale_x_discrete(expand=expansion(0.02,0)) +
    xlab("") + ylab(ylab_txt) + theme_bw() + standard_theme + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),
                                        axis.text.y=element_text(size=15),strip.text=element_text(size=14),legend.position="top",
                                        legend.text=element_text(size=14),axis.title.y=element_text(size=18)) + 
    guides(linetype=guide_legend(nrow=2),color=guide_legend(nrow=2))
  # save
  compar_dirname <- paste0("output/cea_plots/",subfolder_name,"comparisons/")
  if (!dir.exists(compar_dirname)) { dir.create(compar_dirname) }
  if (SAVE_FLAG){
  ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/",ifelse(any(grepl("YLL",sel_vars)),"DALY",
        ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),width=36,height=18,units="cm")
    }
} # end of for-loop for diff vars
