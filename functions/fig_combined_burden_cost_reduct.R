base_prices <- (cea_summary_all %>% group_by(intervention) %>% summarise(price=min(price)))$price
names(base_prices)<-(cea_summary_all %>% group_by(intervention) %>% summarise(price=min(price)))$intervention

plot_list<-list()
for (k_name_categ in 1:length(old_new_names$old)){
  for (k_name in 1:length(old_new_names$old[[k_name_categ]])){
    cea_summary_all <- cea_summary_all %>% 
      mutate(plot_variable=ifelse(variable %in% old_new_names$old[[k_name_categ]][k_name],
                                  old_new_names$new[[k_name_categ]][k_name],plot_variable)) } }

for (k_plot in 1:3) {
  sel_vars <- old_new_names$old[[k_plot]]
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) & # intervention=="mAb" & 
             ((price==base_prices["maternal"] & intervention=="maternal")|(price==base_prices["mAb"] & 
                  intervention=="mAb")) & grepl("new",source)) %>%
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
           burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
           source=ifelse(grepl("projection",as.character(source)),
                         "projection (from [Shi 2017])",as.character(source)),
           vartype=gsub("_averted","",variable)) %>% group_by(vartype) %>% 
    mutate(plot_variable=unique(plot_variable[!grepl("averted",variable)])) %>% ungroup() %>% 
    mutate(vec=as.character( round((10^floor(log10(median/norm_median)))*round(
      median/norm_median/(10^floor(log10(median/norm_median))),3)) ) )
  if (any(grepl("e\\+",df_plot$vec))){
    df_plot$vec <- gsub("e\\+07",paste0(rep("0",7),collapse=""),df_plot$vec)
    df_plot$vec <- gsub("e\\+06",paste0(rep("0",6),collapse=""),df_plot$vec) }
  if (any(grepl("\\.",df_plot$vec)) ){ 
    subs_vec<-df_plot$vec[grepl("\\.",df_plot$vec)]
    df_plot$vec[grepl("\\.",df_plot$vec)] <- substr(gsub("\\.","",subs_vec),1,nchar(gsub("\\.","",subs_vec))-1) }
  # rounded values
  df_plot$orig_burden_round=NA; df_plot$orig_burden_round[!is.na(df_plot$vec)]=gsub("^\\.","",
           sapply(df_plot$vec[!is.na(df_plot$vec)], 
          function(x) paste0(substring(x,first=c(1,seq(nchar(x)-floor(nchar(x)/3)*3,
        nchar(x)-1,by=3)+1),last=c(seq(nchar(x)-floor(nchar(x)/3)*3,nchar(x)-1,by=3),nchar(x))),collapse=",")))
  df_plot <- df_plot %>% mutate(orig_burden_round=ifelse(as.numeric(vec)<1e4,
                                                         gsub("\\,","",orig_burden_round),orig_burden_round),
                                plot_variable=ifelse(grepl("cost",plot_variable),
                                                     paste0(plot_variable," (US$)"),plot_variable))
  # for burden plot, order factors
  if (any(df_plot$vartype %in% "total_DALY")) {
    df_plot$vartype=factor(df_plot$vartype,
                           levels=unique(df_plot$vartype)[c(which(!grepl("total_DALY",unique(df_plot$vartype))),
                                                            which(grepl("total_DALY",unique(df_plot$vartype))))])
    df_plot$plot_variable=factor(df_plot$plot_variable,
                                 levels=c("YLD non-hospitalised cases","YLD medically attended cases",
                                          "YLD hospitalised cases","total YLD","total YLL","total DALY")) 
  }
  # plot for 2 cntrs, 2 intervents
  dodge_val=1; round_val=2; caption_txt=paste0("Numbers above medians are pre-intervention",
                ifelse(any(grepl("YLL",sel_vars))," DALYs",
                ifelse(any(grepl("death",sel_vars))," case/death numbers"," costs (USD)") ),
                ifelse(any(grepl("YLL",sel_vars)),". YLD=years lived with disability. YLL=years of life lost",""))
  # labels for y-axis
  ylab_txt=paste0("% reduction in ",ifelse(any(grepl("YLL",sel_vars)),"DALYs",
                                           ifelse(any(grepl("death",sel_vars)),"cases/deaths","cost")))
  ### ### ### ### ### ### ### ### ### ### ### ###
  if (k_plot==1) {
    df_combined_burden_cost_reduct <- df_plot } else {
    df_combined_burden_cost_reduct <- bind_rows(df_combined_burden_cost_reduct,df_plot %>% 
                                                  mutate(across(where(is.numeric),round,3)))
    }

  ### ### ### ### ### ### ### ### ### ### ### ### 
  # PLOT with cntr on x-axis, MV/mAb as colors
  plot_list[[k_plot]] <- ggplot(df_plot %>% filter(grepl("averted",burden_interv)) %>% 
                                  mutate(cnt_int=paste0(country_iso,"\n(",intervention,")")) ) +
    geom_boxplot(aes(x=cnt_int,color=intervention,middle=norm_median*1e2,
                     ymin=norm_CI95_low*1e2,ymax=norm_CI95_high*1e2,lower=norm_CI50_low*1e2,upper=norm_CI50_high*1e2),
                 position=position_dodge(width=dodge_val),stat="identity") +
    scale_color_manual(values=c("red","blue")) +
    facet_wrap(~plot_variable) + geom_vline(xintercept=2.5,linetype="dashed",size=0.3) + # vartype
    theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
    geom_text(aes(x=as.numeric(interaction(intervention,country_iso))+1/2,
                  y=max(norm_CI95_high*1e2)*1.05,group=intervention,
                  label=ifelse(intervention!="MV",gsub("^,","",paste0(orig_burden_round,
                               ifelse(grepl("YLD non",plot_variable),"*",""))),"")),
              position=position_dodge(width=dodge_val),size=ifelse(exists("geom_text_font_size"),geom_text_font_size,5)) + 
    labs(color="",linetype="",
         caption=ifelse(any(grepl("costs",df_plot$plot_variable)),"*pre-intervention median value","")) + 
    scale_y_continuous(breaks=(0:15)*10) + # scale_x_discrete(expand=expansion(0.1,0)) + 
    theme(axis.text.x=element_text(angle=0,vjust=1/2,size=17),axis.text.y=element_text(size=17),
          strip.text=element_text(size=18),legend.position=ifelse(save_flag,"top",
                                  ifelse(k_plot<3,"none","bottom")),legend.text=element_text(size=18),
          axis.title.y=element_text(size=16.5,margin=margin(t=0,r=12,b=0,l=0)),
          legend.spacing.x=unit(0.7,'cm'),legend.key.size=unit(3,"line"),plot.caption=element_text(size=14))
  if (k_plot<3) { plot_list[[k_plot]] <- plot_list[[k_plot]] + theme(axis.text.x=element_blank()) }
  # save
  if (save_flag) {
    plot_list[[k_plot]]
    ggsave(paste0("output/cea_plots/",subfolder_name,ifelse(any(grepl("YLL",sel_vars)),"DALY",
                               ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),
           width=36,height=18,units="cm") } 
} ##### end of loop