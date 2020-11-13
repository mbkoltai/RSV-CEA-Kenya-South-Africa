fcn_calc_quantls_minmax=function(burden_mcmarcel_owndata_comp,quantl_vals){
  quantls_min_max=burden_mcmarcel_owndata_comp %>% group_by(variable,source) %>% 
    summarise(value=quantile(value,quantl_vals),q=quantl_vals) %>% group_by(variable,q) %>% 
    summarise(min_q=min(value),max_q=max(value)) %>% group_by(variable) %>% summarise(min_q=min(min_q),max_q=max(max_q))
}

# manually modify x axis limits
fcn_set_xlims<-function(sel_cntr,burden_mcmarcel_owndata_comp,quantls_min_max,xmin_adj_vars,xmin_adj_values,
                        xmax_adj_vars,icer_maxval,scale_max_val,man_adj_flag){ # xlimvals_cols,xlimvals
  if (nchar(man_adj_flag)>0){
  xlimvals=list(); xlimvals_cols=list()
  xlimvals_cols[[sel_cntr]][[1]]=xmin_adj_vars; xlimvals[[sel_cntr]][[1]]=xmin_adj_values
  minvals=quantls_min_max$min_q[quantls_min_max$variable %in% xmin_adj_vars]
  overshoot_xmin_vals = xmin_adj_values > minvals # overshoot_xmin_vals[length(overshoot_xmin_vals)]=F
  xlimvals[[sel_cntr]][[1]][overshoot_xmin_vals]=minvals[overshoot_xmin_vals] # -1.5*abs()
  # if the minimum stretches the axis too much (too low value)
  kxmin=10; neg_overshoot=(xmin_adj_values<(-kxmin*(abs(minvals-1)))); neg_overshoot[length(neg_overshoot)]=F
  # print(rbind(xmin_adj_values,minvals,-kxmin*(abs(minvals-1)),neg_overshoot))
  xlimvals[[sel_cntr]][[1]][neg_overshoot]=-kxmin*abs(minvals[neg_overshoot])
  # set xlim maximum
  xlimvals_cols[[sel_cntr]][[2]]=xmax_adj_vars
  xlimvals[[sel_cntr]][[2]]=array(t(scale_max_val*quantls_min_max[quantls_min_max$variable %in% xmax_adj_vars,3]))
  xlimvals[[sel_cntr]][[2]][xlimvals_cols[[sel_cntr]][[2]] %in% "net_cost/DALY_averted" ]=icer_maxval}
  col_inds=c(2,3) # minvalues (positive), maxvalues
  if (any(sapply(xlimvals_cols,nchar)>0)){ # min values positive, minvalues negative, maxvalues
    for (k in 1:length(xlimvals_cols[[sel_cntr]])){ 
      if (any(!sapply(xlimvals_cols[[sel_cntr]][[k]],is.na))){ # any(sapply(xlimvals_cols[[sel_cntr]][k],nchar)>0)
        quantls_min_max[quantls_min_max$variable %in% unlist(xlimvals_cols[[sel_cntr]][[k]]),col_inds[k]]=
          unlist(xlimvals[[sel_cntr]][[k]]) }   }   }
  scales_list=list(); for (k in 1:length(cols_burden_sel)) {
    scales_list[[k]]=scale_override(k,scale_x_continuous(limits=c(quantls_min_max[k,]$min_q,quantls_min_max[k,]$max_q))) }
  scales_list
}

### process burden calc output
fcn_process_burden_output<-function(sim_output_flexible,sim_output,sel_cntr,cols_burden_sel){
  icercolname='net_cost/DALY_averted'; icercols=c('incremental_cost_0to1y','total_DALY_0to1y_averted')
  sim_output[,icercolname]=sim_output[,icercols[1]]/sim_output[,icercols[2]]
  sim_output_flexible[,icercolname]=sim_output_flexible[,icercols[1]]/sim_output_flexible[,icercols[2]]
  # outputs to display
  # cols_burden_sel=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y',
  #                   'cost_rsv_hosp_0to1y','total_medical_cost_0to1y','incremental_cost_0to1y','total_medical_cost_averted',
  #                   "hosp_cases_averted", "rsv_deaths_averted", "total_DALY_0to1y_averted",icercolname) 
  # "total_YLD_1to5y_averted","total_YLL_1to5y_averted",
  # histograms: mcmarcel vs kemri
  burden_mcmarcel_owndata_comp=melt(
    rbind( cbind(sim_output[,cols_burden_sel],
                data.frame(source=paste0(sel_cntr,' (mcmarcel)'),iter=1:nrow(sim_output))),
           cbind(sim_output_flexible[,cols_burden_sel],
                data.frame(source=paste0(sel_cntr,' (own)'),iter=1:nrow(sim_output_flexible) ) ) ),
    id.vars=c('iter','source'))
  burden_mcmarcel_owndata_comp
}

### calculate mean values as intercepts --------------
fcn_mean_intercepts=function(burden_mcmarcel_owndata_comp,color_vals){
mean_intercepts=burden_mcmarcel_owndata_comp %>% group_by(source,variable) %>% summarize(int = mean(value))
mean_intercepts[,'colorval']=color_vals[1]; mean_intercepts$colorval[grepl('own',mean_intercepts$source)]=color_vals[2]
mean_intercepts
}

### plot distributions ---------
fcn_plot_cea_distribs <- function(burden_mcmarcel_owndata_comp,mean_intercepts,scales_list,
                                standard_theme,n_interv,sel_cntr_fullname,folderpath,save_flag){
  if (n_interv==1) {interv_tag='_Mat_Vacc' } else {interv_tag='_monocl_Ab'}
p<-ggplot(burden_mcmarcel_owndata_comp,aes(x=value,group=source,color=source)) + geom_freqpoly(size=1.2) + 
  geom_vline(data=mean_intercepts,aes(xintercept=int,color=source),size=0.8,linetype='dashed',show.legend = FALSE) + 
  theme_bw() + standard_theme + # scale_linetype_manual(values=c('solid','dotdash')) +
  geom_rect(data=subset(burden_mcmarcel_owndata_comp, variable %in% "net_cost/DALY_averted"),fill=NA,colour="blue",
            size=2,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
  ggtitle(paste('RSV burden & intervention estimates',gsub('_','',interv_tag),": ",sel_cntr_fullname)) +
  labs(color='data source') # + guides(xintercept=FALSE,linetype=guide_legend(ncol=2)) # ,linetype='data source'
if (length(scales_list)>1){
p <- p + facet_wrap_custom(~variable,scales='free',scale_overrides=scales_list,labeller=label_wrap_gen(width=10))} else {
p <- p + facet_wrap(~variable,scales = "free") }
# SAVE
if (nchar(save_flag)>1){ print('saving plot')
  cea_plot_filename=paste(folderpath,gsub(" ","_",sel_cntr_fullname),"_burden_estimates_n",n_iter,interv_tag,".png",sep="")
  ggsave(cea_plot_filename,width=30,height=18,units="cm") }
p
}

###
# calculate burden median, mean, stdev for subsah afr projection
fcn_burden_median_mean=function(burden_mcmarcel_owndata){
burden_comparison_median=burden_mcmarcel_owndata %>% group_by(variable,source) %>%
  summarise(meanvals=mean(value),medianval=median(value),stdevvals=sd(value))
burden_comparison_quants=burden_mcmarcel_owndata %>% group_by(variable,source) %>% 
  summarise(value=quantile(value,probs=quantl_vals_90),q=quantl_vals_90) %>% group_by(variable,source,q) %>% 
  summarise(min_q=min(value),max_q=max(value)) %>% group_by(variable,source) %>% summarise(min_q=min(min_q),max_q=max(max_q))
burden_comparison_median=merge(burden_comparison_median,burden_comparison_quants,by=c("variable","source"))
burden_comparison_median[,"country"]=sapply(strsplit(as.character(burden_comparison_median$source)," "),"[[",1)
burden_comparison_median$source=gsub("\\(|\\)","",sapply(strsplit(as.character(burden_comparison_median$source)," "),"[[",2))
burden_comparison_median
# list_all_subsah_afr[[k]]=burden_comparison_median
}

fcn_process_median_burden <- function(list_all_subsah_afr,interv_tag){
burden_comparison_all_subsah_afr=data.frame(do.call(rbind,list_all_subsah_afr))
# source of the age distribution
burden_comparison_all_subsah_afr[,'age_distrib']=age_distrib_used; burden_comparison_all_subsah_afr[,"intervention"]=interv_tag
colclass=sapply(1:ncol(burden_comparison_all_subsah_afr), function(x) {is.numeric(burden_comparison_all_subsah_afr[,x])})
# no point in storing decimals beyond the first (numbers are 10^2 --> 10^7), so round to 1st decimal
burden_comparison_all_subsah_afr[,colclass]=round(burden_comparison_all_subsah_afr[,colclass],1)
burden_comparison_all_subsah_afr
}
