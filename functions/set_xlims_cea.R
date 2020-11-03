# minvalues 0
fcn_set_xlims<-function(f_country_iso,burden_mcmarcel_owndata_comp,quantl_vals,xlimvals_cols,xlimvals){
  
  # quantl_vals=c(0.01,0.99) 
  quantls_min_max=burden_mcmarcel_owndata_comp %>% group_by(variable,source) %>% 
    summarise(value=quantile(value,quantl_vals),q=quantl_vals) %>% group_by(variable,q) %>% 
    summarise(min_q=min(value),max_q=max(value)) %>% group_by(variable) %>% summarise(min_q=min(min_q),max_q=max(max_q))
  
  if (!(is.na(xlimvals_cols) | xlimvals_cols=='')){
quantls_min_max[quantls_min_max$variable %in% unlist(xlimvals_cols[[f_country_iso]][1]),2]=unlist(xlimvals[[f_country_iso]][1])
# minvalues negative
quantls_min_max[quantls_min_max$variable %in% unlist(xlimvals_cols[[f_country_iso]][2]),2]=unlist(xlimvals[[f_country_iso]][2])
# maxvalues
quantls_min_max[quantls_min_max$variable %in% unlist(xlimvals_cols[[f_country_iso]][3]),3]=unlist(xlimvals[[f_country_iso]][3])
}
scales_list=list(); for (k in 1:length(cols_burden_sel)) {
  scales_list[[k]]=scale_override(k,scale_x_continuous(limits=c(quantls_min_max[k,]$min_q,quantls_min_max[k,]$max_q))) }
scales_list
}

### process burden calc output
fcn_process_burden_output<-function(sim_output_flexible,sim_output,sel_interv,cols_burden_sel){
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
                data.frame(source=paste0(sel_interv$country_iso,' (mcmarcel)'),iter=1:nrow(sim_output))),
           cbind(sim_output_flexible[,cols_burden_sel],
                data.frame(source=paste0(sel_interv$country_iso,' (own)'),iter=1:nrow(sim_output_flexible) ) ) ),
    id.vars=c('iter','source'))
  burden_mcmarcel_owndata_comp
}

### calculate mean values as intercepts
fcn_mean_intercepts=function(burden_mcmarcel_owndata_comp,color_vals){
mean_intercepts=burden_mcmarcel_owndata_comp %>% group_by(source,variable) %>% summarize(int = mean(value))
mean_intercepts[,'colorval']=color_vals[1]; mean_intercepts$colorval[grepl('own',mean_intercepts$source)]=color_vals[2]
mean_intercepts
}