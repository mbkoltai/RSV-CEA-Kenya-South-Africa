k_par=1
n_cntr_output=par_table$n_cntr_output[k_par]; n_interv=par_table$n_interv[k_par]
# intervention config table
sel_interv=sim_config_matrix[which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])[n_interv],]
if (cntrs_cea[n_cntr_output]=="ZAF"){sel_interv$country_iso=cntrs_cea[n_cntr_output]}
k_price=1
doseprice=c("mat_vacc"=ifelse(n_interv==1,pricelist$mat_vacc[k_price],pricelist$mat_vacc[1]),
            "mAb"=ifelse(n_interv==2,pricelist$mAb[k_price],pricelist$mAb[1]))
# calculation with data from mcmarcel (community-based)
sim_output=get_burden_flexible(sel_interv,NA,NA,exp_wane=exp_wane_val,doseprice)
