#' ---
#' title: RSV burden&cost calculations with incidence data from Kenya and South Africa 
#' author: Mihaly Koltai
#' date: 29/October 2020
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#'      fig_width: 15
#'      fig_height: 10
#' ---
# Libraries
# to render as html: rmarkdown::render("RSV_kenya_SA_calculation.R",output_dir='output/cea_plots/')
library(tidyverse); library(reshape2); library(matrixStats); library(rstudioapi); # library(fitdistrplus)
# from: http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/GammaParmsFromQuantiles.html
# sessionInfo()
# path
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# base MCMARCEL functions
source('functions/RSV_load_all.R')
# custom made functions
source('functions/GammaParmsFromQuantiles.R'); source('functions/get_burden_flexible.R');source('functions/facet_wrap_fcns.R')
# clear vars: rm(list=ls())
# plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
      plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
      legend.title=element_text(size=14),legend.text=element_text(size=12),
      axis.title=element_text(size=14), text=element_text(family="Calibri"))
# run tag
run_tag  <- 'RSV_gavi72_basecase'  # 72 Gavi countries (basecase) 
# number of stochastic samples in the probabilistic sensitivity analysis (PSA).
num_sim <- 5000
##### SELECT COUNTRY
cntrs_cea=c('KEN','ZAF')
n_cntr_output=2; cntr_sel=cntrs_cea[n_cntr_output]
# load config params
source('functions/load_config_pars.R')
# if running for many cntrs and we want to parallelise
# start_parallel_workers()
### Demographic data -------------------------
#' ## Load demographic data
#### Load life tables and incidence used by MCMARCEL (this is not needed for CEA below)
source('functions/load_incidence_lifetables.R')
# plot lifetable
fcn_plot_lifetable('')
# plot incidence for cntrs (Set argument to 'save' to save)
fcn_plot_lmic_incidence('')
### OWN BURDEN DATA --------------------------------------------------
#' ## Incidence data from partners
source('functions/load_own_data.R')
### Kenya incidence data -------------------------
#' ### Kenya
kenya_data_path='../path_rsv_data/SARI_Rates_2010_2018/SARI_Rates_2010_2018_tidydata_cleaned.csv'
n_iter=5e3; age_maxval=60; CI95_const=1.96; CI_intervals=c(0.025,0.975); randsampl_distrib_type='gamma'
kenya_hosp_rate_stdev=0.05656925 # value from: unique(apply(config$hosp_prob,1,sd))
# a list of incidence and hospit rate matrices: list(kemri_incid_rate_matrix,kemri_hosp_rate_matrix)
kenya_incid_hosp_rate=fcn_load_kenya(kenya_data_path,n_iter,age_maxval,
                                     CI95_const,CI_intervals,kenya_hosp_rate_stdev,randsampl_distrib_type)
#' ### South Africa
safr_data_path='../path_rsv_data/s_afr_incidence_data.csv'; s_afr_province_means_filepath='input/safr_hosp_rate_province_means.csv'
s_afr_incid_hosp_rate=fcn_load_s_afr(safr_data_path,n_iter,age_maxval,CI95_const,CI_intervals,
                        kenya_hosp_rate_stdev,randsampl_distrib_type)
####
# put all matrices into 1 list; 1 matrix: burden_list_own_data[[2]][[2]]
burden_list_own_data=list(list(kenya_incid_hosp_rate[[1]],kenya_incid_hosp_rate[[2]]),
                          list(s_afr_incid_hosp_rate[[1]],s_afr_incid_hosp_rate[[2]]))
### BURDEN & CEA CALCULATION --------------------------------------------------
#' ## BURDEN CALCULATION with OWN DATA
# South Africa was not in the original analysis, so we cannot compare to default mcmarcel output
# cntrs_cea=c('KEN','ZAF'); n_cntr_output=2
# read in config parameters for burden calculation
burden_cntr_ind=which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])
# MV=sim_config_matrix[burden_cntr_ind[1],]; mAb=sim_config_matrix[burden_cntr_ind[2],]
# SELECT INTERVENTION: 1=MatVacc, 2=monocl Abs
n_interv=1; sel_interv=sim_config_matrix[burden_cntr_ind[n_interv],]; # config <- get_rsv_ce_config(sel_interv)
# if (n_cntr_output==2){sel_interv$country_iso=cntrs_cea[n_cntr_output]}
# with original MCMARCEL: 
# sim_output=get_burden(sel_interv)
####
# own data: need to provide incidence matrix [60*5e3] and hospitalisation matrix [60*5e3] as arguments
sim_output_flexible=get_burden_flexible(sel_interv,burden_list_own_data[[n_cntr_output]][[1]],
                                        burden_list_own_data[[n_cntr_output]][[2]])
# for original
sim_output=get_burden_flexible(sel_interv,NA,NA)
### Plot results -------------------------
#' ## Process & plot results
# outputs to display 
icercolname="net_cost/DALY_averted"
cols_burden_sel=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y',
                  'cost_rsv_hosp_0to1y','total_medical_cost_0to1y','incremental_cost_0to1y','total_medical_cost_averted',
                  "hosp_cases_averted", "rsv_deaths_averted", "total_DALY_0to1y_averted",icercolname) 
# "total_YLD_1to5y_averted","total_YLL_1to5y_averted",
# source('functions/set_xlims_cea.R')
burden_mcmarcel_owndata=fcn_process_burden_output(sim_output_flexible,sim_output,sel_interv$country_iso,cols_burden_sel)
# we'll show mean values by vertical lines
color_vals=c('red','green'); mean_intercepts=fcn_mean_intercepts(burden_mcmarcel_owndata,color_vals)
# set xlimits by quantiles
source('functions/set_xlims_cea.R')
quantl_vals=c(0.01,0.99); xlimvals=list(); xlimvals_cols=list()
# define xlims (leave xlimvals_cols, xlimvals)
quantls_min_max=fcn_calc_quantls_minmax(burden_mcmarcel_owndata,quantl_vals)
# leave last argument empty for automatic x limits
xmin_adj_vars=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y','cost_rsv_hosp_0to1y','total_medical_cost_0to1y',
                'total_medical_cost_averted','hosp_cases_averted','rsv_deaths_averted','total_DALY_0to1y_averted','net_cost/DALY_averted')
xmax_adj_vars=c('rsv_cases','incremental_cost_0to1y','net_cost/DALY_averted')
xmin_adj_values=c(0,0,-1e2,-1e3,-1e6,-1e6,-1e6,-1e3,-10,-2e3,-1e3); icer_maxval=3e3; scale_max_val=2
scales_list=fcn_set_xlims(sel_interv$country_iso,burden_mcmarcel_owndata_comp,quantls_min_max,xmin_adj_vars,xmin_adj_values,
                          xmax_adj_vars,icer_maxval,scale_max_val,'man_adj_flag')
### Plot CEA histograms ------------------
if (n_interv==1) {interv_tag='_Mat_Vacc' } else {interv_tag='_monocl_Ab'}
ggplot(burden_mcmarcel_owndata,aes(x=value,group=source)) + geom_freqpoly(aes(color=source),size=1.2) + 
  geom_vline(data=mean_intercepts,aes(xintercept=int,linetype=source),size=0.8) + 
  # facet_wrap(~variable,scales='free',labeller=label_wrap_gen(width=10)) + # ,ncol=3
  facet_wrap_custom(~variable,scales='free',scale_overrides=scales_list,labeller=label_wrap_gen(width=10)) +
  theme_bw() + standard_theme + scale_linetype_manual(values=c('solid','dotdash')) +
  geom_rect(data=subset(burden_mcmarcel_owndata, variable %in% icercolname),fill=NA,colour="blue",
            size=2,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
  ggtitle(paste(sel_interv$country_iso,'RSV burden & intervention estimates:',gsub('_','',interv_tag))) +
  labs(color='data source',linetype='mean') + guides(xintercept=FALSE,linetype=guide_legend(ncol=2))
# save plot
cea_plot_filename=paste("output/cea_plots/",sel_interv$country_iso,"_mcmarcel_burden_estimates_n",n_iter,
                        interv_tag,'_',randsampl_distrib_type,'samples',".png",sep="")
# if (!dir.exists('output/cea_plots')) {dir.create('output/cea_plots')}
# SAVE
# ggsave(cea_plot_filename,width=30,height=18,units="cm")
### Net cost/DALY averted (mean, median, stdev) -------------------------
burden_mcmarcel_owndata_comp[burden_mcmarcel_owndata_comp$variable %in% 'net_cost/DALY_averted',] %>% 
  group_by(source) %>% summarise(meanvals=mean(value),medianval=median(value),stdevvals=sd(value))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Extrapolation to other Sub-Sah Afr countries -------------------------

# fraction of cases by age groups
kenya_agegroup_fractions=rowMeans(kenya_incid_hosp_rate[[1]]/colSums(kenya_incid_hosp_rate[[1]]))
if (!sum(kenya_agegroup_fractions)==1){kenya_agegroup_fractions=kenya_agegroup_fractions/sum(kenya_agegroup_fractions)}
# S Afr
s_afr_agegroup_fractions=rowMeans(s_afr_incid_hosp_rate[[1]]/colSums(s_afr_incid_hosp_rate[[1]]))
if (!sum(s_afr_agegroup_fractions)==1){s_afr_agegroup_fractions=s_afr_agegroup_fractions/sum(s_afr_agegroup_fractions)}
# SELECT AGE DISTRIBUTION
age_distrib_used=c('KEN','ZAF')[1]
# SELECT FOLDER TO SAVE/WRITE FILES TO
if (age_distrib_used %in% 'ZAF') {
  agegroup_fractions=s_afr_agegroup_fractions; folder_projection='output/sub_sah_afr_proj/s_afr_age_distrib/' } else {
    agegroup_fractions=kenya_agegroup_fractions; folder_projection='output/sub_sah_afr_proj/kenya_age_distrib/' }
# mean incidence
RSV_burden_Shi_2017=read_csv('input/RSV_burden_Shi_2017.csv')
# sub-Saharan Afr countries
sub_sah_afr_cntrs=c("Angola","Benin","Botswana","Burkina Faso","Burundi","Cameroon","Cape Verde","Central African Republic",
"Chad","Comoros","Congo, Rep.","Congo, Dem. Rep.","C_te d'Ivoire","Djibouti","Equatorial Guinea","Eritrea","Ethiopia","Gabon",
"Gambia, The","Ghana","Guinea","Guinea-Bissau","Kenya","Lesotho","Liberia","Madagascar","Malawi","Mali","Mauritania","Mauritius",
"Mozambique","Namibia","Niger","Nigeria","Rwanda","S_o Tom_ and Principe","Senegal","Seychelles","Sierra Leone","Somalia",
"South Africa","Sudan","South Sudan","Swaziland","Tanzania","Togo","Uganda","Zambia","Zimbabwe")
# convert to ISO3
sub_sah_afr_cntrs_iso3=RSV_burden_Shi_2017$country_iso[RSV_burden_Shi_2017$location_name %in% sub_sah_afr_cntrs]
sub_sah_afr_cntrs_iso3=sub_sah_afr_cntrs_iso3[sub_sah_afr_cntrs_iso3 %in% sim_config_matrix$country_iso]
# which intervention?
n_interv=1
### Calculate CEA  -------------------------
set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))
start_parallel_workers(); all_packages=append(all_packages,'readr')
### 
# this is time consuming, ~30mins depending on # CPUs
# ptm <- proc.time()
foreach(k=1:length(sub_sah_afr_cntrs_iso3),.combine='rbind',.packages=all_packages,.verbose=FALSE) %dopar% {
  # LOOP through SubSah cntrs
  cntr_name=sub_sah_afr_cntrs_iso3[k]; df_country <- get_incidence(cntr_name,output_dir)
  # Set config values
  config=list(); config$rsv_rate <- df_country[[1]]; config$hosp_prob <- df_country[[2]]; config$hosp_CFR <- df_country[[3]]
  rsv_rate_fractions=sapply(1:ncol(config$rsv_rate),function(x) {config$rsv_rate[,x]/sum(config$rsv_rate[,x])})
  subsah_incid_redist=matrix(unlist(lapply(colSums(config$rsv_rate),function(x) {agegroup_fractions*x})),ncol=n_iter)
  # can generate values x>1 or x<0: 
# sapply(1:n_iter, function(x){ (rsv_rate_fractions[,x]/rowMeans(rsv_rate_fractions))*agegroup_fractions*sum(config$rsv_rate[,x]) })
# generates low noise: matrix(unlist(lapply(colSums(config$rsv_rate),function(x) {agegroup_fractions*x})),ncol=n_iter)
  sel_interv=sim_config_matrix[which(sim_config_matrix$country_iso %in% cntr_name)[n_interv],] 
  if (!is.na(sel_interv$rng_seed)){ write_csv(get_burden_flexible(sel_interv,subsah_incid_redist,config$hosp_prob),
            paste0(folder_projection,"indiv_countries/",cntr_name,"_own_data.csv"))
  print(paste("done with",cntr_name,"own data",sep=" ")) # mcmarcel_outputs=get_burden_flexible(sel_interv,NA,NA)}
  write_csv(get_burden_flexible(sel_interv,NA,NA),paste0(folder_projection,"indiv_countries/",cntr_name,"_mcmarcel.csv")) }
}
# proc.time() - ptm

#### Process extrapolation results  -------------------------
# select datafiles
cea_file_list=list.files(paste0(folder_projection,'indiv_countries')); cea_file_list=cea_file_list[grepl('.csv',cea_file_list)]
cea_file_id=sapply(strsplit(cea_file_list,"_"),'[[',1); burden_comparison_all_subsah_afr=list()
# x axis limits
icercolname="net_cost/DALY_averted"
# cols_burden_sel=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y',
#                   'cost_rsv_hosp_0to1y','total_medical_cost_0to1y','incremental_cost_0to1y','total_medical_cost_averted',
#                   "hosp_cases_averted", "rsv_deaths_averted", "total_DALY_0to1y_averted",icercolname) 
# xmin_adj_vars=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y','cost_rsv_hosp_0to1y','total_medical_cost_0to1y',
#   'total_medical_cost_averted','hosp_cases_averted','rsv_deaths_averted','total_DALY_0to1y_averted','net_cost/DALY_averted')
# xmax_adj_vars=c('rsv_cases','incremental_cost_0to1y','net_cost/DALY_averted')
# PLOT CONFIG
# quantile to draw x-axis limits
quantl_vals=c(0.01,0.99); color_vals=c('red','green')
# quantiles to store in csv file
quantl_vals_90=c(0.05,0.9)
# minvals for x-axis, maxval for ICER, scaling of maxvals of x-axis
xmin_adj_values=c(0,0,-1e2,-1e3,-1e6,-1e6,-1e6,-1e3,-10,-2e3,-1e3); icer_maxval=3e3; scale_max_val=2
save_flag='save'; list_all_subsah_afr=list()
for (k in 1:length(sub_sah_afr_cntrs_iso3)){
  sel_cntr=sub_sah_afr_cntrs_iso3[k]; CEA_files_indcntr=cea_file_list[cea_file_id %in% sel_cntr]
  sel_cntr_fullname=RSV_burden_Shi_2017$location_name[RSV_burden_Shi_2017$country_iso %in% sel_cntr]
  # own data will be always the 2nd of the two files for each country
  burden_mcmarcel_owndata=fcn_process_burden_output(
        read_csv(paste0(folder_projection,'indiv_countries/',CEA_files_indcntr[[2]])),
        read_csv(paste0(folder_projection,'indiv_countries/',CEA_files_indcntr[[1]])),sel_cntr,cols_burden_sel)
  # PLOT distribution (density plot)
  data_to_plot=burden_mcmarcel_owndata; data_to_plot$source=gsub(' ','',gsub("\\(|\\)","",gsub(sel_cntr,"",data_to_plot$source)))
  ## we'll show mean values by vertical lines
  mean_intercepts=fcn_mean_intercepts(data_to_plot,color_vals)
  # set xlimits by quantiles
  quantls_min_max=fcn_calc_quantls_minmax(burden_mcmarcel_owndata,quantl_vals)
  # leave last argument empty for automatic x limits
  scales_list=fcn_set_xlims(sel_cntr,burden_mcmarcel_owndata_comp,quantls_min_max,xmin_adj_vars,xmin_adj_values,
               xmax_adj_vars,icer_maxval,scale_max_val,'man_adj_flag')
  # leave 'scales_list' empty to have default x axes
  # PLOT distributions of variables for selected country
  # fcn_plot_cea_distribs(data_to_plot,mean_intercepts,scales_list, standard_theme,n_interv, # scales_list
  #                       sel_cntr_fullname,paste0(folder_projection,'plots/'),save_flag) # save_flag
  ## WRITE TABLE WITH mean/median/stdev + QUINTILES
  burden_comparison_median=burden_mcmarcel_owndata %>% group_by(variable,source) %>%
    summarise(meanvals=mean(value),medianval=median(value),stdevvals=sd(value))
  burden_comparison_quants=burden_mcmarcel_owndata %>% group_by(variable,source) %>% 
    summarise(value=quantile(value,probs=quantl_vals_90),q=quantl_vals_90) %>% group_by(variable,source,q) %>% 
    summarise(min_q=min(value),max_q=max(value)) %>% group_by(variable,source) %>% summarise(min_q=min(min_q),max_q=max(max_q))
  burden_comparison_median=merge(burden_comparison_median,burden_comparison_quants,by=c("variable","source"))
  burden_comparison_median[,"country"]=sapply(strsplit(as.character(burden_comparison_median$source)," "),"[[",1)
 burden_comparison_median$source=gsub("\\(|\\)","",sapply(strsplit(as.character(burden_comparison_median$source)," "),"[[",2))
  list_all_subsah_afr[[k]]=burden_comparison_median
}

burden_comparison_all_subsah_afr=data.frame(do.call(rbind,list_all_subsah_afr))
# source of the age distribution
burden_comparison_all_subsah_afr[,'age_distrib']=age_distrib_used
#''; burden_comparison_all_subsah_afr[burden_comparison_all_subsah_afr$source %in% "own",'age_distrib']=
colclass=sapply(1:ncol(burden_comparison_all_subsah_afr), function(x) {is.numeric(burden_comparison_all_subsah_afr[,x])})
# no point in storing decimals beyond the first (numbers are 10^2 --> 10^7), so round to 1st decimal
burden_comparison_all_subsah_afr[,colclass]=round(burden_comparison_all_subsah_afr[,colclass],1)
write_csv(burden_comparison_all_subsah_afr,paste0(folder_projection,'burden_comparison_all_subsah_afr.csv'))

### combine two tables with diff age distributions --------------------
burden_comparison_all_subsah_afr=rbind(read_csv('output/sub_sah_afr_proj/s_afr_age_distrib/burden_comparison_all_subsah_afr.csv'),
      read_csv('output/sub_sah_afr_proj/kenya_age_distrib/burden_comparison_all_subsah_afr.csv'))
write_csv(burden_comparison_all_subsah_afr,'output/sub_sah_afr_proj/burden_comparison_all_subsah_afr.csv')

# read in results using both age distribs
burden_comparison_all_subsah_afr=read_csv('output/sub_sah_afr_proj/burden_comparison_all_subsah_afr.csv')
# the results from two cntrs with mcmarcel age distrib have to be the same, removing one of them
mcmarcel_agedistrib_dupl_rows=burden_comparison_all_subsah_afr$source %in% 'mcmarcel' & 
  burden_comparison_all_subsah_afr$age_distrib %in% "KEN"
burden_comparison_all_subsah_afr=burden_comparison_all_subsah_afr[!mcmarcel_agedistrib_dupl_rows,]
burden_comparison_all_subsah_afr$age_distrib[burden_comparison_all_subsah_afr$source %in% 'mcmarcel']='Shi2017'
# set order by factoring
burden_comparison_all_subsah_afr$variable=factor(burden_comparison_all_subsah_afr$variable,levels = cols_burden_sel)
burden_comparison_all_subsah_afr$age_distrib=factor(burden_comparison_all_subsah_afr$age_distrib,levels=c("Shi2017","ZAF","KEN"))
# filter out cntrs below X million
# using "wpp2017" and "countrycode"
sub_sah_afr_cntrs_popul2020=popproj[popproj$country_code %in% countrycode(sub_sah_afr_cntrs_iso3,origin='iso3c',destination='iso3n'),"2020"]
cntr_list=popproj[popproj$country_code %in% countrycode(sub_sah_afr_cntrs_iso3,origin='iso3c',destination='iso3n'),"name"]
# data.frame(cntrs=cntr_list,popul=sub_sah_afr_cntrs_popul2020)
# cntrs with popul larger than 5 million
cntr_list_popul_lt=countrycode(cntr_list,origin='country.name',destination='iso3c')[sub_sah_afr_cntrs_popul2020>5e3]

# plot all countries together
burden_vars=cols_burden_sel[1:6]; interv_vars=cols_burden_sel[7:length(cols_burden_sel)]
# select vars to plot
varselect=burden_comparison_all_subsah_afr$country %in% cntr_list_popul_lt & 
  burden_comparison_all_subsah_afr$variable %in% burden_vars
#
ggplot(burden_comparison_all_subsah_afr[varselect,],aes(x=country,y=medianval,ymin=min_q,ymax=max_q)) + # medianval-stdevvals
 geom_pointrange(aes(color=age_distrib,shape=age_distrib,linetype=age_distrib),
                 fatten=1.8,size=0.8,position=position_dodge(width=0.6)) +
  facet_wrap(~variable,scales='free',ncol=2) + scale_color_manual(values=c('blue','red2','red3')) +
  scale_linetype_manual(values=c('twodash','solid','solid')) + theme_bw() + standard_theme + theme(legend.key.size = grid::unit(3, "lines")) + 
  labs(color='data source',shape='data source',linetype='data source') + # guides(linetype=FALSE) +
  xlab('') + ylab('') + ggtitle('Effect of age distribution on projected CEA')
# guide_legend(override.aes=list(shape=NA)),color=FALSE
####
# save
if (any(burden_comparison_all_subsah_afr$variable[varselect] %in% interv_vars)) {
    ggsave('output/sub_sah_afr_proj/plots/intervention_variable_faceted.png',width=50,height=30,units="cm")} else {
    ggsave('output/sub_sah_afr_proj/plots/burden_variable_faceted.png',width=50,height=30,units="cm")
}
