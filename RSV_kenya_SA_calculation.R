#' ---
#' title: RSV burden&cost calculations with incidence data from Kenya and South Africa 
#' author: Mihaly Koltai
#' date: October 2020
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
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path)
# path
setwd(currentdir_path)
# clear vars: rm(list=ls())
# plotting theme I reuse
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
    plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
    axis.title=element_text(size=14), text=element_text(family="Calibri"))
# run tag
run_tag  <- 'RSV_gavi72_basecase'  # 72 Gavi countries (basecase) 
# number of stochastic samples in the probabilistic sensitivity analysis (PSA)
num_sim <- 1000
# random number generater seed
rng_seed <- gsub('-','',Sys.Date()) # 20190118
# option to create geographical and country-specific plots # note: this might require substantial processing time
boolean_country_plots <- FALSE; boolean_map_plots     <- FALSE # TRUE
source('functions/RSV_load_all.R')
# output directory postfix
output_dir_postfix <- paste0(run_tag,'_n',num_sim)
# add timestap to output directory name
output_dir <- paste0('output/',format(Sys.time(),'%m%d%H%M%S_'),output_dir_postfix)
# set seed
set.seed(rng_seed)
# config filename
config_filename <- paste0('./config/',run_tag,'.csv'); time_stamp_main <- Sys.time()
# always clear temporary results
cli_print('Clear all temporary output'); unlink(file.path(get_temp_output_folder(output_dir)),recursive = T)
# start parallel workers: if running for many cntrs and we want to parallelise
# start_parallel_workers()
# 144x12 table
sim_config_matrix <- read.table(config_filename,sep=',', dec='.',stringsAsFactors = F,header = T)
# set output file name prefix
sim_output_filename  <- file.path(output_dir,run_tag)
# add simulation details
sim_config_matrix$num_sim <- num_sim; sim_config_matrix$scenario_id <- 1:nrow(sim_config_matrix)
sim_config_matrix$rng_seed <- rng_seed; sim_config_matrix$outputFileDir <- get_output_folder(output_dir)
######################################################
#' ## Load demographic data
cntrs_cea=c('KEN','ZAF'); 
##### SELECT COUNTRY
n_cntr_output=2; cntr_sel=cntrs_cea[n_cntr_output]
######
# S Afr is not in the original study so needs to be appended
if (!cntr_sel %in% sim_config_matrix$country_iso) {
  df_append=sim_config_matrix[(nrow(sim_config_matrix)-1):nrow(sim_config_matrix),]
  df_append$country_iso=cntr_sel; sim_config_matrix=rbind(sim_config_matrix,df_append) }
# create UN country data (calls wpp package)
create_UN_country_database(output_dir)
# pre-process WPP2017 data
load_wpp2017_databases(output_dir)
# loads sex ratio and mortality by age groups; dataframes: sexratio, mxM, mxF
n_cntr=length(unique(sim_config_matrix$country_iso))
country_period_opt=data.frame(country_iso=rep(unique(sim_config_matrix$country_iso),2),
                            year=unlist(lapply(c(2020,2015),function(x){rep(x,n_cntr)})))
country_period_opt <- cbind(as.character(country_period_opt$country_iso), t(sapply(country_period_opt$year,get_year_category)))
# generate life table
# select country
i_life=which(country_period_opt[,1] %in% cntr_sel)[1]; f_country_iso=country_period_opt[i_life,1]
f_year=country_period_opt[i_life,3]; f_outputFileDir=output_dir
# with or without discounting
for (f_disc_rate_qaly in c(0,0.03)) {
# this saves into an Rdata file
generate_life_table(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
# 'life_table_KEN_2020_2025_disc0' ,'p03','.RData' # , 'life_table_KEN_2020_2025_disc0.RData'
rdata_prefix=paste0('life_table_',f_country_iso,'_2020_2025_disc0')
if (f_disc_rate_qaly>0) {rdata_filename=paste(rdata_prefix,'p03','.RData',sep='')} else {rdata_filename=paste(rdata_prefix,'.RData',sep='')}
# load life_table
load(paste(f_outputFileDir,"temp",rdata_filename,sep='/'))
# life_table contains: "age","lx": # left alive, "life_expectancy", "nMx": mortality, "lx_rate": # left alive/(init popul)
# "life_expectancy_disc": disc life yrs
if (f_disc_rate_qaly==0){
life_table_tidy=data.frame(melt(life_table,id.vars='age'),disc_rate=f_disc_rate_qaly)} else {
life_table_tidy=rbind(life_table_tidy,data.frame(melt(life_table,id.vars='age'),disc_rate=f_disc_rate_qaly)) }
}
####
# plot life tables
lifetable_varlist=list(as.character(unique(life_table_tidy$variable)),
                       c('alive/100.000 births','life expectancy','mortality','alive/birth','discounted life expectancy'))
life_table_tidy$variable_expl=as.character(life_table_tidy$variable)
for (k in unique(life_table_tidy$variable)){
  life_table_tidy$variable_expl[life_table_tidy$variable %in% k]=lifetable_varlist[[2]][lifetable_varlist[[1]] %in% k]
  }
ggplot(life_table_tidy,aes(x=age,y=value,group=disc_rate,color=factor(disc_rate),linetype=factor(disc_rate))) + 
  geom_line(size=1.2) + facet_wrap(~variable_expl,scales='free') + 
  theme_bw() + standard_theme + xlab('age (in months)') + ylab('') + 
  ggtitle(paste(cntr_sel,'demographics')) + labs(color='discount rate',linetype='discount rate')
# ggsave(paste0("output/life_table_",f_country_iso,"_2020_2025_disc0p03.png"),width=30,height=18,units="cm")
######################################################
# get unique country codes
country_opt <- data.frame(country_iso = unique(sim_config_matrix$country_iso)); 
burden_cntr_ind=which(country_opt$country_iso %in% cntr_sel)
incidence_one_table=get_incidence(country_opt$country_iso[burden_cntr_ind],output_dir)
# incidence data for all LMICs
# RSV_burden_Shi_2017=read_csv('input/RSV_burden_Shi_2017.csv')
# histogram: ggplot(RSV_burden_Shi_2017,aes(x=incidence_RSV_associated_ALRI)) + geom_histogram(binwidth=2)
# lineplot: national averages with CIs
RSV_burden_Shi_2017_tidy=read_csv('input/RSV_burden_Shi_2017_tidy.csv')
RSV_burden_Shi_2017_tidy$variable[grepl('incidence',RSV_burden_Shi_2017_tidy$variable)]='incidence of RSV-assoc. ALRI per 1000'
RSV_burden_Shi_2017_tidy$variable[grepl('episode',RSV_burden_Shi_2017_tidy$variable)]='number of episodes'
standard_theme_mod=standard_theme; standard_theme_mod$axis.text.x$size=4.5; standard_theme_mod$axis.text.x$hjust=0.99
standard_theme_mod$axis.text.x$vjust=0.5
non_na_vals=!(is.na(RSV_burden_Shi_2017_tidy$value) | is.na(RSV_burden_Shi_2017_tidy$lower_CI))
ggplot(RSV_burden_Shi_2017_tidy[non_na_vals,],aes(x=location_name,y=value,ymin=lower_CI,ymax=upper_CI,group=1)) + 
  geom_pointrange(aes(color=variable),size=0.6,fatten=0.1,fill="white",shape=22) + scale_color_manual(values=c('blue','red')) +
  facet_wrap(~variable,nrow=2,scales='free') + theme_bw() + standard_theme_mod + 
  scale_y_continuous(trans='log10') + xlab('') + ylab('') + guides(color=FALSE)
# ggsave('output/incid_per_cntr.png',width=30,height=18,units="cm")
# theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5))
# coord_cartesian(ylim=c(0, 7))
######################################################
#' ## Incidence data from partners
#' ### Kenya
kenya_data_path='../path_rsv_data/SARI_Rates_2010_2018/SARI_Rates_2010_2018_tidydata_cleaned.csv'
SARI_Rates_2010_2018_tidydata=read_csv(kenya_data_path)
# remove summary age groups
nonsumm_truthvals=!grepl("<|24-59|12-23",SARI_Rates_2010_2018_tidydata$age_in_months) | 
  grepl("<1$",SARI_Rates_2010_2018_tidydata$age_in_months)
KEMRI_kenya_rsv_incidence_ageinf=SARI_Rates_2010_2018_tidydata[nonsumm_truthvals & # don't include summary variables
                                          SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & # national average
                                          SARI_Rates_2010_2018_tidydata$RSV_association,] # RSV-associated
KEMRI_kenya_rsv_incidence_ageinf$age_in_months=factor(KEMRI_kenya_rsv_incidence_ageinf$age_in_months,
                                                      levels=unique(KEMRI_kenya_rsv_incidence_ageinf$age_in_months))
KEMRI_kenya_rsv_incidence_ageinf$period=factor(KEMRI_kenya_rsv_incidence_ageinf$period,
                                               levels=unique(KEMRI_kenya_rsv_incidence_ageinf$period))
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.character(KEMRI_kenya_rsv_incidence_ageinf$age_in_months)
KEMRI_kenya_rsv_incidence_ageinf$age_inf[KEMRI_kenya_rsv_incidence_ageinf$age_inf %in% "<1"]='0'
KEMRI_kenya_rsv_incidence_ageinf$freq=1; 
KEMRI_kenya_rsv_incidence_ageinf$freq[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',
                                  KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'))),ncol=2,byrow=T))+1
KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  as.numeric(sapply(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'),'[[',1))
KEMRI_kenya_rsv_incidence_ageinf=KEMRI_kenya_rsv_incidence_ageinf %>% uncount(weights=freq, .id="n",.remove=F)
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.numeric(KEMRI_kenya_rsv_incidence_ageinf$age_inf)+(KEMRI_kenya_rsv_incidence_ageinf$n-1)
# create per capita matrix
kemri_rsv_incidence_per_capita=array(t(KEMRI_kenya_rsv_incidence_ageinf[KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'&
                                       KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,'rate']))
popul_denom=1e5; CI95_const=1.96
if (max(kemri_rsv_incidence_per_capita>1)){ kemri_rsv_incidence_per_capita=kemri_rsv_incidence_per_capita/popul_denom }
kemri_rsv_incidence_CIs=KEMRI_kenya_rsv_incidence_ageinf[KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'&
                                  KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,c('CI_lower','CI_upper')]/popul_denom
# CI_lower = mu - 1.96*stdev/sqrt(sample size)
stdev_est=cbind((kemri_rsv_incidence_per_capita-kemri_rsv_incidence_CIs[,1])/CI95_const,
                (kemri_rsv_incidence_CIs[,2] - kemri_rsv_incidence_per_capita)/CI95_const)
kemri_rsv_incidence_stdevs=rowMeans(stdev_est)
# create matrix with 5000 iterations; max age in months is 60 (5yrs)
n_iter=5e3; age_maxval=60
# what distribution are we assuming? McMarcel data is gamma-distributed
kemri_incid_rate_matrix=sapply(1:n_iter, function(iters) {
  sapply(1:age_maxval, function(x) {rnorm(1,mean=kemri_rsv_incidence_per_capita[x],sd=kemri_rsv_incidence_stdevs[x] )})})
# write_csv(kemri_incid_rate_matrix,'input/kemri_incid_rate_matrix.csv')
####
# kemri hospitalisation rate
kemri_hosp_rate=KEMRI_kenya_rsv_incidence_ageinf %>% group_by(age_in_months,period) %>% 
  summarise(hosp_rate=rate[hospitalisation==TRUE]/(rate[hospitalisation==FALSE]+rate[hospitalisation==TRUE]))
kemri_hosp_val=as.numeric(unique(round(array(kemri_hosp_rate[kemri_hosp_rate$period %in% '2010-2018','hosp_rate']),2)))
# how much variation around this value? we dont know, lets use values from mcmarcel (there mean is 9%, ours is 24%)
kemri_hosp_rate_stdev=0.05656925 # value from: unique(apply(config$hosp_prob,1,sd))
# generate a hosp matrix with 5000 samples from our data
kemri_hosp_rate_matrix=sapply(1:n_iter, function(x) {rep(rnorm(1,mean=kemri_hosp_val,
                                                    sd=unique(round(rep(kemri_hosp_rate_stdev,age_maxval),4))),age_maxval)})
# write_csv(kemri_hosp_rate_matrix,'input/kemri_hosp_rate_matrix.csv')
######################################################
#' ### South Africa
s_afr_incidence_data=read_csv('../path_rsv_data/s_afr_incidence_data.csv')
s_afr_nat_average_totalcases=s_afr_incidence_data[s_afr_incidence_data$Province %in% 'South Africa' 
                                                  & s_afr_incidence_data$data_type %in% 'total',]
safr_hosp_rate_province_means=read_csv('input/safr_hosp_rate_province_means.csv')
# the incidence rate in data file is per 100K, normalize to per capita
if (max(s_afr_nat_average_totalcases$rate)>1){ s_afr_nat_average_totalcases[,c("rate","rate_CI_lower","rate_CI_upper")]=
  s_afr_nat_average_totalcases[,c("rate","rate_CI_lower","rate_CI_upper")]/popul_denom}
# st devs from CI95s
s_afr_nat_average_totalcases[,"stdev"]=rowMeans(cbind((s_afr_nat_average_totalcases$rate-s_afr_nat_average_totalcases$rate_CI_lower)/CI95_const,
                                     (s_afr_nat_average_totalcases$rate_CI_upper - s_afr_nat_average_totalcases$rate)/CI95_const))
# we need a 60x5000 matrix of per capita RSV cases - need to pick a distribution to generate random samples
s_afr_incid_rate_matrix=sapply(1:n_iter, function(iters) {
  sapply(1:age_maxval, function(x) {rnorm(1,mean=s_afr_nat_average_totalcases$rate[x],sd=s_afr_nat_average_totalcases$stdev[x])})})
# hosp rate
# we don't have a stdev value, i'll use MCMARCEL again
safr_hosp_rate_stdev=kemri_hosp_rate_stdev
s_afr_hosp_rate_matrix=sapply(1:n_iter, function(x) {
  rep(rnorm(1,mean=safr_hosp_rate_province_means$mean_hosp_rate[safr_hosp_rate_province_means$Province %in% 'South Africa'],
      sd=safr_hosp_rate_stdev),age_maxval)})
####
# put all matrices into 1 list; 1 matrix: burden_list_own_data[[2]][[2]]
burden_list_own_data=list(list(kemri_incid_rate_matrix,kemri_hosp_rate_matrix),
                          list(s_afr_incid_rate_matrix,s_afr_hosp_rate_matrix))
######################################################
######################################################
#' ## BURDEN CALCULATION with OWN DATA
# South Africa was not in the original analysis, so we cannot compare to default mcmarcel output
# cntrs_cea=c('KEN','ZAF'); n_cntr_output=2
# read in config parameters for burden calculation
burden_cntr_ind=which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])
# MV=sim_config_matrix[burden_cntr_ind[1],]; mAb=sim_config_matrix[burden_cntr_ind[2],]
# which intervention? 1=MatVacc, 2=monocl Abs
n_interv=2; sel_interv=sim_config_matrix[burden_cntr_ind[n_interv],]; config <- get_rsv_ce_config(sel_interv)
# if (n_cntr_output==2){sel_interv$country_iso=cntrs_cea[n_cntr_output]}
# with original MCMARCEL: 
# sim_output=get_burden(sel_interv)
####
# need to provide config$rsv_rate [60*5e3], config$hosp_prob [60*5e3] as arguments to get_burden_flexible
source('functions/get_burden_flexible.R') # source('functions/RSV_get_cost_data.R')
# for own data
sim_output_flexible=get_burden_flexible(sel_interv,burden_list_own_data[[n_cntr_output]][[1]],
                                        burden_list_own_data[[n_cntr_output]][[2]])
# for original
sim_output=get_burden_flexible(sel_interv,NA,NA)
######################################################
#' ## Process & plot results
# calculate ICER
icercolname='net_cost/DALY_averted'; icercols=c('incremental_cost_0to1y','total_DALY_0to1y_averted')
sim_output[,icercolname]=sim_output[,icercols[1]]/sim_output[,icercols[2]]
sim_output_flexible[,icercolname]=sim_output_flexible[,icercols[1]]/sim_output_flexible[,icercols[2]]
# outputs to display
cols_burden_sel=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y',
          "cost_rsv_hosp_0to1y",'total_medical_cost_0to1y','incremental_cost_0to1y','total_medical_cost_averted',
          "hosp_cases_averted", "rsv_deaths_averted", "total_DALY_0to1y_averted",icercolname) 
# "total_YLD_1to5y_averted","total_YLL_1to5y_averted",
# histograms: mcmarcel vs kemri
burden_mcmarcel_owndata_comp=melt(rbind(cbind(sim_output[,cols_burden_sel],data.frame(source='mcmarcel',iter=1:nrow(sim_output))),
      cbind(sim_output_flexible[,cols_burden_sel],data.frame(source=paste0(sel_interv$country_iso,' (own)')),iter=1:nrow(sim_output_flexible) ) ),
      id.vars=c('iter','source'))
# mean values
mean_intercepts=burden_mcmarcel_owndata_comp %>%  group_by(source,variable) %>% summarize(int = mean(value))
mean_intercepts[,'colorval']='red'; mean_intercepts$colorval[mean_intercepts$source %in% 'kemri']='green'
# PLOT histograms
if (n_interv==1) {interv_tag='_Mat_Vacc' } else {interv_tag='_monocl_Ab'}
ggplot(burden_mcmarcel_owndata_comp,aes(x=value,group=source)) +
  # geom_histogram(aes(y=..density..,fill=source),color="NA",size=0.4) + 
  geom_freqpoly(aes(color=source),size=1.2) + 
  geom_vline(data=mean_intercepts,aes(xintercept=int,linetype=source),size=0.8) + # color='black'
  facet_wrap(~variable,scales='free',labeller=label_wrap_gen(width=10)) + # ,ncol=3
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=7),
                     axis.title=element_text(size=14),text=element_text(family="Calibri")) +
  scale_linetype_manual(values=c('solid','dotdash')) +
  geom_rect(data=subset(burden_mcmarcel_owndata_comp, variable %in% icercolname),fill=NA,colour="blue",size=2,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
  ggtitle(paste(sel_interv$country_iso,'RSV burden & intervention estimates:',gsub('_','',interv_tag))) +
  labs(color='data source',linetype='mean') + guides(xintercept=FALSE,linetype=guide_legend(ncol=2)) # xlab('')+ylab('')
######
# save plot
# cea_plot_filename=paste("output/cea_plots/",sel_interv$country_iso,"_mcmarcel_burden_estimates_1000samples",interv_tag,".png",sep="")
# if (!dir.exists('output/cea_plots')) {dir.create('output/cea_plots')}
# ggsave(cea_plot_filename,width=30,height=18,units="cm")