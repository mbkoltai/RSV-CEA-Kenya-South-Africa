# MCMARCEL predictions
# for all LMICs
# incidence_lmic_ts_n5000=as.data.frame(read_csv("./input/incidence_lmic_ts_n5000.csv"))
# cfr_lmic_ts_n5000=as.data.frame(read_csv("./input/cfr_lmic_ts_n5000.csv"))
# hosp_prob_ts_n5000=as.data.frame(read_csv("input/hosp_prob_ts_n5000.csv"))
# # mean by month
# age_maxval=nrow(incidence_month_means)
# incidence_month_means=incidence_lmic_ts_n5000 %>% group_by(mos) %>% 
#   summarise(meanpred=mean(pred),minpred=mean(pred)-sd(pred),maxpred=mean(pred)+sd(pred))
# incidence_month_means$type='cases_per_1000'
# cfr_month_means=cfr_lmic_ts_n5000 %>% 
#   group_by(mos) %>% summarise(meanpred=mean(pred),minpred=mean(pred)-sd(pred),maxpred=mean(pred)+sd(pred))
# cfr_month_means$type='hospit_cfr'
# hosp_prob_means = hosp_prob_ts_n5000 %>% group_by(mos) %>% 
#   summarise(meanpred=mean(pred),minpred=mean(pred)-sd(pred),maxpred=mean(pred)+sd(pred))
# hosp_prob_means$type='probab_hospit'
# incidence_cfr_hospprob_means=rbind(incidence_month_means,cfr_month_means,hosp_prob_means)
# # INCIDENCE by age for ALL LMICs
# ggplot(incidence_cfr_hospprob_means,aes(x=mos,y=meanpred)) + geom_line() + geom_point() +
#   geom_ribbon(aes(ymin=minpred,ymax=maxpred),alpha=0.3,colour=NA,fill="red") + facet_wrap(~type,scales='free',nrow=2) + 
#   theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),text=element_text(family="Calibri"),
#         plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=11,angle=90,vjust=0.5)) + 
#   scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
#   xlab('month') + ylab('incidence') + ggtitle('LMIC RSV incidence, CFR, hospitalisation (mean [-/+ stdev])')
######
# path
library(rstudioapi); currentfile_path=rstudioapi::getActiveDocumentContext()$path
currentfile_path=paste0(unlist(strsplit(
  currentfile_path,"\\/"))[1:(length(unlist(strsplit(currentfile_path,"\\/")))-1)],collapse="/")
setwd(currentfile_path)
# libraries
library(tidyverse); library(reshape2)
# run tag
run_tag  <- 'RSV_gavi72_basecase'  # 72 Gavi countries (basecase) 
# number of stochastic samples in the probabilistic sensitivity analysis (PSA)
num_sim <- 1000
# random number generater seed
rng_seed <- gsub('-','',Sys.Date()) # 20190118
# option to create geographical and country-specific plots
# note: this might require substantial processing time
boolean_country_plots <- FALSE # TRUE
boolean_map_plots     <- FALSE # TRUE
source('functions/RSV_load_all.R')
# output directory postfix
output_dir_postfix <- paste0(run_tag,'_n',num_sim)
# add timestap to output directory name
output_dir <- paste0('output/',format(Sys.time(),'%m%d%H%M%S_'),output_dir_postfix)
# set seed
set.seed(rng_seed)
# config filename
config_filename <- paste0('./config/',run_tag,'.csv')
time_stamp_main <- Sys.time()
# always clear temporary results
cli_print('Clear all temporary output'); unlink(file.path(get_temp_output_folder(output_dir)),recursive = T)
# start parallel workers: if running for many cntrs and we want to parallelise
start_parallel_workers()
# 144x12 table
sim_config_matrix <- read.table(config_filename,sep=',', dec='.',stringsAsFactors = F,header = T)
# set output file name prefix
sim_output_filename  <- file.path(output_dir,run_tag)
# add simulation details
sim_config_matrix$num_sim <- num_sim; sim_config_matrix$scenario_id <- 1:nrow(sim_config_matrix)
sim_config_matrix$rng_seed <- rng_seed; sim_config_matrix$outputFileDir <- get_output_folder(output_dir)
# create UN country data
create_UN_country_database(output_dir)
# pre-process WPP2017 data
load_wpp2017_databases(output_dir)
# loads sex ratio and mortality by age groups; dataframes: sexratio, mxM, mxF
country_year_opt <- sim_config_matrix[,c('country_iso','year')]
country_year_opt <- rbind(country_year_opt, cbind(country_iso=country_year_opt$country_iso, year=2015))
# get unique combinations: this creates a table with 72 cntrs with dates 2015 and 2020
country_year_opt <- unique(country_year_opt)
# make summary matrix, including the 5-year period notation
# 3 columns: countries, '2020-25', 2020
country_period_opt    <- cbind(country_year_opt$country_iso, t(sapply(country_year_opt$year,get_year_category)))



ken_inds=which(country_opt$country_iso %in% 'KEN')
incidence_one_table=get_incidence(country_opt$country_iso[ken_inds],output_dir)
# list of 3 elements
# this contains: RSV_rate [60x5000 dataframe], hosp_prob [60x5000 dataframe], hCFR_prob [60x1000 dataframe]
# incidence data is from: this is the incidence (+CIs) by country
RSV_burden_Shi_2017=read_csv('input/RSV_burden_Shi_2017.csv')
RSV_burden_Shi_2017_tidy=read_csv('output/RSV_burden_Shi_2017_tidy.csv')
# NATIONAL AVERAGES (incidence)
ggplot(RSV_burden_Shi_2017_tidy,aes(x=location_name,y=value,group=1)) + geom_line() + geom_point(size=1) +
  geom_ribbon(aes(ymin=lower_CI,ymax=upper_CI),alpha=0.3,colour=NA,fill="red") + facet_wrap(~variable,nrow=2,scales='free') +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5)) + scale_y_continuous(trans='log10')
#####################
# Kenya RSV data
# LMIC incidence by age scaled by data from 'RSV_burden_Shi_2017'
kenya_rate_files=c('input/country_rsv_rate.csv', 'input/hosp_prob_mat.csv', 'input/cfr_mat.csv')
# country_rsv_rate=read_csv(kenya_rate_files[1]);
hosp_prob_mat=read_csv(kenya_rate_files[2]); cfr_mat=read_csv(kenya_rate_files[3])
f_country_iso='KEN'; f_outputFileDir="output/1005164300_RSV_gavi72_basecase_n100" 
burden_country <- read.csv('input/RSV_burden_Shi_2017.csv',stringsAsFactors = FALSE);
burden_country_reference <- burden_country$incidence_RSV_associated_ALRI[burden_country$country_iso==f_country_iso]
# Life table
life_table <- get_life_table(f_country_iso,2015,f_outputFileDir,0); under5_pop <- life_table$lx[1:60]
# incidence matrix for all LMICs
incidence_mat <- convert_pred_into_model_input("./input/incidence_lmic_ts_n5000.csv") # spline_datafiles$incidence
# cases = incidence per 1000 * population by age
country_rsv_cases <- (incidence_mat/1000)*under5_pop
# rsv rate = cases / total cases (per column)
# take # of cases in an age group and divide by TOTAL number of cases 
# --> this gives the FRACTION of RSV cases per age group out of ALL RSV cases
country_rsv_fraction <- country_rsv_cases / rep(colSums(country_rsv_cases),each=dim(country_rsv_cases)[1])
# country rate = (rsv rate) * (country_reference per 1000) * (country population)
# (burden_country_reference / 1000) is predicted burden per 1 million -->  x population = predicted total burden
# country_rsv_rate is the distribution of cases by age groups (as fractions of total) so product distributes total burden
# this is cases per age group, but # of ppl in age groups not the same!
country_rsv_pred_cases <- country_rsv_fraction * (burden_country_reference / 1000) * sum(under5_pop)
# so the colsum is constant, being the total burden
# rsv rate = rsv cases / population (cases/capita)
country_rsv_rate=data.frame(country_rsv_pred_cases/under5_pop); age_maxval=nrow(country_rsv_rate)
# collate the tables
kenya_burden_abs=data.frame(age=1:age_maxval,mean=rowMeans(country_rsv_pred_cases),
                            sd=apply(country_rsv_pred_cases,1,sd),type='number_of_cases');
kenya_burden_percap=data.frame(age=1:age_maxval,mean=rowMeans(country_rsv_rate),sd=apply(country_rsv_rate,1,sd),type='case_per_capita')
hosp_prob_mat_tidy=data.frame(age=1:age_maxval,mean=rowMeans(hosp_prob_mat),
                              sd=apply(hosp_prob_mat,1,sd),type='hosp_admissions_per_capita');
cfr_mat_tidy=data.frame(age=1:age_maxval,mean=rowMeans(cfr_mat),sd=apply(cfr_mat,1,sd),type='cfr_hosp_admissions')
kenya_burden_all=rbind(kenya_burden_abs,kenya_burden_percap,hosp_prob_mat_tidy,cfr_mat_tidy)
# KENYA RSV incidence rates
ggplot(kenya_burden_all,aes(x=age,y=mean)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.3,colour=NA,fill="red") +
  facet_wrap(~type,scales='free',nrow=2) + scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
        axis.text.x=element_text(size=11,angle=90,vjust=0.5),axis.text.y=element_text(size=11),
        axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  xlab('Age (month)') + ylab('Total burden') + labs(color="Samples") + ggtitle('Kenya RSV burden (mean +/- stdev)')
# ggsave("output/RSV_burden_kenya_randomsamples.png",width=30,height=18,units="cm") 

#####
# normalise Kenya data to /1000 ppl
lmic_kenya_burden=kenya_burden_all; case_cols=lmic_kenya_burden$type %in% 'number_of_cases'
lmic_kenya_burden=lmic_kenya_burden[!case_cols,]
# scale to 1000
if (max(lmic_kenya_burden$mean)<1){
  # to have incidence per 1000, multiply per capita rate by 1e3
  lmic_kenya_burden[c('mean','sd')]=lmic_kenya_burden[,c('mean','sd')]*1e3
  # to have incidence per 1000, divide case number (which was for 100e3) by 100
  # lmic_kenya_burden[case_cols,c('mean','SD')]=lmic_kenya_burden[case_cols,c('mean','sd')]/1e2
  lmic_kenya_burden$type=str_replace_all(str_replace_all(lmic_kenya_burden$type,'cfr','deaths_per_1000'),'per_capita','per_1000')
  lmic_kenya_burden$type=str_replace_all(lmic_kenya_burden$type,'case','cases')
}
# collate with all LMIC data
lmic_kenya_burden$country='Kenya'
lmic_incidence_means=data.frame(age=1:age_maxval,mean=rowMeans(incidence_mat),sd=apply(incidence_mat,1,sd),type='cases_per_1000')
lmic_incidence_means$country='lmic'
# rbind two dfs
if (length(unique(lmic_kenya_burden$country))==1){ lmic_kenya_burden=rbind(lmic_kenya_burden,lmic_incidence_means) }

### plot
ggplot(lmic_kenya_burden,aes(x=age,y=mean,color=country)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=country),alpha=0.3,colour=NA) +
  facet_wrap(~type,scales='free',nrow=2)+ scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=11,angle=90,vjust=0.5),axis.text.y=element_text(size=11),
                     axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  xlab('Age (month)') + ylab('incidence') + labs(color="Samples") + ggtitle('Kenya RSV burden (mean +/- stdev)')
ggsave("output/RSV_incidence_kenya_lmic_comparison.png",width=30,height=18,units="cm") 

#############################################################################
# our own data from KEMRI
kenya_data_path='../path_rsv_data/SARI_Rates_2010_2018/SARI_Rates_2010_2018_tidydata_cleaned.csv'
SARI_Rates_2010_2018_tidydata=read_csv(kenya_data_path)
nonsumm_truthvals=!grepl("<|24-59|12-23",SARI_Rates_2010_2018_tidydata$age_in_months) | 
  grepl("<1$",SARI_Rates_2010_2018_tidydata$age_in_months)
kenya_rsv_data_timemeans=SARI_Rates_2010_2018_tidydata[nonsumm_truthvals & # don't include summary variables that subsume others
                                                         SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & 
                                                         SARI_Rates_2010_2018_tidydata$RSV_association & 
                                                         SARI_Rates_2010_2018_tidydata$period %in% '2010-2018',]
kenya_rsv_data_timemeans$age_inf=kenya_rsv_data_timemeans$age_in_months;
kenya_rsv_data_timemeans$age_inf[kenya_rsv_data_timemeans$age_inf %in% '<1']=0
kenya_rsv_data_timemeans$age_inf[grepl('-',kenya_rsv_data_timemeans$age_inf)]=rowMeans(matrix(as.numeric(unlist(
  sapply(kenya_rsv_data_timemeans$age_inf[grepl('-',kenya_rsv_data_timemeans$age_inf)], function(x) {strsplit(x,'-')}))),
  ncol=2,byrow=T)); kenya_rsv_data_timemeans$age_inf=as.numeric(kenya_rsv_data_timemeans$age_inf)
# rates are per 100k popul
# predicted values: kenya_burden_all
kenya_pred_data_comb = kenya_rsv_data_timemeans[,c('age_inf','rate','CI_lower','CI_upper','hospitalisation')];
datasum=kenya_pred_data_comb %>% group_by(age_inf) %>% 
  summarize(rate=sum(rate),CI_lower=sum(CI_lower),CI_upper=sum(CI_upper),hospitalisation='RSV_LRTI')
if (length(unique(kenya_pred_data_comb$hospitalisation))==2){
  kenya_pred_data_comb$hospitalisation[kenya_pred_data_comb$hospitalisation]='hosp_RSV_LRTI'
  kenya_pred_data_comb$hospitalisation[kenya_pred_data_comb$hospitalisation %in% 'FALSE']='nonhosp'
  kenya_pred_data_comb=rbind(kenya_pred_data_comb,datasum) 
}
kenya_pred_data_comb$datasource='KEMRI_data'
##########################
# MCMARCEL predictions
kenya_pred=kenya_burden_all[kenya_burden_all$type %in% 'cases_per_cap',c("age","mean","sd")]; 
kenya_pred$CI_lower=kenya_pred$mean-kenya_pred$sd; kenya_pred$CI_upper=kenya_pred$mean+kenya_pred$sd
kenya_pred=kenya_pred[,c('age','mean','CI_lower','CI_upper')]; kenya_pred$hospitalisation='RSV_LRTI'
kenya_pred$datasource='MCMARCEL_pred'
colnames(kenya_pred)=colnames(kenya_pred_data_comb); 
# normalize to cases / 1000
popul_denom=1000; popul_denom_kemri=1e5;  
if (length(unique(kenya_pred_data_comb$datasource))==1){
  kenya_pred[,c("rate","CI_lower","CI_upper")]=kenya_pred[,c("rate","CI_lower","CI_upper")]*popul_denom
  kenya_pred_data_comb[,c("rate","CI_lower","CI_upper")]=
    kenya_pred_data_comb[,c("rate","CI_lower","CI_upper")]*(popul_denom/popul_denom_kemri)
  kenya_pred_data_comb=rbind(kenya_pred_data_comb,kenya_pred)
}
kenya_pred_data_comb[,'infection_stage']='unknown'

# data from Nokes 2008 article (https://doi.org/10.1086/524019)
incidence_nokes2008=read_csv('input/incidence_nokes2008.csv');incidence_nokes2008=incidence_nokes2008[,!grepl('perc',colnames(incidence_nokes2008))]
incidence_nokes2008[,'datasource']='nokes2008'; incidence_nokes2008_tidy=melt(incidence_nokes2008,id.vars=c('age','CYO','datasource','infection_stage'))
# normalisation by CYO is probably wrong, lets work with the absolute numbers
incidence_nokes2008_tidy[,'rate']=incidence_nokes2008_tidy$value # /incidence_nokes2008_tidy$CYO)*1000
incidence_nokes2008_tidy[,'hospitalisation']='RSV_LRTI' # "hosp","nonhosp","both"   
incidence_nokes2008_tidy$hospitalisation[incidence_nokes2008_tidy$variable %in% 'hospitalisation']='hosp_RSV_LRTI'
# 'nmbr_RSV_infections','nmbr_RSV_LRTI','RSV_severe_LRT'
incidence_nokes2008_tidy$hospitalisation[incidence_nokes2008_tidy$variable %in% 'nmbr_RSV_infections']='all_infections'
incidence_nokes2008_tidy$hospitalisation[incidence_nokes2008_tidy$variable %in% 'nmbr_RSV_LRTI']='RSV_LRTI'
incidence_nokes2008_tidy$hospitalisation[incidence_nokes2008_tidy$variable %in% 'RSV_severe_LRTI']='severe_RSV_LRTI'
incidence_nokes2008_tidy[,'datasource']='nokes2008data'
incidence_nokes2008_tidy[,'age_inf']=rowMeans(matrix(as.numeric(unlist(strsplit(incidence_nokes2008_tidy$age,"-"))),ncol=2,byrow=T))
incidence_nokes2008_tidy$rate=
  incidence_nokes2008_tidy$rate/array(diff(t(matrix(as.numeric(unlist(strsplit(incidence_nokes2008_tidy$age,"-"))),ncol=2,byrow=T))))
# scale to 1000 ppl
cohort_size=635; incidence_nokes2008_tidy$rate=incidence_nokes2008_tidy$rate*(1e3/cohort_size)
incidence_nokes2008_tidy[,c("CI_lower","CI_upper")]=incidence_nokes2008_tidy$rate
# create a sum of both primary and reinfection
# incidence_nokes2008_tidy %>% group_by(age) %>% summarise(rate=sum(rate),infection_stage='both')

###
# nokes 2009 article
incidence_nokes2009=read_csv('input/incidence_nokes2009.csv')
incidence_nokes2009[,'rate']=(incidence_nokes2009$positiveresult/incidence_nokes2009$no_tested)*1e3
incidence_nokes2009[,'age_inf']=rowMeans(matrix(as.numeric(unlist(strsplit(incidence_nokes2009$age,"-"))),ncol=2,byrow=T))
incidence_nokes2009[,c("CI_lower","CI_upper")]=incidence_nokes2009$rate; incidence_nokes2009[,'datasource']='nokes2009'
incidence_nokes2009$infection_stage='unknown'; incidence_nokes2009$hospitalisation=paste('hosp',incidence_nokes2009$severity,sep='_')
incidence_nokes2009[,'age_inf']=rowMeans(matrix(as.numeric(unlist(strsplit(incidence_nokes2009$age,"-"))),ncol=2,byrow=T))
incidence_nokes2009[,colnames(kenya_pred_data_comb)]

# rbind different data sources
if (!(any(grepl('nokes2008data',kenya_pred_data_comb$datasource)) | any(grepl('nokes2009data',kenya_pred_data_comb$datasource)))) {
  kenya_pred_data_comb=rbind(kenya_pred_data_comb,incidence_nokes2008_tidy[,colnames(kenya_pred_data_comb)],
                             incidence_nokes2009[,colnames(kenya_pred_data_comb)])
}

# create column with data datasource
kenya_pred_data_comb$type='per 1000 population'
kenya_pred_data_comb$type[kenya_pred_data_comb$datasource %in% "nokes2008data"] = 'per 1000 population (birth cohort, n=635)'
kenya_pred_data_comb$type[kenya_pred_data_comb$datasource %in% "nokes2009"] = 'per 1000 hospital admissions'
kenya_pred_data_comb$type=factor(kenya_pred_data_comb$type)

####
# PLOT
data_plot=kenya_pred_data_comb[!(kenya_pred_data_comb$hospitalisation %in% c('nonhosp')),]
data_plot=data_plot[!((data_plot$datasource %in% "KEMRI_data" | data_plot$datasource %in% "MCMARCEL_pred") & data_plot$age_inf>42),]
y_step=1e2; ymaxval=ceiling(max(data_plot$rate)/y_step)*y_step; ymin=0; xlimval=30
####
ggplot(data_plot,aes(x=age_inf,y=rate,color=hospitalisation,shape=datasource)) + # linetype=infection_stage,
  geom_line(size=1) + geom_point(size=2.5) + # aes(linetype=infection_stage)
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=hospitalisation),alpha=0.3,linetype=0) +
  facet_wrap(~type+infection_stage,scales='free',ncol=2) + # infection_stage,labeller=label_wrap_gen(multi_line = T)
  scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  # scale_y_continuous(labels=as.character(seq(ymin,ymaxval,y_step)),breaks=seq(ymin,ymaxval,y_step)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=13,angle=90,vjust=0.5),axis.text.y=element_text(size=13),
                     axis.title=element_text(size=16),text=element_text(family="Calibri"),
                     legend.title=element_text(size=14),legend.text=element_text(size=10),
                     strip.text=element_text(size=12)) + 
  guides(colour=guide_legend(override.aes=list(shape=NA))) + #coord_cartesian(xlim=c(0.5,46)) 
  xlab('Age (month)') + ylab('incidence') + ggtitle('Kenya RSV incidence under 5yr') 
# ylim=c(-0.5,ymaxval),expand=F
# labs(color="Data source",fill="Data source")
#######
ggsave("output/RSV_burden_kenya_data_pred_compare_datasources.png",width=32,height=20,units="cm") 
# ggsave("output/RSV_burden_kenya_data_pred_compare_nokes2008_2009.pdf",width=30,height=24,units="cm",device=grDevices::pdf) 

#############################################################################
# calculate disease burden
# first with MCMARCEL data
sim_config_matrix <- read.table("./config/RSV_gavi72_basecase.csv",sep=',',dec='.',stringsAsFactors=F,header=T)
sim_config_matrix$num_sim=num_sim; sim_config_matrix$scenario_id=1:nrow(sim_config_matrix)
output_dir= "output/1005164300_RSV_gavi72_basecase_n100"
sim_config_matrix$rng_seed=rng_seed;sim_config_matrix$outputFileDir=get_output_folder(output_dir)
# sim_config_matrix has the parameters:
# "config_tag", "scenario", "intervention", "efficacy_maternal", "efficacy_infant"
# [6] "dur_protection_maternal", "dur_protection_infant", "price_mAb_sens_factor", "year", "country_iso"
# [11] "coverage_maternal", "coverage_infant", "num_sim", "scenario_id", "rng_seed"
# [16] "outputFileDir"
# total number of cases under 5yr
sum(kenya_burden_all[kenya_burden_all$type %in% 'number_of_cases',]$mean)

# kenya intervention indices
ken_inds=which(sim_config_matrix$country_iso %in% 'KEN')
# first is maternal coverage, second infant vaccination
# both assuming 70% efficacy. maternalVacc: 5 months, infant Ab: 6m protection. 86% coverage. 
# incidence and hospitalisation rate should be normalised to per capita
# lets take the 2010-2018 year average of incidence
# popul_denom=1e5
kemri_rsv_incidence_per_capita=array(t(KEMRI_kenya_rsv_incidence_ageinf[KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'&
                                                                  KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,'rate']))
if (max(kemri_rsv_incidence_per_capita>1)){ kemri_rsv_incidence_per_capita=kemri_rsv_incidence_per_capita/popul_denom }
kemri_rsv_incidence_CIs=KEMRI_kenya_rsv_incidence_ageinf[KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'&
                                  KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,c('CI_lower','CI_upper')]/popul_denom
# CI_lower = mu - 1.96*
stdev_est=cbind((kemri_rsv_incidence_per_capita-kemri_rsv_incidence_CIs[,1])/1.96,
      (kemri_rsv_incidence_CIs[,2] - kemri_rsv_incidence_per_capita)/1.96)
kemri_rsv_incidence_stdevs=rowMeans(stdev_est)
# create matrix with 5000 iterations
kemri_incid_rate_matrix=sapply(1:n_iter, function(iters) {
  sapply(1:age_maxval, function(x) {rnorm(1,mean=kemri_rsv_incidence_per_capita[x],sd=kemri_rsv_incidence_stdevs[x] )})})
####
# kemri hospitalisation rate
# age_maxval=60
kemri_hosp_rate=KEMRI_kenya_rsv_incidence %>% group_by(age_in_months,period) %>% 
  summarise(hosp_rate=rate[hospitalisation==TRUE]/(rate[hospitalisation==FALSE]+rate[hospitalisation==TRUE]))
kemri_hosp_val=as.numeric(unique(round(array(kemri_hosp_rate[kemri_hosp_rate$period %in% '2010-2018','hosp_rate']),2)))
# kemri_hosp_rate_average=rep(kemri_hosp_val,age_maxval)
# how much variation around this value? we dont know, lets use values from mcmarcel (there the mean is 9%, ours is 24%)
hosp_rate_stdev=apply(config$hosp_prob,1,sd) # it's a single value
# generate a hosp matrix with 5000 samples from our data
n_iter=5e3
kemri_hosp_rate_matrix=sapply(1:n_iter, function(x) {rep(rnorm(1,mean=kemri_hosp_val,sd=unique(round(hosp_rate_stdev,4))),age_maxval)})
####################################################################################################
####################################################################################################
####################################################################################################
#####
# BURDEN CALCULATION with OWN DATA
ken_inds = which(sim_config_matrix$country_iso %in% 'KEN')
MV=sim_config_matrix[ken_inds[1],]; mAb=sim_config_matrix[ken_inds[2],]
n_interv=1
sel_interv=sim_config_matrix[ken_inds[n_interv],]
sim_output=get_burden(sel_interv)
####
# with own data
# we need to provide config$rsv_rate [60*5e3], config$hosp_prob [60*5e3], they are:
# kemri_hosp_rate_matrix, kemri_incid_rate_matrix
# source('functions/get_burden_flexible.R')
sim_output_flex=get_burden_flexible(sel_interv,kemri_incid_rate_matrix,kemri_hosp_rate_matrix) # get_burden_flexible(MV,NA,NA)
# how many of the outputs are unique? # sum(sapply(1:ncol(sim_output_flex),function(x) {length(unique(sim_output_flex[,x]))})>1)
### PROCESS RESULTS
# add net cost per DALY averted
icercolname='net_cost/DALY_averted'
sim_output[,icercolname]=sim_output$incremental_cost_0to1y/sim_output$total_DALY_0to1y_averted
sim_output_flex[,icercolname]=sim_output_flex$incremental_cost_0to1y/sim_output_flex$total_DALY_0to1y_averted
# outputs to display
cols_burden_sel=c('rsv_cases','hosp_cases','rsv_deaths','total_DALY_0to1y',
                  "cost_rsv_hosp_0to1y",'total_medical_cost_0to1y','incremental_cost_0to1y','total_medical_cost_averted', # 'intervention_cost'
                  "hosp_cases_averted", "rsv_deaths_averted", "total_DALY_0to1y_averted",icercolname) 
# "total_YLD_1to5y_averted","total_YLL_1to5y_averted",
# histograms: mcmarcel vs kemri
burden_mcmarcel_kemri_comp=melt(rbind(cbind(sim_output[,cols_burden_sel],data.frame(source='mcmarcel',iter=1:nrow(sim_output))),
      cbind(sim_output_flex[,cols_burden_sel],data.frame(source='kemri'),iter=1:nrow(sim_output_flex) ) ),id.vars=c('iter','source'))
# mean values
mean_intercepts <- burden_mcmarcel_kemri_comp %>%  group_by(source,variable) %>% summarize(int = mean(value))
mean_intercepts[,'colorval']='red'; mean_intercepts$colorval[mean_intercepts$source %in% 'kemri']='green'
# PLOT histograms
if (n_interv==1) {   interv_tag='_mat_vacc' } else {interv_tag='_monocl_Ab'}
ggplot(burden_mcmarcel_kemri_comp,aes(x=value,group=source)) + 
  geom_histogram(aes(y=..density..,fill=source),color="black",size=0.4) + 
  geom_vline(data=mean_intercepts, aes(xintercept=int,linetype=source,color=source),size=1.25) +
  # geom_density(aes(color=source),alpha=0.2) + # binwidth=1e-4
  facet_wrap(~variable,scales='free',labeller = label_wrap_gen(width=10)) + 
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                                axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=7),
                                axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  geom_rect(data=subset(burden_mcmarcel_kemri_comp, variable %in% icercolname),fill=NA,colour="blue",size=2,
      xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + ggtitle(paste('RSV burden and intervention estimates:',gsub('_','',interv_tag))) + 
  labs(fill='data source') + guides(xintercept=FALSE,linetype=FALSE,color=FALSE) # xlab('') + ylab('') 
######
# save
ggsave(paste("output/kemri_mcmarcel_burden_estimates_1000samples",interv_tag,".png",sep=""),width=30,height=18,units="cm")

# sample gamma
# sample_rgamma <- function(mean,stdev,num_sim){  alpha <- (mean^2) / (stdev^2) # shape parameter
#   beta  <- (stdev^2) / (mean) # rate parameter (beta=1/theta)
#   sample <- rgamma(num_sim,shape=alpha,rate = 1/beta) }

# daly for severe and nonsevere
# daly_vals=melt(data.frame(cbind(config$severe_rsv_DALYloss,config$non_severe_rsv_DALYloss),stringsAsFactors=F))
# daly_vals$variable=as.character(daly_vals$variable); daly_vals$variable[daly_vals$variable %in% 'X1']='severe RSV'; daly_vals$variable[daly_vals$variable %in% 'X2']='nonsevere RSV'
#
# ggplot(daly_vals,aes(x=value)) + geom_histogram(aes(y=..density..),binwidth=1e-4,color="black",fill="white") + 
#   geom_density(alpha=0.2,fill="#FF6666") + facet_wrap(~variable,scales='free_x')
#   geom_vline(aes(xintercept=mean(value)),color="blue", linetype="dashed", size=1) + 

# CFR
# ggplot(data.frame(rowMeans(config$hosp_CFR)),aes(x=1:60,y=rowMeans(config$hosp_CFR)*100)) + geom_line()

# our own RSV-SARI incidence data
nonsumm_truthvals=!grepl("<|24-59|12-23",SARI_Rates_2010_2018_tidydata$age_in_months) | 
  grepl("<1$",SARI_Rates_2010_2018_tidydata$age_in_months)
KEMRI_kenya_rsv_incidence=SARI_Rates_2010_2018_tidydata[nonsumm_truthvals & # don't include summary variables that subsume others
                                                         SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & # national average
                                                         SARI_Rates_2010_2018_tidydata$RSV_association,] # RSV-associated
KEMRI_kenya_rsv_incidence$age_in_months=factor(KEMRI_kenya_rsv_incidence$age_in_months,levels=unique(KEMRI_kenya_rsv_incidence$age_in_months))
KEMRI_kenya_rsv_incidence$period=factor(KEMRI_kenya_rsv_incidence$period,levels=unique(KEMRI_kenya_rsv_incidence$period))
# this is incidence per 100e5!!!
# ggplot(KEMRI_kenya_rsv_incidence,aes(x=age_in_months,y=rate,color=hospitalisation,group=hospitalisation)) + geom_line() +
#   geom_point() + geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=hospitalisation),alpha=0.3,colour=NA) + 
#   facet_wrap(~period,ncol=3) +
#   theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
#         axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
#         axis.title=element_text(size=14), text=element_text(family="Calibri")) + ggtitle('Kenya RSV incidence per 100k popul')
# ggsave("output/KEMRI_kenya_rsv_incidence.png",width=30,height=18,units="cm") 

# hospitalisation rate
kemri_hosp_rate=KEMRI_kenya_rsv_incidence %>% group_by(age_in_months,period) %>% 
  summarise(hosp_rate=rate[hospitalisation==TRUE]/(rate[hospitalisation==FALSE]+rate[hospitalisation==TRUE]))

# have all age groups 1 t0 60 by assuming months within larger age bands have the same incidence
KEMRI_kenya_rsv_incidence_ageinf=KEMRI_kenya_rsv_incidence
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.character(KEMRI_kenya_rsv_incidence_ageinf$age_in_months)
KEMRI_kenya_rsv_incidence_ageinf$age_inf[KEMRI_kenya_rsv_incidence_ageinf$age_inf %in% "<1"]='0'
# library(matrixStats)
KEMRI_kenya_rsv_incidence_ageinf$freq=1; 
KEMRI_kenya_rsv_incidence_ageinf$freq[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'))),ncol=2,byrow=T))+1
KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  as.numeric(sapply(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'),'[[',1))
KEMRI_kenya_rsv_incidence_ageinf=KEMRI_kenya_rsv_incidence_ageinf %>% uncount(weights=freq, .id="n",.remove=F)
# KEMRI_kenya_rsv_incidence_ageinf=as.data.frame(lapply(KEMRI_kenya_rsv_incidence_ageinf,rep,KEMRI_kenya_rsv_incidence_ageinf$ntimes))
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.numeric(KEMRI_kenya_rsv_incidence_ageinf$age_inf)+(KEMRI_kenya_rsv_incidence_ageinf$n-1)

# plot KEMRI data with inferred datapoints within multi-month age bands
ggplot(KEMRI_kenya_rsv_incidence_ageinf,aes(x=age_inf,y=rate,color=hospitalisation,group=hospitalisation)) + geom_line() + 
  geom_point(size=0.5) + geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=hospitalisation),alpha=0.3,colour=NA) + 
  facet_wrap(~period,ncol=3) + scale_x_continuous(labels=as.character(seq(0,age_maxval,5)),breaks=seq(0,age_maxval,5)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                axis.title=element_text(size=14), text=element_text(family="Calibri")) + ggtitle('Kenya RSV incidence per 100k popul')
# ggsave("output/KEMRI_kenya_rsv_incidence_infdatapoints.png",width=30,height=18,units="cm") 

#####
# we want to input our own incidence and hosp data to get_burden_flexible()
# currently it is
# config$rsv_rate [60*5000]: this is a RATE, so per capita
# config$hosp_prob [60*5000]
config_hosp_rate=config$rsv_rate*config$hosp_prob

# Kenya RSV rates (per capita) is in 'input/kenya_rsv_rate.csv

# compare to our data: KEMRI_kenya_rsv_incidence_ageinf[]
# KEMRI_kenya_rsv_incidence_ageinf$period
popul_denom=1e5
mcmarcel_incid_100k=data.frame(age_inf=1:60, period=2015, hospitalisation=FALSE, rate=rowMeans(config$rsv_rate)*popul_denom, 
                               CI_lower=(rowMeans(config$rsv_rate)-apply(config$rsv_rate,1,sd))*popul_denom, 
                               CI_upper=(rowMeans(config$rsv_rate)+apply(config$rsv_rate,1,sd))*popul_denom, datatype='mcmarcel')
# hospitalised
z=mcmarcel_incid_100k; z$hospitalisation=TRUE; 
z[,c('rate',"CI_lower","CI_upper")]=data.frame(rate=rowMeans(config_hosp_rate)*popul_denom,
                                               CI_lower=(rowMeans(config_hosp_rate)-apply(config_hosp_rate,1,sd))*popul_denom,
                                               CI_upper=(rowMeans(config_hosp_rate)+apply(config_hosp_rate,1,sd))*popul_denom)
if (sum(mcmarcel_incid_100k$hospitalisation)==0){ mcmarcel_incid_100k=rbind(mcmarcel_incid_100k,z) }
kemri_mcmarcel_compar=rbind(mcmarcel_incid_100k,
  cbind(KEMRI_kenya_rsv_incidence_ageinf[,c('age_inf','period',"hospitalisation","rate","CI_lower","CI_upper")],
        data.frame(datatype='kemri'))  )
kemri_mcmarcel_compar$period=factor(kemri_mcmarcel_compar$period,levels=unique(KEMRI_kenya_rsv_incidence$period))
kemri_mcmarcel_compar$hospitalisation[kemri_mcmarcel_compar$hospitalisation]='hospitalised'
kemri_mcmarcel_compar$hospitalisation[kemri_mcmarcel_compar$hospitalisation==FALSE]='non-hospitalised'
# linewidth
kemri_mcmarcel_compar$timespan='single year'; 
kemri_mcmarcel_compar$timespan[kemri_mcmarcel_compar$period %in% '2010-2018']='multiyear average 2010-18'
multiyear_color='blue'; kemri_indivyear='gray55'
#
ggplot(kemri_mcmarcel_compar,aes(x=age_inf,y=rate,color=interaction(datatype,timespan),group=interaction(datatype,period))) + 
  geom_line(aes(size=timespan,linetype=interaction(datatype,timespan))) + 
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=interaction(datatype,timespan)),alpha=0.3,colour=NA) + 
  scale_size_manual(values = c(2,1)) + scale_color_manual(values=c(multiyear_color,'red',kemri_indivyear)) + 
  scale_fill_manual(values=c(multiyear_color,'red',kemri_indivyear)) +
  scale_linetype_manual(values=c("solid",'solid',"blank")) +
  facet_wrap(~hospitalisation,nrow=2,scales='free') + 
  scale_x_continuous(labels=as.character(seq(0,age_maxval,5)),breaks=seq(0,age_maxval,5)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                axis.title=element_text(size=14), text=element_text(family="Calibri")) + xlab('age (months)') + ylab('incidence') +
  ggtitle('Kenya RSV incidence/100k: MCMARCEL prediction vs KEMRI data') + 
  labs(color='data source',fill='data source') + guides(size=FALSE,linetype=FALSE,color=FALSE)
# 
ggsave("output/kemri_mcmarcel_compar.png",width=24,height=18,units="cm") 
