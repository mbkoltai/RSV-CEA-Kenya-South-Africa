##############################################################################
# This file is part of the RSV modelling project McMarcel.
# 
# => MAIN SCRIPT TO RUN THE MODEL FOR 72 GAVI COUNTRIES
#
# Multi-Country Model Application for RSV Cost-Effectiveness poLicy (McMarcel) 
#
#  Copyright 2020, CHERMID, UNIVERSITY OF ANTWERP, BELGIUM
##############################################################################
# The objective of this modelling project is to evaluate the impact and cost-
# effectiveness of potential maternal and neonatal RSV immunisation strategies 
# in 72 Gavi countries. See our README file for more info.
#
# Citation: Li, Willem, Antillon, Bilcke, Jit, Beutels. Health and economic 
# burden of Respiratory Syncytial Virus (RSV) disease and the cost-
# effectiveness of potential interventions against RSV among children under 
# 5 years in 72 Gavi-eligible countries. BMC Medicine. (2020)
##############################################################################

# clear workspace
rm(list=ls())

## set working directory (or open RStudio with this script)
# setwd("C:/User/path/to/the/rcode/folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# path
library(rstudioapi); currentfile_path=rstudioapi::getActiveDocumentContext()$path
currentfile_path=paste0(unlist(strsplit(
  currentfile_path,"\\/"))[1:(length(unlist(strsplit(currentfile_path,"\\/")))-1)],collapse="/")
setwd(currentfile_path)

library(tidyverse); library(reshape2)

#######################
## SETTINGS          ##
#######################

# select the model configuration
# => set 'run_tag' to find the config file at ./config/<run_tag>.csv 
run_tag  <- 'RSV_gavi72_basecase'  # 72 Gavi countries (basecase) 
#run_tag <- 'RSV_gavi72_all'       # 72 Gavi countries (all scenarios) 
#run_tag <- 'RSV_gavi72_efficacy'  # 72 Gavi countries (severity-specific efficacy) 

# number of stochastic samples in the probabilistic sensitivity analysis (PSA)
num_sim <- 100

# random number generater seed
rng_seed <- gsub('-','',Sys.Date()) # 20190118

# option to create geographical and country-specific plots
# note: this might require substantial processing time
boolean_country_plots <- FALSE # TRUE
boolean_map_plots     <- FALSE # TRUE

#######################
## MODEL SETUP       ##
#######################

# (re)load packages and functions
source('functions/RSV_load_all.R')

# output directory postfix
output_dir_postfix <- paste0(run_tag,'_n',num_sim)

# add timestap to output directory name
output_dir <- paste0('output/',format(Sys.time(),'%m%d%H%M%S_'),output_dir_postfix)

# set seed
set.seed(rng_seed)

# config filename
config_filename <- paste0('./config/',run_tag,'.csv')

# log timings
time_stamp_main <- Sys.time()

# always clear temporary results
cli_print('Clear all temporary output')
unlink(file.path(get_temp_output_folder(output_dir)),recursive = T)

# start parallel workers
start_parallel_workers()

cli_print("****** START MC MARCEL ******")
cli_print("WORK DIR:",system('pwd',intern = T))
cli_print("OUTPUT DIR:",output_dir)

#######################
## LOAD CONFIG       ##
#######################

# load config file in csv format
sim_config_matrix <- read.table(config_filename,sep=',', dec='.',stringsAsFactors = F,header = T)
# this is a list of cntrs with the MV and mAb scenarios. it assumes MV protections lasts 5 months, mAb for 6.
# estimates for coverage: for kenya both are 86%

# set output file name prefix
sim_output_filename  <- file.path(output_dir,run_tag)

# add simulation details
sim_config_matrix$num_sim          <- num_sim
sim_config_matrix$scenario_id      <- 1:nrow(sim_config_matrix)
sim_config_matrix$rng_seed         <- rng_seed
sim_config_matrix$outputFileDir    <- get_output_folder(output_dir)

# Count number of scenarios
num_scen <- length(sim_config_matrix$scenario_id)

###############################################
## PRE-PROCESSING: country databases         ##
###############################################
cli_print('START PRE-PROCESSING',run_tag); time_stamp <- Sys.time()

# create UN country data
create_UN_country_database(output_dir)
# this loads UN world popul database
# data("UNlocations",package='wpp2017')
# creates a table 'UNcountries' with cntr names and iso3 codes
###
# list variables (but not functions): setdiff(ls(), lsf.str())

# pre-process WPP2017 data
load_wpp2017_databases(output_dir)
# loads sex ratio and mortality by age groups
# dataframes: sexratio, mxM, mxF
# mort_all = mxF*(1/(1+sexratio) + mxM*sexratio/(1+sexratio)
#######################################
## PRE-PROCESSING: life tables       ##
#######################################
## note: for a sequential run, replace %dopar% by %do%
cli_print('PRE-PROCESSING LIFE TABLES...'); time_stamp <- Sys.time()

# construct matrix with all unique [country, year] combinations
country_year_opt <- sim_config_matrix[,c('country_iso','year')]

# add [country, 2015] to derive the reference incidence, based on Shi et al (2017)
country_year_opt <- rbind(country_year_opt, cbind(country_iso=country_year_opt$country_iso, year=2015))
# get unique combinations: this creates a table with 72 cntrs with dates 2015 and 2020
country_year_opt <- unique(country_year_opt)

# make summary matrix, including the 5-year period notation
# 3 columns: countries, '2020-25', 2020
country_period_opt    <- cbind(country_year_opt$country_iso,
                                  t(sapply(country_year_opt$year,get_year_category)))
######################
# get life table for each [country, period] combination
# this needs to be run within 10 mins of starting parallel workers
par_out <- foreach(i_life=1:nrow(country_period_opt), .combine='rbind', .packages=all_packages,.verbose=FALSE) %dopar%
{  # print progress
  cli_progress(i_life,nrow(country_period_opt),time_stamp)
  # life table with discounting
  generate_life_table(country_period_opt[i_life,1],
                 country_period_opt[i_life,3],
                 output_dir,
                 0.03)
  # life table without discounting
  generate_life_table(country_period_opt[i_life,1],
                 country_period_opt[i_life,3],
                 output_dir,
                 0)
  # dummy return, the results are printed to a file
  return(0)
}

####
# generate_life_table creates a table with mortality, life expectancy and saves it to Rdata file
# the cntr specific infos here are: mortality, sex ratio
# it gives life expectancy and discounted life years remaining
i_life=which(country_period_opt[,1] %in% 'KEN')[1]
f_country_iso=country_period_opt[i_life,1]; f_year=country_period_opt[i_life,3]; 
f_outputFileDir=output_dir; f_disc_rate_qaly=0.03
# this saves into an Rdata file
generate_life_table(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
load(paste(f_outputFileDir,"temp/life_table_KEN_2020_2025_disc0p03.RData",sep='/'))
# for one country
life_table_tidy=melt(life_table,id.vars='age')
ggplot(life_table_tidy,aes(x=age,y=value)) + geom_line() + facet_wrap(~variable,scales='free')
ggsave("output/life_table_KEN_2020_2025_disc0p03.png",width=30,height=18,units="cm") # ,device=grDevices::pdf
# how foreach works
# df_foreach=foreach(i_life=1:5,.combine='rbind',.packages=all_packages,.verbose=FALSE) %do% {runif(5)}
#            foreach(i_life=1:5,.combine='rbind',.packages=all_packages,.verbose=FALSE) %dopar% {runif(5)}

#######################################
## PRE-PROCESSING: incidence         ##
#######################################
## note: for a sequential run, replace %dopar% by %do%
cli_print('PRE-PROCESSING INCIDENCE...' ); time_stamp <- Sys.time()

# get unique country codes
country_opt <- data.frame(country_iso = unique(sim_config_matrix$country_iso))

# preprocess incidence data for each country
# only for (eg) Kenya

par_out <- foreach(i_country=1:nrow(country_opt),.combine='rbind', .packages=all_packages, .verbose=FALSE) %dopar% 
{ # print progress
  cli_progress(i_country,nrow(country_opt),time_stamp)
  # get country-specific incidence data
  get_incidence(country_opt$country_iso[i_country],output_dir)
  # dummy return, the results are printed to a file
  return(0)
}

# for one country
ken_inds=which(country_opt$country_iso %in% 'KEN')
z=get_incidence(country_opt$country_iso[ken_inds],output_dir)
# this contains: RSV_rate [60x5000 dataframe], hosp_prob [60x5000 dataframe], hCFR_prob [60x1000 dataframe]
# incidence data is from
RSV_burden_Shi_2017=read_csv('input/RSV_burden_Shi_2017.csv')
# 'nrb episodes' is number of RSV episodes, incidence is incidence under 5yrs older/1000popul , bc eg for india
# incidence is 56.7, nrb_episodes=7.013.468, popul=(nrb_episodes/incidence)*1000=123.694.300, which is correct

# histogram of incidences. median=52.8. 
ggplot(RSV_burden_Shi_2017,aes(x=incidence_RSV_associated_ALRI)) + geom_histogram(binwidth=2)
# lineplot
id_vars=c("country_iso","location_name","lower_Incidence", "upper_Incidence", "lower_episodes","uper_episodes","reference")
RSV_burden_Shi_2017_tidy=melt(RSV_burden_Shi_2017,id.vars=id_vars)
RSV_burden_Shi_2017_tidy[grepl('incidence',RSV_burden_Shi_2017_tidy$variable),
                         grepl('episode',colnames(RSV_burden_Shi_2017_tidy))]=NA
RSV_burden_Shi_2017_tidy[grepl('nrb',RSV_burden_Shi_2017_tidy$variable),
                         grepl('Incidence',colnames(RSV_burden_Shi_2017_tidy))]=NA
RSV_burden_Shi_2017_tidy=RSV_burden_Shi_2017_tidy %>% 
  unite(lower_CI, c(lower_Incidence,lower_episodes), sep="",remove=TRUE,na.rm=TRUE) %>% 
  unite(upper_CI,c(upper_Incidence,uper_episodes),sep="",remove=TRUE,na.rm=TRUE)
RSV_burden_Shi_2017_tidy[,c('lower_CI','upper_CI')]=sapply(RSV_burden_Shi_2017_tidy[,c('lower_CI','upper_CI')],as.numeric)
write_csv(RSV_burden_Shi_2017_tidy,path='output/RSV_burden_Shi_2017_tidy.csv')
ggplot(RSV_burden_Shi_2017_tidy,aes(x=location_name,y=value,group=1)) + geom_line() + geom_point(size=1) +
  geom_ribbon(aes(ymin=lower_CI,ymax=upper_CI),alpha=0.3,colour=NA,fill="red") + facet_wrap(~variable,nrow=2,scales='free') +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5)) + scale_y_continuous(trans='log10')
ggsave("output/RSV_burden_Shi_2017.png",width=30,height=18,units="cm")
# intermediate steps: spline_datafiles contains the 3 input files
# ./input/incidence_lmic_ts_n5000.csv [300.000x4], ./input/hosp_prob_ts_n5000.csv [300.000x4], 
# ./input/cfr_lmic_ts_n5000.csv [60.000x4]

# it is the fcn 'convert_pred_into_model_input' that takes the input file
# 'incidence_lmic_ts_n5000.csv', 'hosp_prob_ts_n5000.csv', 'cfr_lmic_ts_n5000.csv'
# and generates 60x5000 matrix (cfr_mat 60x1000)

###
# how the total burden is distributed across age groups
samples_burden_distrib=melt(as.matrix(country_rsv_pred_cases[,1:25])); samples_burden_distrib$type='burden_abs'
samples_burden_permillion=melt(as.matrix(country_rsv_rate[,1:25])); samples_burden_permillion$type='burden_per_popul'
hosp_prob_mat_tidy=melt(as.matrix(hosp_prob_mat[,1:25])); hosp_prob_mat_tidy$type='hosp'
cfr_mat_tidy=melt(as.matrix(cfr_mat[,1:25])); cfr_mat_tidy$type='cfr'
samples_burden_all=rbind(samples_burden_permillion,samples_burden_distrib,hosp_prob_mat_tidy,cfr_mat_tidy); 
age_maxval=max(samples_burden_all$Var1)
ggplot(samples_burden_all,aes(x=Var1,y=value,group=Var2,color=Var2)) + geom_line() + 
  facet_wrap(~type,scales='free',nrow=2) + theme_bw() + 
  scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
  axis.text.x=element_text(size=11,angle=90,vjust=0.5),axis.text.y=element_text(size=11),
  axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  xlab('Age (month)') + ylab('Total burden') + labs(color="Samples") + ggtitle('Kenya burden random samples')
ggsave("output/RSV_burden_kenya_randomsamples.png",width=30,height=18,units="cm")

#####
# input file: incidence_lmic_ts_n5000
incidence_sample=incidence_lmic_ts_n5000[incidence_lmic_ts_n5000$iter<=25,]; age_maxval=max(incidence_sample$mos)
incidence_sample$type='incidence'
cfr_lmic_ts_n5000=as.data.frame(read_csv("./input/cfr_lmic_ts_n5000.csv")); cfr_sample=cfr_lmic_ts_n5000[cfr_lmic_ts_n5000$iter<=25,]
cfr_sample$type='cfr'; colnames(cfr_sample)[1]='X'; 
hosp_prob_ts_n5000=as.data.frame(read_csv("input/hosp_prob_ts_n5000.csv")); hosp_prob_sample=hosp_prob_ts_n5000[hosp_prob_ts_n5000$iter<=25,]
hosp_prob_sample$type='hosp'; colnames(hosp_prob_sample)[1]='X'; 
lmic_ts_all=rbind(incidence_sample,cfr_sample,hosp_prob_sample)
ggplot(lmic_ts_all,aes(x=mos,y=pred,group=iter,color=as.factor(iter))) + geom_line() + 
  facet_wrap(~type,scales='free',nrow=2) + 
  scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) + theme_bw() + 
  theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
        axis.text.x=element_text(size=11,angle=90,vjust=0.5),axis.text.y=element_text(size=11),
        axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  xlab('Age (month)') + ylab('Total burden') + labs(color="Samples") + ggtitle('Kenya incidence (random samples)')
ggsave("output/RSV_burden_kenya_randomsamples_prob_ts_n5000.png",width=30,height=18,units="cm")
#######################################
## PROCESSING: burden                ##
#######################################
## note: for a sequential run, replace %dopar% by %do%
cli_print('PROCESSING BURDEN [FOREACH]:',run_tag); time_stamp <- Sys.time()

# loop over each configuration
# start_parallel_workers()
ken_inds = which(sim_config_matrix$country_iso %in% 'KEN')
sim_output <- foreach(i_scen=ken_inds,.combine='rbind',.packages=all_packages,.verbose=FALSE) %dopar% # 1:num_scen
{  # print progress
  cli_progress(i_scen,num_scen,time_stamp)
  # run burden function
  run_output_long <- get_burden(sim_config_matrix[i_scen,])
  # write results to file
  save(run_output_long,file=file.path(get_temp_output_folder(
    sim_config_matrix$outputFileDir[i_scen],'burden'),
    paste0('run_output_long_',i_scen,'.Rdata')))
  # dummy return, the results are printed to a fle
  return(0)
}

###############################
## POST-PROCESSING           ##
###############################
cli_print('COLLECT BURDEN OUTPUT [FOREACH]:',run_tag); time_stamp <- Sys.time()

# check parallel workers
check_parallel_workers()

# loop over each scenario
sim_output <- foreach(i_scen=ken_inds, .combine = 'rbind', .verbose = FALSE) %dopar% # 1:num_scen
{ # print progress
  cli_progress(i_scen,num_scen,time_stamp)
  # get output name
  var_name <- load(file.path(get_temp_output_folder(sim_config_matrix$outputFileDir[i_scen],
                   'burden'),paste0('run_output_long_',i_scen,'.Rdata')))
  # load data and add scenario id
  run_output_long             <- get(var_name)
  run_output_long$scenario_id <- i_scen
  # return
  return(run_output_long)
}

# add config details to output
sim_output <- merge(sim_config_matrix,sim_output,all=T)

# sort output on scenario_id and save as RData file
sim_output <- sim_output[order(sim_output$scenario_id),]
save(sim_output,sim_output_filename,file=paste0(sim_output_filename,'.RData'))

#######################
## PLOT RESULTS      ##
#######################

# plot CEAF table overview
plot_CEAF_table(sim_output_filename)

# plot CEA results by country
if(boolean_country_plots) { plot_CEA_country_results(sim_output_filename) }

# get geographic figures (optional)
if(boolean_map_plots) { plot_maps(sim_output_filename) }

# get aggregated global and country tables
write_global_summary_tables(sim_output_filename)

# stop parallel workers
stop_parallel_workers()

###
# when output contains only 1 cntr (Kenya)
sim_output_kenya = sim_output[sim_output$country_iso %in% 'KEN',]
df_datatypes=as.data.frame(sapply(1:ncol(sim_output_kenya), function(x){length(unique(sim_output_kenya[,x]))}))
rownames(df_datatypes)=colnames(sim_output_kenya); colnames(df_datatypes)='unique_vals'
View(df_datatypes)

# histogram of RSV cases
ggplot(sim_output_kenya,aes(x=rsv_cases)) + geom_histogram(binwidth=20) + 
  geom_vline(aes(xintercept=mean(rsv_cases)),col='red',size=2)
# 3.8e5 cases according to RSV
mean(sim_output_kenya$rsv_cases)

# mean rate for <5 yr olds from our own data
N_under5_kenya=7e6
kenya_data_path='../path_rsv_data/SARI_Rates_2010_2018/SARI_Rates_2010_2018_tidydata_cleaned.csv'
SARI_Rates_2010_2018_tidydata=read_csv(kenya_data_path)
kenya_rsv_data_means=SARI_Rates_2010_2018_tidydata[SARI_Rates_2010_2018_tidydata$age_in_months %in% "<60" & 
                  SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & SARI_Rates_2010_2018_tidydata$RSV_association & 
                  SARI_Rates_2010_2018_tidydata$period %in% '2010-2018',]
# 6.6e4 cases according to our kenya data
sum(kenya_rsv_data_means$rate)*N_under5_kenya/1e5

####
library(outbreaks)
