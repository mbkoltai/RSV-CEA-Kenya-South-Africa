### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Script to reproduce analysis and figures in the article at [..TBC..]
rm(list=ls())
package_names <- c("tidyverse","rstudioapi","fitdistrplus","rstudioapi","matrixStats","ungeviz",
  "stringi","rriskDistributions","cowplot","here","conflicted")
lapply(package_names, function(x) if (!any(row.names(installed.packages()) %in% x)) {install.packages(x)})
lapply(package_names,library,character.only=TRUE)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
conflict_prefer("select", "dplyr"); conflict_prefer("filter", "dplyr")
num_sim <- 5000; source(here::here('functions/RSV_load_all.R'))
# LOAD FUNCTIONs required for analysis 
source(here::here('functions/load_config_pars.R'))
lapply(here::here(c("functions/set_xlims_cea.R","functions/get_burden_flexible.R",
         "functions/get_burden_flexible_ari_sari.R",'functions/GammaParmsFromQuantiles.R',
         "functions/load_own_data.R")),function(x) {source(x)})
# conflict_prefer("select","dplyr"); # conflict_prefer("filter","dplyr")
# load data
### Kenya incidence data -------------------------
# hosp rate (p): p/(1-p) ~ norm
kenya_data_file_path<-"custom_input/Kenya_ARI_SARI_Rates_2010_2018_tidydata_updated_2021_08.csv"
# "custom_input/ARI_SARI_Rates_2010_2018_tidydata.csv"
### PLOT Kenya incidence data with error bars
bind_rows(fcn_load_kenya(kenya_data_path=kenya_data_file_path,sel_disease="ARI")$rsv_incidence_ageinf,
          fcn_load_kenya(kenya_data_path=kenya_data_file_path,sel_disease="SARI")$rsv_incidence_ageinf) %>% 
  mutate(disease_type_medic_status=paste(disease_type, 
                                         ifelse(medically_attended,"medically attended","not attended")) ) %>%
  mutate(disease_type_medic_status=gsub("SARI medically attended","SARI hospitalised",
                                        disease_type_medic_status)) %>%
  mutate(disease_type_medic_status=gsub("SARI not attended","SARI non-hospitalised",
                                        disease_type_medic_status)) %>%
  dplyr::select(!c(RSV_assoc,freq,n)) %>% relocate(age_inf,.before=variable) %>% 
  relocate(metric_per_popul,.after=value) %>% 
  group_by(age_in_months,disease_type_medic_status) %>% summarise(value=unique(value),
    CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),
    metric_per_popul=unique(metric_per_popul)) %>%
# PLOT
ggplot(aes(x=age_in_months)) + 
  geom_bar(aes(y=value/metric_per_popul,fill=disease_type_medic_status),position="stack",stat="identity") +
  geom_errorbar(aes(ymin=CI_95_lower/metric_per_popul,ymax=CI_95_upper/metric_per_popul),size=0.4) +
  facet_wrap(~disease_type_medic_status,nrow=2,scales = "free") + # 
  xlab("age (months)") + ylab("cases per person year") + labs(fill="")+
  scale_y_continuous(expand=expansion(0.01,0)) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=13),
                                      axis.text.y=element_text(size=13),legend.text=element_text(size=13),
                                      legend.background=element_rect(fill=NA),legend.position=c(0.92,0.925),
                                      axis.title.x=element_text(size=17),axis.title.y=element_text(size=17),
                                      strip.text=element_text(size=14))
# geom_text(data=hosp_rate_kenya,aes(x=45,y=c(0.17,0.2),
# label=paste0(status,"=",round(mean_medic_attended*1e2,1),"%")),size=6) +
# scale_x_continuous(breaks=(0:30)*2,expand=expansion(0.01,0))
# SAVE
# ggsave("output/ari_sari_burden/kenya_ari_sari_burden_errorbars_grouped_updated_jul2021.png",
#   width=35,height=22,units="cm")
ggsave("output/ari_sari_burden/kenya_ari_sari_burden_errorbars_grouped_updated_2021_08.png",
       width=42,height=22,units="cm")

# LOAD data, fit (gamma) distrib to CI95 values, generate matrix with 5e3 columns, age groups from 0 to 59mts
# Kenya
ci50_range <- c(25,75)/1e2; ci95_range <- c(2.5,97.5)/1e2
kenya_nonhosp_hosp_incid_ari_sari <- lapply(c("ARI","SARI"), function(x)
  fcn_gen_nonhosp_hosp_incid_samples_kenya(kenya_data_file_path,sel_disease=x,n_iter=5e3,age_maxval=60,
            CI_intervals=ci95_range,randsampl_distrib_type="gamma"))
names(kenya_nonhosp_hosp_incid_ari_sari)=c("ARI","SARI")
### deaths
# Kenya deaths
deaths_kenya <- read_csv("custom_input/deaths_kenya_tidy.csv") %>% 
  filter(variable=="rate" & !age_in_months %in% c("<12","12-23","<24","24-59","<60")) %>% 
  mutate(age_in_months=ifelse(age_in_months=="<1","0",age_in_months),freq=1) %>% 
  mutate(freq=ifelse(grepl('-',age_in_months),as.numeric(sapply(age_in_months, 
                              function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,freq)) %>% 
  mutate(age_in_months_orig=factor(age_in_months,levels=unique(age_in_months)),
         age_in_months=ifelse(grepl('-',age_in_months), 
                              sapply(strsplit(age_in_months,'-'),'[[',1),age_in_months)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% mutate(age_inf=as.numeric(age_in_months)+(n-1)) %>% 
  dplyr::select(!c(n,freq,age_in_months)) %>% relocate(age_inf,.before=value)
# fit distributions to CI95
kenya_deaths_distrib_params <- bind_rows(lapply(c("yes","no"), function(y_no) data.frame(age_inf=0:59, 
  t(sapply(1:(nrow(deaths_kenya)/2), function(x) gamma.parms.from.quantiles(p=c(2.5,97.5)/100, 
  q=as.numeric((deaths_kenya %>% mutate(CI_95_lower=ifelse(CI_95_lower==0,0.1,CI_95_lower)) %>% 
                  filter(in_hospital==y_no) %>% 
  dplyr::select(c(CI_95_lower,CI_95_upper)))[x,]))[c("shape","rate")])),in_hospital=y_no))) %>% 
  mutate(shape=unlist(shape),rate=unlist(rate))
# generate samples from fitted distribs
kenya_deaths_incid <- lapply(c("yes","no"), function(y_no) t(sapply(0:59, function(x) 
    rgamma(5e3,shape=(kenya_deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$shape,
          rate=(kenya_deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$rate)))/1e5)
names(kenya_deaths_incid)=c("hosp","nonhosp")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# SA deaths
deaths_SA <- read_csv("custom_input/mortality South Africa updated17_01_2022_TIDY.csv") %>% 
  rename(age_in_months=`age in months`,CI_95_lower=LCI,CI_95_upper=UCI,in_hospital=`medically attended`) %>% 
  filter(!in_hospital %in% "total") %>%
  mutate(in_hospital=ifelse(grepl("in-hospital",in_hospital),"yes","no")) %>%
  filter(variable %in% "rate" & !age_in_months %in% c("<1 year","<5 years")) %>% 
  mutate(age_in_months=ifelse(age_in_months=="<1","0",age_in_months),freq=1) %>% 
  mutate(freq=ifelse(grepl('-',age_in_months),as.numeric(sapply(age_in_months, 
                                      function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,freq)) %>% # 
  mutate(age_in_months_orig=factor(age_in_months,unique(age_in_months)),
         age_in_months=ifelse(grepl('-',age_in_months), 
                              sapply(strsplit(age_in_months,'-'),'[[',1),age_in_months)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% mutate(age_inf=as.numeric(age_in_months)+(n-1)) %>% 
  dplyr::select(!c(n,freq,age_in_months)) %>% relocate(age_inf,.before=value)
# plot
# ggplot(deaths_SA,aes(x=age_inf,y=value,color=in_hospital)) + geom_line() + geom_point() +
#   geom_ribbon(aes(ymin=CI_95_lower,ymax=CI_95_upper,fill=in_hospital),alpha=0.2) + theme_bw() +
#   scale_x_continuous(breaks=(0:30)*2,expand=expansion(0.01,0)) + scale_y_log10(breaks=2^(-3:7))
deaths_data <- bind_rows(deaths_SA %>% select(!age_inf) %>% group_by(age_in_months_orig,in_hospital) %>%
            summarise(value=unique(value),CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),
            country="South Africa") %>% ungroup() %>% 
            mutate(age_in_months_orig=as.character(age_in_months_orig)),
        deaths_kenya %>% group_by(age_in_months_orig,in_hospital) %>%
            summarise(value=unique(value),CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),
            country="Kenya")) %>%
        mutate(age_in_months_orig=factor(age_in_months_orig,levels=unique(age_in_months_orig))) %>% 
        group_by(country,age_in_months_orig) %>% 
        mutate(CI_95_lower_sum=sum(CI_95_lower),CI_95_upper_sum=sum(CI_95_upper)) %>% 
  mutate(in_hospital=ifelse(in_hospital=="yes","in-hospital","out-of-hospital"))
# 
p <- ggplot(deaths_data,aes(x=age_in_months_orig)) + 
  geom_bar(aes(y=value,fill=in_hospital),stat="identity",position="dodge") +
  geom_errorbar(aes(ymin=CI_95_lower,ymax=CI_95_upper,group=in_hospital),size=0.2,position="dodge") +
  facet_wrap(~country,scales="free_y",nrow=2) + # 
  xlab("age (months)") + ylab("deaths per 100,000 person year") + labs(fill="") +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:12)*25) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(vjust=0.5,size=13),axis.text.y=element_text(size=13),
    legend.text=element_text(size=13),legend.background=element_rect(fill=NA),strip.text=element_text(size=14),
    legend.position=c(0.88,0.925),axis.title.x=element_text(size=17),axis.title.y=element_text(size=17)); p
# save
ggsave("output/cea_plots/ALL_deaths_data_dodged_2rows.png",width=32,height=18,units="cm") #  # _yfixed
# fit distributions to CI95
SA_deaths_distrib_params <- bind_rows(lapply(c("yes","no"), 
        function(y_no) data.frame(age_inf=0:(nrow(deaths_SA)/2-1), 
        t(sapply(1:(nrow(deaths_SA)/2), function(x) gamma.parms.from.quantiles(p=c(2.5,97.5)/100, 
        q=as.numeric((deaths_SA %>% mutate(CI_95_lower=ifelse(CI_95_lower==0,0.1,CI_95_lower)) %>% 
        filter(in_hospital==y_no) %>% 
        dplyr::select(c(CI_95_lower,CI_95_upper)))[x,]))[c("shape","rate")])),in_hospital=y_no))) %>% 
  mutate(shape=unlist(shape),rate=unlist(rate))
# generate samples from fitted distribs
popul_denom<-unique(deaths_SA$metric_per_popul)
SA_deaths_incid <- lapply(c("yes","no"), function(y_no) t(sapply(0:(nrow(deaths_SA)/2-1), function(x) 
  rgamma(5e3,shape=(SA_deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$shape,
         rate=(SA_deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$rate)))/popul_denom)
names(SA_deaths_incid)=c("hosp","nonhosp")
### ### ### ### ### ### ### ### ### ### ### ### ### ###
# check if fits match the data
# ggplot(bind_rows(data.frame(age_inf=0:59,value=rowMeans(SA_deaths_incid$hosp),
#             t(sapply(1:nrow(SA_deaths_incid$hosp),
#     function(x) quantile(SA_deaths_incid$hosp[x,],probs=c(2.5,97.5)/100))),source="fit" ) %>%
#       rename(CI_95_lower=`X2.5.`,CI_95_upper=`X97.5.`),
#   deaths_SA %>% filter(in_hospital %in% "yes") %>%
#       select(age_inf,value,CI_95_lower,CI_95_upper) %>%
#     mutate(source="data",value=value/1e5,CI_95_lower=CI_95_lower/1e5,CI_95_upper=CI_95_upper/1e5) ),
#       aes(x=age_inf,color=source,fill=source)) + geom_point(aes(y=value)) + geom_line(aes(y=value)) +
#   geom_ribbon(aes(ymin=CI_95_lower,ymax=CI_95_upper),alpha=1/5) + scale_y_log10() +
#   theme_bw() + standard_theme

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SA
SA_SARI_data <- fcn_load_s_afr(safr_data_path = "custom_input/s_afr_incidence_data_rate.csv") %>%
  mutate(disease_type_medic_status=paste(disease_type,
        ifelse(hospitalisation,"hospitalised","not hospitalised")) ) %>%
  rename(agegroup_mts=age) %>% 
  relocate(age_inf,.before=Province) %>% relocate(Province,.after=disease_type_medic_status) %>% 
  relocate(year,.before=Province)
# SA ILI data
SA_ILI_data <- read_csv("custom_input/s_afr_ILI_incidence_rate.csv") %>% 
  filter(!(grepl("<",agegroup) | agegroup %in% c("0-5m","6-11m","12-23m","24-59m","<5y"))) %>%
  mutate(agegroup_mts=agegroup,agegroup=gsub("m","",agegroup), freq=1) %>% 
  mutate(freq=ifelse(grepl('-',agegroup),as.numeric(sapply(agegroup,
                function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,freq)) %>% 
  mutate(agegroup=ifelse(grepl('-',agegroup), sapply(strsplit(agegroup,'-'),'[[',1),agegroup)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% 
  mutate(age_inf=as.numeric(agegroup)+(n-1)) %>% dplyr::select(!c(n,freq,agegroup)) %>% # 
  relocate(age_inf,.before=rate) %>% relocate(disease_type,.before=hospitalisation) %>%
  mutate(disease_type=ifelse(disease_type=="ILI","ARI",""),
         disease_type_medic_status=ifelse(hospitalisation,paste0("medically attended ",disease_type),
                paste0("non medically attended ",disease_type))) %>% 
  relocate(disease_type_medic_status,.before=Province)
# concatenate
if (!exists("SA_data")){ SA_data=bind_rows(SA_ILI_data,SA_SARI_data) }
## estimate from Kenya on % of ARI cases with fever (=ILI) -> take this percentage to expand ILI to ARI
ILI_adjust_SA=TRUE
if (ILI_adjust_SA) { # 
  print("divide by fever proportion")
  SA_data[SA_data$disease_type=="ARI",c("rate","rate_CI_lower","rate_CI_upper")]=
  SA_data[SA_data$disease_type=="ARI",c("rate","rate_CI_lower","rate_CI_upper")]/mean(c(0.333,0.205,0.16))
  divided_fever=TRUE }
# CI lower limit should not be 0 -> setting it to 1 
# (since average value of mean rate is > 1000, this does not change results more than 0.1%)
SA_data$rate_CI_lower[SA_data$age_inf==0 & SA_data$disease_type=="ARI"]=1
### generate 5e3 sample paths for CEA
sa_nonhosp_hosp_incid_ari_sari=lapply(c("ARI","SARI"), 
      function(x) fcn_gen_nonhosp_hosp_incid_samples_SA(SA_data,diseasetype=x,
      n_iter=5e3, age_maxval=60,CI_intervals=ci95_range,randsampl_distrib_type="gamma"))
names(sa_nonhosp_hosp_incid_ari_sari)=c("ARI","SARI")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# South Africa cost estimates
# inpatient
s_afr_inpatient_cost <- read_csv("custom_input/s_afr_PDE_calcs.csv") %>% mutate(freq=ifelse(grepl('-',age),
  as.numeric(sapply(age, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,1)) %>% 
  mutate(age=ifelse(grepl('-',age), sapply(strsplit(age,'-'),'[[',1),age)) %>% 
  uncount(weights=freq, .id="n",.remove=F) %>%
  mutate(age=as.numeric(age)+(n-1)) %>% dplyr::select(!c(n,freq)) %>% rename(mean=`Mean cost per illness episode (USD)`)
  
# outpatient
s_afr_outpatient_cost <- cbind(data.frame(age="all",mean=25,LCI=18.3,UCI=31.8),
        data.frame(t(unlist(gamma.parms.from.quantiles(q=c(18.3,31.8),p=c(2.5,97.5)/100)[c("shape","rate")]))) )
if (!any(grepl("shape",colnames(s_afr_inpatient_cost)))){
 s_afr_inpatient_cost=cbind(s_afr_inpatient_cost, t(sapply(1:nrow(s_afr_inpatient_cost), function(x) 
  unlist(gamma.parms.from.quantiles(q=c(s_afr_inpatient_cost$LCI[x],s_afr_inpatient_cost$UCI[x]),
                                    p=ci95_range)[c("shape","rate")]))) ) }
list_SA_costs <- list("inpatient"=s_afr_inpatient_cost,"outpatient"=s_afr_outpatient_cost)
### ### ### ### ### ### ### ### ### ### ### ### ###
# KENYA costs
kenya_costs <- read_csv("custom_input/kenya_costing_tables_tidy.csv")
# using inpatient/outpatient ratio in South Africa
SA_total_av_cost<-median(s_afr_inpatient_cost$mean[s_afr_inpatient_cost$name %in% "total"])+s_afr_outpatient_cost$mean
SA_inpatient_cost_share <- (SA_total_av_cost-s_afr_outpatient_cost$mean)/SA_total_av_cost; rm(SA_total_av_cost)
# assemble list
inpat_rows <- kenya_costs %>% filter(grepl("Siaya",site) & grepl("Total patient",variable))
list_KEN_costs <- list(inpatient_household=bind_cols(
  kenya_costs %>% filter(grepl("Siaya",site) & grepl("Total patient",variable)), 
  t(sapply(1:nrow(inpat_rows), function(x)
  unlist(get.gamma.par(q=c(inpat_rows$ci95_low[x],inpat_rows$median[x],inpat_rows$ci95_up[x])/100,
                                    p=c(2.5,50,97.5)/100,plot=F)[c("shape","rate")]))),scaling=100),
  inpatient_healthcare_system=c(mean=SA_inpatient_cost_share*(kenya_costs %>% 
          filter(grepl("Siaya",site) & grepl("LRTI",variable)))$mean),
  outpatient_cost=c(mean=(1-SA_inpatient_cost_share)*(kenya_costs %>% 
          filter(grepl("Siaya",site) & grepl("LRTI",variable)))$mean))
# gamma can fit median and CIs well, but not so much the mean
# get.gamma.par(q=c(inpat_rows$ci95_low[1],inpat_rows$median[1],inpat_rows$ci95_up[1])/100,
#              p=c(2.5,50,97.5)/100,plot=F)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cntrs_cea=c("KEN","ZAF")
# efficacy figures for vaccine for RESVAX (Novavax trial)
# for mAb from NIRSEVIMAB (from https://www.nejm.org/doi/full/10.1056/nejmoa1913556)
# if this flag is set to TRUE, then using published efficacy data. if FALSE --> interim results
flag_publ_effic <- T
if (flag_publ_effic){
efficacy_figures <- list(mat_vacc=list(sympt_disease=c(mean=0.394,CI95_low=0.053,CI95_high=0.612),
                            hospit=c(mean=0.444,CI95_low=0.196,CI95_high=0.615),
                            severe=c(mean=0.483,CI95_low=-8.2/100,CI95_high=0.753), # c(mean=0.483,-8.2/100,0.753)
                            half_life=36.5/30,duration=3),
                      monocl_ab=list(sympt_disease=c(mean=0.701,CI95_low=0.523,CI95_high=0.812),
                                     hospit=c(mean=0.784,CI95_low=0.519,CI95_high=0.903),
                                     half_life=59.3/30,duration=5)) } else {
  # new data (from conference)
efficacy_figures <- list(mat_vacc=list(sympt_disease=c(mean=0.847,CI95_low=0.216,CI95_high=0.976),
                              hospit=c(mean=0.847,CI95_low=0.216,CI95_high=0.976), # no sep data
                              severe=c(mean=0.915,CI95_low=-5.6/100,CI95_high=0.998),
                              half_life=36.5/30,duration=3),
                         monocl_ab=list(sympt_disease=c(mean=0.745,CI95_low=0.496,CI95_high=0.871),
                              hospit=c(mean=0.621,CI95_low=-8.6/100,CI95_high=0.868),
                                       half_life=59.3/30,duration=5)) 
}

# fitting efficacy figures with a beta distribution
source("functions/fit_efficacy.R")
g(list_effic_betafit,allfits) %=% fcn_betafit_efficacy(effic_figs=efficacy_figures,
                                   scan_range_resol_nsample=c(min=-2,max=2,by=1/100,n_sample=2e4),
                                   optim_range_res=c(min=-1,max=2,by=0.04),optim_initguess=c(-0.05,1))
# fit as: beta_fit <- rbeta(n=1e4,shape1=alphaval,shape2=betaval)*scale_val + shift_val

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# parameters for exponential waning model
g(list_exp_waning_param,df_exp_waning_param) %=% fcn_exp_waning_rate(efficacy_figures,n_row=60)

# prices for doses
pricelist=list("mat_vacc"=c(3,10,30),"mAb"=c(6,20,60))
# loop thru: cntrs * interventions * dose prices
# n_cntr_output=1:length(cntrs_cea); n_interv=1:2
par_table <- expand_grid(n_cntr_output=1:length(cntrs_cea),n_interv=1:2); read_calc_flag=c("calc","read")[1]
kenya_deaths_input=TRUE; SA_deaths_input=TRUE
# exponential waning model used for efficacy
exp_wane_val=TRUE
# distribution used to fit efficacy figures
effic_dist_fit <- "beta"
# lower_cov=FALSE; lower_cov_val=0.7; 
subfolder_name <- paste0("new_price_efficacy_",ifelse(kenya_deaths_input,"KENdeaths",""),
                       ifelse(SA_deaths_input,"_SAdeaths",""),
        "_CIs_SA_ILI_",ifelse(ILI_adjust_SA,"broader","narrow"),ifelse(exp_wane_val,"_expwaning",""),
        ifelse(lower_cov,paste0("_coverage",lower_cov_val),""),
        ifelse(grepl("gamma",effic_dist_fit),"","_effic_betafit"),
        ifelse(flag_publ_effic,"","_interim"),"/") 
# ifelse(min(unlist(efficacy_figures))<=0,"_nonposit_effic","")
### before starting loop need to create temp folder
source("init_cea_calc_parallel.R")
###
# outputs to display
all_cols=c("non_hosp_cases","hosp_cases",
           "rsv_deaths","rsv_deaths_disc",
           "total_DALY","total_DALY_disc",
           "cost_rsv_hosp",
           "total_medical_cost","total_medical_cost_disc",
           "incremental_cost","incremental_cost_disc",
           "total_medical_cost_averted","total_medical_cost_disc_averted",
           "non_hosp_cases_averted","hosp_cases_averted",
           "rsv_deaths_averted", 
           "total_DALY_averted","total_DALY_disc_averted",
           "SARI_averted",
           "incremental_cost/DALY_averted",
           "total_YLD","total_YLD_disc",
           "total_YLL","total_YLL_disc",
           "total_YLD_disc_averted","total_YLL_disc_averted",
           "hosp_SARI","non_hosp_SARI","ARI_averted") 
####
burden_cols <- all_cols[!grepl("cost|averted",all_cols)]; cost_cols <- all_cols[grepl("cost",all_cols)]
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# cl=parallel::makeCluster(8); registerDoParallel(cl)
# foreach (k_par=1:nrow(par_table),.packages=c("dplyr","ggplot2","tidyr","readr","rriskDistributions")) %dopar% 
# parallelisation might crash R
ci50_range <- c(25,75)/1e2; ci95_range <- c(2.5,97.5)/1e2

for (k_par in 1:nrow(par_table)) {
    n_cntr_output <- par_table$n_cntr_output[k_par]; n_interv <- par_table$n_interv[k_par]
    # intervention config table
    sel_interv <- sim_config_matrix[which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])[n_interv],]
    if (cntrs_cea[n_cntr_output]=="ZAF"){
      sel_interv$country_iso=cntrs_cea[n_cntr_output]}
    # modify duration
    if (n_interv==1) { sel_interv$dur_protection_maternal=3/12 
    } else { sel_interv$dur_protection_infant<-5/12 }
    # half-life
    half_life <- ifelse(n_interv==1,36.5,59.3)
    # lower coverage
    if (lower_cov) {if (n_interv==1) { sel_interv$coverage_maternal <- lower_cov_val } else {
      sel_interv$coverage_infant=lower_cov_val} }
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    # calculate
    if (grepl("calc",read_calc_flag)) {
      # loop through PRICE levels
      for (k_price in 1:length(pricelist[[n_interv]])) {
        doseprice=c("mat_vacc"=ifelse(n_interv==1,pricelist$mat_vacc[k_price],pricelist$mat_vacc[1]),
                    "mAb"=ifelse(n_interv==2,pricelist$mAb[k_price],pricelist$mAb[1]))
        # calculation with data from mcmarcel (community-based)
        sim_output=get_burden_flexible(sel_interv,NA,NA,exp_wane=exp_wane_val,doseprice)
        # SARI, ARI distinguished, hosp/nonhosp distinguished
        if (cntrs_cea[n_cntr_output]=="ZAF") {
          cost_input <- list("inpatient"=list_SA_costs$inpatient %>% filter(name %in% "total"),
                             "outpatient"=list_SA_costs$outpatient) } else {
            cost_input<-list_KEN_costs }
        # input death data if available
      if (kenya_deaths_input) { print("using propr death data for Kenya")
        kenya_nonhosp_hosp_incid_ari_sari$deaths$hosp=kenya_deaths_incid$hosp
        kenya_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=kenya_deaths_incid$nonhosp } else {
          kenya_nonhosp_hosp_incid_ari_sari$deaths=NULL}
        if (SA_deaths_input) { print("using proprietary death data for SA")
          sa_nonhosp_hosp_incid_ari_sari$deaths$hosp=SA_deaths_incid$hosp
          sa_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=SA_deaths_incid$nonhosp } else {
            sa_nonhosp_hosp_incid_ari_sari$deaths=NULL}
    # RUN CALCULATIONS
        sim_output_user_input_ari_sari <- get_burden_flexible_ari_sari(sel_interv,
          list(kenya_nonhosp_hosp_incid_ari_sari,sa_nonhosp_hosp_incid_ari_sari)[[n_cntr_output]],
          efficacy_figures,effic_prob=T,effic_distr=effic_dist_fit,list_effic_fit=list_effic_betafit,
          exp_wane=exp_wane_val,list_exp_waning_param,dose_price=doseprice,cost_data=cost_input)
  ### processing samples -> mean, median, CIs
  # ICER with discounted DALYs
    with_discounted_DALYs <- fcn_process_burden_output(
    user_output=sim_output_user_input_ari_sari,default_output=sim_output,
    sel_cntr=sel_interv$country_iso,cols_burden_sel="",
    plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",
                  own="new data (hospital-based + HUS)"),
    icercolname="incremental_cost/DALY_averted",icercols=c("incremental_cost","total_DALY_disc_averted")) %>%
    filter(variable %in% "incremental_cost/DALY_averted") %>% 
      mutate(variable="incremental_cost/DALY_disc_averted")
  # every other outcome (appending `with_discounted_DALYs` to it)
    burden_mcmarcel_owndata <- bind_rows(fcn_process_burden_output(
    user_output=sim_output_user_input_ari_sari,default_output=sim_output,
    sel_cntr=sel_interv$country_iso,cols_burden_sel="",
    plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",
                  own="new data (hospital-based + HUS)"),
    icercolname="incremental_cost/DALY_averted",
    icercols=c("incremental_cost","total_DALY_averted")),
    with_discounted_DALYs) %>%
      mutate(name_root=gsub("_averted","",variable)) %>% group_by(source,iter,name_root) %>% 
      mutate(value_norm=ifelse(grepl("avert",variable)&!grepl("incremental",variable),
            value[grepl("avert",variable)]/value[!grepl("avert",variable)],NA)) %>% 
      ungroup() %>% dplyr::select(!name_root)
    
    # calculate mean, median, CI95
    x=burden_mcmarcel_owndata %>% group_by(source,variable) %>% 
      summarise(mean=mean(value,na.rm=T),median=median(value,na.rm=T),
          CI50_low=quantile(value,probs=ci50_range,na.rm=T)[1],
          CI50_high=quantile(value,probs=ci50_range,na.rm=T)[2],
          CI95_low=quantile(value,probs=ci95_range,na.rm=T)[1],
          CI95_high=quantile(value,probs=ci95_range,na.rm=T)[2],
          norm_mean=mean(value_norm,na.rm=T),norm_median=median(value_norm,na.rm=T),
          norm_CI50_low=quantile(value_norm,probs=ci50_range,na.rm=T)[1],
          norm_CI50_high=quantile(value_norm,probs=ci50_range,na.rm=T)[2],
          norm_CI95_low=quantile(value_norm,probs=ci95_range,na.rm=T)[1],
          norm_CI95_high=quantile(value_norm,probs=ci95_range,na.rm=T)[2]) %>%
          mutate(price=doseprice[n_interv],source_num=c(0.2,2)[as.numeric(source)],
                 country_iso=sel_interv$country_iso,intervention=sel_interv$intervention)
    # collect all outputs
        if (k_price==1) {cea_summary=x} else {cea_summary=bind_rows(cea_summary,x)}  } # END price scan
      # save CEA summary table
      folder_path=paste0("output/cea_plots/",subfolder_name)
      if (!dir.exists(folder_path)) {dir.create(folder_path)}
      write_csv(cea_summary,paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,
            "_cea_summary_mean_CI_",gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv") ) } 
    # LOAD EXISTING DATA IF already available
    else { cea_summary <- read_csv(paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,
                  "_cea_summary_mean_CI_",gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv")) }
} # loop country
# stop cluster
# stopCluster(cl)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# save in one file
filenames <- list.files(paste0("output/cea_plots/",subfolder_name),pattern=".csv")[
  !(list.files(paste0("output/cea_plots/",subfolder_name),pattern=".csv") %in% "cea_summary_all.csv")]
filenames <- filenames[grepl("_cea_summary_mean_",filenames)]
for (k_filename in 1:length(filenames)) {
  x=read_csv(paste0("output/cea_plots/",subfolder_name,filenames[k_filename]))
  if (k_filename==1){ cea_summary_all=x } else {cea_summary_all=bind_rows(cea_summary_all,x)} 
}
# save
write_csv(cea_summary_all,paste0("output/cea_plots/",subfolder_name,"cea_summary_all.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Make summary plots (after running calcs by "run_cea_calc_parallel.R" above)

### READ IN results of simulations if already available 
# cea_summary_all <- read_csv(paste0("output/cea_plots/",subfolder_name,"cea_summary_all.csv"))
# composition of total burden
# total_DALY <- total_YLD + total_YLL
# total_YLL <- rsv_deaths*config$hosp_CFR_DALYloss
# total_YLD <- hosp_med_att_YLD + non_hosp_YLD
# hosp_med_att_YLD <- hosp_SARI * config$severe_rsv_DALYloss + med_att_ARI*config$non_severe_rsv_DALYloss
# hosp_YLD <- hosp_SARI* config$severe_rsv_DALYloss
# non_hosp_YLD <- non_hosp_SARI*config$severe_rsv_DALYloss + non_med_att_ARI*config$non_severe_rsv_DALYloss
#
# replace names
old_new_names <- list(
  "old"=list(c("total_YLD","total_YLL","hosp_YLD","hosp_med_att_YLD","non_hosp_YLD","total_DALY"),
            c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),  # ,"ARI_YLD","SARI_YLD"
            c("admin_cost","cost_rsv_hosp","cost_rsv_outpatient","total_medical_cost")),
  "new"=list(c("total YLD","total YLL","YLD hospitalised cases","YLD medically attended cases",
                "YLD non-hospitalised cases","total DALY"),
              c("deaths","hospitalised SARI","non-hospitalised SARI","hospitalised cases","non-hospitalised cases"),
              c("admin. costs","hospitalisation costs","outpatient costs","total medical cost")))
# add column of country
cea_summary_all <- cea_summary_all %>% mutate(plot_variable=NA,
                                              country_plot=ifelse(country_iso=="KEN","Kenya","South Africa"))

#######################
plot_list<-list()
for (k_name_categ in 1:length(old_new_names$old)){
  for (k_name in 1:length(old_new_names$old[[k_name_categ]])){
  cea_summary_all <- cea_summary_all %>% 
    mutate(plot_variable=ifelse(variable %in% old_new_names$old[[k_name_categ]][k_name],
                              old_new_names$new[[k_name_categ]][k_name],plot_variable)) } }

save_flag=FALSE # TRUE
for (k_plot in 1:3) {
  sel_vars <- old_new_names$old[[k_plot]]
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) & # intervention=="mAb" & 
              ((price==3&intervention=="maternal")|(price==6&intervention=="mAb")) & grepl("new",source)) %>%
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
  if (k_plot==1) {df_combined_fig3_4_5 <- df_plot} else {
  df_combined_fig3_4_5 <- bind_rows(df_combined_fig3_4_5,df_plot %>% mutate(across(where(is.numeric),round,3))) }
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
                  position=position_dodge(width=dodge_val),size=5.5) + 
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
       width=36,height=18,units="cm")} 
} ##### end of loop
# combine figures
plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],nrow=3,rel_heights=c(1.5,1,1.1),
          labels="auto",label_size=19)
ggsave(paste0("output/cea_plots/",subfolder_name,"combined_fig3_4_5.png"),width=35,height=40,units="cm")
# save table
write_csv(df_combined_fig3_4_5,paste0("output/cea_plots/",subfolder_name,"combined_fig3_4_5.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables (discounted DALYs or not)

for (k_daly in c("incremental_cost/DALY_averted","incremental_cost/DALY_disc_averted")) {
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c("intervention_cost","incremental_cost",k_daly) & grepl("new",source)) %>%
        mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
        vec=as.character(abs(round((10^floor(log10(abs(median))))*
                                     round(median/(10^floor(log10(abs(median)))),3)))),
    price_interv=factor(paste0(price,"$ (",intervention,")"),
                            levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")"))),
    variable=factor(gsub("disc","(disc)",gsub("_"," ",variable)),levels=c("intervention cost","incremental cost",
                                 gsub("disc","(disc)",gsub("_"," ",k_daly))))) %>%
  select(!(contains("norm")|plot_variable)) %>%
  mutate(cnt_int=paste0(country_iso,"\n$",price,"\n",intervention), 
         mean=round(ifelse(!grepl("DALY",variable),mean/1e6,mean),2),
         median=round(ifelse(!grepl("DALY",variable),median/1e6,median),2),
         CI50_low=round(ifelse(!grepl("DALY",variable),CI50_low/1e6,CI50_low),2),
         CI50_high=round(ifelse(!grepl("DALY",variable),CI50_high/1e6,CI50_high),2),
         CI95_low=round(ifelse(!grepl("DALY",variable),CI95_low/1e6,CI95_low),2),
         CI95_high=round(ifelse(!grepl("DALY",variable),CI95_high/1e6,CI95_high),2),
         variable=ifelse(!grepl("DALY",variable),
                         paste0(variable," (million US$)"),paste0(variable, " (US$)"))) %>%
  mutate(cnt_int=factor(cnt_int,levels=unique(cnt_int[order(country_iso,intervention,price)]),ordered=T),
         variable=factor(variable,levels=c("intervention cost (million US$)",
                      "incremental cost (million US$)",
                      ifelse(grepl("disc",k_daly),"incremental cost/DALY (disc) averted (US$)",
                             "incremental cost/DALY averted (US$)"))),
        CI95_low=ifelse(grepl("incremental cost/DALY",variable) & (CI95_high>5e3|abs(CI95_low)>2e3),NA,CI95_low),
        CI95_high=ifelse(grepl("incremental cost/DALY",variable) & (CI95_high>5e3|abs(CI95_low)>2e3),NA,CI95_high))

# save table
if (k_daly %in% "incremental_cost/DALY_averted"){ df_interv_incremcosts_icer <- df_plot } else {
df_interv_incremcosts_icer <- bind_rows(df_plot,df_interv_incremcosts_icer) %>% 
  mutate(across(where(is.numeric),round,3)) }
# plot
ylab_txt<-"cost in USD (median, CI95)"
ggplot(df_plot) + geom_boxplot(aes(x=cnt_int,middle=median, color=price_interv,
                      lower=CI50_low,upper=CI50_high,ymin=CI95_low,ymax=CI95_high),
                      position=position_dodge(width=dodge_val),stat="identity",width=0.85,size=1.1) +
  facet_wrap(~variable,scales="free_y",nrow=3) +
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),
                              colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=c(3.5,6.5,9.5),linetype="dashed",size=0.3) + 
  geom_vline(xintercept=c(6.5),size=1/2) + geom_hline(yintercept=0,linetype="dashed",size=1/2) +
  xlab("") + ylab(ylab_txt) + labs(color="",linetype="") +
  scale_x_discrete(expand=expansion(0.05,0)) + guides(color=guide_legend(ncol=2)) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=18),
                                      axis.text.y=element_text(size=19),strip.text=element_text(size=14),
                                      legend.text=element_text(size=20),legend.position="top",
                                      axis.title.y=element_text(size=20),strip.text.x=element_text(size=18)) 
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"interv_incremcosts_icer",ifelse(grepl("disc",k_daly),"_disc",""),
       "_KEN_ZAF_3rows",ifelse(any(is.na(df_plot$CI95_high)),"_ci95_removed",""),".png"),
       width=25,height=40,units="cm")
}

write_csv(df_interv_incremcosts_icer,paste0("output/cea_plots/",subfolder_name,"df_interv_incremcosts_icer.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plots COMPARING CEA results with projected data vs new data
# subfolder_name<-"new_price_efficacy_kenyadeaths_CIs/"

for (k_plot in 1:3) {
  sel_vars <- list(c("total_YLD","total_YLL","hosp_YLD","non_hosp_YLD","total_DALY"),
                   c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),
                   c("admin_cost","cost_rsv_hosp","cost_rsv_outpatient", # ,"hosp_cost"
                      "outpatient_cost","total_medical_cost"))[[k_plot]]
  df_plot <- cea_summary_all %>% 
    filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) &
            ((price==3&intervention=="maternal")|(price==6&intervention=="mAb"))) %>%
    mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
       burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
       source=ifelse(grepl("projection",as.character(source)),
                     "projection (from [Shi 2017])",as.character(source)),
       vartype=gsub("hosp","hospitalised",
                    gsub("rsv","RSV",gsub("_"," ",gsub("_averted","",variable)))) ) %>% 
    mutate(vec=(10^floor(log10(median/norm_median)))*round(
      median/norm_median/(10^floor(log10(median/norm_median))),3)) %>%
    mutate(orig_burden_round=ifelse(vec>1e6,paste0(round(vec/1e6,1),"e6"),
                                    ifelse(vec>1e3,paste0(round(vec/1e3,1),"e3"),vec))) %>% 
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
  ggplot(df_plot) + geom_hpline(aes(x=country_iso,y=norm_median*1e2,group=interaction(source,intervention),
                    color=intervention,linetype=source),position=position_dodge(width=dodge_val),
                width=0.22,size=1.2) + 
    geom_linerange(aes(x=country_iso,ymin=norm_CI95_low*1e2,ymax=norm_CI95_high*1e2,
                       group=interaction(source,intervention),color=intervention),
                   alpha=0.35,position=position_dodge(width=dodge_val),size=17,show.legend=F) +
    facet_wrap(~vartype) + scale_y_continuous(breaks=(0:10)*10) + scale_color_manual(values=c("red","blue")) +
    geom_vline(xintercept=1.5,size=0.3) + geom_vline(xintercept=c(1,2),linetype="dashed",size=0.3) +
    geom_text(aes(x=country_iso,y=norm_CI95_high*1e2+2,group=interaction(source,intervention),
                label=orig_burden_round ), # ifelse(intervention!="MV",orig_burden_round,"")
                position=position_dodge(width=dodge_val)) +
    labs(color="",linetype="",caption=caption_txt) + scale_x_discrete(expand=expansion(0.02,0)) +
    xlab("") + ylab(ylab_txt) + 
    theme_bw() + standard_theme + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),
          axis.text.y=element_text(size=15),strip.text=element_text(size=14),legend.position="top",
          legend.text=element_text(size=14),axis.title.y=element_text(size=18)) + 
    guides(linetype=guide_legend(nrow=2),color=guide_legend(nrow=2))
  # save
  compar_dirname<-paste0("output/cea_plots/",subfolder_name,"comparisons/")
  if (!dir.exists(compar_dirname)) {dir.create(compar_dirname)}
  ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/",ifelse(any(grepl("YLL",sel_vars)),"DALY",
              ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),
        width=36,height=18,units="cm")
} # end of for-loop for diff vars

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables
price_interv_vals <- gsub("maternal","MV",unique(paste0(unique(cea_summary_all$price),
                        "$ (",unique(cea_summary_all$intervention),")")))
# process dataframe for plotting
df_total_DALY_medcost_averted_KEN_ZAF <- cea_summary_all %>% select(!(contains("norm")|plot_variable)) %>%
  filter(variable %in% c("total_DALY_averted","total_DALY_disc_averted",
                         "total_medical_cost_averted","incremental_cost")) %>% #
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention), # grepl("new",source)
         vec=sign(median)*abs(round((10^floor(log10(abs(median))))*round(
           median/(10^floor(log10(abs(median)))),3))),
         price_interv=paste0(price,"$ (",intervention,")"), # factor(,levels=price_interv_vals), 
         variable=factor(gsub("_"," ",variable),
                         levels=gsub("_"," ",c("total_DALY_averted","total_DALY_disc_averted",
                                  "total_medical_cost_averted","incremental_cost")) )) %>%
  mutate(orig_burden_round=ifelse(abs(vec)>1e6,paste0(round(vec/1e6,1),"e6"),
                                  ifelse(abs(vec)>1e3,paste0(round(vec/1e3,1),"e3"),vec)),
         price_interv=factor(price_interv,levels=unique(price_interv))) %>%
  filter(!((grepl("total_medical",variable)|grepl("total_DALY",variable)) & price>6)) %>%
  pivot_longer(!c(price,source_num,country_iso,intervention,source,variable,
                  country_plot,vec,price_interv,orig_burden_round)) %>%
  group_by(variable) %>% mutate(money_unit=ifelse(max(value)>1e6,"million USD","USD")) %>% ungroup() %>%
  mutate(value=ifelse(grepl("mill",money_unit),value/1e6,value/1e3),
    variable=ifelse(grepl("mill",money_unit),paste0(as.character(variable)," (million US$)"),
        paste0(as.character(variable)," (thousand US$)"))) %>% 
  mutate(variable=factor(variable,levels=unique(variable)[c(2,3,4,1)])) %>%
  pivot_wider(names_from=name,values_from=value) %>%
  mutate(orig_burden_round=ifelse(grepl("mill",money_unit),paste0(round(median,1),"m"),paste0(round(median,1),"k")),
         source=gsub(", community-based","",source)) %>%
  mutate(across(where(is.numeric),round,3))
ylab_txt="cost in USD (median, CI50, CI95)"
# plot total DALY AVERTED, total medical cost averted
ggplot(df_total_DALY_medcost_averted_KEN_ZAF %>% filter(!grepl("incremental",variable))) + 
  geom_hpline(aes(x=country_iso,y=median,group=interaction(source,intervention),
            linetype=source),position=position_dodge(width=dodge_val),width=0.23,size=2) + # color=price_interv,
  geom_linerange(aes(x=country_iso,ymin=CI95_low,ymax=CI95_high,group=interaction(source,intervention),
           color=intervention),alpha=0.25,position=position_dodge(width=dodge_val),size=42) +
  # geom_linerange(aes(x=country_iso,ymin=CI50_low,ymax=CI50_high,group=interaction(source,intervention),
  #          color=intervention),alpha=0.25,position=position_dodge(width=dodge_val),size=34,show.legend=F) +
  facet_wrap(~variable,scales="free_y",ncol=1) + # ,nrow=2
  scale_color_manual(values=c("red","blue")) +
  geom_vline(xintercept=1.5,size=0.5) + geom_vline(xintercept=c(1,2),linetype="dashed",size=0.75) +
  geom_vline(xintercept=c(0.75,1.25,1.75,2.25),linetype="dashed",size=0.25) + 
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=CI95_high+2,
                group=interaction(source,intervention),label=orig_burden_round),
             position=position_dodge(width=dodge_val),size=8) +
  scale_x_discrete(expand=expansion(0.02,0)) + labs(color="",linetype="") +#,caption="Numbers show median values."
  guides(color=guide_legend(ncol=2,override.aes=list(size=5)),
         linetype=guide_legend(nrow=2,override.aes=list(size=1.25))) +
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=18),axis.text.y=element_text(size=18),
        strip.text=element_text(size=22),legend.position="top",legend.text=element_text(size=21),
        plot.caption=element_text(10),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/total_DALY_medcost_averted_KEN_ZAF.png"),
       width=30,height=30,units="cm")
write_csv(df_total_DALY_medcost_averted_KEN_ZAF,
          paste0("output/cea_plots/",subfolder_name,"comparisons/total_DALY_medcost_averted_KEN_ZAF.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot incremental costs

ggplot(df_total_DALY_medcost_averted_KEN_ZAF %>% filter(grepl("incremental",variable))) + 
  geom_hpline(aes(x=country_iso,y=median,group=interaction(source,intervention),
                  linetype=source),position=position_dodge(width=dodge_val),width=0.23,size=2) + 
  geom_linerange(aes(x=country_iso,ymin=CI95_low,ymax=CI95_high,group=interaction(source,intervention),
              color=price_interv),alpha=0.4,position=position_dodge(width=dodge_val),size=44) +
  facet_wrap(~variable,scales="free_y",ncol=1) + # ,nrow=2
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),
                              colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=1.5,size=0.5) + geom_vline(xintercept=c(1,2),linetype="dashed",size=1) +
  geom_vline(xintercept=c(0.75,1.25,1.75,2.25),linetype="dashed",size=0.25) + 
  geom_hline(yintercept=0,linetype="dashed",size=0.5) +
  theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=ifelse(median>=0,CI95_high+2,CI95_low-2),
                group=interaction(source,intervention),label=orig_burden_round),
            position=position_dodge(width=dodge_val),size=7) +
  scale_x_discrete(expand=expansion(0.02,0)) + labs(color="",linetype="") +#,caption="Numbers show median values."
  guides(color=guide_legend(ncol=2,override.aes=list(size=4)),
         linetype=guide_legend(nrow=2,override.aes=list(size=1.4))) + #
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=20),axis.text.y=element_text(size=20),
        strip.text=element_text(size=22),legend.position="top",legend.text=element_text(size=17),
        axis.title.x=element_text(21),axis.title.y=element_text(size=18),
        plot.caption=element_text(size=10))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/incremental cost_KEN_ZAF.png"),
       width=35,height=25,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ICER plot

ylab_txt<-"cost in USD (mean, CI95)"
price_interv_vals <- gsub("maternal","MV",
                        unique(paste0(cea_summary_all$price,"$ (",cea_summary_all$intervention,")")))
df_plot_icer_comp <- cea_summary_all %>% select(!(contains("norm")|plot_variable)) %>% 
  filter(grepl("incremental_cost/DALY",variable)) %>%
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
         vec=sign(median)*abs(round((10^floor(log10(abs(median))))*
                                      round(median/(10^floor(log10(abs(median)))),3))),
         price_interv=factor(paste0(price,"$ (",intervention,")"),levels=price_interv_vals),
         variable=gsub("disc","(discounted)",gsub("_"," ",variable)) ) %>%
  mutate(orig_burden_round=paste0(round(median/1e3,1)),source=gsub(", community-based","",source)) %>%
  mutate(across(where(is.numeric),round,1))
# PLOT
ggplot(df_plot_icer_comp) + 
  geom_hpline(aes(x=country_iso,y=median/1e3,group=interaction(price,source),linetype=source),
                position=position_dodge(width=dodge_val),width=0.14,size=1.2) +
  geom_linerange(aes(x=country_iso,ymin=CI50_low/1e3,ymax=CI50_high/1e3,
                     group=interaction(price,source),color=price_interv),
                 alpha=0.35,position=position_dodge(width=dodge_val),size=13) +
  facet_grid(intervention~variable,scales="free") + 
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),
                              colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=c(1.5),size=0.5) + geom_vline(xintercept=c(1,2),size=0.75,linetype="dashed") +
  geom_vline(xintercept=setdiff(seq(2/3,4,by=1/6),c(1,1.5,2)),linetype="dashed",size=0.2) +
  geom_hline(yintercept=0,size=0.25) +
  xlab("") + ylab("incremental cost (thousand US$) per DALY averted") + 
  labs(color="",linetype="",caption="*cost-saving") + # 
  geom_text(aes(x=country_iso,y=(CI50_high+6e3)/1e3,group=interaction(price,source),
                  label=ifelse(median>=0,orig_burden_round,"*")),
                  size=6,position=position_dodge(width=dodge_val)) + 
  scale_x_discrete(expand=expansion(0.2,0)) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=20),
      axis.text.y=element_text(size=18),axis.title.y=element_text(size=19),
      strip.text=element_text(size=21),plot.caption=element_text(size=15),
      legend.position="top",legend.text=element_text(size=16)) + 
  guides(color=guide_legend(ncol=2),linetype=guide_legend(nrow=2))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"comparisons/ICER_KEN_ZAF.png"),width=38,height=32,units="cm")
#
write_csv(df_plot_icer_comp,paste0("output/cea_plots/",subfolder_name,"comparisons/ICER_KEN_ZAF.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# selvars=c("total_medical_cost","total_medical_cost_disc", # "intervention_cost"
#           "total_DALY_disc","total_DALY_averted","total_DALY_disc_averted",
#           "total_medical_cost_averted","total_medical_cost_disc_averted",
#           "incremental_cost","incremental_cost_disc",
#           "incremental_cost/DALY_averted","incremental_cost/DALY_disc_averted",
#           "hosp_cases_averted","SARI_averted")
###################
# MV effic for lmic: 
# list(sympt_disease=c(mean=0.405,CI95_low=-0.031,CI95_high=0.657),
# hospit=c(mean=0.542,CI95_low=0.295,CI95_high=0.702),
# severe=c(mean=0.493,CI95_low=-0.271,CI95_high=0.798) )
# resvax_trial_figures <- read_csv("custom_input/resvax_trial_figures.csv")
# resvax_efficacy_figures <- resvax_trial_figures %>% 
#   mutate(vacc_status=ifelse(vacc_status=="vacc",TRUE,FALSE)) %>% group_by(RSV_LRTI_categ,period,geogr) %>% 
#   summarise(VE=round((1-rate[vacc_status]/rate[!vacc_status])*1e2,1),
#       RR_ci95_l=round(exp(log(rate[vacc_status]/rate[!vacc_status])-
#                       1.96*sqrt( (1-rate[vacc_status])/n_cases[vacc_status] + 
#                         (1-rate[!vacc_status])/n_cases[!vacc_status])),3),
#       RR_ci95_u=round(exp(log(rate[vacc_status]/rate[!vacc_status]) +
#                       1.96*sqrt( (1-rate[vacc_status])/n_cases[vacc_status] + 
#                         (1-rate[!vacc_status])/n_cases[!vacc_status])),3) ) %>%
#   mutate(VE_ci95_l=round((1-RR_ci95_u)*1e2,1),VE_ci95_u=round((1-RR_ci95_l)*1e2,1))
# write_csv(resvax_efficacy_figures,file="custom_input/resvax_efficacy_figures.csv")
###################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# fitting CI95 with a beta distrib
# efficacy_figures$monocl_ab$hospit[c("CI95_low","mean","CI95_high")]
# comparing beta and gamma distrib
# beta_pars <- get.beta.par(q=efficacy_figures$monocl_ab$hospit[c("CI95_low","mean","CI95_high")],
#              p=c(2.5,50,97.5)/100,show.output=F,plot=F)
# gamma_pars <- get.gamma.par(q=c(0.519,0.784,0.903),
#               p=c(2.5,50,97.5)/100,show.output=F,plot=F)[c("shape","rate")]
# gamma_sim <- rgamma(n=1e4,shape=gamma_pars[1],rate=gamma_pars[2]) 
# beta_sim <- rbeta(n=1e4,shape1=beta_pars[1],shape2=beta_pars[2]); 
# efficacy_figures$monocl_ab$hospit
# c(mean=mean(beta_sim),ci95=quantile(beta_sim,probs=c(2.5,97.5)/100))
# c(mean=mean(gamma_sim),ci95=quantile(gamma_sim,probs=c(2.5,97.5)/100))
# #
# norm_pars<-get.norm.par(q=c(-0.01,0.784,0.903),p=c(2.5,50,97.5)/100,show.output=F,plot=F)
# norm_sim <- rnorm(n=1e4,mean=norm_pars[1],sd=norm_pars[2]) 
# c(mean=0.394,CI95_low=0.053,CI95_high=0.612) -> c(5.3,39.4,61.2)
# severe dis: # c(mean=0.483,-8.2/100,0.753)
# beta_pars<-get.beta.par(p=c(2.5,50,97.5)/100,
#                         q=c(0,48.3,75.3)/100,show.output=F,plot=F)[c("shape1","shape2")]
# beta_sim<-rbeta(n=1e6,shape1=beta_pars[1],shape2=beta_pars[2])
# c(mean(beta_sim),quantile(beta_sim,probs=c(2.5,97.5)/100))
# hist(beta_sim,breaks=100); abline(v=mean(beta_sim),col="red")
# #
# gamma_pars<-get.gamma.par(p=c(2.5,50,97.5)/100,q=c(5.3,39.4,61.2)/100,show.output=F,plot=F)
# gamma_sim<-rgamma(n=1e6,shape=gamma_pars[1],rate = gamma_pars[2])
# c(mean(gamma_sim),quantile(gamma_sim,probs=c(2.5,97.5)/100))
# hist(gamma_sim,breaks=100)
# lnorm_sim<-rlnorm(n=1e4,meanlog=lnorm_par["meanlog"],sdlog=lnorm_par["sdlog"])
# c(mean(lnorm_sim),quantile(lnorm_sim,probs=c(2.5,97.5)/100))
