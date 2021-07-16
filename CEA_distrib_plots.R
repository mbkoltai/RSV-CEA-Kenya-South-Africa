# clean up workspace
rm(list=ls())
lapply(c("tidyverse","rstudioapi","fitdistrplus","rstudioapi","matrixStats","ungeviz","stringi","rriskDistributions"),library,character.only=TRUE)
# sessionInfo()
# path
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
num_sim <- 5000; source('functions/RSV_load_all.R'); source('functions/load_config_pars.R')
source("functions/set_xlims_cea.R") # creates 'sim_config_matrix'
# load functions
lapply(c("functions/get_burden_flexible.R","functions/get_burden_flexible_ari_sari.R"),function(x) {source(x)})
# load data
### Kenya incidence data -------------------------
# hosp rate (p): p/(1-p) ~ norm
# load functions
source("functions/load_own_data.R"); source('functions/GammaParmsFromQuantiles.R')
# load & process Kenya SARI data, create 5e3 sample paths with CI95 corresponding to those in data
kenya_data_file_path="../path_rsv_data/SARI_Rates_2010_2018_updated/ARI_SARI_Rates_2010_2018_tidydata.csv"
# hospitalisation rate in mcmarcel
# kenya_hosp_rates_mcmarcel=as.numeric(get_rsv_ce_config(subset(sim_config_matrix,
  # country_iso %in% "KEN"&intervention=="maternal"))$hosp_prob[1,])
# kenya_incid_hosp_rate=fcn_gen_incid_hospprob_samples_kenya(kenya_data_path=kenya_data_file_path,sel_disease="SARI",n_iter=5e3,
#   age_maxval=60,CI_intervals=c(2.5,97.5)/1e2,logit_hosp_rate_stdev=sd(log(kenya_hosp_rates_mcmarcel/(1-kenya_hosp_rates_mcmarcel))),
#   n_length=2e2,x_prec=0.01,randsampl_distrib_type='gamma')
# create incid matrices for nonhosp and hosp cases directly from data
kenya_nonhosp_hosp_incid_ari_sari=lapply(c("ARI","SARI"), function(x)
    fcn_gen_nonhosp_hosp_incid_samples_kenya(kenya_data_file_path,sel_disease=x,n_iter=5e3,age_maxval=60,
    CI_intervals=c(2.5,97.5)/1e2,randsampl_distrib_type="gamma")); names(kenya_nonhosp_hosp_incid_ari_sari)=c("ARI","SARI")
# ggplot(data.frame(age=0:59,nonhosp=rowMeans(kenya_nonhosp_hosp_incid_ari_sari$ARI$nonhosp_incid),
#   hosp=rowMeans(kenya_nonhosp_hosp_incid_ari_sari$ARI$hosp_incid) ) %>% mutate(age=factor(age)) %>% pivot_longer(!age)) +
#   geom_bar(aes(x=age,y=value,group=name,fill=name),position="stack",stat="identity") + theme_bw()
### ### ### ### ### ### ### ### ### ### ### ### ###
# analyse Kenya data
kenya_data=bind_rows(fcn_load_kenya(kenya_data_path=kenya_data_file_path,sel_disease="ARI")$rsv_incidence_ageinf,
                     fcn_load_kenya(kenya_data_path=kenya_data_file_path,sel_disease="SARI")$rsv_incidence_ageinf) %>% 
  mutate(disease_type_medic_status=paste(disease_type, ifelse(medically_attended,"medically attended","not attended")) ) %>%
  mutate(disease_type_medic_status=gsub("SARI medically attended","SARI hospitalised",disease_type_medic_status)) %>%
  mutate(disease_type_medic_status=gsub("SARI not attended","SARI non-hospitalised",disease_type_medic_status)) %>%
  dplyr::select(!c(RSV_assoc,freq,n)) %>% relocate(age_inf,.before=variable) %>% relocate(metric_per_popul,.after=value)
# ratio of medically-attended uniform across age groups (65% for ARI, 24% SARI)
hosp_rate_kenya=kenya_data %>% group_by(disease_type,age_inf) %>% 
  summarise(medic_attended=unique(value[medically_attended]/sum(value)),maxval=unique(max(value)/metric_per_popul),
  status=unique(disease_type_medic_status)) %>% group_by(disease_type) %>% summarise(mean_medic_attended=mean(medic_attended),
  sd=sd(medic_attended),maxval=max(maxval),status=unique(status)[!grepl("not|non",unique(status))])
# PLOT burden ARI + SARI
p <- ggplot(kenya_data,aes(x=age_inf)) + 
  geom_bar(aes(y=value/metric_per_popul,fill=disease_type_medic_status),position="stack",stat="identity") +
  theme_bw() + standard_theme + xlab("age (months)") + ylab("cases per person year") + theme(axis.text.x=element_text(vjust=0.5),
  legend.text=element_text(size=12),legend.background=element_rect(fill=NA),legend.position=c(0.85,0.925)) + # 
  geom_text(data=hosp_rate_kenya,aes(x=45,y=c(0.17,0.2),label=paste0(status,"=",round(mean_medic_attended*1e2,1),"%")),size=6) +
  scale_x_continuous(breaks=0:59,expand=expansion(0.01,0)) + labs(fill="") + facet_wrap(~disease_type,nrow=2,scales="free") + 
  scale_y_continuous(expand=expansion(0.01,0)) + ggtitle("Kenya ARI and SARI burden"); p
# save
ggsave(paste0("output/ari_sari_burden/kenya_ari_sari_burden","_sep"[length(p$facet$params)>0],".png"),width=30,height=22,units="cm")
# with error bars
df_plot <- kenya_data %>% group_by(age_in_months,disease_type_medic_status) %>% 
  summarise(value=unique(value),CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),metric_per_popul=unique(metric_per_popul))
ggplot(df_plot,aes(x=age_in_months)) +
  geom_bar(aes(y=value/metric_per_popul,fill=disease_type_medic_status),position="stack",stat="identity") +
  geom_errorbar(aes(ymin=CI_95_lower/metric_per_popul,ymax=CI_95_upper/metric_per_popul),size=0.4) +
  facet_wrap(~disease_type_medic_status,nrow=2,scales = "free") + # 
  theme_bw() + standard_theme + xlab("age (months)") + ylab("cases per person year") + theme(axis.text.x=element_text(vjust=0.5,size=13),
   axis.text.y=element_text(size=13),legend.text=element_text(size=13),legend.background=element_rect(fill=NA),legend.position=c(0.92,0.925),
   axis.title.x=element_text(size=17),axis.title.y=element_text(size=17),strip.text=element_text(size=14)) + 
  # geom_text(data=hosp_rate_kenya,aes(x=45,y=c(0.17,0.2),label=paste0(status,"=",round(mean_medic_attended*1e2,1),"%")),size=6) +
  scale_y_continuous(expand=expansion(0.01,0)) + labs(fill="") # scale_x_continuous(breaks=(0:30)*2,expand=expansion(0.01,0))
  # ggtitle("Kenya ARI and SARI burden")
# save
ggsave(paste0("output/ari_sari_burden/kenya_ari_sari_burden_errorbars_grouped.png"),width=35,height=22,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ###
# ### ARI/SARI ratio per age group
# kenya_ari_sari_ratio = kenya_data %>% group_by(disease_type,age_inf) %>% summarise(value=sum(value)/unique(metric_per_popul)) %>% 
#   group_by(age_inf) %>% summarise(value=value[disease_type=="ARI"]/value[disease_type=="SARI"])
# # plot ARI/SARI RATIO
# ggplot(kenya_ari_sari_ratio,aes(y=value)) + geom_segment(aes(x=age_inf-0.45,xend=age_inf+0.45,yend=value),size=1.05) +
#   scale_x_continuous(breaks=0:60,expand=expansion(0.01,0)) + scale_y_log10(breaks=round(10^((-4:10)/4),2) ) +
#   geom_segment(aes(x=age_inf,xend=age_inf,y=0,yend=value),size=0.5,linetype="dashed",color="darkgrey") + 
#   theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) + xlab("age (months)") + ylab("ARI/SARI ratio")
# # SAVE
# ggsave(paste0("output/ari_sari_burden/kenya_ari_sari_ratio.png"),width=30,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ###
# Kenya deaths
deaths_kenya <- read_csv("../path_rsv_data/SARI_Rates_2010_2018_updated/deaths_kenya_tidy.csv") %>% 
  filter(variable=="rate" & !age_in_months %in% c("<12","12-23","<24","24-59","<60")) %>% 
  mutate(age_in_months=ifelse(age_in_months=="<1","0",age_in_months),freq=1) %>% 
  mutate(freq=ifelse(grepl('-',age_in_months),
                     as.numeric(sapply(age_in_months, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,freq)) %>% 
  mutate(age_in_months=ifelse(grepl('-',age_in_months), sapply(strsplit(age_in_months,'-'),'[[',1),age_in_months)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% mutate(age_inf=as.numeric(age_in_months)+(n-1)) %>% 
  dplyr::select(!c(n,freq,age_in_months)) %>% relocate(age_inf,.before=value)
# generate samples
deaths_distrib_params = bind_rows(lapply(c("yes","no"), function(y_no) data.frame(age_inf=0:59, 
    t(sapply(1:60, function(x) gamma.parms.from.quantiles(p=c(2.5,97.5)/100, 
    q=as.numeric((deaths_kenya %>% mutate(CI_95_lower=ifelse(CI_95_lower==0,0.1,CI_95_lower)) %>%
      filter(in_hospital==y_no) %>% dplyr::select(c(CI_95_lower,CI_95_upper)))[x,]))[c("shape","rate")])),in_hospital=y_no))) %>%
      mutate(shape=unlist(shape),rate=unlist(rate))
kenya_deaths_incid = lapply(c("yes","no"), function(y_no)
  t(sapply(0:59, function(x) rgamma(5e3,shape=(deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$shape,
      rate=(deaths_distrib_params %>% filter(age_inf==x&in_hospital==y_no))$rate)))/1e5); names(kenya_deaths_incid)=c("hosp","nonhosp")
### ### ### ### ### ### ### ###
# plot CFR
ggplot(deaths_kenya,aes(x=age_inf,y=value,color=in_hospital,fill=in_hospital)) + geom_point() +
  geom_line() + geom_ribbon(aes(ymin=CI_95_lower,ymax=CI_95_upper),alpha=0.2) + theme_bw() + standard_theme +
  scale_x_continuous(expand=expansion(0.01,0),breaks=(0:30)*2) + scale_y_continuous(expand=expansion(0.01,0)) +
  theme(axis.text.x=element_text(vjust=0.5)) + xlab("age in months") + ylab("RSV-associated deaths/100.000 population")
# calculate in-hosp and out-hosp CFR 
cfr_kenya=bind_rows(data.frame(age_inf=0:59,rate=deaths_kenya$value[deaths_kenya$in_hospital=="yes"]/(kenya_data %>% 
        filter(medically_attended==TRUE & disease_type=="SARI"))$value,medically_attended=TRUE),
      data.frame(age_inf=0:59,rate=deaths_kenya$value[deaths_kenya$in_hospital!="yes"]/(kenya_data %>% 
        filter(medically_attended!=TRUE & disease_type=="SARI"))$value,medically_attended=FALSE))
# plot
ggplot(cfr_kenya,aes(x=age_inf,y=rate*100,color=medically_attended)) + geom_line() + geom_point() + theme_bw() + standard_theme +
  ylab("% CFR") + scale_x_continuous(expand=expansion(0.01,0),breaks=(0:30)*2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### South Africa
SA_SARI_data <- fcn_load_s_afr(safr_data_path = "../path_rsv_data/s_afr_incidence_data_rate.csv") %>%
  mutate(disease_type_medic_status=paste(disease_type,ifelse(hospitalisation,"hospitalised","not hospitalised")) ) %>% 
  rename(agegroup_mts=age) %>% relocate(age_inf,.before=Province) %>% relocate(Province,.after=disease_type_medic_status) %>% 
  relocate(year,.before=Province)
# SA ILI data
SA_ILI_data <- read_csv("../path_rsv_data/s_afr_ILI_incidence_rate.csv") %>% 
      filter(!(grepl("<",agegroup) | agegroup %in% c("0-5m","6-11m","12-23m","24-59m","<5y"))) %>%
      mutate(agegroup_mts=agegroup,agegroup=gsub("m","",agegroup), freq=1) %>% mutate(freq=ifelse(grepl('-',agegroup),
        as.numeric(sapply(agegroup, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,freq)) %>% 
  mutate(agegroup=ifelse(grepl('-',agegroup), sapply(strsplit(agegroup,'-'),'[[',1),agegroup)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% mutate(age_inf=as.numeric(agegroup)+(n-1)) %>% dplyr::select(!c(n,freq,agegroup)) %>% # 
  relocate(age_inf,.before=rate) %>% relocate(disease_type,.before=hospitalisation) %>%
  mutate(disease_type=ifelse(disease_type=="ILI","ARI",""),
    disease_type_medic_status=ifelse(hospitalisation,paste0("medically attended ",disease_type),
      paste0("non medically attended ",disease_type))) %>% relocate(disease_type_medic_status,.before=Province)
# concatenate
if (!exists("SA_data")){ SA_data=bind_rows(SA_ILI_data,SA_SARI_data) }
# subtract SARIs from ARIs? don't do this - orig data mutually exclusive!!
##
hosp_rate_sa = SA_data %>% group_by(disease_type,age_inf) %>% summarise(medic_attended=unique(rate[hospitalisation]/sum(rate,na.rm=T)),
  maxval=unique(max(rate)/popul_denom),status=unique(disease_type_medic_status)) %>% group_by(disease_type) %>% 
  summarise(mean_medic_attended=mean(medic_attended,na.rm=T),sd=sd(medic_attended,na.rm=T),maxval=max(maxval,na.rm=T),
  status=unique(status)[!grepl("not|non",unique(status))]) %>% mutate(status=gsub("^ARI hospitalised","ARI medic. attended",status))
# plot
p <- ggplot(SA_data,aes(x=age_inf)) + # %>% mutate(disease_type=gsub("^ARI","ARI (ILI)",disease_type))
  geom_bar(aes(y=rate/popul_denom,fill=disease_type_medic_status),position="stack",stat="identity") +
  facet_wrap(~disease_type,nrow=2,scales="free") + 
  theme_bw() + standard_theme + xlab("age (months)") + ylab("cases per person year") +
  theme(axis.text.x=element_text(vjust=0.5),legend.text=element_text(size=12),legend.background=element_rect(fill=NA),
        legend.position=c(0.85,0.38)) + geom_text(data=hosp_rate_sa,aes(x=33,y=list(c(0.08,0.1),c(0.55,0.6))[[1]],
  label=paste0(status,"=",round(mean_medic_attended*1e2,1),"%")),size=6) + scale_x_continuous(breaks=0:59,expand=expansion(0.01,0)) +
  scale_y_continuous(expand=expansion(0.02,0)) + labs(fill=""); p
# SAVE
ggsave(paste0("output/ari_sari_burden/s_afr_ari_sari_burden","_sep"[length(p$facet$params)>0],"_ili.png"),width=30,height=22,units="cm")
# PLOT with error bars
df_plot <- SA_data %>% mutate(agegroup_mts=factor(gsub("m","",agegroup_mts),levels=unique(gsub("m","",agegroup_mts)))) %>%
  group_by(agegroup_mts,disease_type_medic_status) %>% 
  summarise(rate=unique(rate),rate_CI_lower=unique(rate_CI_lower),rate_CI_upper=unique(rate_CI_upper),popul_denom=unique(popul_denom))
ggplot(df_plot,aes(x=agegroup_mts)) + # %>% mutate(disease_type=gsub("^ARI","ARI (ILI)",disease_type))
  geom_bar(aes(y=rate/popul_denom,fill=disease_type_medic_status),position="stack",stat="identity") +
  geom_errorbar(aes(ymin=rate_CI_lower/popul_denom, ymax=rate_CI_upper/popul_denom),size=0.4) +
  facet_wrap(~disease_type_medic_status,scales="free") + # ,nrow=2
  theme_bw() + standard_theme + xlab("age (months)") + ylab("cases per person year") +
  theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=13),legend.text=element_text(size=13),
        legend.background=element_rect(fill=NA),legend.position=c(0.85,0.38),strip.text=element_text(size=15),
        axis.title.x=element_text(size=17),axis.title.y=element_text(size=17)) + 
  scale_y_continuous(expand=expansion(0.02,0)) + labs(fill="") # scale_x_continuous(breaks=(0:30)*2,expand=expansion(0.01,0)) + 
# save
ggsave(paste0("output/ari_sari_burden/s_afr_ari_sari_burden_ili_errorbars_grouped.png"),width=35,height=22,units="cm")
### adjustment: project ILIs to all cases (not just those with fever) by dividing rate by proportion of ARI cases with fever (in Kenya)
# % of cases with fever: 2015: 33.3%; 2016 20.5%; 2017, 16%
ILI_adjust_SA=TRUE
if (ILI_adjust_SA & ifelse(!exists("divided_fever"),TRUE,!divided_fever)) { 
  print("divide by fever proportion")  
  SA_data[SA_data$disease_type=="ARI",c("rate","rate_CI_lower","rate_CI_upper")]=
    SA_data[SA_data$disease_type=="ARI",c("rate","rate_CI_lower","rate_CI_upper")]/mean(c(33.3,20.5,16)/100); divided_fever=TRUE }
# CI lower limit should not be 0! 
SA_data$rate_CI_lower[SA_data$age_inf==0 & SA_data$disease_type=="ARI"]=1
### generate 5e3 sample paths for CEA
sa_nonhosp_hosp_incid_ari_sari=lapply(c("ARI","SARI"), function(x) fcn_gen_nonhosp_hosp_incid_samples_SA(SA_data,diseasetype=x,
    n_iter=5e3, age_maxval=60,CI_intervals=c(2.5,97.5)/1e2,randsampl_distrib_type="gamma"))
names(sa_nonhosp_hosp_incid_ari_sari)=c("ARI","SARI")
# check plot
# ggplot(data.frame(age=0:59,nonhosp=rowMeans(sa_nonhosp_hosp_incid_ari_sari$ARI$nonhosp_incid),
#   hosp=rowMeans(sa_nonhosp_hosp_incid_ari_sari$ARI$hosp_incid) ) %>% mutate(age=factor(age)) %>% pivot_longer(!age)) + 
#   facet_wrap(~name,nrow=2) + geom_bar(aes(x=age,y=value,group=name,fill=name),position="stack",stat="identity") + theme_bw()
# SA costing data # for clinics the cost is always 25USD
s_afr_inpatient_cost <- read_csv("../path_rsv_data/s_afr_PDE_calcs.csv") %>% mutate(freq=ifelse(grepl('-',age),
  as.numeric(sapply(age, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,1)) %>% 
  mutate(age=ifelse(grepl('-',age), sapply(strsplit(age,'-'),'[[',1),age)) %>% uncount(weights=freq, .id="n",.remove=F) %>%
  mutate(age=as.numeric(age)+(n-1)) %>% dplyr::select(!c(n,freq))
s_afr_outpatient_cost <- cbind(data.frame(age="all",mean=25,LCI=18.3,UCI=31.8),
  data.frame(t(unlist(gamma.parms.from.quantiles(q=c(18.3,31.8),p=c(2.5,97.5)/100)[c("shape","rate")]))) )
# plot SA cost
ggplot(s_afr_inpatient_cost %>% filter(name %in% "total"),
       aes(x=age+0.42,y=`Mean cost per illness episode (USD)`,group=1)) + 
  geom_segment(aes(xend=age-0.42,yend=`Mean cost per illness episode (USD)`),size=1.05) + # geom_bar(position="stack",stat="identity")
  geom_linerange(aes(x=age,ymin=LCI,ymax=UCI),color="red",alpha=0.3,size=3) +
  geom_vline(aes(xintercept=age+0.5),size=0.2,linetype="dashed",color="blue") +
  theme_bw() + standard_theme + scale_x_continuous(breaks=0:59,expand=expansion(0,0)) + scale_y_continuous(breaks=(0:11)*100) + 
  theme(axis.text.x=element_text(vjust=0.5)) + xlab("age (month)") + ggtitle("South Africa RSV inpatient cost")
# save
ggsave(paste0("output/s_afr_costing.png"),width=30,height=18,units="cm")
# generate distributions from mean and CIs
if (!any(grepl("shape",colnames(s_afr_inpatient_cost)))){
s_afr_inpatient_cost = cbind(s_afr_inpatient_cost, t(sapply(1:nrow(s_afr_inpatient_cost), function(x) 
 unlist(gamma.parms.from.quantiles(q=c(s_afr_inpatient_cost$LCI[x],s_afr_inpatient_cost$UCI[x]),p=c(2.5,97.5)/1e2)[c("shape","rate")])))) }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cntrs_cea=c("KEN","ZAF"); subfolder_name="new_price_efficacy/"
# efficacy figures for vaccine for RESVAX (Novavax trial); for mAb from NIRSEVIMAB (AstraZ)
efficacy_figures=list(mat_vacc=list(sympt_disease=c(mean=0.394,CI95_low=5.3/100,CI95_high=0.612),
                      hospit=c(mean=0.444,CI95_low=5.3/100,CI95_high=0.612),severe=c(mean=0.483,CI95_low=-8.2/100,CI95_high=0.753)),
                      monocl_ab=list(sympt_disease=c(mean=0.726,CI95_low=0.565,CI95_high=0.831),
                                     hospit=c(mean=0.8,CI95_low=0.55,CI95_high=0.911)))
# outputs to display
all_cols=c("non_hosp_cases","hosp_cases","rsv_deaths","total_DALY_disc",
           "cost_rsv_hosp","total_medical_cost","incremental_cost","total_medical_cost_averted",
           "non_hosp_cases_averted","hosp_cases_averted", "rsv_deaths_averted", "total_DALY_disc_averted","SARI_averted",
           "incremental_cost/DALY_averted","total_YLD","total_YLL","hosp_SARI","non_hosp_SARI","ARI_averted") 
selvars=c("total_medical_cost","total_DALY_disc", # "total_DALY_disc","intervention_cost"
          "total_DALY_averted","total_medical_cost_averted","incremental_cost","incremental_cost/DALY_averted",
          "hosp_cases_averted","SARI_averted")
burden_cols <- all_cols[!grepl("cost|averted",all_cols)]; cost_cols <- all_cols[grepl("cost",all_cols)]
# prices for doses
pricelist=list("mat_vacc"=c(3,10,30),"mAb"=c(6,20,60))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# CEA calculation for 2 cntrs and 2 interv types
read_calc_flag=c("calc","read")[1]; kenya_deaths_input=TRUE
# loop through: cntrs * interventions * dose prices
for (n_cntr_output in 1:length(cntrs_cea)){
  for (n_interv in 1:2){ # SELECT INTERVENTION: 1=MatVacc, 2=monocl Abs
# intervention config table
sel_interv=sim_config_matrix[which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])[n_interv],]
if (cntrs_cea[n_cntr_output]=="ZAF"){sel_interv$country_iso=cntrs_cea[n_cntr_output]}
######
# calculate
if (grepl("calc",read_calc_flag)) {
for (k_price in 1:length(pricelist[[n_interv]])) {
doseprice=c("mat_vacc"=ifelse(n_interv==1,pricelist$mat_vacc[k_price],pricelist$mat_vacc[1]),
            "mAb"=ifelse(n_interv==2,pricelist$mAb[k_price],pricelist$mAb[1]))
# calculation with data from mcmarcel (community-based)
sim_output=get_burden_flexible(sel_interv,NA,NA,exp_wane=TRUE,doseprice)
# SARI, ARI distinguished, hosp/nonhosp distinguished
if (cntrs_cea[n_cntr_output]=="ZAF") {
  cost_input <- list("inpatient"=s_afr_inpatient_cost %>% filter(name %in% "total"),"outpatient"=s_afr_outpatient_cost)} else {
    cost_input<-NA}
# input death data if available
if (kenya_deaths_input) { print("using propr death data"); kenya_nonhosp_hosp_incid_ari_sari$deaths$non_hosp=kenya_deaths_incid$nonhosp
  kenya_nonhosp_hosp_incid_ari_sari$deaths$hosp=kenya_deaths_incid$hosp }
# if "effic_prob" is set to TRUE than the CIs for efficacy will be used, otherwise only the means
sim_output_user_input_ari_sari <- get_burden_flexible_ari_sari(sel_interv,
    list(kenya_nonhosp_hosp_incid_ari_sari,sa_nonhosp_hosp_incid_ari_sari)[[n_cntr_output]],
    efficacy_figures,effic_prob=TRUE,exp_wane=TRUE,doseprice,cost_data=cost_input)
### Plot results -------------------------
# are nonhosp SARIs accounted as SARIs? (or as nonhosp = mild cases)
burden_mcmarcel_owndata=fcn_process_burden_output(user_output=sim_output_user_input_ari_sari,default_output=sim_output,
 sel_cntr=sel_interv$country_iso,cols_burden_sel="",
 plot_labels=c(mcmarcel="projection ([Li 2020] from [Shi 2017] meta-analysis, community-based)",own="new data (hospital-based + HUS)"),
 icercolname="incremental_cost/DALY_averted",icercols=c("incremental_cost_disc","total_DALY_disc_averted"))
# calculate mean, median, CI95
x=burden_mcmarcel_owndata %>% group_by(source,variable) %>% summarise(mean=mean(value,na.rm=T),median=median(value,na.rm=T),
        CI50_low=quantile(value,probs=c(25,75)/1e2,na.rm=T)[1],CI50_high=quantile(value,probs=c(25,75)/1e2,na.rm=T)[2],
        CI95_low=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[1],CI95_high=quantile(value,probs=c(2.5,97.5)/1e2,na.rm=T)[2]) %>%
    mutate(price=doseprice[n_interv],source_num=c(0.2,2)[as.numeric(source)],
         country_iso=sel_interv$country_iso,intervention=sel_interv$intervention)
if (k_price==1) {cea_summary=x} else {cea_summary=bind_rows(cea_summary,x)}  } # END price scan
# save CEA summary table
write_csv(cea_summary,paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                                               gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv") )} else { 
      cea_summary <- read_csv(paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                               gsub("maternal","Mat_Vacc",sel_interv$intervention),".csv")) 
      if (n_cntr_output==1 & n_interv==1) {cea_summary_all=cea_summary} else {cea_summary_all = bind_rows(cea_summary_all,cea_summary)}
      }
# plot
for (k_col in 1:4){
df_plot=subset(cea_summary, variable %in% list(selvars,all_cols,burden_cols,cost_cols)[[k_col]])
if (k_col==3) {df_plot=subset(df_plot,price==min(price) ); x_dodge_val=0.1 } else {x_dodge_val=0.35}
# PLOT
p <- ggplot(df_plot,aes(x=factor(price))) + 
  # geom_linerange(aes(ymin=CI95_low,ymax=CI95_high,group=interaction(source,price),color=factor(source)),
  # alpha=0.3,position=position_dodge(width=x_dodge_val),size=3) + 
  geom_point(aes(y=mean,group=interaction(source,price)),pch="-",size=ifelse(k_col==2,12,14),position=position_dodge(width=x_dodge_val))+
  geom_linerange(aes(ymin=CI50_low,ymax=CI50_high,group=interaction(source,price),color=factor(source)),alpha=0.5,
                 position=position_dodge(width=x_dodge_val),size=ifelse(k_col==2,3,4)) + facet_wrap(~variable,scales="free") +
  geom_rect(data=subset(df_plot,grepl("cost",variable)),fill=NA,colour="blue",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
  # geom_rect(data=subset(df_plot,!grepl("cost|averted",variable)),fill=NA,colour="red",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
  geom_rect(data=subset(df_plot,grepl("averted",variable)),fill=NA,colour="green",size=1.4,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
  theme_bw() + standard_theme + labs(color="data source") + scale_x_discrete(expand=c(0,ifelse(k_col!=3,0.4,0.1))) + 
  theme(legend.position="bottom",axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=ifelse(k_col==2,8,11)),
        strip.text=element_text(size=ifelse(k_col==2,8,12)),legend.text=element_text(size=11)) + 
  xlab("dose price") + ylab("mean (CI50, CI95)") + geom_text(aes(x=factor(price),y=mean,group=interaction(source,price),
    label=ifelse(mean>0,paste0(round(mean/(10^floor(log10(mean))),1),"e",floor(log10(mean))),
    paste0(round(mean/(10^floor(log10(abs(mean)))),1),"e",floor(log10(abs(mean))))) ),size=ifelse(k_col==2,3,4.5),#check_overlap=T,
    position=position_dodge(width=ifelse(k_col==3,0.13,ifelse(k_col!=2,0.9,1.05))),angle=ifelse(k_col!=2,90,90),show.legend=F) +
  ggtitle(paste0(cntrs_cea[n_cntr_output]," cost effectiveness for: ",gsub("maternal","maternal vaccination",sel_interv$intervention)))
if (min(df_plot$mean,na.rm=T)>0) {p=p+scale_y_log10(expand=expansion(0.19,0));p} else {p=p+scale_y_continuous(expand=expansion(0.3,0));p}
if (length(unique(df_plot$price))==1) {p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + xlab("")+
  ggtitle(paste0(cntrs_cea[n_cntr_output]," disease burden")); p}
# filename
cea_summ_plot_filename=paste("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_",
                             gsub("maternal","Mat_Vacc",sel_interv$intervention),".png",sep="") # calc_tag
if (identical(unique(df_plot$variable), sort(selvars))) {cea_summ_plot_filename=gsub(".png","_SELVARS.png",cea_summ_plot_filename)}
if (identical(unique(df_plot$variable), sort(burden_cols))){
  cea_summ_plot_filename=paste0("output/cea_plots/",subfolder_name,sel_interv$country_iso,"_cea_summary_mean_CI_burden.png")}
if (identical(unique(df_plot$variable), sort(cost_cols))){cea_summ_plot_filename=gsub(".png","_COSTS.png",cea_summ_plot_filename)}
# SAVE
ggsave(cea_summ_plot_filename,width=36,height=18,units="cm")
    } # loop for plot
  } # loop intervention type
} # loop country

# filenames=list.files(paste0("output/cea_plots/",subfolder_name),pattern=".csv");filenames=filenames[!filenames %in% "cea_summary_all.csv"]
# for (k_filename in 1:4) {
#   x=read_csv(paste0("output/cea_plots/",subfolder_name,filenames[k_filename]))
#   if (k_filename==1){ cea_summary_all=x } else {cea_summary_all=bind_rows(cea_summary_all,x)} }
# write_csv(cea_summary_all,paste0("output/cea_plots/",subfolder_name,"cea_summary_all.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Make summary plots
#
# composition of total burden
# total_DALY <- total_YLD + total_YLL
# total_YLL <- rsv_deaths*config$hosp_CFR_DALYloss
# total_YLD <- hosp_med_att_YLD + non_hosp_YLD
# hosp_med_att_YLD <- hosp_SARI * config$severe_rsv_DALYloss + med_att_ARI*config$non_severe_rsv_DALYloss
# hosp_YLD <- hosp_SARI* config$severe_rsv_DALYloss
# non_hosp_YLD <- non_hosp_SARI*config$severe_rsv_DALYloss + non_med_att_ARI*config$non_severe_rsv_DALYloss
for (k_plot in 1:3) {
sel_vars <- list(c("total_YLD","total_YLL","hosp_YLD","hosp_med_att_YLD","non_hosp_YLD","total_DALY"), # ,"ARI_YLD","SARI_YLD"
            c("rsv_deaths","hosp_SARI","non_hosp_SARI","hosp_cases","non_hosp_cases"),
            c("admin_cost","cost_rsv_hosp","hosp_cost","cost_rsv_outpatient","outpatient_cost","total_medical_cost"))[[k_plot]]
df_plot <- cea_summary_all %>% filter(variable %in% c(sel_vars,paste0(sel_vars,"_averted")) & intervention=="mAb" & 
      ((price==3&intervention=="maternal")|(price==6&intervention=="mAb")) & grepl("new",source)) %>%
  mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
        burden_interv=ifelse(grepl("averted",variable),"averted burden","burden"), 
        source=ifelse(grepl("projection",as.character(source)),"projection (from [Shi 2017])",as.character(source)),
        vartype=gsub("_averted","",variable)) %>% group_by(source,vartype,price,country_iso,intervention) %>% 
        mutate(norm_mean=mean/mean[!grepl("averted",variable)], norm_median=median/median[!grepl("averted",variable)],
            norm_CI50_low=CI50_low/mean[!grepl("averted",variable)],norm_CI50_high=CI50_high/mean[!grepl("averted",variable)],
  norm_CI95_low=CI95_low/mean[!grepl("averted",variable)],norm_CI95_high=CI95_high/mean[!grepl("averted",variable)]) %>% ungroup() %>% 
      mutate(vec=as.character((10^floor(log10(median/norm_median)))*round(median/norm_median/(10^floor(log10(median/norm_median))),3)))
df_plot$orig_burden_round=gsub("^\\.","",
        sapply(df_plot$vec, function(vec) paste0(substring(vec,first=c(1,seq(nchar(vec)-floor(nchar(vec)/3)*3,
        nchar(vec)-1,by=3)+1),last=c(seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by =3),nchar(vec))),collapse=".")))
df_plot <- df_plot %>% mutate(orig_burden_round=ifelse(as.numeric(vec)<1e4,gsub("\\.","",orig_burden_round),orig_burden_round))
if (any(df_plot$vartype %in% "total_DALY")) {
  df_plot$vartype=factor(df_plot$vartype,levels=
        unique(df_plot$vartype)[c(which(!grepl("total_DALY",unique(df_plot$vartype))),which(grepl("total_DALY",unique(df_plot$vartype))))])}
# plot for 2 cntrs, 2 intervents
dodge_val=1; round_val=2; caption_txt=paste0("Numbers above medians are pre-intervention",ifelse(any(grepl("YLL",sel_vars))," DALYs",
                            ifelse(any(grepl("death",sel_vars))," case/death numbers"," costs (USD)") ),
                            ifelse(any(grepl("YLL",sel_vars)),". YLD=years lived with disability. YLL=years of life lost",""))
ylab_txt=paste0("% reduction in ",ifelse(any(grepl("YLL",sel_vars)),"DALYs",
            ifelse(any(grepl("death",sel_vars)),"cases/deaths","cost"))," (mean, CI50)")
# plot with cntr on x-axis, MV/mAb as colors
ggplot(df_plot %>% filter(grepl("averted",burden_interv)) ) +
  geom_hpline(aes(x=country_iso,y=norm_median*1e2,group=intervention,color=intervention), # ,linetype=source
              position=position_dodge(width=dodge_val),width=0.42,size=1) + # scale_linetype_manual(values=c("solid","longdash"))+
  geom_linerange(aes(x=country_iso,ymin=norm_CI50_low*1e2,ymax=norm_CI50_high*1e2,group=,color=intervention),
                 alpha=0.35,position=position_dodge(width=dodge_val),size=30,show.legend=F) +
  facet_wrap(~vartype) + scale_color_manual(values=c("red","blue")) +
  geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=norm_CI50_high*1e2+2,group=intervention,label=ifelse(intervention!="MV",orig_burden_round,"")),
   position=position_dodge(width=dodge_val),size=5) + labs(color="",linetype="",caption=caption_txt) +
  scale_x_discrete(expand=expansion(0.02,0)) + # scale_y_continuous(breaks=(0:10)*10) + 
  theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=15),
        strip.text=element_text(size=14),legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=18))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,ifelse(any(grepl("YLL",sel_vars)),"DALY",
    ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),width=36,height=18,units="cm")
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables  # "total_medical_cost_averted"
df_plot <- cea_summary_all %>% filter(variable %in% c("intervention_cost","incremental_cost", "incremental_cost/DALY_averted") &
                    grepl("new",source)) %>% mutate(intervention=ifelse(intervention=="maternal","MV",intervention),
        vec=as.character(abs(round((10^floor(log10(abs(median))))*round(median/(10^floor(log10(abs(median)))),3)))),
        price_interv=factor(paste0(price,"$ (",intervention,")"),levels=unique(paste0(df_plot$price,"$ (",df_plot$intervention,")"))),
 variable=factor(variable,levels=c("intervention_cost","incremental_cost", "incremental_cost/DALY_averted")))
df_plot$orig_burden_round=df_plot$vec; df_plot$orig_burden_round[as.numeric(df_plot$vec)>1e4]=sapply(df_plot$vec[as.numeric(df_plot$vec)>1e4], 
  function(vec) paste0(substring(vec, first=c(1,seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by=3)+1),
                                      last=c(seq(nchar(vec)-floor(nchar(vec)/3)*3,nchar(vec)-1,by=3),nchar(vec))),collapse="."))
df_plot <- df_plot %>% dplyr::select(!vec) %>% mutate(orig_burden_round=gsub("^\\.","",orig_burden_round)) %>% 
  mutate(orig_burden_round=ifelse(median<0,paste0("-",orig_burden_round),orig_burden_round))
# plot
ylab_txt="cost in USD (mean, CI50)"
ggplot(df_plot) + geom_hpline(aes(x=country_iso,y=median,group=intervention,color=price_interv), # 
              position=position_dodge(width=dodge_val),width=0.43,size=1) + # scale_linetype_manual(values=c("solid","longdash"))+
  geom_linerange(aes(x=country_iso,ymin=CI50_low,ymax=CI50_high,group=intervention,color=price_interv),
                 alpha=0.35,position=position_dodge(width=dodge_val),size=28,show.legend=F) +
  facet_wrap(~variable,scales = "free") +
  scale_color_manual(values=c(colorRampPalette(colors=c("rosybrown","red"))(3),colorRampPalette(colors=c("blue","blueviolet"))(3))) +
  geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
  geom_text(aes(x=country_iso,y=ifelse(abs(CI50_high)>1e6,CI50_high+2e6,CI50_high+80),group=intervention,label=orig_burden_round),
            position=position_dodge(width=dodge_val)) + labs(color="",linetype="",caption="Numbers show median values.") + 
  scale_x_discrete(expand=expansion(0.02,0)) + theme(axis.text.x=element_text(angle=0,vjust=1/2,size=14),axis.text.y=element_text(size=15),
      strip.text=element_text(size=14),legend.position="top",legend.text=element_text(size=14),axis.title.y=element_text(size=18)) + 
  guides(color=guide_legend(ncol=2))
# save
ggsave(paste0("output/cea_plots/",subfolder_name,"incremental_costs_KEN_ZAF.png"),width=36,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# calculate with modified function: this function takes user inputs but assumes incidence matrix = all cases (ARI+SARI)
# sim_output_user_input=get_burden_flexible(sel_interv,incid_matrix=burden_list_own_data[[cntrs_cea[n_cntr_output]]][["incid"]],
#                                           hosp_prob_matrix=burden_list_own_data[[cntrs_cea[n_cntr_output]]][["hosp"]])
# this function takes user inputs and takes the incidence matrix as SARIs only (hosp and nonhosp)
# sim_output_user_input_nonhosp_sari=get_burden_flexible_nonhosp_sari(sel_interv,
#                       incid_matrix=burden_list_own_data[[cntrs_cea[n_cntr_output]]]$incid,
#                       hosp_prob_matrix=burden_list_own_data[[cntrs_cea[n_cntr_output]]]$hosp)
# ARIs and SARIs separately accounted for, hosp/nonhosp distinguished
# exp waning
# eff_sympt_disease_fit <- get.gamma.par(q=efficacy_figures$mat_vacc$sympt_disease[c("CI95_low","mean","CI95_high")],
#                                        p=c(2.5,50,97.5)/100,show.output=F,plot=F)[c("shape","rate")]
# effic_sympt_disease <- rgamma(5e3,shape=eff_sympt_disease_fit["shape"],rate=eff_sympt_disease_fit["rate"])
# x <- matrix(0,nrow=60,ncol=5e3); x[1:nrow(x),] <- effic_sympt_disease
# x <- x*exp(-(1/5)*log(2)*(0:59))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot with intervntn type on x-axis
# ggplot(df_plot %>% filter(grepl("averted",burden_interv)) ) +
#   geom_hpline(aes(x=intervention,y=norm_median*1e2,group=country_iso,color=country_iso), # ,linetype=source
#               position=position_dodge(width=dodge_val),width=0.42,size=1) + # scale_linetype_manual(values=c("solid","longdash"))+
#   geom_linerange(aes(x=intervention,ymin=norm_CI50_low*1e2,ymax=norm_CI50_high*1e2,group=country_iso,color=country_iso),
#                  alpha=0.5,position=position_dodge(width=dodge_val),size=30,show.legend=F) +
#   facet_wrap(~vartype) + scale_y_continuous(breaks=(0:10)*10) +
#   geom_vline(xintercept=1.5,linetype="dashed",size=0.3) + theme_bw() + standard_theme + xlab("") + ylab(ylab_txt) + 
# geom_text(aes(x=intervention,y=norm_CI50_high*1e2+2,group=country_iso,label=ifelse(intervention=="MV",orig_burden_round,"")),
# position=position_dodge(width=dodge_val)) + labs(color="",linetype="",caption=caption_txt)+scale_x_discrete(expand=expansion(0.02,0))+
#   theme(axis.text.x=element_text(angle=0,vjust=1/2,size=12),legend.position="top",legend.text=element_text(size=14)) 
# # save
# ggsave(paste0("output/cea_plots/",subfolder_name,ifelse(any(grepl("YLL",sel_vars)),"DALY",
#               ifelse(any(grepl("death",sel_vars)),"case_death","cost")),"_reductions_KEN_ZAF.png"),width=36,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# subtract SARIs from ARIs (don't do this - they're already excluded)
# ARI_SARI_Rates_2010_2018_tidydata = read_csv(kenya_data_file_path)
# tr_vals=ARI_SARI_Rates_2010_2018_tidydata$disease_type == "ARI" & ARI_SARI_Rates_2010_2018_tidydata$RSV_assoc =="yes"  
# ARI_SARI_Rates_2010_2018_tidydata[tr_vals,c("value","CI_95_lower","CI_95_upper")] = 
# ARI_SARI_Rates_2010_2018_tidydata %>% filter(disease_type=="ARI" & RSV_assoc=="yes") %>% 
#     dplyr::select(c(value,CI_95_lower,CI_95_upper)) -
# (ARI_SARI_Rates_2010_2018_tidydata %>% filter(disease_type=="SARI" & RSV_assoc=="yes") %>% 
#     dplyr::select(c(value,CI_95_lower,CI_95_upper)))/100
# ARI_SARI_Rates_2010_2018_tidydata=ARI_SARI_Rates_2010_2018_tidydata %>% mutate(value=ifelse(value<0,0,value),
#                                     CI_95_lower=ifelse(CI_95_lower<0,0,CI_95_lower),CI_95_upper=ifelse(CI_95_upper<0,0,CI_95_upper))
# write_csv(ARI_SARI_Rates_2010_2018_tidydata,gsub(".csv","_ARI_adjusted.csv",kenya_data_file_path))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Compare to age-specific incidence curve generated by mcmarcel
# sel_disease="SARI"
# config=get_rsv_ce_config(subset(sim_config_matrix,country_iso %in% "KEN" & intervention=="maternal"))
# ci95_mcmarcel=t(sapply(1:nrow(config$rsv_rate),function(x)
#   quantile( (config$rsv_rate*ifelse(sel_disease=="SARI",config$hosp_prob,1-config$hosp_prob))[x,],probs = c(2.5,97.5)/100)))
# mcmarcel_agedep_incid=data.frame(age=1:60,disease_type=sel_disease,
#     value=rowMeans(config$rsv_rate*ifelse(sel_disease=="SARI",config$hosp_prob,1-config$hosp_prob) ),ci95_mcmarcel,
#     source="mcmarcel",disease_type_medic_status="all") %>% mutate(`X2.5.`=as.numeric(X2.5.),`X97.5.`=as.numeric(X97.5.))
# colnames(mcmarcel_agedep_incid) = c("age_inf","disease_type","value","CI_95_lower","CI_95_upper","source","disease_type_medic_status")
# mcmarcel_kemri_compare=bind_rows(mcmarcel_agedep_incid,
#   subset(kenya_data,disease_type==sel_disease) %>% group_by(age_inf) %>% summarise(age_inf=unique(age_inf),
#   disease_type=unique(disease_type),value=sum(value/metric_per_popul),CI_95_lower=sum(CI_95_lower/metric_per_popul),
#             CI_95_upper=sum(CI_95_upper/metric_per_popul)) %>% mutate(source="KEMRI",disease_type_medic_status="all"),
#   subset(kenya_data,disease_type==sel_disease & medically_attended==TRUE) %>% mutate(value=value/metric_per_popul,
#     CI_95_lower=CI_95_lower/metric_per_popul,CI_95_upper=CI_95_upper/metric_per_popul,source="KEMRI") %>%
#     dplyr::select(age_inf,disease_type,value,CI_95_lower,CI_95_upper,source,disease_type_medic_status) )
# # plot
# ggplot(mcmarcel_kemri_compare,aes(x=age_inf,y=value))+geom_line(aes(group=disease_type_medic_status,color=disease_type_medic_status))+
#   geom_point(aes(group=disease_type_medic_status,color=disease_type_medic_status)) +
#   geom_ribbon(aes(ymin=CI_95_lower,ymax=CI_95_upper,group=disease_type_medic_status,fill=disease_type_medic_status),alpha=0.3) +
#   facet_wrap(~source,nrow = 2) + theme_bw() + standard_theme + scale_x_continuous(breaks=0:60,expand=expansion(0.02,0)) +
#   scale_y_continuous(breaks=(0:50)/100,expand=expansion(0.02,0)) + xlab("age (months)") + ylab("disease episodes/person-years") +
#   labs(color="data source",fill="data source") + ggtitle(paste0("RSV-associated ",sel_disease," incidence"))
# # save
# ggsave(paste0("output/",sel_disease,"_incidence_comparison_KEN_faceted.png"),width=30,height=18,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Fitting efficacy figures with distributions
# gamma_estims=gamma.parms.from.quantiles(c(efficacy_figures$mat_vacc$sympt_disease["CI95_low"],
#                                           efficacy_figures$mat_vacc$sympt_disease["CI95_high"]))
# severe_gamma_dist=rgamma(5e3,shape=gamma_estims$shape,rate=gamma_estims$rate)
# 
# norm_fits_mat_vacc=lapply(efficacy_figures$mat_vacc, function(x) get.norm.par(q=x[c("CI95_low","CI95_high")],p=c(2.5,97.5)/100))
# norm_fits_mAb=lapply(efficacy_figures$monocl_ab, function(x) get.norm.par(q=x[c("CI95_low","CI95_high")],p=c(2.5,97.5)/100))
# pred_effic_distr_mat_vacc=data.frame(t(data.frame(lapply(norm_fits_mat_vacc, function(x) c("mean"=mean(rnorm(1e4,mean=x["mean"],
#   sd=x["sd"])),quantile(rnorm(1e4,mean=x["mean"],sd=x["sd"]),probs=c(2.5,97.5)/1e2)) )))) %>% mutate(type="pred")
# data.frame(t(data.frame(efficacy_figures$mat_vacc)))
# 
# # quantile(rnorm(1e4,mean=fitnrm["mean"],sd=fitnrm["sd"]),probs=c(2.5,97.5)/1e2)
# 
# fit_gamma=get.gamma.par(q=efficacy_figures$mat_vacc$hospit[c("CI95_low","CI95_high")],p=c(2.5,97.5)/1e2)
# # get.gamma.par(q=efficacy_figures$mat_vacc$severe[c("CI95_low","CI95_high")],p=c(2.5,97.5)/1e2) ---> INVALID (negative)
# get.norm.par(q=efficacy_figures$mat_vacc$severe[c("CI95_low","CI95_high")],p=c(2.5,97.5)/1e2)
# 
# # try different distribs
# fit.perc(p=c(0.025,0.975), q=efficacy_figures$mat_vacc$severe[c("CI95_low","CI95_high")] )
# # cauchy or logistic could work
# fit_norm=get.norm.par(q=efficacy_figures$mat_vacc$severe[c("CI95_low","CI95_high")],p=c(2.5,97.5)/1e2)
# c(mean(rnorm(1e4,mean=fit_norm["mean"],sd=fit_norm["sd"])),
#   quantile(rnorm(1e4,mean=fit_norm["mean"],sd=fit_norm["sd"]),probs=c(2.5,97.5)/1e2))
# # cauchy worse then normal (too high mean)
# # fit_cauch=get.cauchy.par(q=efficacy_figures$mat_vacc$severe[c("CI95_low","CI95_high")],p=c(2.5,97.5)/1e2)
# # cauch_sampls=rcauchy(1e4,location=fit_cauch["location"],scale=fit_cauch["scale"])
# # c(mean(cauch_sampls),quantile(cauch_sampls,probs=c(2.5,97.5)/1e2))
#   efficacy_figures$mat_vacc$severe