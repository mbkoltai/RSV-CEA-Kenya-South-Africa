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
kenya_data_file_path <- "custom_input/Kenya_ARI_SARI_Rates_2010_2018_tidydata_updated_2021_08.csv"
# "custom_input/ARI_SARI_Rates_2010_2018_tidydata.csv"
### PLOT Kenya incidence data with error bars
kenya_ari_sari_incidence <- bind_rows(fcn_load_kenya(kenya_data_path=kenya_data_file_path,sel_disease="ARI")$rsv_incidence_ageinf,
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
    metric_per_popul=unique(metric_per_popul))
# PLOT
plot_flag=F
if (plot_flag){
ggplot(kenya_ari_sari_incidence, aes(x=age_in_months)) + 
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
ggsave("output/cea_plots/kenya_ari_sari_burden_errorbars_grouped_updated_2021_08.png",
        width=42,height=22,units="cm")
}
# LOAD data, fit (gamma) distrib to CI95 values, generate matrix with 5e3 columns, age groups from 0 to 59mts
# Kenya
ci50_range <- c(25,75)/1e2; ci95_range <- c(2.5,97.5)/1e2
kenya_nonhosp_hosp_incid_ari_sari <- lapply(c("ARI","SARI"), function(x)
  fcn_gen_nonhosp_hosp_incid_samples_kenya(kenya_data_file_path,sel_disease=x,n_iter=5e3,age_maxval=60,
            CI_intervals=ci95_range,randsampl_distrib_type="gamma"))
names(kenya_nonhosp_hosp_incid_ari_sari)=c("ARI","SARI")
### deaths
# Kenya deaths
deaths_kenya <- read_csv("custom_input/deaths_kenya_tidy_adjusted_02_2022.csv") %>% 
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
# PLOT DEATHS
deaths_data <- bind_rows(deaths_SA %>% select(!age_inf) %>% group_by(age_in_months_orig,in_hospital) %>%
            summarise(value=unique(value),CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),
            country="South Africa") %>% ungroup() %>% 
            mutate(age_in_months_orig=as.character(age_in_months_orig)),
        deaths_kenya %>% group_by(age_in_months_orig,in_hospital) %>% summarise(value=unique(value),
            CI_95_lower=unique(CI_95_lower),CI_95_upper=unique(CI_95_upper),country="Kenya")) %>%
        mutate(age_in_months_orig=factor(age_in_months_orig,levels=unique(age_in_months_orig))) %>% 
        group_by(country,age_in_months_orig) %>% 
        mutate(CI_95_lower_sum=sum(CI_95_lower),CI_95_upper_sum=sum(CI_95_upper)) %>% 
  mutate(in_hospital=ifelse(in_hospital=="yes","in-hospital","out-of-hospital"))
# 
plot_flag=F
if (plot_flag){
p <- ggplot(deaths_data,aes(x=age_in_months_orig)) + 
  geom_bar(aes(y=value,fill=in_hospital),stat="identity") + # ,position="dodge"
  geom_errorbar(aes(ymin=CI_95_lower_sum,ymax=CI_95_upper_sum,group=in_hospital),size=0.2) + # ,position="dodge"
  facet_wrap(~country,scales="free_y") + # ,nrow=2
  xlab("age (months)") + ylab("deaths per 100,000 person year") + labs(fill="") +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:12)*25) + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(vjust=0.5,size=13),axis.text.y=element_text(size=13),
    legend.text=element_text(size=13),legend.background=element_rect(fill=NA),strip.text=element_text(size=14),
    legend.position=c(0.88,0.925),axis.title.x=element_text(size=17),axis.title.y=element_text(size=17)); p
# save
# ggsave("output/cea_plots/ALL_deaths_data_dodged_2rows.png",width=32,height=18,units="cm") #  # _yfixed
# ggsave("output/cea_plots/ALL_deaths_data_stacked_2rows.png",width=32,height=18,units="cm") #  # _yfixed
# ggsave("output/cea_plots/ALL_deaths_data_stacked.png",width=32,height=18,units="cm") #  # _yfixed
}
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
SA_ILI_data <- read_csv("custom_input/s_afr_ILI_incidence_rate_160921.csv") %>% 
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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# PLOT South Africa incidence
SA_ILI_SARI_rawdata <- bind_rows( 
  distinct(fcn_load_s_afr(safr_data_path="custom_input/s_afr_incidence_data_rate.csv") %>%
  mutate(disease_type_medic_status=paste(disease_type,
                     ifelse(hospitalisation,"hospitalised","not hospitalised"))) %>% select(!age_inf)),
read_csv("custom_input/s_afr_ILI_incidence_rate_160921.csv") %>% 
  filter(!(grepl("<",agegroup) | agegroup %in% c("0-5m","6-11m","12-23m","24-59m","<5y"))) %>%
  rename(age=agegroup) %>% mutate(disease_type=ifelse(disease_type=="ILI","ARI",""),
         disease_type_medic_status=ifelse(hospitalisation,paste0("medically attended ",disease_type),
                                          paste0("non medically attended ",disease_type))) ) %>%
  mutate(age=gsub("m","",age),age=factor(age,levels=unique(age)))
if (plot_flag){
ggplot(SA_ILI_SARI_rawdata,aes(x=age)) + 
  geom_bar(aes(y=rate/popul_denom,fill=disease_type_medic_status),position="stack",stat="identity") +
  geom_errorbar(aes(ymin=rate_CI_lower/popul_denom,ymax=rate_CI_upper/popul_denom),size=0.4) +
  facet_wrap(~disease_type_medic_status,nrow=2,scales = "free") + # 
  xlab("age (months)") + ylab("cases per person year") + labs(fill="")+
  scale_y_continuous(expand=expansion(0.01,0)) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=13),
                                      axis.text.y=element_text(size=13),legend.text=element_text(size=13),
                                      legend.background=element_rect(fill=NA),legend.position=c(0.9,0.925),
                                      axis.title.x=element_text(size=17),axis.title.y=element_text(size=17),
                                      strip.text=element_text(size=14))
ggsave("output/cea_plots/SA_ari_sari_burden_errorbars_grouped_ILI_160921.png",
       width=42,height=22,units="cm")
}
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

# inflation
inflation_data <- (read_csv("custom_input/inflation_data.csv") %>% filter(grepl("South Africa",`Country Name`)))
inflation_rate <- 1+as.numeric(inflation_data[,(ncol(inflation_data)-5):ncol(inflation_data)])/100
# ratio of exchange rates 2014 vs 2021
exch_rate_adj <- 11.2/14.77
# adjustment of 2014 prices
hist_adj <- prod(inflation_rate)*exch_rate_adj
# outpatient
s_afr_outpatient_cost <- bind_rows(data.frame(age="0-59",mean=25,LCI=18.3,UCI=31.8,
                                          cost_type="healthcare",disease="outpatient",name="total"),
      read_csv("custom_input/s_afr_PDE_calcs.csv") %>% filter(grepl("outpatient",disease)) %>%  
  rename(mean=`Mean cost per illness episode (USD)`)) %>% 
  # adjustment for inflation and exch rate change
  mutate(mean=mean*hist_adj,LCI=LCI*hist_adj,UCI=UCI*hist_adj)
# fitting with gamma distribution
if (!any(grepl("shape",colnames(s_afr_outpatient_cost)))) {
  # colbind to original dataframe
for (k_row in 1:nrow(s_afr_outpatient_cost)) {
  gamma_fit <- c(get.gamma.par(p=ci95_range,q=c(s_afr_outpatient_cost$LCI[k_row],s_afr_outpatient_cost$UCI[k_row]),
                             show.output=F,plot=F),scaling=1)
  if (any(is.na(gamma_fit))) { gamma_fit <- c(get.gamma.par(p=ci95_range,
                      q=c(s_afr_outpatient_cost$LCI[k_row],s_afr_outpatient_cost$UCI[k_row])/100,
                                                show.output=F,plot=F),scaling=100) }
  sim_gamma <- rgamma(n=1e4,shape=gamma_fit["shape"],rate=gamma_fit["rate"])*gamma_fit["scaling"]
  gamma_fit <- c(gamma_fit,sim_mean=mean(sim_gamma),
                     sim_ci95_low=as.numeric(quantile(sim_gamma,probs=ci95_range[1])),
                     sim_ci95_up=as.numeric(quantile(sim_gamma,probs=ci95_range[2])))
  if (k_row==1) { SA_outpatient_fit_pars <- gamma_fit } else {
    SA_outpatient_fit_pars <- bind_rows(SA_outpatient_fit_pars,gamma_fit)} 
}
# bind together
  s_afr_outpatient_cost <- bind_cols(s_afr_outpatient_cost,SA_outpatient_fit_pars) %>% 
  mutate(freq=ifelse(grepl('-',age),
                     as.numeric(sapply(age, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,1)) %>%
  mutate(age=ifelse(grepl('-',age), sapply(strsplit(age,'-'),'[[',1),age)) %>%
  uncount(weights=freq, .id="n",.remove=F) %>% # dplyr::
  mutate(age=as.numeric(age)+(n-1)) %>% select(!c(n,freq)) 
}
# types: s_afr_outpatient_cost %>% select(c(name,cost_type,disease)) %>% distinct()
# ggplot(s_afr_outpatient_cost) + geom_point(aes(x=mean,y=sim_mean,color=name,fill=cost_type),shape=21,size=3) +
#   theme_bw()+ scale_x_log10() + scale_y_log10()

# fit inpatient costs by gamma distributions
# inpatient
s_afr_inpatient_cost <- read_csv("custom_input/s_afr_PDE_calcs.csv") %>% 
  filter(!name %in% "PDE" & grepl("inpatient",disease)) %>% rename(mean=`Mean cost per illness episode (USD)`) %>%
  mutate(freq=ifelse(grepl('-',age),
                     as.numeric(sapply(age, function(x) diff(as.numeric(unlist(strsplit(x,"-"))))))+1,1)) %>% 
  mutate(age=ifelse(grepl('-',age), sapply(strsplit(age,'-'),'[[',1),age)) %>% 
  uncount(weights=freq, .id="n",.remove=F) %>% # dplyr::
  mutate(age=as.numeric(age)+(n-1)) %>% select(!c(n,freq)) %>%
 mutate(mean=mean*hist_adj,LCI=LCI*hist_adj,UCI=UCI*hist_adj) # adjustment for inflation and exch rate change
###
if (!any(grepl("shape",colnames(s_afr_inpatient_cost)))){
  sa_costs_unique <- s_afr_inpatient_cost %>% select(c(name,mean,LCI,UCI,cost_type,disease)) %>% distinct()
  for (k_row in 1:nrow(sa_costs_unique)) {
    g_fit_data <- c(ifelse(sa_costs_unique$LCI[k_row]==0,0.01,sa_costs_unique$LCI[k_row]),
                    sa_costs_unique$UCI[k_row])
    fit_gamma <- c(get.gamma.par(q=g_fit_data,p=ci95_range,show.output=F,plot=F)[c("shape","rate")],scaling=1)
    if (any(is.na(fit_gamma))) {
      fit_gamma <- c(get.gamma.par(q=g_fit_data/100,p=ci95_range,show.output=F,plot=F)[c("shape","rate")],scaling=100)
      }
    sim_gam<-rgamma(n=1e4,rate=fit_gamma["rate"],shape=fit_gamma["shape"])*fit_gamma["scaling"]
    gamma_out<-c(n=k_row,fit_gamma,sim_mean=mean(sim_gam),
                 sim_ci95_low=as.numeric(quantile(sim_gam,probs=ci95_range[1])),
                 sim_ci95_up=as.numeric(quantile(sim_gam,probs=ci95_range[2])))
    if (k_row==1) { gamma_fits <- gamma_out  } else { gamma_fits <- bind_rows(gamma_fits,gamma_out)}
  }
  s_afr_inpatient_cost <- left_join(s_afr_inpatient_cost,
            left_join(sa_costs_unique %>% mutate(n=row_number()),gamma_fits,by="n") %>% select(!n),
            by=c("name","mean","LCI","UCI","cost_type","disease"))
}

list_SA_costs <- list("inpatient"=s_afr_inpatient_cost,"outpatient"=s_afr_outpatient_cost)
# cost data types:
# list_SA_costs$inpatient %>% select(c(name,disease,cost_type)) %>% distinct()
# list_SA_costs$outpatient %>% select(c(name,disease,cost_type)) %>% distinct()

# check fit
# ggplot(s_afr_inpatient_cost,aes(x=age)) + geom_line(aes(y=mean)) + geom_line(aes(y=sim_mean),color="red") + 
#   geom_point(aes(y=mean),shape=21) + geom_point(aes(y=sim_mean),shape=21,color="red") + 
#   geom_ribbon(aes(ymin=LCI,ymax=UCI),alpha=0.2) + 
#   geom_ribbon(aes(ymin=sim_ci95_low,ymax=sim_ci95_up),fill="red",alpha=0.2) + 
#   facet_grid(disease~name~cost_type,scales="free_y") + theme_bw() + standard_theme
### ### ### ### ### ### ### ### ### ### ### ### ###
# KENYA costs
kenya_costs <- read_csv("custom_input/kenya_costing_tables_tidy.csv")
# using inpatient/outpatient ratio in South Africa OR study from Malawi
in_outpatient_cost_ratio <- 45.37/9.26
  # median(s_afr_inpatient_cost$mean[s_afr_inpatient_cost$name %in% "total"])/
  #   median(s_afr_outpatient_cost$mean[s_afr_outpatient_cost$cost_type %in% "healthcare"])
# we assume that inpatient/outpatient cost ratio is the same as in SA, and total cost
# (# inpatients)*outpatient_cost*SA_in_outpatient_ratio + (# outpatients)*outpatient_cost = total cost
# outpatient_cost = (total cost)/[(# inpatients)*SA_in_outpatient_ratio + (# outpatients)]
KEN_outpatient_cost <- kenya_costs$mean[kenya_costs$variable %in% "Total healthcare cost"]/
  (kenya_costs$mean[grepl("Total number of inpatients",kenya_costs$variable)]*in_outpatient_cost_ratio + 
     kenya_costs$mean[grepl("Total number of outpatients",kenya_costs$variable)])
# assemble list
ken_inpatient_rows <- kenya_costs %>% filter(grepl("Siaya",site) & grepl("Total patient",variable))
list_KEN_costs <- list(inpatient_household=bind_cols(
  kenya_costs %>% filter(grepl("Siaya",site) & grepl("Total patient",variable)), 
  t(sapply(1:nrow(ken_inpatient_rows), function(x)
  unlist(get.gamma.par(q=c(ken_inpatient_rows$ci95_low[x],ken_inpatient_rows$median[x],
            ken_inpatient_rows$ci95_up[x])/100,p=c(2.5,50,97.5)/100,plot=F)[c("shape","rate")]))),scaling=100),
  inpatient_healthcare_system=KEN_outpatient_cost*in_outpatient_cost_ratio,outpatient_cost=KEN_outpatient_cost)
# gamma can fit median and CIs well, but not so much the mean
# gamma_pars_cost <- get.gamma.par(q=c(ken_inpatient_rows$ci95_low,ken_inpatient_rows$median,
#   ken_inpatient_rows$ci95_up)/100,p=c(2.5,50,97.5)/100,plot=F)
# gamma_costs_sim <- rgamma(n=1e4,shape=gamma_pars_cost["shape"],rate=gamma_pars_cost["rate"])*100
# cbind(source=c("fit","data"),
#   bind_rows(c(mean=mean(gamma_costs_sim),ci95_low=as.numeric(quantile(gamma_costs_sim,probs=c(2.5)/100)),
#   ci95_up=as.numeric(quantile(gamma_costs_sim,probs=c(97.5)/100))),
#   ken_inpatient_rows[,c("mean","ci95_low","ci95_up")]))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cntrs_cea=c("KEN","ZAF")
# efficacy figures for vaccine for RESVAX (Novavax trial)
# for mAb from NIRSEVIMAB (from https://www.nejm.org/doi/full/10.1056/nejmoa1913556)
#
# if this flag is set to TRUE, then using published efficacy data. if FALSE --> interim results
flag_publ_effic <- FALSE
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
                         monocl_ab=list(sympt_disease=c(mean=0.701,CI95_low=0.523,CI95_high=0.812),
                                        hospit=c(mean=0.784,CI95_low=0.519,CI95_high=0.903),
                                        half_life=59.3/30,duration=5)) 
# list(sympt_disease=c(mean=0.745,CI95_low=0.496,CI95_high=0.871),
# hospit=c(mean=0.621,CI95_low=-8.6/100,CI95_high=0.868),
# half_life=59.3/30,duration=5)
}

# fitting efficacy figures with a beta distribution
source("functions/fit_efficacy.R")
g(list_effic_betafit,allfits) %=% fcn_betafit_efficacy(effic_figs=efficacy_figures,
                                   scan_range_resol_nsample=c(min=-2,max=2,by=1/100,n_sample=2e4),
                                   optim_range_res=c(min=-1,max=2,by=0.04),optim_initguess=c(-0.05,1))
# fit as: beta_fit <- rbeta(n=1e4,shape1=alphaval,shape2=betaval)*scale_val + shift_val
#
# if this is set to true, then the fit is shifted to be more aligned with the mean efficacy rather than
# the CI interval (beta distribution often cannot do both)
flag_adjust_to_mean <- FALSE
if (flag_adjust_to_mean){
mv_severe_err <- c(t(allfits %>% filter(interv=="mat_vacc" & disease=="severe") %>% 
                       select(c(data_mean,data_CI95_low,data_CI95_high)))) -
  c(t(allfits %>% filter(interv=="mat_vacc"&disease=="severe") %>% select(c(mean,CI95_low,CI95_high))))
# adjust
list_effic_betafit$mat_vacc$severe["shift_fit"] <- list_effic_betafit$mat_vacc$severe["shift_fit"] + mv_severe_err[1]
# try how mean changed
sim_beta <- (rbeta(n=1e4,shape1=list_effic_betafit$mat_vacc$severe["shape1"],
              shape2=list_effic_betafit$mat_vacc$severe["shape2"])*list_effic_betafit$mat_vacc$severe["scale_fit"])+
  list_effic_betafit$mat_vacc$severe["shift_fit"]
c(mean(sim_beta),quantile(sim_beta,probs = ci95_range))
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# prices for doses
pricelist=list("mat_vacc"=c(3,10,30),"mAb"=c(6,20,60))
# loop thru: cntrs * interventions * dose prices
# n_cntr_output=1:length(cntrs_cea); n_interv=1:2
par_table <- expand_grid(n_cntr_output=1:length(cntrs_cea),n_interv=1:2); read_calc_flag=c("calc","read")[1]
kenya_deaths_input=TRUE; SA_deaths_input=TRUE
# exponential waning model used for efficacy
exp_wane_val <- FALSE
# parameters for exponential waning model
g(list_exp_waning_param,df_exp_waning_param) %=% fcn_exp_waning_rate(efficacy_figures,n_row=60)
# distribution used to fit efficacy figures
effic_dist_fit <- "beta" # 
lower_cov=FALSE # ; lower_cov_val=0.7
subfolder_name <- paste0("new_price_efficacy_",ifelse(kenya_deaths_input,"KENdeaths",""),
                       ifelse(SA_deaths_input,"_SAdeaths",""),
        "_CIs_SA_ILI_",ifelse(ILI_adjust_SA,"broader","narrow"),ifelse(exp_wane_val,"_expwaning",""),
        # ifelse(FALSE,paste0("_coverage",lower_cov_val),""),
        ifelse(flag_adjust_to_mean,"_adj_mv_sev_eff_mean",""),
        ifelse(grepl("gamma",effic_dist_fit),"","_effic_betafit"),
        ifelse(flag_publ_effic,"","_interim"),"/") 
# ifelse(min(unlist(efficacy_figures))<=0,"_nonposit_effic","")
### before starting loop need to create temp folder
# source("init_cea_calc_parallel.R")
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
        if (cntrs_cea[n_cntr_output]=="ZAF") { cost_input <- list_SA_costs } else { cost_input<-list_KEN_costs }
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
    icercolname="incremental_cost/DALY_averted",
    icercols=c("incremental_cost","total_DALY_disc_averted")) %>%
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
# add column of country
cea_summary_all <- cea_summary_all %>% 
                  mutate(plot_variable=NA,country_plot=ifelse(country_iso=="KEN","Kenya","South Africa"))
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Figure combining burden and costs and their relative reduction
save_flag=FALSE # TRUE
source("functions/fig_3_4_5.R")
# combine figures
plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],
          nrow=3,rel_heights=c(1.5,1,1.1),labels="auto",label_size=19)
# SAVE
ggsave(paste0("output/cea_plots/",subfolder_name,"combined_fig3_4_5.png"),width=35,height=40,units="cm")
# save table
df_combined_fig3_4_5[,c("mean","median","CI50_low","CI50_high","CI95_low","CI95_high")] <- 
  round(df_combined_fig3_4_5 %>% select(c(mean,median,CI50_low,CI50_high,CI95_low,CI95_high)))
write_csv(df_combined_fig3_4_5 %>% mutate(across(where(is.numeric),round,3)),
          paste0("output/cea_plots/",subfolder_name,"combined_fig3_4_5.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables (discounted DALYs or not)

source("functions/icer_plots.R")
write_csv(df_interv_incremcosts_icer %>% mutate(across(where(is.numeric),round,1)),paste0("output/cea_plots/",
                                            subfolder_name,"df_interv_incremcosts_icer.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plots COMPARING CEA results with projected data vs new data
# subfolder_name<-"new_price_efficacy_kenyadeaths_CIs/"

source("functions/compar_plots_loop.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot of cost-effectiveness variables

source("functions/compar_plots_cost_red.R")
write_csv(df_total_DALY_medcost_averted_KEN_ZAF,
          paste0("output/cea_plots/",subfolder_name,"comparisons/total_DALY_medcost_averted_KEN_ZAF.csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# plot incremental costs

source("functions/compar_plots_increm_cost.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ICER plot

source("functions/icer_comp_plot.R"); p
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
