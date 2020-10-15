standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))
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
####
### PROCESS RESULTS
# add net cost per DALY averted
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
  theme_bw() + standard_theme + xlab('age (months)') + ylab('incidence') +
  ggtitle('Kenya RSV incidence/100k: MCMARCEL prediction vs KEMRI data') + 
  labs(color='data source',fill='data source') + guides(size=FALSE,linetype=FALSE,color=FALSE)
# 
ggsave("output/kemri_mcmarcel_compar.png",width=24,height=18,units="cm") 

#############################
# S Afr data
s_afr_incidence_data=read_csv('../path_rsv_data/SRI_RSV_PATHproject09102020_tidy.csv')
summary_categs=c("0-6m","6-11m","12-23m","24-59","<1y","<5y")
s_afr_incidence_data=s_afr_incidence_data[!s_afr_incidence_data$age %in% summary_categs,
                                          !grepl('number',colnames(s_afr_incidence_data))]
s_afr_incidence_data$age=factor(s_afr_incidence_data$age,levels=unique(s_afr_incidence_data$age))
s_afr_incidence_data$data_type=str_replace_all(s_afr_incidence_data$data_type,' \\(non medically','\n\\(non medically')
# remove nonhosp, we have totals
s_afr_incidence_data=s_afr_incidence_data[!grepl('Non-hospitalised',s_afr_incidence_data$data_type),]
# PLOT
ggplot(safr_plot_data,aes(x=age,y=rate,group=data_type,color=data_type)) + geom_line() + 
 geom_ribbon(aes(ymin=rate_CI_lower,ymax=rate_CI_upper,fill=data_type),alpha=0.3,colour=NA) + geom_point(color='black',size=0.2) +
 facet_wrap(~Province) + theme_bw() + standard_theme + xlab('age (months)') + ylab('incidence') +
  ggtitle('South Africa RSV incidence per 100.000 (2011-2016)') + geom_rect(data=subset(safr_plot_data, Province %in% 'South Africa'),fill=NA,
            colour="blue",size=1.25,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
# ggsave("output/south_africa_rsv_incidence.png",width=28,height=18,units="cm") 

# hospit rate
s_afr_incidence_data[,'hospitalisation']=TRUE;s_afr_incidence_data$hospitalisation[s_afr_incidence_data$data_type %in% "total"]=NaN
safr_hosp_rate=s_afr_incidence_data %>% group_by(Province,age,year) %>% summarise(hosp_rate=rate[hospitalisation==1]/rate[is.nan(hospitalisation)])
safr_hosp_rate=safr_hosp_rate[!is.na(safr_hosp_rate$hosp_rate),]
safr_hosp_rate_province_means=safr_hosp_rate %>% group_by(Province) %>% 
  summarise(mean_hosp_rate=mean(hosp_rate),stdev_hosp_rate=sd(hosp_rate))
write_csv(safr_hosp_rate_province_means,'safr_hosp_rate_province_means.csv')
####
# create all age groups from 1 to 60
# s_afr_incidence_data[,'hosp_rate']=NaN; s_afr_incidence_data$hosp_rate=
s_afr_incidence_data$age_inf=gsub('m','',s_afr_incidence_data$age)x
s_afr_incidence_data$freq=1; s_afr_incidence_data$freq[grepl('-',s_afr_incidence_data$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(
    s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)],'-'))),ncol=2,byrow=T))+1
s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)]=
  as.numeric(sapply(strsplit(s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)],'-'),'[[',1))
s_afr_incidence_data=s_afr_incidence_data %>% uncount(weights=freq, .id="n",.remove=F)
s_afr_incidence_data$age_inf=as.numeric(s_afr_incidence_data$age_inf)+(s_afr_incidence_data$n-1)
s_afr_incidence_data=s_afr_incidence_data[,!colnames(s_afr_incidence_data) %in% c("freq","n")]
write_csv(s_afr_incidence_data,'../path_rsv_data/s_afr_incidence_data.csv')
#####
age_maxval=max(s_afr_incidence_data$age_inf)
##
# PLOT
ggplot(s_afr_incidence_data,aes(x=age_inf,y=rate,group=data_type,color=data_type)) + geom_line() + 
  geom_ribbon(aes(ymin=rate_CI_lower,ymax=rate_CI_upper,fill=data_type),alpha=0.3,colour=NA) + 
  geom_point(color='black',size=0.1) + scale_y_log10() + facet_wrap(~Province) +
  geom_text(data=data.frame(safr_hosp_rate_province_means,data_type=unique(s_afr_incidence_data$data_type)[1]),
          aes(x=5,y=1.9e4,label=paste0('hosp. rate=',round(mean_hosp_rate,2))),color='black',hjust=0,size=3.5,show.legend=FALSE) +
  theme_bw() + standard_theme + xlab('age (months)') + ylab('incidence') +
  scale_x_continuous(labels=as.character(seq(0,age_maxval,4)),breaks=seq(0,age_maxval,4)) +
  ggtitle('South Africa RSV incidence per 100.000 (2011-2016)')+labs(color='hospitalisation status',fill='hospitalisation status')+
  geom_rect(data=subset(s_afr_incidence_data, Province %in% 'South Africa'),fill=NA,
            colour="blue",size=1.25,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)
# SAVE
ggsave("output/south_africa_rsv_incidence_age_inf.png",width=28,height=18,units="cm")

# national average, with regions around it
s_afr_incidence_data$nat_aver='provinces'
s_afr_incidence_data$nat_aver[s_afr_incidence_data$Province %in% 'South Africa']='national average'
# plot
yscale_min_max=c(10^floor(log10(min(s_afr_incidence_data$rate))),10^ceiling(log10(max(s_afr_incidence_data$rate))))
yscale_breaks=10^(seq(log10(yscale_min_max)[1],log10(yscale_min_max)[2],0.5))
# ceiling(10^(log10(max(s_afr_incidence_data$rate))%%1))*10^floor(log10(max(s_afr_incidence_data$rate)))
ggplot(s_afr_incidence_data,aes(x=age_inf,y=rate,group=Province,color=nat_aver)) + 
  geom_line(aes(linetype=nat_aver,size=as.numeric(factor(nat_aver)))) + 
  geom_ribbon(aes(ymin=rate_CI_lower,ymax=rate_CI_upper,fill=nat_aver),alpha=0.3,colour=NA) + 
  facet_wrap(~data_type) + # ,scales='free' geom_point(color='black',size=0.1) + 
  theme_bw() + standard_theme + xlab('age (months)') + ylab('incidence') +
  scale_x_continuous(labels=as.character(seq(0,age_maxval,4)),breaks=seq(0,age_maxval,4)) + 
  # scale_y_log10(breaks=yscale_breaks,labels=paste0('1e',log10(yscale_breaks))) +
  scale_linetype_manual(values=c("solid","dashed"))+ scale_fill_manual(values=c("red",NA)) + 
  scale_size(range=c(1,0.5),guide=FALSE) + scale_color_manual(values=c('red','blue')) + 
  ggtitle('South Africa RSV incidence per 100.000 (2011-2016)') +
  labs(fill='CI95 national',color='mean incidence',size='mean incidence',linetype='mean incidence')
# breaks=scales::trans_breaks("log10", function(x) 10^x),
# labels=scales::trans_format("log10", scales::math_format(10^.x))
# ggsave("output/south_africa_rsv_incidence_age_inf_NAT_AVER.png",width=28,height=14,units="cm")
ggsave("output/south_africa_rsv_incidence_age_inf_NAT_AVER_linscale.png",width=28,height=14,units="cm")

#####
# birth rates and other data in: 
all_country_data <- read.table('./input/country_details_gavi72.csv',sep=',',header=T,stringsAsFactors=F)
# for south africa we don't have the birth rates. load total fertility rate
data("tfrprojMed"); data("popproj")
# birth rate 
# https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2019_POP_F03_RATE_OF_NATURAL_INCREASE.xlsx
crude_birthrates=read_csv('input/crude_birth_rate_kenya_safr.csv')
s_afr_crude_birthrates=crude_birthrates[crude_birthrates$`Region, subregion, country or area` %in% 'South Africa',
                 which(colnames(crude_birthrates) %in% '2020-2025'):which(colnames(crude_birthrates) %in% '2045-2050')]
s_afr_crude_birthrates_inf=approx(seq(2020,2050,6),as.numeric(s_afr_crude_birthrates),n=31)[[2]]
s_afr_popproj=popproj[popproj$name %in% 'South Africa',which(colnames(popproj) %in% '2020'):which(colnames(popproj) %in% '2050')]
s_afr_popproj_inf=approx(colnames(s_afr_popproj),as.numeric(s_afr_popproj),n=length(2020:2050))
# we want birthrates for 2020-2050 (every year)
s_afr_births=s_afr_popproj_inf[[2]]*s_afr_crude_birthrates_inf[[2]]
# stillbirths
stillbirth_rate_data=read_csv('input/stillbirth_rate_data.csv')
stillbirth_rate_s_afr=as.numeric(stillbirth_rate_data[stillbirth_rate_data$Country %in% 'South Africa',2])/1000
# lets see if we get the same result for Kenya: 1485098 (2020); we get 1432695, 3% less
s_afr_lifebirths=s_afr_births*(1-stillbirth_rate_s_afr)

# c("country","country_iso3","target_population","year","stillbirth_rate","income_region", "incomplete_maternal_transfer_rate")
s_afr_targetpop_df=data.frame(cbind('South Africa','ZAF',s_afr_lifebirths,2020:2050,stillbirth_rate_s_afr,'LMIC',0))
colnames(s_afr_targetpop_df)=colnames(all_country_data)
write_csv(rbind(all_country_data,s_afr_targetpop_df),'input/country_details_gavi72_expanded.csv')