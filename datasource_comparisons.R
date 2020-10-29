# !diagnostics off
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))

####### country-specific RSV incidence from Shi2017 ------------------------------------------------
# LMIC incidence by age scaled by data from 'RSV_burden_Shi_2017'

# steps of calculation
# f_country_iso='ZAF'; f_outputFileDir="output/1027174825_RSV_gavi72_basecase_n1000" 
# burden_country <- read.csv('input/RSV_burden_Shi_2017.csv',stringsAsFactors = FALSE);
# burden_country_reference <- burden_country$incidence_RSV_associated_ALRI[burden_country$country_iso==f_country_iso]
# # Life table
# disc_rate=0.03; life_table <- get_life_table(f_country_iso,2015,f_outputFileDir,disc_rate); under5_pop <- life_table$lx[1:60]
# # incidence matrix for all LMICs
# incidence_mat <- convert_pred_into_model_input("./input/incidence_lmic_ts_n5000.csv") # spline_datafiles$incidence
# # cases = incidence per 1000 * population by age
# country_rsv_cases <- (incidence_mat/1000)*under5_pop
# # rsv rate = cases / total cases (per column)
# # take # of cases in an age group and divide by TOTAL number of cases 
# # --> this gives the FRACTION of RSV cases per age group out of ALL RSV cases
# country_rsv_fraction <- country_rsv_cases / rep(colSums(country_rsv_cases),each=dim(country_rsv_cases)[1])
# # country rate = (rsv rate) * (country_reference per 1000) * (country population)
# # (burden_country_reference / 1000) is predicted burden per 1 million -->  x population = predicted total burden
# # country_rsv_rate is the distribution of cases by age groups (as fractions of total) so product distributes total burden
# # this is cases per age group, but # of ppl in age groups not the same!
# country_rsv_pred_cases <- country_rsv_fraction * (burden_country_reference / 1000) * sum(under5_pop)
# # so the colsum is constant, being the total burden
# # rsv rate = rsv cases / population (cases/capita)
# country_rsv_rate=data.frame(country_rsv_pred_cases/under5_pop); age_maxval=nrow(country_rsv_rate)
# collate the tables

# get rates from 'config' list
sim_config_matrix=read.table("./config/RSV_gavi72_basecase.csv",sep=',', dec='.',stringsAsFactors=F,header=T)
sim_config_matrix$num_sim<-num_sim; sim_config_matrix$scenario_id<-1:nrow(sim_config_matrix)
sim_config_matrix$rng_seed<-rng_seed; sim_config_matrix$outputFileDir<-get_output_folder(output_dir)
cntrs_cea=c('KEN','ZAF'); n_cntr_output=2; cntr_sel=cntrs_cea[n_cntr_output]
# append S Afr
if (!cntr_sel %in% sim_config_matrix$country_iso) { 
  df_append=sim_config_matrix[(nrow(sim_config_matrix)-1):nrow(sim_config_matrix),]
  df_append$country_iso=cntr_sel; sim_config_matrix=rbind(sim_config_matrix,df_append) }
burden_cntr_ind=which(sim_config_matrix$country_iso %in% cntrs_cea[n_cntr_output])
# MV=sim_config_matrix[burden_cntr_ind[1],]; mAb=sim_config_matrix[burden_cntr_ind[2],]
# which intervention? 1=MatVacc, 2=monocl Abs
n_interv=1; sel_interv=sim_config_matrix[burden_cntr_ind[n_interv],]; config <- get_rsv_ce_config(sel_interv)
# get life_table by (run in other file):
# life_table_ZAF_2020_2025_disc0p03.RData
# load(paste(f_outputFileDir,"temp",rdata_filename,sep='/'))

# number of cases
cntr_burden_abs=data.frame(age=1:age_maxval,mean=apply(config$rsv_rate,1,mean)*life_table$lx[1:60],
                           sd=apply(config$rsv_rate*life_table$lx[1:60],1,sd),type='number_of_cases')
# no. cases per capita
cntr_burden_percap=data.frame(age=1:age_maxval,mean=apply(config$rsv_rate,1,mean),
                              sd=apply(config$rsv_rate,1,sd),type='case_per_capita')
# hospitalisation rate
hosp_prob_mat_tidy=data.frame(age=1:age_maxval,mean=rowMeans(config$hosp_prob),
                              sd=apply(config$hosp_prob,1,sd),type='hosp_admissions_per_capita')
# hospitalisations per capita
# hosp_percap_mat_tidy=data.frame(age=1:age_maxval,mean=rowMeans(config$hosp_prob)*apply(config$rsv_rate,1,mean),
#                           sd=apply(config$hosp_prob*config$rsv_rate,1,sd),type='hosp_admissions_per_capita')
# cfr in hospital
cfr_mat_tidy=data.frame(age=1:age_maxval,mean=rowMeans(config$hosp_CFR),
                        sd=apply(config$hosp_CFR,1,sd),type='cfr_hosp_admissions')
cntr_burden_all=rbind(cntr_burden_abs,cntr_burden_percap,hosp_prob_mat_tidy,cfr_mat_tidy)
# cntr RSV incidence rates
ggplot(cntr_burden_all,aes(x=age,y=mean)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.3,colour=NA,fill="red") +
  facet_wrap(~type,scales='free',nrow=2) + scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=11,angle=90,vjust=0.5),axis.text.y=element_text(size=11),
                     axis.title=element_text(size=14), text=element_text(family="Calibri")) +
  xlab('Age (month)') + ylab('Total burden') + labs(color="Samples") + ggtitle(paste(f_country_iso,'RSV burden (mean +/- stdev)'))
# ggsave(paste0("output/RSV_incidence_",f_country_iso,"_lmic_comparison.png"),width=30,height=18,units="cm") 

### ### ### ### ### ### ### ### ###
# compare to own data (s_afr_incid_rate_matrix from other file)
own_data_matrix=kemri_incid_rate_matrix # s_afr_incid_rate_matrix
incid_data_compar=rbind(
  data.frame(age=1:60,mean=apply(own_data_matrix,1,mean),sd=apply(own_data_matrix,1,sd),source='own'),
  data.frame(cntr_burden_percap[,1:3],source='mcmarcel'))
ggplot(incid_data_compar,aes(x=age,y=mean,group=source,color=source)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=source),alpha=0.3,colour=NA) +
  scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme_bw() + standard_theme + xlab('Age (month)') + ylab('incidence per capita') + 
  ggtitle(paste0(f_country_iso,' RSV incidence per capita'))
# ggsave(paste0('output/',f_country_iso,'_incidence_mcmarcel_owndata.png'),width=30,height=18,units="cm")

### ### ### ### ### ### ### ### ###
# normalise cntr data to /1000 ppl
lmic_cntr_burden=cntr_burden_all; case_cols=lmic_cntr_burden$type %in% 'number_of_cases'
lmic_cntr_burden=lmic_cntr_burden[!case_cols,]
# scale to 1000
if (max(lmic_cntr_burden$mean)<1){
  # to have incidence per 1000, multiply per capita rate by 1e3
  lmic_cntr_burden[c('mean','sd')]=lmic_cntr_burden[,c('mean','sd')]*1e3
  # to have incidence per 1000, divide case number (which was for 100e3) by 100
  # lmic_cntr_burden[case_cols,c('mean','SD')]=lmic_cntr_burden[case_cols,c('mean','sd')]/1e2
  lmic_cntr_burden$type=str_replace_all(str_replace_all(lmic_cntr_burden$type,'cfr','deaths_per_1000'),'per_capita','per_1000')
  lmic_cntr_burden$type=str_replace_all(lmic_cntr_burden$type,'case','cases')
}
# collate with all LMIC data
lmic_cntr_burden$country=f_country_iso
lmic_incidence_means=data.frame(age=1:age_maxval,mean=rowMeans(incidence_mat),sd=apply(incidence_mat,1,sd),type='cases_per_1000')
lmic_incidence_means$country='lmic'
# rbind two dfs
if (length(unique(lmic_cntr_burden$country))==1){ lmic_cntr_burden=rbind(lmic_cntr_burden,lmic_incidence_means) }

### plot
lmic_cntr_burden$country=factor(lmic_cntr_burden$country,levels=unique(lmic_cntr_burden$country))
ggplot(lmic_cntr_burden,aes(x=age,y=mean,color=country)) + geom_line() + geom_point() +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,fill=country),alpha=0.3,colour=NA) +
  facet_wrap(~type,scales='free',nrow=2) + # scale_color_manual(values=)
  scale_x_continuous(labels=as.character(seq(0,age_maxval,2)),breaks=seq(0,age_maxval,2)) +
  theme_bw() + standard_theme + xlab('Age (month)') + ylab('incidence') + labs(color="Samples") + 
  ggtitle(paste0(f_country_iso,' RSV burden (mean +/- stdev)'))
ggsave(paste0("output/RSV_incidence_",f_country_iso,"_lmic_comparison.png"),width=30,height=18,units="cm") 

####### our own data from KEMRI ------------------------------------------------
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
### MCMARCEL predictions --------------
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
### ### SAVE 
ggsave("output/RSV_burden_kenya_data_pred_compare_datasources.png",width=32,height=20,units="cm") 
# ggsave("output/RSV_burden_kenya_data_pred_compare_nokes2008_2009.pdf",width=30,height=24,units="cm",device=grDevices::pdf) 

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

####### own RSV-SARI incidence data ------------------------------------------------
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
KEMRI_kenya_rsv_incidence_ageinf[,'freq']=1; 
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

### Compare own incid data w MCMARCEL ------
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

####
# kenya incidence matrix
ggplot(melt(kemri_incid_rate_matrix[,1:60]),aes(x=Var1,y=value,group=Var2,color=Var2)) + geom_line()
####### What distribution was used in mcmarcel to generate random samples? ------------------------------------------------
# histogram
# mcmarcel data
# config <- get_rsv_ce_config(sel_interv)
library(fitdistrplus)
safr_rsv_incid_mcmarcel_kemri=melt(data.frame(rbind(data.frame(config$rsv_rate),data.frame(s_afr_incid_rate_matrix)),
          source=array(sapply(c('mcmarcel','kemri'),function(x) {rep(x,age_maxval)})),age_mts=rep(1:age_maxval,2)),
          id.vars=c('age_mts','source'))
# histograms
ggplot(safr_rsv_incid_mcmarcel_kemri[safr_rsv_incid_mcmarcel_kemri$age_mts<=18,],aes(x=value,color=source,group=source)) + 
  geom_freqpoly() + facet_wrap(~age_mts) + theme_bw() + standard_theme + 
  ggtitle('incidence distribution for Shi_et_al vs own data')
# ggsave(paste0('output/',cntr_sel,'_incidence_distrib.png'),width=30,height=18,units="cm")
for (sel_age in 1:age_maxval){
  incid_params=fitdist(as.numeric(config$rsv_rate[sel_age,]),distr='gamma')
  if (sel_age==1) {incid_simul=data.frame(); shape_rate_param_est=data.frame()}
  incid_simul=rbind(incid_simul,data.frame(age_mts=sel_age,source='simul',
                                           value=rgamma(5000,shape=incid_params$estimate['shape'],incid_params$estimate['rate'])))
  shapeval=as.numeric(incid_params$estimate['shape']); rateval=as.numeric(incid_params$estimate['rate'])
  shape_rate_param_est=rbind(shape_rate_param_est, c(sel_age,shapeval,rateval)) }
colnames(shape_rate_param_est)=c('age_mts','shape','rate')
#### plot predicted (gamma distrib) and actual distributions
agelim=24; rowsel=safr_rsv_incid_mcmarcel_kemri$age_mts<=agelim & safr_rsv_incid_mcmarcel_kemri$source %in% 'mcmarcel'
dataplot=rbind(safr_rsv_incid_mcmarcel_kemri[rowsel,c('age_mts','source','value')],incid_simul[incid_simul$age_mts<=agelim,])
ggplot(dataplot,aes(x=value,color=source,linetype=source)) + geom_freqpoly(size=1.2,bins=30) + facet_wrap(~age_mts) + 
  theme_bw() + standard_theme
# ggsave(paste0('output/',f_country_iso,'_incidence_distrib_gamma_fit.png'),width=30,height=18,units="cm")

# what distrib is it in original data? mean vs stdev plot
means_stdevs_incidence=rbind(data.frame(mean=apply(config$rsv_rate,1,mean),
                                        stdev=apply(config$rsv_rate,1,sd),coeffvar=apply(config$rsv_rate,1,sd)/apply(config$rsv_rate,1,mean),age=1:60,source='mcmarcel'),
                             data.frame(mean=apply(s_afr_incid_rate_matrix,1,mean),
                                        stdev=apply(s_afr_incid_rate_matrix,1,sd),coeffvar=apply(config$rsv_rate,1,sd)/apply(config$rsv_rate,1,mean),age=1:60,source='own'))
ggplot(means_stdevs_incidence, aes(x=mean,y=stdev,color=age)) + facet_wrap(~source) + geom_path(size=1.2) + 
  theme_bw() + standard_theme + ggtitle('incidence random samples')
# ggsave(paste0('output/',f_country_iso,'_incid_randomsamples_mean_stdev.png'),width=30,height=18,units="cm")
## mean-stdev plot
ggplot(melt(means_stdevs_incidence,id.vars=c('age','source')),aes(x=age,y=value,group=variable,color=variable)) + 
  facet_wrap(~source) + geom_line() + theme_bw() + standard_theme + ggtitle() 
# ggsave(paste0('output/',f_country_iso,'_incid_randomsamples_mean_stdev_coeffvar.png'),width=30,height=18,units="cm")

# compare actual means & stdevs w those predicted from gamma distrib
means_stdevs_incidence_gammafit=rbind(
  means_stdevs_incidence[means_stdevs_incidence$source %in% 'mcmarcel',c("age","mean","stdev","source")],
  data.frame(age=shape_rate_param_est$age_mts,mean=shape_rate_param_est$shape/shape_rate_param_est$rate,
             stdev=sqrt(shape_rate_param_est$shape)/shape_rate_param_est$rate,source='pred_gamma_distrib'))
# PLOT: predicted and actual means and stdevs
ggplot(melt(means_stdevs_incidence_gammafit,id.vars=c('age','source')),aes(x=age,y=value,color=source,linetype=source)) + 
  geom_line(size=1.2) + facet_wrap(~variable,nrow=2,scales='free') + # geom_point(aes(shape=source),fill=NA) + 
  scale_x_continuous(labels=as.character(seq(0,age_maxval,4)),breaks=seq(0,age_maxval,4)) +
  theme_bw() + standard_theme + xlab('age (months)') + ylab('')
# ggsave(paste0('output/',f_country_iso,'_incid_randomsamples_mean_stdev_fit_gammadistrib.png'),width=30,height=18,units="cm")
### gamma distrib fits the mcmarcel data very well!

# PLOT
ggplot(melt(shape_rate_param_est,id.vars='age_mts'), aes(x=age_mts,y=value,color=variable)) + geom_line() + geom_point() +
  scale_x_continuous(labels=as.character(seq(0,age_maxval,4)),breaks=seq(0,age_maxval,4)) +
  theme_bw() + standard_theme + xlab('age (months)') + ylab('')
# ggsave(paste0('output/',f_country_iso,'_gamma_shape_rate_param.png'),width=30,height=18,units="cm")
# scatterplot
ggplot(shape_rate_param_est, aes(x=shape,y=rate,color=age_mts)) + geom_point() +
  # scale_x_continuous(labels=as.character(seq(0,age_maxval,4)),breaks=seq(0,age_maxval,4)) +
  theme_bw() + standard_theme + xlab('shape') + ylab('rate')
# ggsave(paste0('output/',f_country_iso,'_gamma_shape_rate_param_scatter.png'),width=30,height=18,units="cm")

# can we infer gamma distrib parameters from CI95 values? use 'gamma.parms.from.quantiles()' function
# kemri_rsv_incidence_per_capita; kemri_rsv_incidence_CIs
# s_afr_incid_rate_matrix=sapply(1:n_iter, function(iters) {sapply(1:age_maxval, 
#          function(x) {rnorm(1,mean=s_afr_nat_average_totalcases$rate[x],sd=s_afr_nat_average_totalcases$stdev[x])})})
s_afr_incid_rate_matrix_gammadistr=matrix(NA,nrow=age_maxval,ncol=5e3)
for (k in 1:nrow(s_afr_nat_average_totalcases)) {
  if (k==1) {s_afr_gamma_estim=data.frame()}
  gammavals=gamma.parms.from.quantiles(q=c(s_afr_nat_average_totalcases$rate_CI_lower[k],
                                           s_afr_nat_average_totalcases$rate_CI_upper[k]),p=c(0.025,0.975))
  s_afr_gamma_estim=rbind(s_afr_gamma_estim,c(gammavals$shape,gammavals$rate))
  s_afr_incid_rate_matrix_gammadistr[k,]=rgamma(5e3,shape=gammavals$shape,rate=gammavals$rate)  }
colnames(s_afr_gamma_estim)=c('rate','shape')

# histograms
safr_rsv_incid_mcmarcel_normal_gamma=melt(data.frame(rbind(data.frame(s_afr_incid_rate_matrix_gammadistr),
                                                           data.frame(s_afr_incid_rate_matrix)), 
                                                     source=array(sapply(c('gamma','normal'),function(x) {rep(x,age_maxval)})),
                                                     age_mts=rep(1:age_maxval,2)),id.vars=c('age_mts','source'))
# histograms
ggplot(safr_rsv_incid_mcmarcel_normal_gamma[safr_rsv_incid_mcmarcel_normal_gamma$age_mts<=18,],
       aes(x=value,color=source,group=source,linetype=source)) + geom_freqpoly(size=1.2) + facet_wrap(~age_mts,scales='free') + 
  theme_bw() + standard_theme + ggtitle('incidence distribution normal vs gamma')
# rm(safr_rsv_incid_mcmarcel_normal_gamma)
ggsave(paste0('output/',f_country_iso,'_gamma_normal_simul_histograms.png'),width=30,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ###
# S Afr data -----------------------------
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
s_afr_incidence_data[,'freq']=1; s_afr_incidence_data$freq[grepl('-',s_afr_incidence_data$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(
    s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)],'-'))),ncol=2,byrow=T))+1
s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)]=
  as.numeric(sapply(strsplit(s_afr_incidence_data$age_inf[grepl('-',s_afr_incidence_data$age_inf)],'-'),'[[',1))
s_afr_incidence_data=s_afr_incidence_data %>% uncount(weights=freq, .id="n",.remove=F)
s_afr_incidence_data$age_inf=as.numeric(s_afr_incidence_data$age_inf)+(s_afr_incidence_data$n-1)
s_afr_incidence_data=s_afr_incidence_data[,!colnames(s_afr_incidence_data) %in% c("freq","n")]
write_csv(s_afr_incidence_data,'../path_rsv_data/s_afr_incidence_data.csv')
### ###
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

### ###
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
s_afr_targetpop_df=data.frame(cbind('South Africa','ZAF',s_afr_lifebirths,2020:2050,stillbirth_rate_s_afr,'LMIC',0),stringsAsFactors=F)
colnames(s_afr_targetpop_df)=colnames(all_country_data); s_afr_targetpop_df$target_population=round(as.numeric(s_afr_targetpop_df$target_population))
write_csv(rbind(all_country_data,s_afr_targetpop_df),'input/country_details_gavi72_expanded.csv')

### ### ###
# we don't have inpatient/outpatient cost for south africa
filename_cost_outpatient='./input/cost_data_outpatient.csv'; filename_cost_inpatient='./input/cost_data_inpatient.csv'
# config$sample_outpatient_cost <- get_cost_data(configList$country_iso,config$num_sim, filename_cost_outpatient)
# config$sample_inpatient_cost  <- get_cost_data(configList$country_iso,config$num_sim, filename_cost_inpatient)
# head(read.table(filename_cost_outpatient,sep=','))[,1:11]
cost_data_outpatient=read.table(filename_cost_outpatient,sep=','); cost_data_inpatient=read.table(filename_cost_inpatient,sep=',')
cost_treatment_gavi=data.frame(round(cbind(rowMeans(cost_data_outpatient),rowMeans(cost_data_inpatient), 
                                           apply(cost_data_outpatient,1,sd),apply(cost_data_inpatient,1,sd) ),2))
cost_treatment_gavi=cbind(iso3c=rownames(cost_treatment_gavi),cost_treatment_gavi); rownames(cost_treatment_gavi)=c()
colnames(cost_treatment_gavi)[2:ncol(cost_treatment_gavi)]=c('mean_outpat','mean_hosp','stdev_outpat','stdev_hosp')
# package for econ data
install.packages('WDI')
gdp_per_cap = WDI(indicator='NY.GDP.PCAP.KD', country=c(unique(all_country_data$country_iso3),'ZAF'),start=2010,end=2019,extra=TRUE)
gdp_per_cap_aver_2010_2020=gdp_per_cap[,c('iso3c','country',"NY.GDP.PCAP.KD")] %>% group_by(iso3c) %>% 
  summarise(gdp_per_cap_time_aver=mean(NY.GDP.PCAP.KD))
gdp_per_cap_aver_2010_2020=gdp_per_cap_aver_2010_2020[!is.na(gdp_per_cap_aver_2010_2020$gdp_per_cap_time_aver),]
cost_treatment_gdp=left_join(gdp_per_cap_aver_2010_2020,cost_treatment_gavi,by='iso3c')
# cost_treatment_gdp=melt(cost_treatment_gdp,id.vars='iso3c')
ggplot(cost_treatment_gdp) + geom_point(aes(x=gdp_per_cap_time_aver,y=mean_outpat),color='blue') + 
  geom_point(aes(x=gdp_per_cap_time_aver,y=mean_hosp),color='red') + xlim(c(0,6e3)) + ylim(c(0,650)) # + labs(color='zz')
# 65 and 80% correlation
# linear regression: x <- rnorm(15); y <- x + rnorm(15); z=predict(lm(y ~ x))
predicted_outpat_cost=predict(lm(cost_treatment_gdp$mean_outpat ~ cost_treatment_gdp$gdp_per_cap_time_aver),
                              newdata=cost_treatment_gdp[,c('gdp_per_cap_time_aver','mean_outpat')])
predicted_hosp_cost=predict(lm(cost_treatment_gdp$mean_hosp ~ cost_treatment_gdp$gdp_per_cap_time_aver),
                            newdata=cost_treatment_gdp[,c('gdp_per_cap_time_aver','mean_hosp')])
s_afr_mean_outpat_hosp_pred=c(predicted_outpat_cost[gdp_per_cap_aver_2010_2020$iso3c %in% 'ZAF'],
                              predicted_hosp_cost[gdp_per_cap_aver_2010_2020$iso3c %in% 'ZAF'])
s_afr_std_outpat_hosp_pred=s_afr_mean_outpat_hosp_pred/2
alpha_outpat_hosp_cost_gammadistrib=(s_afr_mean_outpat_hosp_pred/s_afr_std_outpat_hosp_pred)^2 # shape param
beta_outpat_hosp_cost_gammadistrib=s_afr_mean_outpat_hosp_pred/s_afr_std_outpat_hosp_pred^2 # rate param
# generate gamma distrib random samples
outpat_hosp_rand_samples=sapply(1:2, function(k) {rgamma(5e3,shape=alpha_outpat_hosp_cost_gammadistrib[k],
                                                         rate=beta_outpat_hosp_cost_gammadistrib[k])})
# outpat
cost_data_outpatient_expanded=rbind(cost_data_outpatient,outpat_hosp_rand_samples[,1]); 
rownames(cost_data_outpatient_expanded)[nrow(cost_data_outpatient_expanded)]='ZAF'
# inpat
cost_data_inpatient_expanded=rbind(cost_data_inpatient,outpat_hosp_rand_samples[,2]) 
rownames(cost_data_inpatient_expanded)[nrow(cost_data_inpatient_expanded)]='ZAF'
# save
write.table(cost_data_outpatient_expanded,'./input/cost_data_outpatient_expanded.csv',sep = ',')
write.table(cost_data_inpatient_expanded,'./input/cost_data_inpatient_expanded.csv',sep = ',')

### Nyiro 2018 (Nokes) hospitalisation data ---------------------------
nyiro2018_figure4=read_csv("../path_rsv_data/nyiro2018/figure4.csv")
nyiro2018_figure4_tidy=melt(nyiro2018_figure4,id.vars=c("id","month","year","age_months","Agegrp_mths"))
nyiro2018_figure4_tidy=nyiro2018_figure4_tidy[!is.na(nyiro2018_figure4_tidy$value),
                                              !colnames(nyiro2018_figure4_tidy) %in% 'value']
nyiro2018_figure4_tidy$Agegrp_mths=factor(nyiro2018_figure4_tidy$Agegrp_mths,
  levels=c("0-5M","6-11M","12-23M","24-35M","36-59M","5-9Y","10-19Y",">=20Y"))
# cases by age group and virus
nyiro2018_figure4_rates=nyiro2018_figure4_tidy %>% group_by(Agegrp_mths,variable) %>% tally()
nyiro2018_figure4_rates=dcast(nyiro2018_figure4_rates,Agegrp_mths ~ variable)
nyiro2018_figure4_rates[,'rsv_total']=rowSums(nyiro2018_figure4_rates[,c('rsva_positive','rsvb_positive')],na.rm=T)
nyiro2018_figure4_rates[,'n_agegroup']=rowSums(nyiro2018_figure4_rates[,2:ncol(nyiro2018_figure4_rates)],na.rm=T)
standard_theme_mod=standard_theme; standard_theme_mod$axis.text.x$size=16; standard_theme_mod$axis.text.y$size=16
ggplot(nyiro2018_figure4_rates,aes(x=Agegrp_mths,y=rsv_total/n_agegroup,group=1)) + geom_line(colour='blue') + geom_point() +
  theme_bw() + standard_theme_mod + xlab('age group') + ylab('RSV A+B incidence per capita') + 
  ggtitle('Kenya RSV incidence (Nyiro 2018)')
ggsave("output/kenya_nyiro2018_figure4_RSV_incidence.png",width=20,height=15,units="cm") 

### figure 5
nyiro2018_figure5=read_csv("../path_rsv_data/nyiro2018/figure5.csv")
nyiro2018_figure5_tidy=melt(nyiro2018_figure5,id.vars=c("id","month","year","sample_tested","sample_positives" ))
nyiro2018_figure5_tidy=nyiro2018_figure5_tidy[!is.na(nyiro2018_figure5_tidy$value),
                      !colnames(nyiro2018_figure5_tidy) %in% c('value',"sample_tested","sample_positives")]
nyiro2018_figure5_rates=nyiro2018_figure5_tidy %>% group_by(month,variable) %>% tally()
nyiro2018_figure5_rates=dcast(nyiro2018_figure5_rates,month ~ variable)
nyiro2018_figure5_rates[,'rsv_total']=rowSums(nyiro2018_figure5_rates[,c('rsva_positive','rsvb_positive')],na.rm=T)
nyiro2018_figure5_rates[,'n_agegroup']=rowSums(nyiro2018_figure5_rates[,2:ncol(nyiro2018_figure5_rates)],na.rm=T)
### ###
ggplot(nyiro2018_figure5_rates,aes(x=month,y=rsv_total/n_agegroup,group=1)) + geom_line(colour='blue') + geom_point() +
  theme_bw() + standard_theme_mod + xlab('age (months)') + ylab('RSV A+B incidence per capita') + 
  scale_x_continuous(breaks=1:12) + scale_y_continuous(breaks=(0:6)/20) + ggtitle('Kenya RSV incidence (Nyiro 2018)')
ggsave("output/kenya_nyiro2018_figure5_RSV_incidence.png",width=24,height=15,units="cm") 
