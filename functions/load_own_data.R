# load own data
### Kenya ------------------
fcn_load_kenya <- function(kenya_data_path,sel_disease){ # ,n_iter,age_maxval
SARI_Rates_2010_2018_tidydata=read_csv(kenya_data_path)
SARI_Rates_2010_2018_tidydata=subset(SARI_Rates_2010_2018_tidydata,disease_type %in% sel_disease & variable %in% "rate")
# remove summary age groups
nonsumm_truthvals=!grepl("<|24-59|12-23",SARI_Rates_2010_2018_tidydata$age_in_months) | 
  grepl("<1$",SARI_Rates_2010_2018_tidydata$age_in_months)
rsv_incidence_ageinf=SARI_Rates_2010_2018_tidydata[nonsumm_truthvals & # don't include summary variables
                                              SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & # national average
                                              SARI_Rates_2010_2018_tidydata$RSV_assoc %in% "yes",] # RSV-associated
rsv_incidence_ageinf$age_in_months=factor(rsv_incidence_ageinf$age_in_months,
                                                      levels=unique(rsv_incidence_ageinf$age_in_months))
rsv_incidence_ageinf$age_inf=as.character(rsv_incidence_ageinf$age_in_months)
rsv_incidence_ageinf$age_inf[rsv_incidence_ageinf$age_inf %in% "<1"]='0'
rsv_incidence_ageinf[,'freq']=1
rsv_incidence_ageinf$freq[grepl('-',rsv_incidence_ageinf$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(rsv_incidence_ageinf$age_inf[grepl('-',rsv_incidence_ageinf$age_inf)],'-'))),
                  ncol=2,byrow=T))+1
rsv_incidence_ageinf$age_inf[grepl('-',rsv_incidence_ageinf$age_inf)]=
  as.numeric(sapply(strsplit(rsv_incidence_ageinf$age_inf[grepl('-',rsv_incidence_ageinf$age_inf)],'-'),'[[',1))
rsv_incidence_ageinf=rsv_incidence_ageinf %>% uncount(weights=freq, .id="n",.remove=F)
rsv_incidence_ageinf$age_inf=as.numeric(rsv_incidence_ageinf$age_inf)+(rsv_incidence_ageinf$n-1)
# name for not medic attended / hospitalised
medattended_hosp_str=unique(rsv_incidence_ageinf$medically_attended[!grepl("non",rsv_incidence_ageinf$medically_attended)])
# print(medattended_hosp_str)
rsv_incidence_ageinf$medically_attended = rsv_incidence_ageinf$medically_attended %in% medattended_hosp_str
# create per capita array
popul_denom=unique(rsv_incidence_ageinf$metric_per_popul)
rsv_incidence_per_persyear=(rsv_incidence_ageinf %>% group_by(age_inf) %>% summarise(value=sum(value)))$value
kemri_rsv_incidence_CIs=rsv_incidence_ageinf %>% group_by(age_inf) %>% 
  summarise(CI_95_lower=sum(CI_95_lower),CI_95_upper=sum(CI_95_upper)) %>% dplyr::select(CI_95_lower,CI_95_upper)
if (max(rsv_incidence_per_persyear)>1){
  rsv_incidence_per_persyear=rsv_incidence_per_persyear/popul_denom; kemri_rsv_incidence_CIs=kemri_rsv_incidence_CIs/popul_denom
  }

list("rsv_incidence_ageinf"=rsv_incidence_ageinf,"rsv_incidence_per_persyear"=rsv_incidence_per_persyear,
     "kemri_rsv_incidence_CIs"=kemri_rsv_incidence_CIs)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### extract data and generate sample paths ------------------
fcn_gen_incid_hospprob_samples_kenya=function(kenya_data_path,sel_disease,n_iter,age_maxval,CI_intervals,
                        logit_hosp_rate_stdev,n_length,stepsize,x_prec,randsampl_distrib_type){
# load data
rsv_incidence_per_persyear=fcn_load_kenya(kenya_data_path,sel_disease)[["rsv_incidence_per_persyear"]]

if (max(rsv_incidence_per_persyear)>1){rsv_incidence_per_persyear=rsv_incidence_per_persyear/popul_denom }
kemri_rsv_incidence_CIs=fcn_load_kenya(kenya_data_path,sel_disease)[["kemri_rsv_incidence_CIs"]]
# what distribution are we assuming? McMarcel can be fit well by gamma distrib # randsampl_distrib_type='gamma'
if (grepl("gamma",randsampl_distrib_type)){
        kemri_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
        for (k in 1:length(rsv_incidence_per_persyear)) {
          if (k==1) {kenya_gamma_estim=data.frame()}
          gammavals=gamma.parms.from.quantiles(q=c(kemri_rsv_incidence_CIs[k,1],kemri_rsv_incidence_CIs[k,2]),p=CI_intervals)
          kenya_gamma_estim=rbind(kenya_gamma_estim,c(gammavals$shape,gammavals$rate))
          kemri_incid_rate_matrix[k,]=rgamma(n_iter,shape=gammavals$shape,rate=gammavals$rate)} } else {
            print("gamma distribution should be used")}
# load data
rsv_incidence_ageinf=fcn_load_kenya(kenya_data_path,sel_disease)[["rsv_incidence_ageinf"]]
kemri_hosp_rate=rsv_incidence_ageinf %>% group_by(age_in_months,period) %>%
  summarise(hosp_rate=value[medically_attended==TRUE]/(value[medically_attended==FALSE]+value[medically_attended==TRUE]))
kemri_hosp_val=mean(kemri_hosp_rate$hosp_rate) # as.numeric(unique(round(array(kemri_hosp_rate$hosp_rate),2)))
# how much variation around this value? distrib of hosp rate in mcmarcel: logit(p) ~ norm
# we generate logit(p) ~ normal distribution with the mean from our data and same stdev as in mcmarcel
print(paste0("hosp rate in data=",round(kemri_hosp_val,3)))
y=fcn_approx_logit_distrib_hosp_rate(n_length,stepsize,x_prec,des_mean=kemri_hosp_val,
                                     mean_guess=log(kemri_hosp_val/(1-kemri_hosp_val)),logit_hosp_rate_stdev,n_iter)
# generate a hosp matrix with 5000 samples from our data
kemri_hosp_rate_matrix=t(sapply(1:age_maxval, function(x) {1/(1+exp(-rnorm(n_iter,mean=y["mean_logit"],sd=y["sd_logit"]))) }))
print(paste("KEN mean hosp rate (approx)=",round(mean(as.numeric(kemri_hosp_rate_matrix[1,])),3) ) )
list(incid=kemri_incid_rate_matrix,hosp=kemri_hosp_rate_matrix)
}
# write_csv(kemri_hosp_rate_matrix,'input/kemri_hosp_rate_matrix.csv')

###
# generate sample paths for nonhosp and hosp incidence
fcn_gen_nonhosp_hosp_incid_samples_kenya <- function(kenya_data_path,sel_disease,n_iter,age_maxval,CI_intervals,randsampl_distrib_type){
  # load data
  hosp_incid_data=subset(fcn_load_kenya(kenya_data_path,sel_disease)$rsv_incidence_ageinf %>% 
   dplyr::select(age_inf,value,CI_95_lower,CI_95_upper,medically_attended,metric_per_popul,medically_attended),medically_attended==TRUE)
  nonhosp_incid_data=subset(fcn_load_kenya(kenya_data_path,sel_disease)$rsv_incidence_ageinf %>% 
   dplyr::select(age_inf,value,CI_95_lower,CI_95_upper,medically_attended,metric_per_popul,medically_attended),medically_attended==FALSE)
  popul_denom=unique(hosp_incid_data$metric_per_popul)
  # message("read in data")
  
  # McMarcel can be fit well by gamma distrib
  if (grepl("gamma",randsampl_distrib_type)){
    hosp_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter); nonhosp_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
    for (k in 1:length(unique(hosp_incid_data$age_inf))) {
      if (k==1) {hosp_gamma_estim=data.frame(); nonhosp_gamma_estim=data.frame()}
      if (hosp_incid_data$value[k]>0){
      hosp_gammavals=gamma.parms.from.quantiles(q=c(hosp_incid_data$CI_95_lower[k],hosp_incid_data$CI_95_upper[k]),p=CI_intervals)
      hosp_gamma_estim=rbind(hosp_gamma_estim,c(hosp_gammavals$shape,hosp_gammavals$rate))
      hosp_incid_rate_matrix[k,]=rgamma(n_iter,shape=hosp_gammavals$shape,rate=hosp_gammavals$rate)} else {
        hosp_incid_rate_matrix[k,]=rep(0,n_iter)  
      } 
      if (nonhosp_incid_data$value[k]>0){
  # message(nonhosp_incid_data$CI_95_lower[k]); message(nonhosp_incid_data$CI_95_upper[k]); message(nonhosp_incid_data$value[k])
  nonhosp_gammavals=gamma.parms.from.quantiles(q=c(nonhosp_incid_data$CI_95_lower[k],nonhosp_incid_data$CI_95_upper[k]),p=CI_intervals)
      nonhosp_gamma_estim=rbind(nonhosp_gamma_estim,c(nonhosp_gammavals$shape,nonhosp_gammavals$rate))
      nonhosp_incid_rate_matrix[k,]=rgamma(n_iter,shape=nonhosp_gammavals$shape,rate=nonhosp_gammavals$rate)} else {
        nonhosp_incid_rate_matrix[k,]=rep(0,n_iter)
      }
    } 
  } else { print("use gamma distribution.") }
  list(nonhosp_incid=nonhosp_incid_rate_matrix/popul_denom,hosp_incid=hosp_incid_rate_matrix/popul_denom)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### South African data --------------------------------
fcn_load_s_afr <- function(safr_data_path){
  s_afr_incidence_data=read_csv(safr_data_path)
  s_afr_nat_average_totalcases=s_afr_incidence_data[s_afr_incidence_data$Province %in% 'South Africa' ,]
  # & s_afr_incidence_data$data_type %in% 'total'
  # the incidence rate in data file is per 100K, normalize to per capita
  s_afr_nat_average_totalcases
}

### generate nonhosp and hosp sample paths --------------------------------
fcn_gen_nonhosp_hosp_incid_samples_SA <- function(safr_data_path,diseasetype,n_iter,age_maxval,CI_intervals,randsampl_distrib_type){
  # load data
  if (is.character(safr_data_path)) {
  all_incid_data=fcn_load_s_afr(safr_data_path)} else {
    all_incid_data=safr_data_path
  }
  popul_denom=unique(all_incid_data$popul_denom)
  # hosp cases
  hosp_incid_data=subset(all_incid_data %>% 
          dplyr::select(age_inf,rate,rate_CI_lower,rate_CI_upper,disease_type,hospitalisation,popul_denom),
          hospitalisation==TRUE & disease_type==diseasetype) %>% 
    mutate(rate=rate/popul_denom,rate_CI_lower=rate_CI_lower/popul_denom,rate_CI_upper=rate_CI_upper/popul_denom)
  # nonhosp cases
  nonhosp_incid_data=subset(all_incid_data %>% dplyr::select(age_inf,rate,rate_CI_lower,rate_CI_upper,disease_type,hospitalisation,popul_denom),
                            hospitalisation==FALSE & disease_type==diseasetype) %>% 
    mutate(rate=rate/popul_denom,rate_CI_lower=rate_CI_lower/popul_denom,rate_CI_upper=rate_CI_upper/popul_denom)
  # McMarcel data can be fit well by gamma distrib
  if (grepl("gamma",randsampl_distrib_type)){
    print("fitting gamma distrib. to data means and CI95s")
    hosp_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter); nonhosp_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
    for (k in 1:length(unique(hosp_incid_data$age_inf))) {
      if (k==1) {hosp_gamma_estim=data.frame(); nonhosp_gamma_estim=data.frame()}
      if (hosp_incid_data$rate[k]>0){
      hosp_gammavals=gamma.parms.from.quantiles(q=c(hosp_incid_data$rate_CI_lower[k],hosp_incid_data$rate_CI_upper[k]),p=CI_intervals)
      hosp_gamma_estim=rbind(hosp_gamma_estim,c(hosp_gammavals$shape,hosp_gammavals$rate))
      hosp_incid_rate_matrix[k,]=rgamma(n_iter,shape=hosp_gammavals$shape,rate=hosp_gammavals$rate) } else {
        hosp_incid_rate_matrix[k,]=rep(0,n_iter) }
      # non-hosp
      nonhosp_fitvals=c(nonhosp_incid_data$rate_CI_lower[k],nonhosp_incid_data$rate_CI_upper[k]) # print(c(k,nonhosp_fitvals))
      if (nonhosp_incid_data$rate[k]>0){ 
        if (nonhosp_fitvals[1]==0){
        message(k); message(nonhosp_fitvals[1]); message(nonhosp_incid_data$rate[k]); message(nonhosp_fitvals[2])
        nonhosp_fitvals[1]=1e-6 }
      nonhosp_gammavals=gamma.parms.from.quantiles(q=nonhosp_fitvals,p=CI_intervals)
      nonhosp_gamma_estim=rbind(nonhosp_gamma_estim,c(nonhosp_gammavals$shape,nonhosp_gammavals$rate)) 
      nonhosp_incid_rate_matrix[k,]=rgamma(n_iter,shape=nonhosp_gammavals$shape,rate=nonhosp_gammavals$rate) } else {
        nonhosp_incid_rate_matrix[k,]=rep(0,n_iter)  }
    } 
  } else {  print("use gamma distribution.") }
  
  list(nonhosp_incid=nonhosp_incid_rate_matrix,hosp_incid=hosp_incid_rate_matrix) # 
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### load data and generate sample paths -----------------
fcn_gen_incid_hospprob_samples_safr=function(safr_data_path,s_afr_hosp_rates_filepath,n_iter,age_maxval,CI_intervals,
                        logit_hosp_rate_stdev,n_length,stepsize,x_prec,randsampl_distrib_type){

  s_afr_nat_average_totalcases=fcn_load_s_afr(safr_data_path)
  popul_denom=s_afr_nat_average_totalcases$popul_denom  
  if (max(s_afr_nat_average_totalcases$rate)>1) {
    s_afr_nat_average_totalcases[,"rate","rate_CI_lower","rate_CI_upper"]=
      s_afr_nat_average_totalcases[,"rate","rate_CI_lower","rate_CI_upper"]/popul_denom}
  
# we need a 60x5000 matrix of per capita RSV cases - need to pick a distribution to generate random samples
# normal distribution, assuming n=1

if (randsampl_distrib_type %in% 'gamma'){
  # inferring gamma distrib from CIs
  s_afr_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
        for (k in 1:nrow(s_afr_nat_average_totalcases)) {  if (k==1) {s_afr_gamma_estim=data.frame()}
                gammavals=gamma.parms.from.quantiles(q=c(s_afr_nat_average_totalcases$rate_CI_lower[k],
                                s_afr_nat_average_totalcases$rate_CI_upper[k]),p=CI_intervals)
                                s_afr_gamma_estim=rbind(s_afr_gamma_estim,c(gammavals$shape,gammavals$rate))
                                s_afr_incid_rate_matrix[k,]=rgamma(n_iter,shape=gammavals$shape,rate=gammavals$rate)  } } else {
                                  print("use gamma distrib!") }
# hosp rate | we don't have a stdev value, we'll use MCMARCEL again
safr_hosp_rate_province_means=read_csv(s_afr_hosp_rates_filepath)
safr_mean_hosp=safr_hosp_rate_province_means$mean_hosp_rate[safr_hosp_rate_province_means$Province %in% 'South Africa']
print(paste0("SAFR mean hosp rate (data)=",round(safr_mean_hosp,3)))

y=fcn_approx_logit_distrib_hosp_rate(n_length,stepsize,x_prec,des_mean=safr_mean_hosp,
                                     mean_guess=log(safr_mean_hosp/(1-safr_mean_hosp)),logit_hosp_rate_stdev,n_iter)
# generate a hosp matrix with 5000 samples from our data
s_afr_hosp_rate_matrix=t(sapply(1:age_maxval, function(x) {1/(1+exp(-rnorm(n_iter,mean=y["mean_logit"],sd=y["sd_logit"]))) }))
print(paste0("SAFR mean hosp rate (approx)=",round(mean(as.numeric(s_afr_hosp_rate_matrix[1,])),3) ) )

list(incid=s_afr_incid_rate_matrix,hosp=s_afr_hosp_rate_matrix)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### find distrib parameters for log(p_hosp/(1-p_hosp)) ~ normal distrib of hosp probability
fcn_approx_logit_distrib_hosp_rate <- function(n_length,stepsize,x_prec,des_mean,mean_guess,sd_guess,n_iter){
for (k in 1:n_length){
  if (k==1) {mean_guess_logit=mean(1/(1+exp(-rnorm(n_iter,mean=mean_guess,sd=sd_guess)))); 
    if (mean_guess_logit<des_mean) {direction=1} else {direction=-1}; 
  print(paste0(paste0(c("initguess=","target=","direction="),round(c(mean_guess_logit,des_mean,direction),3))),collapse=",") }
  mean_rnorm=mean_guess*(1+k*direction*stepsize)
  mean_approx=mean(1/(1+exp(-rnorm(n_iter,mean=mean_rnorm,sd=sd_guess))))
  # if (k %% 50 == 1) {print(c(k,mean_guess_logit,mean_approx))}
  if (abs(mean_approx-des_mean)/des_mean<x_prec) {
    print(paste0("approximation reached after ",k," steps, ",x_prec*1e2,"% accurate")); break}
  if (k==n_length) {print("satisfact. approximation not reached, change <n_length> or <stepsize>")}
}
  c(mean_logit=mean_rnorm,sd_logit=sd_guess)
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### get objects larger than x Mb (memory) --------------
fcn_objs_mem_use <- function(min_size){
  mem_use_df=round(data.frame(unlist(sapply(ls(envir=.GlobalEnv), function(x) object.size(get(x)), simplify = FALSE)))/1e6,1)
  colnames(mem_use_df)[1]<-"size (Mb)"; mem_use_df[,"objs"]=rownames(mem_use_df)
  mem_use_df<-mem_use_df[order(mem_use_df$size,decreasing=T),]; rownames(mem_use_df)<-c()
  mem_use_df[mem_use_df$size>min_size,c(2,1)]
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# concatenate kenya and SAFR data -------------------
fcn_concat_incid_data <- function(disease_categ,incid_data,cntrs){
  incid_table=bind_rows(lapply(1:length(disease_categ),
       function(x){bind_rows(data.frame(age=1:nrow(incid_data[[1]]),name=disease_categ[x],country=cntrs[1],
            mean=rowMeans(incid_data[[x]]), sd=apply(incid_data[[x]],1,sd)) )}))
  nonhosp_incid=incid_table
  incid_table[incid_table$name %in% disease_categ[2],c("mean","sd")]=incid_table[incid_table$name %in% disease_categ[2],"mean"]*
    incid_table[incid_table$name %in% disease_categ[1],c("mean","sd")]
  # print(1-nonhosp_incid[nonhosp_incid$name %in% disease_categ[2],"mean"])
  nonhosp_incid[nonhosp_incid$name %in% disease_categ[2],c("mean","sd")]=(1-nonhosp_incid[nonhosp_incid$name %in% disease_categ[2],"mean"])*
    nonhosp_incid[nonhosp_incid$name %in% disease_categ[1],c("mean","sd")]
  nonhosp_incid=nonhosp_incid[!grepl("total",nonhosp_incid$name),]; nonhosp_incid$name="nonhospitalised"
  
  x=rbind(incid_table,nonhosp_incid); rownames(x)=c()
  x
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#  plot to compare community-based and hospital-based data -------------------
fcn_plot_compare_comm_hosp_data <- function(df_comm_based,df_hosp_based,grouping_type){
p1 <- ggplot(subset(bind_rows(df_comm_based,df_hosp_based),!name %in% "total"),aes(x=age,y=mean,group=name,fill=name)) + 
  geom_area(color="black",size=0.25,position=position_stack(reverse=T)) + labs(fill="") +
  facet_wrap(country~source,ncol=2) + scale_x_continuous(breaks=(1:30)*2-1,expand=expansion(0.003,0)) + theme_bw() + standard_theme + 
  ylab("cases/person-year") + theme(axis.text.x=element_text(vjust=0.75,angle=90)) + scale_y_continuous(expand=expansion(0.01,0))
# PLOT (countries and med categs separate)
p2 <- ggplot(subset(bind_rows(df_comm_based,df_hosp_based),!name %in% "total"),aes(x=age,y=mean,group=source,color=source)) + 
  geom_line() + geom_point() + facet_grid(name~country,scales = "free") + 
  scale_x_continuous(breaks=(1:30)*2-1,expand=expansion(0.003,0)) + scale_y_continuous(expand=expansion(0.01,0)) + 
  theme_bw() + standard_theme + ylab("cases/person-year") + labs(color="") + 
  theme(axis.text.x=element_text(vjust=0.75,angle=90),legend.position="top")
if (grepl("countr|cntr",grouping_type)) {p1} else {p2}
}

# theme for plots
standard_theme=theme(# panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))