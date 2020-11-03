# load own data
### Kenya ------------------
fcn_load_kenya=function(kenya_data_path,n_iter,age_maxval,CI95_const,CI_intervals,
                        kemri_hosp_rate_stdev,randsampl_distrib_type){
SARI_Rates_2010_2018_tidydata=read_csv(kenya_data_path)
# remove summary age groups
nonsumm_truthvals=!grepl("<|24-59|12-23",SARI_Rates_2010_2018_tidydata$age_in_months) | 
  grepl("<1$",SARI_Rates_2010_2018_tidydata$age_in_months)
KEMRI_kenya_rsv_incidence_ageinf=SARI_Rates_2010_2018_tidydata[nonsumm_truthvals & # don't include summary variables
                                                                 SARI_Rates_2010_2018_tidydata$region %in% 'Kenya' & # national average
                                                                 SARI_Rates_2010_2018_tidydata$RSV_association,] # RSV-associated
KEMRI_kenya_rsv_incidence_ageinf$age_in_months=factor(KEMRI_kenya_rsv_incidence_ageinf$age_in_months,
                                                      levels=unique(KEMRI_kenya_rsv_incidence_ageinf$age_in_months))
KEMRI_kenya_rsv_incidence_ageinf$period=factor(KEMRI_kenya_rsv_incidence_ageinf$period,levels=unique(KEMRI_kenya_rsv_incidence_ageinf$period))
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.character(KEMRI_kenya_rsv_incidence_ageinf$age_in_months)
KEMRI_kenya_rsv_incidence_ageinf$age_inf[KEMRI_kenya_rsv_incidence_ageinf$age_inf %in% "<1"]='0'
KEMRI_kenya_rsv_incidence_ageinf[,'freq']=1
KEMRI_kenya_rsv_incidence_ageinf$freq[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  rowDiffs(matrix(as.numeric(unlist(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',
                                    KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'))),ncol=2,byrow=T))+1
KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)]=
  as.numeric(sapply(strsplit(KEMRI_kenya_rsv_incidence_ageinf$age_inf[grepl('-',KEMRI_kenya_rsv_incidence_ageinf$age_inf)],'-'),'[[',1))
KEMRI_kenya_rsv_incidence_ageinf=KEMRI_kenya_rsv_incidence_ageinf %>% uncount(weights=freq, .id="n",.remove=F)
KEMRI_kenya_rsv_incidence_ageinf$age_inf=as.numeric(KEMRI_kenya_rsv_incidence_ageinf$age_inf)+(KEMRI_kenya_rsv_incidence_ageinf$n-1)
# create per capita matrix
kemri_rsv_incidence_per_capita=array(t(KEMRI_kenya_rsv_incidence_ageinf[
  KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'& KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,'rate']))
popul_denom=1e5
if (max(kemri_rsv_incidence_per_capita>1)){kemri_rsv_incidence_per_capita=kemri_rsv_incidence_per_capita/popul_denom }
kemri_rsv_incidence_CIs=KEMRI_kenya_rsv_incidence_ageinf[KEMRI_kenya_rsv_incidence_ageinf$period %in% '2010-2018'&
              KEMRI_kenya_rsv_incidence_ageinf$hospitalisation==FALSE,c('CI_lower','CI_upper')]/popul_denom
# generate random samples from gamma or normal distribution
# create matrix with 5000 iterations; max age in months is 60 (5yrs)
# n_iter=5e3; age_maxval=60; CI95_const=1.96
# what distribution are we assuming? McMarcel can be fit well by gamma distrib
# randsampl_distrib_type='gamma'
if (randsampl_distrib_type %in% 'normal'){
  # CI_lower = mu - 1.96*stdev/sqrt(sample size)
  stdev_est=cbind((kemri_rsv_incidence_per_capita-kemri_rsv_incidence_CIs[,1])/CI95_const,
                  (kemri_rsv_incidence_CIs[,2] - kemri_rsv_incidence_per_capita)/CI95_const)
  kemri_rsv_incidence_stdevs=rowMeans(stdev_est)
  kemri_incid_rate_matrix=sapply(1:n_iter, function(iters) {
    sapply(1:age_maxval, function(x) {
      rnorm(1,mean=kemri_rsv_incidence_per_capita[x],sd=kemri_rsv_incidence_stdevs[x] )})})} else {
        # inferring gamma distrib from CIs
        kemri_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
        for (k in 1:length(kemri_rsv_incidence_per_capita)) {
          if (k==1) {kenya_gamma_estim=data.frame()}
          gammavals=gamma.parms.from.quantiles(q=c(kemri_rsv_incidence_CIs[k,1],
                                                   kemri_rsv_incidence_CIs[k,2]),p=CI_intervals)
          kenya_gamma_estim=rbind(kenya_gamma_estim,c(gammavals$shape,gammavals$rate))
          kemri_incid_rate_matrix[k,]=rgamma(n_iter,shape=gammavals$shape,rate=gammavals$rate)}
      }
# write_csv(kemri_incid_rate_matrix,'input/kemri_incid_rate_matrix.csv')
####
# kemri hospitalisation rate
kemri_hosp_rate=KEMRI_kenya_rsv_incidence_ageinf %>% group_by(age_in_months,period) %>%
  summarise(hosp_rate=rate[hospitalisation==TRUE]/(rate[hospitalisation==FALSE]+rate[hospitalisation==TRUE]))
kemri_hosp_val=as.numeric(unique(round(array(kemri_hosp_rate[kemri_hosp_rate$period %in% '2010-2018','hosp_rate']),2)))
# how much variation around this value? we dont know, lets use values from mcmarcel (there mean is 9%, ours is 24%)
# kemri_hosp_rate_stdev=0.05656925 # value from: unique(apply(config$hosp_prob,1,sd))
# generate a hosp matrix with 5000 samples from our data
kemri_hosp_rate_matrix=sapply(1:n_iter, function(x) {rep(rnorm(1,mean=kemri_hosp_val,
                sd=unique(round(rep(kemri_hosp_rate_stdev,age_maxval),4))),age_maxval)})

list(kemri_incid_rate_matrix,kemri_hosp_rate_matrix)}
# write_csv(kemri_hosp_rate_matrix,'input/kemri_hosp_rate_matrix.csv')
### ### ### ### ### ### ### ### ### ### ### ### ### ###
### South African data --------------------------------
fcn_load_s_afr=function(safr_data_path,n_iter,age_maxval,CI95_const,CI_intervals,
                        kemri_hosp_rate_stdev,randsampl_distrib_type){
s_afr_incidence_data=read_csv(safr_data_path)
s_afr_nat_average_totalcases=s_afr_incidence_data[s_afr_incidence_data$Province %in% 'South Africa' 
                                                  & s_afr_incidence_data$data_type %in% 'total',]
safr_hosp_rate_province_means=read_csv(s_afr_province_means_filepath)
# the incidence rate in data file is per 100K, normalize to per capita
popul_denom=1e5
if (max(s_afr_nat_average_totalcases$rate)>1){ s_afr_nat_average_totalcases[,c("rate","rate_CI_lower","rate_CI_upper")]=
  s_afr_nat_average_totalcases[,c("rate","rate_CI_lower","rate_CI_upper")]/popul_denom}
# st devs from CI95s
s_afr_nat_average_totalcases[,"stdev"]=rowMeans(cbind((s_afr_nat_average_totalcases$rate-s_afr_nat_average_totalcases$rate_CI_lower)/CI95_const,
                                (s_afr_nat_average_totalcases$rate_CI_upper-s_afr_nat_average_totalcases$rate)/CI95_const))
# we need a 60x5000 matrix of per capita RSV cases - need to pick a distribution to generate random samples
# normal distribution, assuming n=1
randsampl_distrib_type='gamma'
if (randsampl_distrib_type %in% 'normal'){
  s_afr_incid_rate_matrix=sapply(1:n_iter, function(iters) {
    sapply(1:age_maxval, function(x) {rnorm(1,mean=s_afr_nat_average_totalcases$rate[x],
                  sd=s_afr_nat_average_totalcases$stdev[x])})})} else {
                  # inferring gamma distrib from CIs
                  s_afr_incid_rate_matrix=matrix(NA,nrow=age_maxval,ncol=n_iter)
                        for (k in 1:nrow(s_afr_nat_average_totalcases)) {
                              if (k==1) {s_afr_gamma_estim=data.frame()}
                              gammavals=gamma.parms.from.quantiles(q=c(s_afr_nat_average_totalcases$rate_CI_lower[k],
                                s_afr_nat_average_totalcases$rate_CI_upper[k]),p=c(0.025,0.975))
                                s_afr_gamma_estim=rbind(s_afr_gamma_estim,c(gammavals$shape,gammavals$rate))
                                s_afr_incid_rate_matrix[k,]=rgamma(n_iter,shape=gammavals$shape,rate=gammavals$rate)  }
                                            }
# hosp rate
# we don't have a stdev value, we'll use MCMARCEL again
safr_hosp_rate_stdev=kemri_hosp_rate_stdev
s_afr_hosp_rate_matrix=sapply(1:n_iter, function(x) {
  rep(rnorm(1,mean=safr_hosp_rate_province_means$mean_hosp_rate[safr_hosp_rate_province_means$Province %in% 'South Africa'],
            sd=safr_hosp_rate_stdev),age_maxval)})
list(s_afr_incid_rate_matrix,s_afr_hosp_rate_matrix)
}