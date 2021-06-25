hosp_prob=get_rsv_ce_config(subset(sim_config_matrix,country_iso %in% "KEN" &intervention=="maternal"))$hosp_prob
ken_rsv_rate=get_rsv_ce_config(subset(sim_config_matrix,country_iso %in% "KEN" &intervention=="maternal"))$rsv_rate

n_age=50
ggplot(data.frame(age=1:n_age,ken_rsv_rate[1:n_age,]) %>% pivot_longer(!age) %>% mutate(name="kenya")) + 
  geom_density(aes(x=log10(value),color=factor(age),group=age)) + facet_wrap(~age,scales="free") + theme_bw() + standard_theme
# incidence
coeff_var=sd(log(df_hosp_prob$value))/abs(mean(log(df_hosp_prob$value)))
df_hosp_prob=data.frame(age=1,hosp_prob[1,]) %>% pivot_longer(!age) %>% 
  mutate(simulvals=rlnorm(n=5e3,meanlog=mean(log(value)),sdlog=sd(log(value))),
         simulvals_new=rlnorm(n=5e3,meanlog=mean(log(0.24)),sdlog=sd(log(value))),
         simulvals_logit_lambda=rnorm(5e3,mean=mean(log(value/(1-value))),sd=sd(log(value/(1-value)))),
         simulvals_logit=1/(1+exp(-simulvals_logit_lambda)) )

# distribution of hospitalisation rate
ggplot(df_hosp_prob) + 
  # geom_density(aes(x=log10(value/(1-value))),color="blue") + 
  geom_density(aes(x=log(value)),color="blue") + 
  geom_density(aes(x=log(simulvals_logit)),color="red") + # geom_density(aes(x=log10(simulvals_logit)),color="black") + 
  theme_bw() + standard_theme + geom_vline(xintercept=mean(log(df_hosp_prob$value))) + 
  geom_vline(xintercept=mean(log(df_hosp_prob$value))-1/2,linetype="dashed") + 
  geom_vline(xintercept=mean(log(df_hosp_prob$value))+1/2,linetype="dashed")

mean(kenya_incid_hosp_rate[[2]][1,])
ggplot(data.frame(age=1,kenya_incid_hosp_rate[[2]][1,]) %>% pivot_longer(!age) %>% mutate(name="kenya")) + 
  geom_density(aes(x=log10(value)),color="red") + # geom_density(aes(x=value),color="blue") +  
  theme_bw() + standard_theme


# generate logit-normal distrib when we have the desired mean and sd
des_mean=0.24; mean_guess=log(des_mean/(1-des_mean)); sd_guess=sd(log(df_hosp_prob$value/(1-df_hosp_prob$value)))
for (k in 1:1000){
  if (k==1) {if (mean(1/(1+exp(-rnorm(5e3,mean=mean_guess,sd=sd_guess))))<des_mean) {direct=-1} else {direct=1}}
  mean_rnorm=mean_guess+(1-k*direct/100)
  mean_approx=mean(1/(1+exp(-rnorm(5e3,mean=mean_rnorm,sd=sd_guess))))
  if (abs(mean_approx-des_mean)/des_mean<0.01) {break}
  print(c(k,mean_approx,mean_rnorm ))
}

approx_distr=1/(1+exp(-rnorm(5e3,mean=mean_rnorm,sd=sd_guess)))
ggplot(data.frame(mcmarcel=df_hosp_prob$value,value=approx_distr,logit_prob=log(approx_distr/(1-approx_distr)) )) + 
  geom_density(aes(x=logit_prob)) + geom_density(aes(x=log(mcmarcel/(1-mcmarcel))),color="blue") +
  theme_bw() + standard_theme + scale_x_continuous(breaks=(-8:4)/2) # 

### compare means and CI95 of data vs simulated
# compare CI95 values
x=data.frame(sapply(2:3, function(x) {
  fcn_load_kenya("../path_rsv_data/SARI_Rates_2010_2018_updated/ARI_SARI_Rates_2010_2018_tidydata.csv",
                                    sel_disease="SARI",n_iter=5e3,age_maxval=60)[[x]]})); colnames(x)[1]="mean_incidence"
ggplot(bind_rows(cbind(age=1:60,x,type="data"),
      cbind(data.frame(age=1:60,t(sapply(1:60, function(x) {c(mean(kenya_incid_hosp_rate[[1]][x,]),
      quantile(kenya_incid_hosp_rate[[1]][x,],probs=c(2.5,97.5)/1e2))}) )),type="simul") %>% 
  rename(mean_incidence=V1,CI_95_lower=X2.5.,CI_95_upper=X97.5.)),aes(x=age)) + 
  geom_line(aes(y=mean_incidence,group=type,color=type)) + geom_point(aes(y=mean_incidence,group=type,color=type)) +
  geom_ribbon(aes(ymin=CI_95_lower,ymax=CI_95_upper,group=type,color=type,linetype=type),size=1.1,alpha=0.2) + 
  theme_bw() + standard_theme + scale_x_continuous(limits=c(0,24),expand=expansion(0.01,0),breaks=1:60) # 
ggsave("output/sample_paths_test/kenya_incidence.png",width=15,height=10,units="cm")

# kenya sari incidence
ken_data_filepath="../path_rsv_data/SARI_Rates_2010_2018_updated/ARI_SARI_Rates_2010_2018_tidydata.csv"
dim(data.frame(sapply(2:3, function(x) {fcn_load_kenya(kenya_data_path=ken_data_filepath,
                 sel_disease="SARI",n_iter=5e3,age_maxval=60)[[x]]})) %>% dplyr::select(mean=1,everything()))
dim(data.frame(sapply(2:3, function(x) {fcn_load_kenya(ken_data_filepath,
                sel_disease="ARI",n_iter=5e3,age_maxval=60)[[x]]})) %>% dplyr::select(mean=1,everything()))

#########################
ggplot(bind_rows(fcn_load_s_afr(safr_data_path="../path_rsv_data/s_afr_incidence_data.csv",popul_denom=1e5) %>% 
                   dplyr::select(age,rate,rate_CI_lower,rate_CI_upper) %>% mutate(type="data",age=1:60),
                 cbind(data.frame(age=1:60,t(sapply(1:60, function(x) {c(mean(s_afr_incid_hosp_rate[[1]][x,]),
                                      quantile(s_afr_incid_hosp_rate[[1]][x,],probs=c(2.5,97.5)/1e2)) }) )),type="simul") %>% 
                   rename(rate=V1,rate_CI_lower=X2.5.,rate_CI_upper=X97.5.)),aes(x=age)) + 
  geom_line(aes(y=rate,group=type,color=type)) + geom_point(aes(y=rate,group=type,color=type)) +
  geom_ribbon(aes(ymin=rate_CI_lower,ymax=rate_CI_upper,group=type,color=type,linetype=type),size=1.1,alpha=0.2) + 
  theme_bw() + standard_theme + scale_x_continuous(limits = c(1,20),expand = expansion(0.01,0),breaks = 1:20)
ggsave("output/sample_paths_test/s_afr_incidence.png",width=15,height=10,units="cm")

###
# plot distrib of (simulated around our data points) incidence values
ggplot(data.frame(age=1:nrow(s_afr_incid_hosp_rate[[1]]),s_afr_incid_hosp_rate[[1]]) %>% 
         pivot_longer(!age) %>% mutate(logit_val=log(value/(1-value)))) + 
  facet_wrap(~age,scales="free") + geom_freqpoly(aes(x=value,group=age,color=factor(age)),bins=50) + theme_bw() + standard_theme
