# create UN country data (calls wpp package)
create_UN_country_database(output_dir)
# pre-process WPP2017 data
load_wpp2017_databases(output_dir)
# loads sex ratio and mortality by age groups; dataframes: sexratio, mxM, mxF
n_cntr=length(unique(sim_config_matrix$country_iso))
country_period_opt=data.frame(country_iso=rep(unique(sim_config_matrix$country_iso),2),
                              year=unlist(lapply(c(2020,2015),function(x){rep(x,n_cntr)})))
country_period_opt<-cbind(as.character(country_period_opt$country_iso), t(sapply(country_period_opt$year,get_year_category)))
# generate life table
# select country
i_life=which(country_period_opt[,1] %in% cntr_sel)[1]; f_country_iso=country_period_opt[i_life,1]
f_year=country_period_opt[i_life,3]; f_outputFileDir=output_dir
# with or without discounting
for (f_disc_rate_qaly in c(0,0.03)) {
  # this saves into an Rdata file
  generate_life_table(f_country_iso,f_year,f_outputFileDir,f_disc_rate_qaly)
  # 'life_table_KEN_2020_2025_disc0' ,'p03','.RData' # , 'life_table_KEN_2020_2025_disc0.RData'
  rdata_prefix=paste0('life_table_',f_country_iso,'_2020_2025_disc0')
  if (f_disc_rate_qaly>0) {rdata_filename=paste(rdata_prefix,'p03','.RData',sep='')} else {
    rdata_filename=paste(rdata_prefix,'.RData',sep='')}
  # load life_table
  load(paste(f_outputFileDir,"temp",rdata_filename,sep='/'))
  # life_table contains: "age","lx": # left alive, "life_expectancy", "nMx": mortality, "lx_rate": # left alive/(init popul)
  # "life_expectancy_disc": disc life yrs
  if (f_disc_rate_qaly==0){ life_table_tidy=data.frame(melt(life_table,id.vars='age'),disc_rate=f_disc_rate_qaly)} else {
    life_table_tidy=rbind(life_table_tidy,data.frame(melt(life_table,id.vars='age'),disc_rate=f_disc_rate_qaly)) }
}
####
# plot life tables
lifetable_varlist=list(as.character(unique(life_table_tidy$variable)),
                       c('alive/100.000 births','life expectancy','mortality','alive/birth','discounted life expectancy'))
life_table_tidy$variable_expl=as.character(life_table_tidy$variable)
for (k in unique(life_table_tidy$variable)){
  life_table_tidy$variable_expl[life_table_tidy$variable %in% k]=lifetable_varlist[[2]][lifetable_varlist[[1]] %in% k] }
# PLOT LIFETABLE
fcn_plot_lifetable <-function(save_tag){
p<-ggplot(life_table_tidy,aes(x=age,y=value,group=disc_rate,color=factor(disc_rate),linetype=factor(disc_rate))) + 
  geom_line(size=1.2) + facet_wrap(~variable_expl,scales='free') + 
  theme_bw() + standard_theme + xlab('age (in months)') + ylab('') + 
 ggtitle(paste(cntr_sel,'demographics')) + labs(color='discount rate',linetype='discount rate')
  if (save_tag=='save') {ggsave(paste0("output/life_table_",f_country_iso,
                                       "_2020_2025_disc0p03.png"),width=30,height=18,units="cm")}
p
}
### Shi incidence data -------------------------
# get unique country codes
country_opt <- data.frame(country_iso = unique(sim_config_matrix$country_iso))
incidence_one_table=get_incidence(country_opt$country_iso[which(country_opt$country_iso %in% cntr_sel)],output_dir)
# incidence data for all LMICs
# RSV_burden_Shi_2017=read_csv('input/RSV_burden_Shi_2017.csv')
# histogram: ggplot(RSV_burden_Shi_2017,aes(x=incidence_RSV_associated_ALRI)) + geom_histogram(binwidth=2)
# lineplot: national averages with CIs
RSV_burden_Shi_2017_tidy=read_csv('input/RSV_burden_Shi_2017_tidy.csv')
RSV_burden_Shi_2017_tidy$variable[grepl('incidence',RSV_burden_Shi_2017_tidy$variable)]='incidence of RSV-assoc. ALRI per 1000'
RSV_burden_Shi_2017_tidy$variable[grepl('episode',RSV_burden_Shi_2017_tidy$variable)]='number of episodes'
standard_theme_mod=standard_theme; standard_theme_mod$axis.text.x$size=6; standard_theme_mod$axis.text.x$hjust=0.99
standard_theme_mod$axis.text.x$vjust=0.5; non_na_vals=!(is.na(RSV_burden_Shi_2017_tidy$value) | is.na(RSV_burden_Shi_2017_tidy$lower_CI))

fcn_plot_lmic_incidence <-function(save_tag){
p<-ggplot(RSV_burden_Shi_2017_tidy[non_na_vals,],aes(x=location_name,y=value,ymin=lower_CI,ymax=upper_CI,group=1)) + 
  geom_pointrange(aes(color=variable),size=0.6,fatten=0.1,fill="white",shape=22) + scale_color_manual(values=c('blue','red')) +
  facet_wrap(~variable,nrow=2,scales='free') + theme_bw() + standard_theme_mod + 
  scale_y_continuous(trans='log10') + xlab('') + ylab('') + guides(color=FALSE)
  if (save_tag=='save') {  ggsave('output/incid_per_cntr.png',width=30,height=18,units="cm")}
p
  }
