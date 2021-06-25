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
fcn_gen_lifetable <- function(cntrsel){
  # select country
  i_life=which(country_period_opt[,1] %in% cntrsel)[1]; f_country_iso=country_period_opt[i_life,1]
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
    if (f_disc_rate_qaly==0){life_table0=life_table
      life_table_tidy=data.frame(life_table %>% pivot_longer(cols=!age,names_to="variable"),disc_rate=f_disc_rate_qaly)
    } else {  life_table_tidy=rbind(life_table_tidy,
                                    data.frame(life_table %>% pivot_longer(cols=!age,names_to="variable"),disc_rate=f_disc_rate_qaly))}
  }
  ####
  # plot life tables
  lifetable_varlist=list(as.character(unique(life_table_tidy$variable)),
                         c('alive/100.000 births','life expectancy','mortality','alive/birth','discounted life expectancy'))
  life_table_tidy$variable_expl=as.character(life_table_tidy$variable)
  for (k in unique(life_table_tidy$variable)){
    life_table_tidy$variable_expl[life_table_tidy$variable %in% k]=lifetable_varlist[[2]][lifetable_varlist[[1]] %in% k] }
  
    list(life_table0,life_table,life_table_tidy)
}
# fcn plot lifetable

### PLOT LIFETABLE -------------------------
fcn_plot_lifetable <-function(cntrsel,save_tag){
  g(life_table0,life_table,lifetabletidy) %=% fcn_gen_lifetable(cntrsel)
  
  p<-ggplot(lifetabletidy,aes(x=age,y=value,group=disc_rate,color=factor(disc_rate),linetype=factor(disc_rate))) + 
  geom_line(size=1.2) + facet_wrap(~variable_expl,scales='free') + 
  theme_bw() + standard_theme + xlab('age (in months)') + ylab('') + 
 ggtitle(paste(cntrsel,'demographics')) + labs(color='discount rate',linetype='discount rate')
  if (save_tag=='save') {ggsave(paste0("output/life_table_",cntrsel,"_2020_2025_disc0p03.png"),width=30,height=18,units="cm")}
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
fcn_get_RSV_burden_shi2017 <- function(filename){
  RSVburdenShi2017tidy=read_csv(filename)
  RSVburdenShi2017tidy$variable[grepl('incidence',RSVburdenShi2017tidy$variable)]='incidence of RSV-assoc. ALRI per 1000'
  RSVburdenShi2017tidy$variable[grepl('episode',RSVburdenShi2017tidy$variable)]='number of episodes'
  RSVburdenShi2017tidy
}

fcn_plot_lmic_incidence <-function(RSVburden_Shi_2017_tidy,save_tag){
  standard_theme_mod=standard_theme; standard_theme_mod$axis.text.x$size=6; standard_theme_mod$axis.text.x$hjust=0.99
  standard_theme_mod$axis.text.x$vjust=0.5
  non_na_vals=!(is.na(RSVburden_Shi_2017_tidy$value) | is.na(RSVburden_Shi_2017_tidy$lower_CI))
  p <- ggplot(RSVburden_Shi_2017_tidy[non_na_vals,],aes(x=location_name,y=value,ymin=lower_CI,ymax=upper_CI,group=1)) + 
    geom_pointrange(aes(color=variable),size=0.6,fatten=0.1,fill="white",shape=22) + scale_color_manual(values=c('blue','red')) +
    facet_wrap(~variable,nrow=2,scales='free') + theme_bw() + standard_theme_mod + scale_y_continuous(trans='log10') + 
    xlab('') + ylab('') + guides(color=FALSE)
  if (save_tag=='save') {  ggsave('output/incid_per_cntr.png',width=30,height=18,units="cm")}
  p
}

###########
# assing multiple vars
### assign multiple variables -----------------
'%=%' = function(l, r, ...) UseMethod('%=%')
# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}
# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

### standard theme for plotting -------------------
standard_theme=theme(#panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
  plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
  legend.title=element_text(size=14),legend.text=element_text(size=12),
  axis.title=element_text(size=14), text=element_text(family="Calibri"))
