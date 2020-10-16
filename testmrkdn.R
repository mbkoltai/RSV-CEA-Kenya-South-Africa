#+ fig.width=15, fig.height=10
ggplot(burden_mcmarcel_owndata_comp,aes(x=value,group=source)) +
  # geom_histogram(aes(y=..density..,fill=source),color="NA",size=0.4) + 
  geom_freqpoly(aes(color=source),size=1.2) + 
  geom_vline(data=mean_intercepts,aes(xintercept=int,linetype=source),size=0.8) + # color='black'
  facet_wrap(~variable,scales='free',labeller=label_wrap_gen(width=10)) + # ,ncol=3
  theme_bw() + theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),plot.title=element_text(hjust=0.5,size=16),
                     axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=7),
                     axis.title=element_text(size=14),text=element_text(family="Calibri")) +
  scale_linetype_manual(values=c('solid','dotdash')) +
  geom_rect(data=subset(burden_mcmarcel_owndata_comp, variable %in% icercolname),fill=NA,colour="blue",size=2,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) + 
  ggtitle(paste(sel_interv$country_iso,'RSV burden & intervention estimates:',gsub('_','',interv_tag))) +
  labs(color='data source',linetype='mean') + guides(xintercept=FALSE,linetype=guide_legend(ncol=2)) # xlab('')+ylab('')
####
### stafafzda
ggplot(RSV_burden_Shi_2017_tidy,aes(x=location_name,y=value,group=1)) + geom_line() + geom_point(size=1) +
  geom_ribbon(aes(ymin=lower_CI,ymax=upper_CI),alpha=0.3,colour=NA,fill="red") + facet_wrap(~variable,nrow=2,scales='free') +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5)) + scale_y_continuous(trans='log10')
