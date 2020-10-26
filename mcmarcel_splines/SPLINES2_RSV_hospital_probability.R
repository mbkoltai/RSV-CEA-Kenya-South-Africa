#############################################################################
# This file is part of the RSV modelling project.
# 
# => FIT HOSPITAL PROBABILITY
#
#  Copyright 2018, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
# note: code based on work by Marina Antillon
#############################################################################

######################
## Clear workspace
######################

rm(list=ls())

######################
## Useful packages
######################

# load packages
library('Hmisc')
library('mgcv')
library('graphics')
library("mvtnorm")
library('ggplot2')
library("gamm4")

# help functions
logistic = function(x){exp(x)/(1+exp(x))}

logit = function(x){
  x[x<0.001] = 0.001
  x[x>0.999] = 0.999
  log(x/(1-x))
}

############################################
## Settings and output directory
############################################

# set number of samples
num_splines <- 5000

# create plots?
get_figures <- TRUE

# perform sensitivity analysis?
sens_analysis <- FALSE 

# set RNG seed
rng_seed <- 20190617
set.seed(rng_seed)

# set regression type: cr = 'cubic regression spline' and ts = 'thin plate regression spline'
bs_type <- "ts" # options: cs, ts

# fix ouput file tag
fitting_tag <- paste0("hosp_",bs_type,'_n',num_splines,'_rnorm2')

# set plot prefix
plot_prefix = file.path("splines/output",fitting_tag)

# check if the output folder exist, and create the folder if not
if(!dir.exists(plot_prefix)){
  print(paste('create folder:',plot_prefix))
  dir.create(plot_prefix,recursive=T)
}

################################
## Estimate knot positions
## Logistic regression
################################

# load data
hosp = read.csv("./splines/data/hospitalizations.csv")
hosp$est = hosp$hosp/hosp$rsv_cases
hosp$lci = binconf(hosp$hosp, hosp$rsv_cases, method="exact")[,"Lower"]
hosp$uci = binconf(hosp$hosp, hosp$rsv_cases, method="exact")[,"Upper"]
hosp$Title = paste(hosp$Location, "\n(", hosp$Author, ", ", hosp$Year, ")", 
                   "\n", hosp$Study_period, sep="")

################################
## LMIC: NOKES AND HOMAIRA
################################

hosp_lmic = hosp[hosp$citation %in% c("Homaira 2012", "Nokes 2008"),]
hosp_lmic$study_no = factor(as.numeric(factor(as.numeric(hosp_lmic$citation))))
hosp_lmic$dummy=1
rm(list=c('hosp'))

b_gamm = gamm4(cbind(hosp, rsv_cases)~s(log(midpoint), k=-1, bs=bs_type) + 
                 t2(log(midpoint), study_no, k=3, bs=c(bs_type, 're'), by=dummy),
               random=~(1|study_no), data=hosp_lmic, family=binomial(link="logit"), REML=T)

newdata = data.frame(midpoint=seq(0.1, 60, 0.1), 
                     Cases=100,
                     study_no=1, dummy=0)
bpred = data.frame(predict.gam(b_gamm$gam, newdata=newdata, se.fit=T))
bpred$lfit = bpred$fit-1.96*bpred$se.fit
bpred$ufit = bpred$fit+1.96*bpred$se.fit

bpred=cbind(newdata, bpred)

if(get_figures) {
  matplot(bpred$midpoint,logistic(as.matrix(bpred[,c("fit", "lfit", "ufit")])), 
        type="l", lty=1, ylim=c(0,1))
  points(hosp_lmic$midpoint, hosp_lmic$est, pch=20)
}

# this predicts the basis at new values
tmp=predict.gam(b_gamm$gam, newdata=newdata, type="lpmatrix", se.fit=T)
# Under ML or REML smoothing parameter estimation it is 
# possible to correct the covariance matrix Vp for smoothing 
# parameter uncertainty. This is the corrected version.

somebetas = rmvnorm(n=num_splines, mean = coef(b_gamm$gam), sigma = b_gamm$gam$Vp)
someiterates = (tmp %*% t(somebetas))
# MUST DO: should check that it gives me uncertainty 
# equivalent to that of the prediction functions
allpred = logistic(someiterates)

allpred_df_lmic = data.frame(pred=as.vector(allpred[seq(5, 595, 10),]), 
                                 mos=rep(seq(0.5, 59.5, 1), times=dim(allpred)[2]), 
                                 iter=rep(1:dim(allpred)[2], each=length(seq(0.5, 59.5, 1))))


write.csv(allpred_df_lmic, file=file.path(plot_prefix, paste0('hosp_lmic_predictions_',bs_type,'_n',num_splines,'.csv')))
save(allpred_df_lmic, allpred, file=file.path(plot_prefix, paste0('hosp_lmic_predictions_',bs_type,'_n',num_splines,'.Rdata')))

## EXPLORE NUMBER OF SPLINES (optional)
if(sens_analysis){
  
  run_design <- expand.grid(num_splines = seq(1000,10000,1000),
                            run_id      = 1:100)
  library(foreach)
  i_design <- 1
  for_out <- foreach(i_design = 1:nrow(run_design), .combine=rbind) %do%
  {
    somebetas      = rmvnorm(n=run_design$num_splines[i_design], coef(b_gamm$gam), b_gamm$gam$Vp)
    someiterates   = (tmp %*% t(somebetas))
    allpred_single = logistic(someiterates[1,])
    
    quantile(allpred_single,c(0,0.5,0.975,0.99,0.999,1))
  }
  
  dim(for_out)
  
  pdf(file.path(plot_prefix,"splines_stochastic_both.pdf"),10,10)
  ii <- 2
  colnames(for_out)[c(1,ncol(for_out))] <- c('min','max')
  par(mfrow=c(2,3))
  for(ii in 1:ncol(for_out)){
    boxplot(for_out[,ii] ~ run_design$num_splines,
            main=colnames(for_out)[ii],
            ylab='hospital probability',
            xlab='number of splines')
  }
  
  for(ii in 1:ncol(for_out)){
    boxplot(for_out[,ii] ~ run_design$num_splines,
            main=colnames(for_out)[ii],
            ylab='hospital probability',
            xlab='number of splines',
            ylim=0:1)
  }
  dev.off()
}

################################
# LMIC: NOKES ET AL
################################

# select the Nokes data
hosp_Nokes = hosp_lmic[hosp_lmic$citation=="Nokes 2008",]

# fit GAM model with binomial logit link, 3 knots and (penalized) smoothing basis = 'tin plate regression spline':
# - hospital cases
# - rsv cases
# - log(age)
b = gam(cbind(hosp, rsv_cases)~s(log(midpoint), k=3, bs=bs_type),
        data=hosp_Nokes, family=binomial(link="logit"))

# create new data matrix with ages from 0.1 up to 60 months
newdata = data.frame(midpoint=seq(0.1, 60, 0.1),
                     rsv_cases=rep(100, length(seq(0.1, 60, 0.1))))

# get predicted values from the GAM model
bpred      = data.frame(predict.gam(b, newdata=newdata, se.fit=T))

# get lower and upper limit, based on the standard error
bpred$lfit = bpred$fit-1.96*bpred$se.fit*2
bpred$ufit = bpred$fit+1.96*bpred$se.fit*2
bpred      = cbind(newdata, bpred)

diff(range(bpred$fit)) - bpred$fit[1]

# sample from normal distribution using the mean of the fitted values and doubling the standard error
somebetas    <- rnorm(num_splines, mean = mean(bpred$fit), sd = mean(bpred$se.fit)*2)

# replicate the value for all age groups
someiterates <- matrix(rep(somebetas,each=nrow(newdata)),nrow=nrow(newdata),byrow=F)

{# tmp=predict.gam(b, newdata=newdata, type="lpmatrix", se.fit=T)
  # somebetas = rmvnorm(n=num_splines, mean=coef(b), sigma=b$Vp)
  # someiterates = (tmp %*% t(somebetas))
}## Previous code to sample new values

# convert to regular scale
allpred_Nokes = logistic(someiterates)

allpred_df_Nokes = data.frame(pred=as.vector(allpred_Nokes[seq(5, 595, 10),]),
                              mos=rep(seq(0.5, 59.5, 1), times=dim(allpred_Nokes[seq(5, 595, 10),])[2]),
                              iter=rep(1:dim(allpred_Nokes[seq(5, 595, 10),])[2], each=length(seq(0.5, 59.5, 1))))

write.csv(allpred_df_Nokes, file=file.path(plot_prefix, paste0("hosp_NOKES_predictions_",bs_type,"_n",num_splines,".csv")))
save(allpred_df_Nokes, allpred_Nokes, file=file.path(plot_prefix, paste0("hosp_NOKES_predictions_",bs_type,"_n",num_splines,".Rdata")))

if(get_figures){
  matplot(bpred$midpoint,logistic(as.matrix(bpred[,c("fit", "lfit", "ufit")])),
          type="l", lty=1, ylim=c(0,1))
}


## EXPLORE NUMBER OF SPLINES (optional)
if(sens_analysis){
  
  for_out <- foreach(i_design = 1:nrow(run_design), .combine=rbind) %do%
  {
    somebetas       = rmvnorm(n=run_design$num_splines[i_design], coef(b), b$Vp)
    someiterates    = (tmp %*% t(somebetas))
    allpred_single  = logistic(someiterates[1,])
    quantile(allpred_single,c(0,0.5,0.975,0.99,0.999,1))
  }
  
  
  pdf(file.path(plot_prefix,"splines_stochastic_nokes.pdf"),10,10)
  colnames(for_out)[c(1,ncol(for_out))] <- c('min','max')
  par(mfrow=c(2,3))
  for(ii in 1:ncol(for_out)){
    boxplot(for_out[,ii] ~ run_design$num_splines,
            main=colnames(for_out)[ii],
            ylab='hospital probability',
            xlab='number of splines')
  }
  
  for(ii in 1:ncol(for_out)){
    boxplot(for_out[,ii] ~ run_design$num_splines,
            main=colnames(for_out)[ii],
            ylab='hospital probability',
            xlab='number of splines',
            ylim=c(0,0.3))
  }
  dev.off()
}


################################
# PLOTS
################################
if(get_figures)
{  
  themebar2 = theme(axis.text.x = element_text(color="black", size=8, angle=0),
                    axis.title.x = element_text(size = 8, angle = 0, face="bold"),
                    axis.text.y = element_text(color="black", size=8, angle=0),
                    axis.title.y = element_text(size = 8, angle = 90, face="bold"),
                    panel.border = element_rect(linetype = "solid", colour = "black", fill=NA),
                    legend.text = element_text(size = 8, face = "bold", lineheight=0.8),
                    legend.position = "bottom",
                    legend.box = "vertical",
                    legend.background = element_rect(fill=NA, size=0.25, linetype="solid", colour ="black"),
                    legend.key = element_rect(fill = "white"),
                    # legend.title = element_blank(),
                    panel.grid.major = element_line(colour="gray70", linetype = "dotted"),
                    panel.grid.minor = element_line(colour="gray60", linetype = "dotted"),
                    panel.background = element_rect(fill = NA),
                    strip.background = element_rect(fill = NA),
                    strip.text = element_text(size=8, face="bold")) # , strip.text = element_blank()
  
  jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_predictions_',bs_type,'.jpeg')), width = 3, height = 3, units = 'in', res=600)
    print(ggplot(data=hosp_lmic, aes(x=midpoint, y=est)) + themebar2 +
      geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.1) +
      geom_point(size=1) + facet_wrap(~citation, ncol=5) +
      geom_errorbar(data=hosp_lmic, aes(x=midpoint, ymin=lci, ymax=uci), width=.25) +
      xlab("Age in months") + ylab("Hospitalisation probability") +
      scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + scale_y_continuous(limits = c(0, 1)))
  dev.off()
  
  #########################
  # Ribbon plots - Xiao's analysis.
  #########################
  
  # columns: lowci, hici, mos
  hosp_ribbons = data.frame(t(apply(allpred, 1, quantile, c(0.5, 0.025, 0.975))))
  colnames(hosp_ribbons) = c("est", "lowci", "hici")
  hosp_ribbons$mos = seq(0.1, 60, 0.1)
  
  jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_ribbons_lmic_',bs_type,'.jpeg')),
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=hosp_ribbons, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
    theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
    geom_ribbon(alpha=0.5, size=0, fill = "rosybrown1") +
    xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
    scale_y_continuous(breaks=seq(0, 0.5, 0.1), limits = c(0, 0.5)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_ribbons_logit_lmic_',bs_type,'.jpeg')),
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=hosp_ribbons, aes(x=mos, ymin=logit(lowci), ymax=logit(hici))) + themebar2 +
    theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
    theme(panel.grid.minor = element_blank()) +
    geom_ribbon(alpha=0.5, size=0, fill = "rosybrown1") +
    xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients (log-odds scale)") +
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
    scale_y_continuous(breaks=logit(c(0.001, 0.0025, 0.005, 0.025, 0.01, 0.05, 0.1, 0.25, 0.5)),
                       labels = c("0.001", "0.0025", "0.005", "0.025", "0.01","0.05","0.1", "0.25", "0.5"), limits = c(-7, 0)))
  dev.off()

################################
# PLOTS - Nokes by itself
################################

jpeg(filename = file.path(plot_prefix, "plot_hosp_NOKES.jpeg"),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_Nokes, aes(x=midpoint, y=est)) + themebar2 +
  theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
  geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +
  geom_line(data=allpred_df_Nokes, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.05) +
  geom_point(size=1) +
  facet_wrap(~Title) +
  geom_errorbar(data=hosp_Nokes, aes(x=midpoint, ymin=lci, ymax=uci), width=.25) +
  xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
  scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + scale_y_continuous(limits = c(0, 1)))
dev.off()

jpeg(filename = file.path(plot_prefix,"plot_hosp_NOKES_logit.jpeg"),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_Nokes, aes(x=midpoint, y=logit(est+0.001))) + themebar2 +
  theme(strip.text = element_text(size=6, face="bold")) +
  theme(panel.grid.minor = element_blank()) +
  geom_line(data=allpred_df_lmic, aes(x=mos, y=logit(pred+0.001), group=iter), col="rosybrown2", alpha=0.05) +
  geom_line(data=allpred_df_Nokes, aes(x=mos, y=logit(pred+0.001), group=iter), col="lightskyblue", alpha=0.05) +
  geom_point(size=1) + facet_wrap(~Title) +
  geom_errorbar(data=hosp_Nokes, aes(x=midpoint, ymin=logit(lci+0.001), ymax=logit(uci+0.001)), width=.25) +
  xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients (log-odds scale)") +
  scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
  scale_y_continuous(breaks=logit(c(0.001, 0.0025, 0.005, 0.025, 0.01, 0.05, 0.1, 0.25, 0.5)),
                     labels = c("0.001", "0.0025", "0.005", "0.025", "0.01","0.05","0.1", "0.25", "0.5"), limits = c(-7, 2)))
dev.off()

jpeg(filename = file.path(plot_prefix,"plot_hosp_nokes_fit_obs_logit.jpeg"),
     width = 7, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_Nokes, aes(x=midpoint, y=logit(est+0.001))) + themebar2 +
  theme(strip.text = element_text(size=6, face="bold")) +
  theme(panel.grid.minor = element_blank()) +
  geom_line(data=allpred_df_lmic, aes(x=mos, y=logit(pred+0.001), group=iter), col="rosybrown2", alpha=0.05) +
  geom_line(data=allpred_df_Nokes, aes(x=mos, y=logit(pred+0.001), group=iter), col="lightskyblue", alpha=0.05) +
  geom_point(size=1) + facet_wrap(~Title, ncol=3) +
  geom_errorbar(aes(x=midpoint, ymin=logit(lci+0.001), ymax=logit(uci+0.001)), width=.25) +
  xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients (log-odds scale)") +
  scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
  scale_y_continuous(breaks=logit(c(0.001, 0.0025, 0.005, 0.025, 0.01, 0.05, 0.1, 0.25, 0.5)),
                     labels = c("0.001", "0.0025", "0.005", "0.025", "0.01","0.05","0.1", "0.25", "0.5"), limits = c(-7, 2)))
dev.off()


#########################
# Ribbon plots - Nokes vs. Both
#########################

# columns: lowci, hici, mos
hosp_ribbons_both  = data.frame(t(apply(allpred, 1, quantile, c(0.5, 0.025, 0.975))),
                                apply(allpred, 1, mean))
hosp_ribbons_nokes = data.frame(t(apply(allpred_Nokes, 1, quantile, c(0.5, 0.025, 0.975))),
                                apply(allpred_Nokes, 1, mean))
colnames(hosp_ribbons_both)  = c("est", "lowci", "hici","mean")
colnames(hosp_ribbons_nokes) = c("est", "lowci", "hici","mean")
hosp_ribbons_both$mos  = seq(0.1, 60, 0.1)
hosp_ribbons_nokes$mos = seq(0.1, 60, 0.1)


jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_ribbons_95CI_compare_',bs_type,'.jpeg')),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_ribbons_both, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
        theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
        geom_ribbon(data=hosp_ribbons_both, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "lightskyblue") +
        geom_ribbon(data=hosp_ribbons_nokes, aes(x=mos, ymin=lowci, ymax=hici), alpha=1, size=0, fill = "rosybrown2" ) +
        xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
        scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
        scale_y_continuous(breaks=seq(0, 0.5, 0.1), limits = c(0, 0.5)))
dev.off()

jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_ribbons_95CI_mean_compare_',bs_type,'.jpeg')),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_ribbons_both, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
        theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
        geom_ribbon(data=hosp_ribbons_both, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "lightskyblue") +
        geom_ribbon(data=hosp_ribbons_nokes, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.7, size=0, fill = "rosybrown2" ) +
        
        geom_ribbon(data=hosp_ribbons_both, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "blue" ) +
        geom_ribbon(data=hosp_ribbons_nokes, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
        
        xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
        scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
        scale_y_continuous(breaks=seq(0, 0.5, 0.1), limits = c(0, 0.5)))
dev.off()

jpeg(filename = file.path(plot_prefix, paste0('plot_hosp_ribbons_95CI_mean_',bs_type,'.jpeg')),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_ribbons_both, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
        theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
        geom_ribbon(data=hosp_ribbons_nokes, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "rosybrown2" ) +
        
        geom_ribbon(data=hosp_ribbons_nokes, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
        
        xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
        scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
        scale_y_continuous(breaks=seq(0, 0.5, 0.1), limits = c(0, 0.5)))
dev.off()

jpeg(filename = file.path(plot_prefix,paste0('plot_hosp_fit_compare_',bs_type,'.jpeg')),
     width = 2.5, height = 3, units = 'in', res=600)
print(ggplot(data=hosp_Nokes, aes(x=midpoint, y=est)) + themebar2 +
        theme(strip.text = element_text(size=6, face="bold")) + # , strip.text = element_blank()
        geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.25) +
        geom_line(data=allpred_df_Nokes, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.25) +
        xlab("Age in months") + ylab("Probability of Hospitalisation\nAmong All Patients") +
        scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
        scale_y_continuous(limits = c(0, 1)))
dev.off()

}
