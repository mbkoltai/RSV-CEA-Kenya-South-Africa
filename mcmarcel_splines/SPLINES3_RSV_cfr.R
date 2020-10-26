#############################################################################
# This file is part of the RSV modelling project.
# 
# => FIT HOSPILTAL CASE FATALITY RATIO (hCFR) SPLINES
#
#  Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
# based on code from Marina Antillon
#############################################################################

# clear workspace
rm(list=ls())

# load packages
library('Hmisc')
library('mgcv')
library('graphics')
library("mvtnorm")
library('ggplot2')
library("gamm4")
library("plyr")


# help functions
logit = function(x){
  x[x<0.001] = 0.001
  x[x>0.999] = 0.999
  log(x/(1-x))
}
logistic = function(x){exp(x)/(1+exp(x))}
substrRight = function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

############################################
## Settings and output directory
############################################
# set number of splines
num_splines <- 5000

# set RNG seed
rng_seed <- 20190118
set.seed(rng_seed)

# create plots?
get_figures <- TRUE

# set regression type: cr = 'cubic regression spline' and ts = 'thin plate regression spline'
bs_type <- "ts" # options: cs, ts

# fix ouput file tag
fitting_tag <- paste0("cfr_",bs_type,'_n',num_splines)

# set plot prefix
plot_prefix = file.path("splines/output",fitting_tag)

# check if the output folder exist, and create the folder if not
if(!dir.exists(plot_prefix)){
  print(paste('create folder:',plot_prefix))
  dir.create(plot_prefix,recursive=T)
}

############################################
## Call data in; make into long-form
############################################

load("./splines/data/Shi_CFR_workingfile.Rdata")
cfr_shi$study_no = 1:dim(cfr_shi)[1]
cfr_shi = cfr_shi[substrRight(cfr_shi$Study_period, 4)>1999,]

cases_vars = which(substr(colnames(cfr_shi), 1, 5) == "Cases")
deaths_vars = which(substr(colnames(cfr_shi), 1, 6) == "Deaths")
redundant_vars = which(substr(colnames(cfr_shi), 1, 9) == "Redundant")
id_vars = which(colnames(cfr_shi) %in% c("study_no", "Author", "Year", "Study_period", "Location"))

table(cfr_shi$Economic_setting, cfr_shi$Development_status)

# cases
agegroups=c()
for (i in 1:length(colnames(cfr_shi)[cases_vars])){
  agegroups[i] = paste(strsplit(colnames(cfr_shi)[cases_vars], "_")[[i]][2], "_", 
                       strsplit(colnames(cfr_shi)[cases_vars], "_")[[i]][3], sep="")
}

tmp_cases = reshape(cfr_shi[,-c(deaths_vars, redundant_vars)], 
                    idvar = c("study_no", "Author", "Year", "Study_period", "Location"), 
                    varying = colnames(cfr_shi)[cases_vars],
                    v.names="Cases", 
                    timevar="Age_Groups", times=agegroups, 
                    direction="long")
tmp_cases = tmp_cases[!is.na(tmp_cases$Cases),]
rownames(tmp_cases) = NULL

# deaths
agegroups=c()
for (i in 1:length(colnames(cfr_shi)[deaths_vars])){
  agegroups[i] = paste(strsplit(colnames(cfr_shi)[deaths_vars], "_")[[i]][2], "_", 
                       strsplit(colnames(cfr_shi)[deaths_vars], "_")[[i]][3], sep="")
}

tmp_deaths = reshape(cfr_shi[,c(id_vars, deaths_vars)], 
                     idvar = c("study_no", "Author", "Year", "Study_period", "Location"), 
                     varying = colnames(cfr_shi)[deaths_vars],
                     v.names="Deaths", 
                     timevar="Age_Groups", times=agegroups, 
                     direction="long")
tmp_deaths = tmp_deaths[!is.na(tmp_deaths$Deaths),]
rownames(tmp_deaths) = NULL

# redundancy vars
agegroups=c()
for (i in 1:length(colnames(cfr_shi)[redundant_vars])){
  agegroups[i] = paste(strsplit(colnames(cfr_shi)[redundant_vars], "_")[[i]][2], "_", 
                       strsplit(colnames(cfr_shi)[redundant_vars], "_")[[i]][3], sep="")
}

tmp_red = reshape(cfr_shi[,c(id_vars, redundant_vars)], 
                  idvar = c("study_no", "Author", "Year", "Study_period", "Location"), 
                  varying = colnames(cfr_shi)[redundant_vars],
                  v.names="Redundant", 
                  timevar="Age_Groups", times=agegroups, 
                  direction="long")
tmp_red = tmp_red[!is.na(tmp_red$Redundant),]
rownames(tmp_red) = NULL

# join it all together
cfr_long = join(tmp_cases, tmp_deaths)
cfr_long = join(cfr_long, tmp_red, type="left")

sum(cfr_long$Redundant==0, na.rm=T)
sum(tmp_red$Redundant==0, na.rm=T)

######################
## Assign midpoints
######################

agemdpt = data.frame(Age_Groups = c("0d_27d", "28d_3m", "0m_3m", "3m_5m", "1m_5m", "0m_5m", 
                                    "6m_8m", "9m_11m", "6m_11m", "1m_11m", "0m_11m",
                                    "12m_23m", "24m_35m", "0m_23m", "0m_35m", "0m_47m", "12m_59m", 
                                    "24m_59m", "36m_59m", "0m_59m", "1m_59m"),
                     midpoint = c(0.5, 2, 1.5, 4.5, 3.5, 3, 7.5, 10.5, 9, 6.5, 6, 18, 30, 12, 18, 
                                  24, 36, 42, 48, 30, 30.5))

cfr_long = join(cfr_long, agemdpt, type="left")

## Titles for graphs
cfr_long$Title = paste(cfr_long$Location, "\n(", cfr_long$Author, ", ", cfr_long$Year, ")", 
                       "\n", cfr_long$Study_period, sep="")

######################
# round all cases to 
# nearest whole (some are 
# X.9995 for some reason)
######################

cfr_long$Cases = round(cfr_long$Cases)
cfr_long$Deaths = round(cfr_long$Deaths)

cfr_long = cfr_long[cfr_long$Economic_setting!="High income",]
cfr_long$Economic_setting = factor(cfr_long$Economic_setting)
levels(cfr_long$Economic_setting)

############################################
## calculate CFR
############################################

cfr_long$cfr = cfr_long$Deaths/cfr_long$Cases
cfr_long$cfr_lci = binconf(cfr_long$Deaths, cfr_long$Cases, method="exact")[,"Lower"]
cfr_long$cfr_uci = binconf(cfr_long$Deaths, cfr_long$Cases, method="exact")[,"Upper"]

################################
## Logistic regression
## LIC+LMIC - altogether
###############################

cfr_long_min = cfr_long[cfr_long$Redundant!=1 | is.na(cfr_long$Redundant),]
cfr_long_min = cfr_long_min[!(cfr_long_min$Age_Groups %in% c( "0m_47m")),] # "0m_35m",
cfr_long_min = cfr_long_min[cfr_long_min$Economic_setting %in% c("Low income","Lower middle income"),]
cfr_long_min_all = cfr_long_min
cfr_long_min = cfr_long_min[cfr_long_min$study_no %in% as.numeric(names(table(cfr_long_min$study_no)))[table(cfr_long_min$study_no)>2],]
cfr_long_min_all = cfr_long_min_all[cfr_long_min_all$study_no %in%
                                      as.numeric(names(table(cfr_long_min_all$study_no)))[table(cfr_long_min_all$study_no)<3],]

cfr_long_min$study_no = factor(as.numeric(as.factor(cfr_long_min$study_no)))
cfr_long_min$Economic_setting = as.factor(as.character(cfr_long_min$Economic_setting))
cfr_long_min$dummy=1

b_gamm = gamm4(cbind(Deaths, Cases)~s(log(midpoint), k=-1, bs=bs_type) +
                 t2(log(midpoint), study_no, k=3, bs=c(bs_type, 're'), by=dummy),
               random=~(1|study_no), data=cfr_long_min, family=binomial(link="logit"), REML=T)

newdata = data.frame(midpoint=seq(0.1, 60, 0.1),
                     Cases=1000,
                     study_no=1, dummy=0)
bpred = data.frame(predict.gam(b_gamm$gam, newdata=newdata, se.fit=T))
bpred$lfit = bpred$fit-1.96*bpred$se.fit
bpred$ufit = bpred$fit+1.96*bpred$se.fit

bpred=cbind(newdata, bpred)

if(get_figures){
  matplot(bpred$midpoint,logistic(as.matrix(bpred[,c("fit", "lfit", "ufit")])), type="l", lty=1, ylim=c(0,.4))
  points(cfr_long_min$midpoint, cfr_long_min$cfr, pch=20)
}

# this predicts the basis at new values
tmp=predict.gam(b_gamm$gam, newdata=newdata, type="lpmatrix", se.fit=T)
#coef(b_gamm$gam)
#b_gamm$gam$Vp # seems to be the same as # vcov(b$gam,unconditional=TRUE)
# gam.check(b_gamm$gam)
# uncertainty for smoothing parameter:
# vcov(b$gam,unconditional=TRUE)

# b_gamm$gam$Vc
# Under ML or REML smoothing parameter estimation it is
# possible to correct the covariance matrix Vp for smoothing
# parameter uncertainty. This is the corrected version.

somebetas = rmvnorm(n=1000, coef(b_gamm$gam), b_gamm$gam$Vp)
someiterates = (tmp %*% t(somebetas))
# MUST DO: should check that it gives me uncertainty equivalent to that of the prediction functions
allpred = logistic(someiterates)
matplot(allpred, type="l", lty=1, col=rgb(0,0,0,alpha=0.1))

allpred_df_lic_lmic = data.frame(pred=as.vector(allpred[seq(5, 595, 10),]),
                                 mos=rep(seq(0.5, 59.5, 1), times=dim(allpred)[2]),
                                 iter=rep(1:dim(allpred)[2], each=length(seq(0.5, 59.5, 1))))

write.csv(allpred_df_lic_lmic, file=file.path(plot_prefix, paste0("cfr_lic_lmic_predictions_",bs_type,"_n",num_splines,".csv")))
save(allpred_df_lic_lmic, allpred, file=file.path(plot_prefix, paste0("cfr_lic_lmic_predictions_",bs_type,"_n",num_splines,".Rdata")))

################################
## Analysis by economic stratum 
## (LIC vs LMIC only)
################################

v_gamm = gamm4(cbind(Deaths, Cases)~s(log(midpoint), k=-1, bs=bs_type, by=Economic_setting) + 
                 t2(log(midpoint), study_no, k=3, bs=c(bs_type, 're'), by=dummy),
               random=~(1|study_no), data=cfr_long_min, family=binomial(link="logit"), REML=T)

modcomp = anova(b_gamm$mer, v_gamm$mer)
write.csv(modcomp, file=file.path(plot_prefix, paste0("modcomp_",bs_type,"_n",num_splines,".csv")))

tmp = unique(cfr_long_min[,c("Economic_setting")])
newdata = data.frame(midpoint=rep(seq(0.1, 60, 0.1), length(tmp)), 
                     Pop=1000,
                     Economic_setting = rep(rev(tmp), each=length(seq(0.1, 60, 0.1))), 
                     study_no=2, dummy=0)

vpred = data.frame(predict.gam(v_gamm$gam, newdata=newdata, se.fit=T))
vpred$lfit = vpred$fit-1.96*vpred$se.fit
vpred$ufit = vpred$fit+1.96*vpred$se.fit

vpred=cbind(newdata, vpred)

if(get_figures){
  par(mfrow=c(1,2))
  
  matplot(vpred$midpoint[vpred$Economic_setting=="Low income"], 
          exp(as.matrix(vpred[vpred$Economic_setting=="Low income",c("fit", "lfit", "ufit")])), 
          ylim=c(0, 0.4), type="l", lty=1)
  points(cfr_long_min$midpoint[cfr_long_min$Economic_setting=="Low income"], 
         cfr_long_min$cfr[cfr_long_min$Economic_setting=="Low income"], pch=20)
  
  matplot(vpred$midpoint[vpred$Economic_setting=="Lower middle income"], 
          exp(as.matrix(vpred[vpred$Economic_setting=="Lower middle income",c("fit", "lfit", "ufit")])), 
          ylim=c(0, .4), type="l", lty=1)
  points(cfr_long_min$midpoint[cfr_long_min$Economic_setting=="Lower middle income"], 
         cfr_long_min$cfr[cfr_long_min$Economic_setting=="Lower middle income"], pch=20)
  
  par(mfrow=c(1,1))
}

# this predicts the basis at new values
tmp=predict.gam(v_gamm$gam, newdata=newdata, type="lpmatrix", se.fit=T)
# coef(v_gamm$gam)
# v_gamm$gam$Vp # seems to be the same as # vcov(v_gamm$gam,unconditional=TRUE)
# gam.check(v_gamm$gam)

somebetas = rmvnorm(n=1000, coef(v_gamm$gam), v_gamm$gam$Vp)
someiterates = (tmp %*% t(somebetas))
allpred_lic = logistic(someiterates[newdata$Economic_setting=="Low income",])
allpred_lmic = logistic(someiterates[newdata$Economic_setting=="Lower middle income",])

if(get_figures){
  par(mfrow=c(1,2))
  matplot(seq(0.1, 60, 0.1), allpred_lic, type="l", ylim=c(0, .1), lty=1, col=rgb(0,0,0,alpha=0.1))
  matplot(seq(0.1, 60, 0.1), allpred_lmic, type="l", ylim=c(0, .1), lty=1, col=rgb(0,0,0,alpha=0.1))
}

allpred_df_lic = data.frame(pred=as.vector(allpred_lic[seq(5, 595, 10),]), 
                            mos=rep(seq(0.5, 59.5, 1), times=dim(allpred_lic)[2]), 
                            iter=rep(1:dim(allpred_lic)[2], each=length(seq(0.5, 59.5, 1))))
allpred_df_lmic = data.frame(pred=as.vector(allpred_lmic[seq(5, 595, 10),]), 
                             mos=rep(seq(0.5, 59.5, 1), times=dim(allpred_lmic)[2]), 
                             iter=rep(1:dim(allpred_lmic)[2], each=length(seq(0.5, 59.5, 1))))

write.csv(allpred_df_lic, file=file.path(plot_prefix, paste0("cfr_lic_predictions_",bs_type,"_n",num_splines,".csv")))
save(allpred_df_lic, allpred_lic, file=file.path(plot_prefix, paste0("cfr_lic_predictions_",bs_type,"_n",num_splines,".Rdata")))

write.csv(allpred_df_lmic, file=file.path(plot_prefix, paste0("cfr_lmic_predictions_",bs_type,"_n",num_splines,".csv")))
save(allpred_df_lmic, allpred_lmic, file=file.path(plot_prefix, paste0("cfr_lmic_predictions_",bs_type,"_n",num_splines,".Rdata")))

################################
# PLOTS
################################
if(get_figures){
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
  
  
  for(i in 1:5){
    jpeg(filename = file.path(plot_prefix, paste0("cfr_lic_lmic_comb_", LETTERS[i], ".jpeg")), 
         width = 7, height = 5.5, units = 'in', res=600)
    print(ggplot(data=cfr_long_min[cfr_long_min$Title %in% unique(cfr_long_min$Title)[((i-1)*6+1):(i*6)],],
                 aes(x=midpoint, y=cfr)) + themebar2 +
            geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +  
            geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
            geom_errorbar(aes(x=midpoint, ymin=cfr_lci, ymax=cfr_uci), width=.25) + 
            xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients") + 
            scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + scale_y_continuous(limits = c(0, 1)))
    dev.off()
  }
  
  for(i in 1:5){
    jpeg(filename = file.path(plot_prefix, paste0("cfr_lic_lmic_comb_", LETTERS[i], "logit.jpeg")), 
         width = 7, height = 5.5, units = 'in', res=600)
    print(ggplot(data=cfr_long_min[cfr_long_min$Title %in% unique(cfr_long_min$Title)[((i-1)*6+1):(i*6)],],
                 aes(x=midpoint, y=logit(cfr))) + themebar2 +
            geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +  
            geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
            geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
            xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients (log-odds scale)") + 
            scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
            scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995)), 
                               labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                          "0.75", "0.95", "0.99", "0.995"), limits = c(-7, 7))) 
    dev.off()
  }
  
  # Low income
  jpeg(filename = file.path(plot_prefix, "cfr_lic.jpeg"), 
       width = 7, height = 9.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Low income",], 
               aes(x=midpoint, y=cfr)) + themebar2 +
          geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +  
          geom_line(data=allpred_df_lic, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.05) +
          geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
          geom_errorbar(aes(x=midpoint, y=cfr, ymin=cfr_lci, ymax=cfr_uci), width=.25) + 
          xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients") + 
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
          scale_y_continuous(limits = c(0, 1)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lic_logit.jpeg"), 
       width = 7, height = 9.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Low income",], 
               aes(x=midpoint, y=logit(cfr))) + themebar2 +
          geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
          geom_line(data=allpred_df_lic, aes(x=mos, y=logit(pred), group=iter), col="lightskyblue", alpha=0.05) + 
          geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
          geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
          xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients (log-odds scale)") + 
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
          scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995)), 
                             labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                        "0.75", "0.95", "0.99", "0.995"), limits = c(-7, 7)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lic_val.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_long_min_all[cfr_long_min_all$Economic_setting=="Low income",], 
               aes(x=midpoint, y=cfr)) + themebar2 +
          geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +
          geom_line(data=allpred_df_lic, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.05) +
          geom_point(size=1) + facet_wrap(~Title, ncol=1) + 
          geom_errorbar(aes(x=midpoint, ymin=cfr_lci, ymax=cfr_uci), width=.25) + 
          xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients (log-odds scale)") + 
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
          scale_y_continuous(limits = c(0, 1)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lic_logit_val.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_long_min_all[cfr_long_min_all$Economic_setting=="Low income",], 
         aes(x=midpoint, y=logit(cfr))) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lic, aes(x=mos, y=logit(pred), group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=1) + 
    geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients (log-odds scale)") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999)), 
                       labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                  "0.75", "0.95", "0.99", "0.999"), limits = c(-7, 7))) 
  dev.off()
  
  # Lower middle income
  # unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])
  # make one image with 10 (3x4) and another with 9 (3x3)
  # cfr_long_min$Title %in% unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])
  
  jpeg(filename = file.path(plot_prefix, "cfr_lmic_A.jpeg"), 
       width = 7, height = 9.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Lower middle income" & 
                             cfr_long_min$Title %in% unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])[1:10],], 
         aes(x=midpoint, y=cfr)) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
    geom_errorbar(aes(x=midpoint, ymin=cfr_lci, ymax=cfr_uci), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + scale_y_continuous(limits = c(0, 1))) 
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lmic_B.jpeg"), 
       width = 7, height = 7.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Lower middle income" & 
                             cfr_long_min$Title %in% unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])[11:19],], 
         aes(x=midpoint, y=cfr)) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
    geom_errorbar(aes(x=midpoint, ymin=cfr_lci, ymax=cfr_uci), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + scale_y_continuous(limits = c(0, 1))) 
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lmic_Alogit.jpeg"), 
       width = 7, height = 9.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Lower middle income" & 
                             cfr_long_min$Title %in% unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])[1:10],], 
         aes(x=midpoint, y=logit(cfr))) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=logit(pred), group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
    geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients (log-odds scale)") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995)), 
                       labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                  "0.75", "0.95", "0.99", "0.999"), limits = c(-7, 7))) 
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_lmic_Blogit.jpeg"), 
       width = 7, height = 7.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min[cfr_long_min$Economic_setting=="Lower middle income" & 
                             cfr_long_min$Title %in% unique(cfr_long_min$Title[cfr_long_min$Economic_setting=="Lower middle income"])[11:19],], 
         aes(x=midpoint, y=logit(cfr))) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=logit(pred), group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=3) + 
    geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients (log-odds scale)") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995)), 
                       labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                  "0.75", "0.95", "0.99", "0.999"), limits = c(-7, 7))) 
  dev.off()
  
  # Validation
  # unique(cfr_long_min_all$Title[cfr_long_min_all$Economic_setting=="Lower middle income"])
  
  jpeg(filename = file.path(plot_prefix, "cfr_lmic_comb_logit_val.jpeg"),
       width = 5, height = 5.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min_all[cfr_long_min_all$Economic_setting=="Lower middle income",], 
         aes(x=midpoint, y=logit(cfr))) + themebar2 +
    theme(panel.grid.minor =  element_blank()) +
    geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=logit(pred), group=iter), col="lightskyblue", alpha=0.05) +
    geom_point(size=1) + facet_wrap(~Title, ncol=2) + 
    geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) + 
    xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients (log-odds scale)") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999)), 
                       labels = c("0.005", "0.01", "0.05", "0.25", "0.50", 
                                  "0.75", "0.95", "0.99", "0.999"), limits = c(-7, 7))) 
  dev.off()
  
  # Validation, no economic stratum
  jpeg(filename = file.path(plot_prefix, "cfr_lic_lmic_comb_logit_val.jpeg"),
       width = 7, height = 5.5, units = 'in', res=600)
  print(ggplot(data=cfr_long_min_all[cfr_long_min_all$Economic_setting %in% c("Low income", "Lower middle income"),],
               aes(x=midpoint, y=logit(cfr))) + themebar2 +
          geom_line(data=allpred_df_lic_lmic, aes(x=mos, y=logit(pred), group=iter), col="rosybrown2", alpha=0.05) +
          geom_point(size=1) + facet_wrap(~Title, ncol=3) +
          geom_errorbar(aes(x=midpoint, ymin=logit(cfr_lci), ymax=logit(cfr_uci)), width=.25) +
          xlab("Age in months") + ylab("Probability of Death among Hospitalised Patients (log-odds scale)") +
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
          scale_y_continuous(breaks=logit(c(0.005, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995)),
                             labels = c("0.005", "0.01", "0.05", "0.25", "0.50",
                                        "0.75", "0.95", "0.99", "0.995"), limits = c(-7, 7)))
  dev.off()
  
  ########################
  ## Plot with ribbons
  ########################
  
  # columns: lowci, hici, mos, econ
  cfr_ribbons_lic = data.frame(t(apply(allpred_lic, 1, quantile, c(0.5, 0.025, 0.975))),
                               apply(allpred_lic, 1, mean))
  colnames(cfr_ribbons_lic) = c("est", "lowci", "hici","mean")
  cfr_ribbons_lic$Economic_setting = "Low income"
  cfr_ribbons_lic$mos = seq(0.1, 60, 0.1)
  
  cfr_ribbons_lmic = data.frame(t(apply(allpred_lmic, 1, quantile, c(0.5, 0.025, 0.975))),
                    apply(allpred_lmic, 1, mean))
  colnames(cfr_ribbons_lmic) = c("est", "lowci", "hici","mean")
  cfr_ribbons_lmic$Economic_setting = "Lower middle income"
  cfr_ribbons_lmic$mos = seq(0.1, 60, 0.1)
  
  cfr_ribbons = rbind(cfr_ribbons_lic, cfr_ribbons_lmic)
  cfr_ribbons$Economic_setting = factor(cfr_ribbons$Economic_setting)
  
  cfr_ribbons_global = data.frame(t(apply(allpred, 1, quantile, c(0.5, 0.025, 0.975))))
  colnames(cfr_ribbons_global) = c("est", "lowci", "hici")
  cfr_ribbons_global$mos = seq(0.1, 60, 0.1)
  
  jpeg(filename = file.path(plot_prefix, "cfr_ribbons_lmic_95CI_mean.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_ribbons_lmic, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
          geom_ribbon(alpha=0.5, size=0, fill = "rosybrown2") +
          geom_ribbon(data=cfr_ribbons_lmic, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
          xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients") +
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
          scale_y_continuous(breaks=seq(0, 0.05, 0.01), limits = c(0, 0.05)))
  dev.off()
  
  
  jpeg(filename = file.path(plot_prefix, "cfr_ribbons_lic_lmic_comb.jpeg"), 
       width = 5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_ribbons, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
    geom_ribbon(data=cfr_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "lightskyblue") +
    geom_ribbon(data=cfr_ribbons_global, aes(x=mos, ymin=est, ymax=est), alpha=0, size=1, col = "blue") +
    geom_ribbon(alpha=0.5, size=0, fill = "burlywood1") +
    geom_ribbon(data=cfr_ribbons_lmic, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
    geom_ribbon(data=cfr_ribbons_lic,  aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
    facet_wrap(~Economic_setting, ncol=2) +
    xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients") +
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
    scale_y_continuous(breaks=seq(0, 0.05, 0.01), limits = c(0, 0.05)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_ribbons_lic_lmic_comb2.jpeg"), 
       width = 5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_ribbons, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
          geom_ribbon(data=cfr_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "lightskyblue") +
          geom_ribbon(alpha=0.5, size=0, fill = "burlywood1") +
          geom_ribbon(data=cfr_ribbons_lmic, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
          geom_ribbon(data=cfr_ribbons_lic,  aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
          facet_wrap(~Economic_setting, ncol=2) +
          xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients") +
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
          scale_y_continuous(breaks=seq(0, 0.05, 0.01), limits = c(0, 0.05)))
  dev.off()
  
  
  jpeg(filename = file.path(plot_prefix, "cfr_ribbons_lic_lmic_comb_old.jpeg"),
       width = 5, height = 3, units = 'in', res=600)
  print(ggplot(data=cfr_ribbons, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
    geom_ribbon(alpha=0.5, size=0, fill = "lightskyblue") +
    geom_ribbon(data=cfr_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici), alpha=0.5, size=0, fill = "rosybrown1") +
    facet_wrap(~Economic_setting, ncol=2) +
    xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients") +
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
    scale_y_continuous(breaks=seq(0, 0.05, 0.01), limits = c(0, 0.05)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "cfr_ribbons_lic_lmic_comb_logit.jpeg"), width = 5, height =3, units = 'in', res=600)
  print(ggplot(data=cfr_ribbons, aes(x=mos, ymin=logit(lowci), ymax=logit(hici))) + themebar2 + 
    theme(panel.grid.minor =  element_blank()) +
    geom_ribbon(alpha=0.5, size=0, fill = "lightskyblue") +
    geom_ribbon(data=cfr_ribbons_global, aes(x=mos, ymin=logit(lowci), ymax=logit(hici)), alpha=0.5, size=0, fill = "rosybrown1") +
    facet_wrap(~Economic_setting, ncol=2) +
    xlab("Age in months") + ylab("Probability of Death among\nHospitalised Patients (log-odds scale)") +
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) +
    scale_y_continuous(breaks=logit(c(0.0025, 0.005, seq(0.01, 0.05, 0.01))),
                       labels = c("0.0025", "0.005","0.01", "0.02", "0.03", "0.04", "0.05"), limits = c(-6, -3)))
  dev.off()
}
