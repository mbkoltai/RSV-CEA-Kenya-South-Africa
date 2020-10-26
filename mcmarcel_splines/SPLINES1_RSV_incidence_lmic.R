#############################################################################
# This file is part of the RSV modelling project.
# 
# => FIT RSV INCIDENCE SPLINES
#
#  Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
# based on code from Marina Antillon
#############################################################################

######################
## Clear workspace
######################

rm(list=ls())

######################
## Useful packages
######################

library("reshape2")
library("plyr")
library("epitools")
library("gamm4")
library("mvtnorm")
library("openxlsx")
library("ggplot2")

############################################
## Settings and output directory
############################################

# set number of splines
num_splines <- 5000

# set random number generator seed
set.seed(20190118)

# create plots?
get_figures <- TRUE

# set regression type: cr = 'cubic regression spline' and ts = 'thin plate regression spline'
bs_type <- "ts" # options: cs, ts

# set plot prefix
plot_prefix = paste0("./splines/output/inc_",bs_type,'_n',num_splines)

# check if the output folder exist, and create the folder if not
if(!dir.exists(plot_prefix)){
  print(paste('create folder:',plot_prefix))
  dir.create(plot_prefix,recursive = T)
}

#########################################################
## Read in cleaned data from Shi et al (by Marina)
#########################################################

# initialize variable
data_list = list()

# add data from the sheets in the excel file
for (i in 1:16){
  data_list[[i]] = read.xlsx("./splines/data/Inc_data_R_readable.xlsx", i)
}

# rbind the sheets 
data_long = data_list[[7]]
data_long = data_long[is.na(data_long$Cases),] # clear data... use only the columns

# convert to a long format
for (i in 1:16){
  data_long = rbind.fill(data_long, data_list[[i]])
}


################################
## Select LMIC data
################################

# select setting
flag      <- data_long$Economic_Setting == 'Lower middle income'
data_long <- data_long[flag,]


################################
## Select non-missing data
################################

# select studies with RSV incidence data
flag      <- !is.na(data_long$Cases)
data_long <- data_long[flag,]


################################
## Assign midpoints
################################

data_long$Age_Groups = c()
for (i in 1:dim(data_long)[1]){
  data_long$Age_Groups[i] = strsplit(as.character(data_long$Age_cat), "Ages_")[[i]][2]
}
data_long$Age_Groups[data_long$Age_Groups=="27d_3m"] = "28d_3m"

agemdpt = data.frame(Age_Groups = c("0d_27d", "28d_3m", "0m_3m", "3m_5m", "1m_5m", "0m_5m", 
              "6m_8m", "9m_11m", "6m_11m", "1m_11m", "0m_11m",
              "12m_23m", "24m_35m", "36m_59m", "0m_23m", "0m_35m", "0m_47m", "12m_59m", 
              "24m_59m", "0m_59m", "1m_59m"),
              midpoint = c(0.5, 2, 1.5, 4.5, 3.5, 3, 7.5, 10.5, 9, 6.5, 6, 18, 30, 48, 12, 18, 
                           24, 36, 42, 30, 30.5))

data_long = join(data_long, agemdpt, type="left")

## Titles for graphs
data_long$Title = paste(data_long$Location, "\n(", data_long$Author, ", ", data_long$Year, ")", 
                       "\n", data_long$Study_period, sep="")

################################################################
# round all cases to nearest whole to obtain 'count data'
################################################################
data_long$Cases = round(data_long$Cases)
data_long$Pop   = round(data_long$Pop)

################################
## calculate incidence and CIs
################################
data_long$inc = data_long$Cases/data_long$Pop*1000
data_long$inc_lci = NaN
data_long$inc_hci = NaN
data_long[, c("inc_lci", "inc_hci")] = 
  pois.exact(data_long$Cases, data_long$Pop, 0.95)[,c("lower", "upper")]*1000

################################
## Poisson regression
################################

data_long_min           <- data_long[!is.na(data_long$Cases),]          # no missing data
data_long_min$study_no  <- as.numeric(as.factor(data_long_min$Title))   # add study number

# make sure each study has at least 3 observations
tbl_study_no            <- data.frame(table(data_long_min$study_no,dnn = 'study_no'))
data_long_min           <- data_long_min[data_long_min$study_no %in% tbl_study_no$study_no[tbl_study_no$Freq>2],]

# rebase study_no 
data_long_min$study_no  <- factor(as.numeric(as.factor(data_long_min$Title)))
table(data_long_min$study_no) 

# add dummy variable 
data_long_min$dummy=1

# fit a poisson model to the current selection
# cases ~ log(age) with offset(log(sample size))
# k: the dimension of the basis used to represent the smooth term.
# bs: (penalized) smoothing basis with cr = 'cubic regression spline' and ts = 'thin plate regression spline'
# offset: added to the linear predictor
# random effect: using study number (!! HAS TO BE A FACTOR !!)
# REML: restricted maximum likelihood
b_gamm = gamm4(Cases~s(log(midpoint), k=-1, bs=bs_type) + 
                 t2(log(midpoint), study_no, k=3, bs=c(bs_type, 're'), by=dummy)+offset(log(Pop)),
               random=~(1|study_no), data=data_long_min, family=poisson, REML=T)

newpred=data.frame(midpoint=seq(0.1, 60, 0.1), Pop=rep(1000, length(seq(0.1, 60, 0.1))),
                   study_no=7, dummy=0)
bpred = data.frame(predict.gam(b_gamm$gam, newdata=newpred, se.fit=T))
bpred$lfit = bpred$fit-1.96*bpred$se.fit
bpred$ufit = bpred$fit+1.96*bpred$se.fit
bpred=cbind(newpred, bpred)
matplot(bpred$midpoint,exp(as.matrix(bpred[,c("fit", "lfit", "ufit")])), type="l", lty=1, ylim=c(0,1000))
points(data_long_min$midpoint, data_long_min$inc, pch=20)

# this predicts the basis at new values
tmp=predict.gam(b_gamm$gam, newdata=newpred, type="lpmatrix", se.fit=T)
# size: number of studies (res), number of knots, plus intercept
# Under ML or REML smoothing parameter estimation it is 
# possible to correct the covariance matrix Vp for smoothing 
# parameter uncertainty. This is the corrected version.


somebetas = rmvnorm(n=num_splines, coef(b_gamm$gam), b_gamm$gam$Vp)
someiterates = (tmp %*% t(somebetas))
allpred = exp(someiterates)*1000
matplot(allpred, type="l", col=rgb(0,0,0,alpha=.1), lty=1, ylim=c(0, 1000))
points(data_long_min$midpoint*10, data_long_min$inc, pch=20, col="red")

allpred_df_lmic   = data.frame(pred=as.vector(allpred[seq(5, 595, 10),]),
                               mos=rep(seq(0.5, 59.5, 1), times=dim(allpred)[2]),
                               iter=rep(1:dim(allpred)[2], each=length(seq(0.5, 59.5, 1))))

# save to file (csv and Rdata)
write.csv(allpred_df_lmic, file=file.path(plot_prefix,paste0("inc_lmic_predictions_",bs_type,"_n",num_splines,".csv")))
save(allpred, allpred_df_lmic, file=file.path(plot_prefix,paste0("inc_lmic_predictions_",bs_type,"_n",num_splines,".Rdata")))

#######################################
# Plot estimates
#######################################
if(get_figures)
{
  newpred_ages <- unique(allpred_df_lmic$mos)

  # modify titles to prevent title-warnings
  data_long_min$Title <- sub('\r\n','',data_long_min$Title)
  
  pdf(file.path(plot_prefix,"inc_lmic_fit_obs.pdf"),12,8)
  par(mfrow=c(2,4))
  i_title <- data_long_min$Title[31]
  for(i_title in unique(data_long_min$Title)){
    flag <- data_long_min$Title == i_title
    plot(data_long_min$midpoint[flag],data_long_min$inc[flag],pch=20,main=i_title,ylim=c(0,200),xlim=c(0,60))
    arrows(data_long_min$midpoint[flag], data_long_min$inc_lci[flag], data_long_min$midpoint[flag], data_long_min$inc_hci[flag], length=0.05, angle=90, code=3)
    boxplot(pred ~ mos, data=allpred_df_lmic,add=T,outline=F,boxcol=2,at=newpred_ages,col=2,xaxt='n')
    points(data_long_min$midpoint[flag],data_long_min$inc[flag],pch=3,lwd=3,col=4)
  }
  dev.off()
  
  #######################################
  # Plot estimates (ggplot2)
  #######################################
  
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
  
  pdf(file = file.path(plot_prefix, paste0('plot_incidence_ggplot2.pdf')),8,4)
  print(ggplot(data=data_long_min, aes(x=midpoint, y=inc)) + themebar2 +
          geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.1) +  
          geom_point(size=1) + facet_wrap(~Title, ncol=4) + 
          geom_errorbar(data=data_long_min, aes(x=midpoint, ymin=inc_lci, ymax=inc_hci), width=.25) + 
          xlab("Age in months") + ylab("Incidence") + 
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, 'inc_lmic_fit_obs.jpeg'),
       width = 8, height = 5, units = "in", # pointsize = 12,
       quality = 100, bg = "white", res = 600)
  print(ggplot(data=data_long_min, aes(x=midpoint, y=inc)) + themebar2 +
    theme(strip.text = element_text(size=6, face="bold")) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) + 
    geom_point(size=1) + facet_wrap(~Title, ncol=4) + 
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    geom_errorbar(aes(x=midpoint, ymin=inc_lci, ymax=inc_hci), width=.25) + 
    xlab("Age in months") + ylab("Cases per 1,000 person-years") + 
    scale_y_continuous(limits = c(0, 1000)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_fit_obs_log.jpeg"), 
       width = 8, height = 5, units = "in", # pointsize = 12,
       quality = 100, bg = "white", res = 600)
  print(ggplot(data=data_long_min, aes(x=midpoint, y=log10(inc+0.1))) + themebar2 +
    theme(strip.text = element_text(size=6, face="bold")) +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=log10(pred), group=iter), col="rosybrown2", alpha=0.05) + 
    geom_point(size=1) + facet_wrap(~Title, ncol=4) + 
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    geom_errorbar(aes(x=midpoint, ymin=log10(inc_lci+0.1), ymax=log10(inc_hci+0.1)), width=.25) +
    xlab("Age in months") + ylab("Cases per 1,000 person-years (log scale)") + 
    scale_y_continuous(breaks=-1:3, labels = c("0.1","1", "10", "100", "1000"), limits = c(-1, 3)))
  dev.off()
  
  allpred_df_lmic_short = allpred_df_lmic[allpred_df_lmic$mos %in% seq(0.5, 12.5, 0.5), ]
  allpred_df_lmic_short$pred[allpred_df_lmic_short$pred<0.1] = 0.1
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_fit_obs_log_babies.jpeg"), 
       width = 8, height = 5, units = "in", # pointsize = 12,
       quality = 100, bg = "white", res = 600)
  print(ggplot(data=data_long_min, aes(x=midpoint, y=log10(inc+0.1))) + themebar2 +
    theme(strip.text = element_text(size=6, face="bold")) +
    geom_line(data=allpred_df_lmic_short, aes(x=mos, y=log10(pred), group=iter), col="rosybrown2", alpha=0.05) + 
    geom_point(size=1) + facet_wrap(~Title, ncol=4) + 
    scale_x_continuous(limits = c(0, 12.5), breaks=seq(0,12,1)) +
    geom_errorbar(aes(x=midpoint, ymin=log10(inc_lci+0.1), ymax=log10(inc_hci)), width=.25) +
    xlab("Age in months") + ylab("Cases per 1,000 person-years (log scale)") + 
    scale_y_continuous(breaks=-1:3, labels = c("0.1", "1", "10", "100", "1000"), limits = c(-1, 3)))
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_fit_obs_val.jpeg"), 
       width = 8, height = 5, units = "in", # pointsize = 12,
       quality = 100, bg = "white", res = 600)
  print(ggplot(data=data_long_min, aes(x=midpoint, y=inc)) + themebar2 +
    geom_line(data=allpred_df_lmic, aes(x=mos, y=pred, group=iter), col="rosybrown2", alpha=0.05) + 
    geom_point(size=1) + facet_wrap(~Title) + 
    scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
    geom_errorbar(aes(x=midpoint, ymin=inc_lci, ymax=inc_hci), width=.25) + 
    xlab("Age in months") + ylab("Cases per 1,000 person-years") + 
    scale_y_continuous(limits = c(0, 1000)))
  dev.off()
  
  #########################
  # Ribbon plots for LMIC 
  #########################
  
  inc_ribbons_global = data.frame(t(apply(allpred, 1, quantile, c(0.5, 0.025, 0.975))),
                                  apply(allpred, 1, mean))
  colnames(inc_ribbons_global) = c("est", "lowci", "hici","mean")
  # cfr_ribbons_global$Economic_setting = "Global"
  inc_ribbons_global$mos = seq(0.1, 60, 0.1)
  inc_ribbons_global$est[inc_ribbons_global$est<0.1] = 0.1
  inc_ribbons_global$lowci[inc_ribbons_global$lowci<0.1] = 0.1
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_ribbons.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=inc_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
    geom_ribbon(alpha=0.5, size=0, fill = "rosybrown2") +
    xlab("Age in months") + ylab("Cases per 1,000 person-years") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_continuous(breaks=seq(0, 200, 25), limits = c(0, 200))) 
  dev.off()
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_ribbons_95CI_mean.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=inc_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
          geom_ribbon(alpha=0.5, size=0, fill = "rosybrown2") +
          geom_ribbon(data=inc_ribbons_global, aes(x=mos, ymin=mean, ymax=mean), alpha=0, size=1, col = "brown" ) +
          xlab("Age in months") + ylab("Cases per 1,000 person-years") + 
          scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
          scale_y_continuous(breaks=seq(0, 200, 25), limits = c(0, 200))) 
  dev.off()
  
  
  
  inc_ribbons_global$lowci[inc_ribbons_global$lowci<1] = 1
  inc_ribbons_global$hici[inc_ribbons_global$hici<1] = 1
  
  jpeg(filename = file.path(plot_prefix, "inc_lmic_ribbons_log10.jpeg"), 
       width = 2.5, height = 3, units = 'in', res=600)
  print(ggplot(data=inc_ribbons_global, aes(x=mos, ymin=lowci, ymax=hici)) + themebar2 +
    geom_ribbon(alpha=0.5, size=0, fill = "rosybrown1") +
    xlab("Age in months") + ylab("Cases per 1,000 person-years (log-scale)") + 
    scale_x_continuous(breaks=seq(0, 60, 12), limits = c(0, 60)) + 
    scale_y_log10(limits = c(1, 550), breaks=c(1, 10, 100, 500), 
                  labels = c("1", "10", "100", "500"))) 
  dev.off()
}

