# VE=1-IRR, IRR=illness_vacc_group/illness_placebo_group
n=405+103; y=3+5
rate_treat<-3/405; rate_plac<-5/103
vacc_eff<-1-rate_treat/rate_plac
theta_ve<-(1-vacc_eff)/(2-vacc_eff)
# theta is the relative rate of infx in the treatment group 
theta_ve <- rate_treat/(rate_treat+rate_plac)
# prior: assuming 30% VE and uninformative beta=1
efficacy_prior<-3/10; beta_prior<-1
# theta = (1-VE)/(1-2*VE)
theta_prior<-(1-efficacy_prior)/(2-efficacy_prior)
# theta_prior is the mean of beta(alpha_prior,beta_prior=1) -> theta_prior=alpha_prior/(alpha_prior+1)
alpha_prior<-theta_prior/(1-theta_prior)
# prior: beta(alpha_prior,beta_prior)
(theta_prior^(alpha_prior-1))*(1-theta_prior)^(beta_prior-1)
theta_prior^(-0.3) # *k proportionality factor
# LIKELIHOOD: p(y|theta) ~ binom_coeff*(theta^y)*(1-theta)^(n-y) (binomial sampling distrib)
# choose(n,y)*(theta^y)*(1-theta)^(n-y)
# posterior = prior*likelihood= Beta(alpha_prior + y, beta_prior + n-y)
# (theta^(alpha_prior-1+y))*(1-theta)^(n-y+beta_prior-1)
theta_vals<-seq(0.01,0.99,by=1/1e2) # rm(theta)
theta_posterior_pdf<-(theta_vals^(alpha_prior-1+y))*(1-theta_vals)^(n-y+beta_prior-1)
# CI95 for theta
# quantile((1-theta_posterior_pdf)/(1-theta_posterior_pdf),probs=c(2.5,97.5)/100)
theta_ci95<-qbeta(c(0.025, 0.975), shape1=alpha_prior+y, shape2 = beta_prior + n-y)
ve_ci95<-(1-2*theta_ci95[2:1])/(1-theta_ci95[2:1])

# quantiles
quantile(theta_posterior_pdf,probs=c(5,95)/100)
# integral of beta distribution:
(x^alpha)*dhyper(alpha,1-betaval,alpha+1,x)/(alpha*beta(alpha,betaval))
# obj function for integral
solve_beta_distrib_pars <- function(x,p){
  vacc_eff_val <- p[3]
  alphaval<-x; betaval<-alphaval*(1-vacc_eff_val)/vacc_eff_val
  # objective fcn is the integral
  norm_const <- beta(alphaval,betaval)
  integral_CI95_up <- (p[2]^alphaval)*phyper(m=alphaval,n=1-betaval,k=alphaval+1,q=p[2])/(alphaval*norm_const)
  integral_CI95_down <- (p[1]^alphaval)*phyper(m=alphaval,n=1-betaval,k=alphaval+1,q=p[1])/(alphaval*norm_const) 
  print(c(integral_CI95_up,integral_CI95_down))
  return(integral_CI95_up - integral_CI95_down - 0.95)
}
# library(nleqslv) # control=list(btol=0.01),
sol_beta_fit <- nleqslv(x=2,fn=solve_beta_distrib_pars,p=c(0.216,0.976,vacc_eff)) 

# MV efficacy
ci95_vals<-c(0.216,0.976); rate_treat<-3/405; rate_plac<-5/103; vacc_eff<-1-rate_treat/rate_plac
# mAb efficacy
ci95_vals<-c(0.496,0.871); rate_treat<-12/994; rate_plac<-5/103; vacc_eff<-1-rate_treat/rate_plac
#
alpha_inputs <- 10^seq(log10(5),1.5,by=1/50) # c(0.1,0.2,10^seq(-1/2,1/2,by=1/20),4,5,8,10)
betafit_ci95 <- data.frame(t(sapply(alpha_inputs, function(x)
    c(quantile(rbeta(n=1e5,shape1=x,shape2=x*(1-vacc_eff)/vacc_eff),probs=c(2.5,97.5)/100),
    alpha=x,beta=x*(1-vacc_eff)/vacc_eff,mean=mean(rbeta(n=1e6,shape1=x,shape2=x*(1-vacc_eff)/vacc_eff))) ))) %>% 
  rename(ci95_low=`X2.5.`,ci95_up=`X97.5.`) %>% 
  mutate(MSE=((ci95_low-ci95_vals[1])^2+(ci95_up-ci95_vals[2])^2)/2,
         MAE=(abs(ci95_low-ci95_vals[1])+(ci95_up-ci95_vals[2]))/2,
         MPE=(abs(ci95_low-ci95_vals[1])+(ci95_up-ci95_vals[2]))/2 )
# plot
ggplot(betafit_ci95 %>% select(c(alpha,MSE,MAE,MPE)) %>% pivot_longer(!alpha),
       aes(x=alpha,y=value,color=name)) + geom_line() + geom_point() + facet_wrap(~name,scales="free_y",nrow=3) +
  scale_x_log10(breaks=round(alpha_inputs,2)) + ylab("SSE(CI95_modeled,CI95_trial)") + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(hjust=1,vjust=1/2))
# extract best estim
betadist_pars <- as.numeric((betafit_ci95 %>% filter(MSE==min(MSE)))[c("alpha","beta")])
mean(rbeta(n=1e4,shape1=betadist_pars[1],shape2=betadist_pars[2]))
quantile(rbeta(n=1e4,shape1=betadist_pars[1],shape2=betadist_pars[2]),probs=c(2.5,97.5)/100)

# efficacy for severe LRTI has negative value for CI95 lower bound
vacc_eff_sev_LRTI <- 0.915; vacc_eff_sev_LRTI_ci95 <- c(-5.6,99.8)/100

mean(rsn(n=1e3, xi=0.915, omega=1, alpha=0))
quantile(rsn(n=1e3, xi=0.915, omega=1, alpha=0),probs=c(2.5,97.5)/100)

solve_skewnorm_distr_pars <- function(x,p){
  mean_sn <- mean(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]))
  ci95_down <- quantile(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]),probs=c(2.5,97.5)/100)[1]
  ci95_up <- quantile(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]),probs=c(2.5,97.5)/100)[2]
  return(c(p[1]-mean_sn,p[2]-ci95_down,p[3]-ci95_up))
}

# resol<-seq(-4,(-2.5),by=1/10)
par_combs <- expand.grid((0:10)/10,seq(-1,1,by=1/10),(-3:2)*6)
for (k_p in 1:nrow(par_combs)) {
  sol_beta_fit <- nleqslv(x=as.numeric(par_combs[k_p,]),
                          fn=solve_skewnorm_distr_pars,p=c(vacc_eff_sev_LRTI,vacc_eff_sev_LRTI_ci95))
  rsn_out <- rsn(n=1e4, xi=sol_beta_fit$x[1], omega=sol_beta_fit$x[2], alpha=sol_beta_fit$x[3]) # nleq_output
  mean_sn <- mean(rsn_out); ci95_down <- quantile(rsn_out,probs=c(2.5,97.5)/100)[1]
  ci95_up <- quantile(rsn_out,probs=c(2.5,97.5)/100)[2]
  nleq_output <- c(sol_beta_fit$x,mean_sn,ci95_down,ci95_up,sol_beta_fit$termcd,sum(sol_beta_fit$fvec^2))
  if (k_p==1) { sol_df <- data.frame(t(nleq_output)) } else {
    sol_df <- bind_rows(sol_df,data.frame(t(nleq_output)))
  }
  if (k_p %% 100 == 0) {print(k_p)}
}
colnames(sol_df) <- c(paste0("x",1:3),"mean_sn","ci95_down","ci95_up","exit_flag","SSE")
sol_df <- sol_df %>% rownames_to_column(var = "n")

sol_df %>% filter(SSE<0.5) %>% select(c(x1,x2,x3)) %>% pivot_longer(cols = c(x1,x2,x3)) %>%
  ggplot(aes(value)) + geom_histogram() + facet_wrap(~name,scales = "free") + theme_bw()