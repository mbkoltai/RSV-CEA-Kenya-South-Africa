# multiple assignment
# Generic form
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

# Used if LHS is larger than RHS
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

solve_beta_distrib_pars <- function(x,p){
  vacc_eff_mean <- p["mean"]; CI95_down <- p["CI95_low"]; CI95_up <- p["CI95_high"]
  alphaval<-x[1]; shift_val<-x[2]; scale_val<-x[3]; betaval <- alphaval*(1-vacc_eff_mean)/vacc_eff_mean
  y_data <- c(vacc_eff_mean,CI95_down,CI95_up)
  beta_fit <- rbeta(n=1e4,shape1=alphaval,shape2=betaval)*scale_val + shift_val
  y_sol <- as.numeric(abs(c(mean(beta_fit)-vacc_eff_mean,
             quantile(beta_fit,probs=2.5/100)-CI95_down,
             quantile(beta_fit,probs=97.5/100)-CI95_up))) # /abs(y_data)
  if (p["err_type"]>1) {y_sol <- y_sol/abs(y_data)}
  return(sum(y_sol))
}

# fitting effic figures by beta distribution
fcn_betafit_efficacy <- function(effic_figs,scan_range_resol_nsample,optim_range_res,optim_initguess){
  for (interv_type in names(effic_figs)) {
    dis_categ <- intersect(names(effic_figs[[interv_type]]),c("sympt_disease","hospit","severe"))
    for (k_dis in dis_categ) {
      print(c(interv_type,k_dis))
      if (interv_type==names(effic_figs)[1] & k_dis==dis_categ[1]) {
        list_effic_betafit <- list(); best_fit_beta_all <- data.frame() }
      input_var <- effic_figs[[interv_type]][[k_dis]] # ; input_var_orig<-input_var
      # if CI95_low is negative, we need to shift & scale the beta distrib, do this by optim
      if (input_var["CI95_low"]<0) {
        param_estims <- lapply(10^seq(optim_range_res[1],optim_range_res[2],by=optim_range_res[3]), 
                               function(x) optim(par=c(x,optim_initguess[1],optim_initguess[2]),
                                    fn=solve_beta_distrib_pars, p=c(input_var,err_type=1) )$par)
        # means_fit <- lapply(param_estims, rbeta(n=1e4,shape1=x[[1]],shape2=))
        # print(param_estims)
        # vacc_eff_mean <- p["mean"]; CI95_down <- p["CI95_low"]; CI95_up <- p["CI95_high"]; p["err_type"]
        # alphaval<-x[1]; shift_val<-x[2]; scale_val<-x[3]; betaval <- alphaval*(1-vacc_eff_mean)/vacc_eff_mean
  # compare to data
    fitted_distribs <- lapply(param_estims, 
          function(x) rbeta(n=1e5,shape1=x[1],shape2=x[1]*(1-input_var["mean"])/input_var["mean"] )*x[3]+x[2])
        beta_scan_fits <- bind_rows(lapply(fitted_distribs,
        function(x) c(mean=mean(x), CI95_low=quantile(x,probs=2.5/100), CI95_high=quantile(x,probs=97.5/100)))) %>%
          rename(CI95_low=`CI95_low.2.5%`,CI95_high=`CI95_high.97.5%`) %>%
          mutate(alphaval=sapply(param_estims,`[[`,1),betaval=alphaval*(1-input_var["mean"])/input_var["mean"],
             MSE=(mean-input_var["mean"])^2+(CI95_low-input_var["CI95_low"])^2+(CI95_high-input_var["CI95_high"])^2,
             MAE=(abs(mean-input_var["mean"])+abs(CI95_low-input_var["CI95_low"])+
                        abs(CI95_high-input_var["CI95_high"]))/3,interv=interv_type,disease=k_dis,
                 shift_fit=sapply(param_estims,`[[`,2),scale_fit=sapply(param_estims,`[[`,3),)
    # message("alphaval"); print(beta_scan_fits$alphaval)
    # message("betaval"); print(beta_scan_fits$betaval)
    # message("mean"); print(beta_scan_fits$mean)
        # best fit
        best_fit_beta <- beta_scan_fits %>% filter(mean>0 & CI95_high<1) %>% filter(MAE==min(MAE)) %>%
          mutate(data_mean=input_var["mean"],data_CI95_low=input_var["CI95_low"],
                 data_CI95_high=input_var["CI95_high"])
      } else { 
        # scan in parameters to fit
        fitted_beta_distrs <- data.frame(t(
          sapply(10^seq(scan_range_resol_nsample[1],scan_range_resol_nsample[2],by=scan_range_resol_nsample[3]),
                 function(x) c(mean=mean(rbeta(n=scan_range_resol_nsample[4],
                                               shape1=x,shape2=x*(1-input_var["mean"])/input_var["mean"])),
                quantile(rbeta(n=scan_range_resol_nsample[4],
                shape1=x,shape2=x*(1-input_var["mean"])/input_var["mean"]),probs=c(2.5,97.5)/100),
             alphaval=x,betaval=x*(as.numeric((1-input_var["mean"])/input_var["mean"]) ) ) ))) %>%
          rename(CI95_low=`X2.5.`,CI95_high=`X97.5.`) %>%
          mutate(MSE=((input_var["CI95_low"]-CI95_low)^2+(input_var["CI95_high"]-CI95_high)^2+
                        (input_var["mean"]-mean)^2)/3,
                 MAE=(abs(input_var["CI95_low"]-CI95_low)+abs(input_var["CI95_high"]-CI95_high)+
                        abs(input_var["mean"]-mean))/3)
        # best fit
        best_fit_beta <- fitted_beta_distrs %>% filter(MAE==min(MAE)) %>% 
          mutate(shift_fit=0,scale_fit=1,interv=interv_type,disease=k_dis,
                 data_mean=input_var["mean"],data_CI95_low=input_var["CI95_low"],
                 data_CI95_high=input_var["CI95_high"])
        }
      # collect param estimates
      if (is.null(nrow(best_fit_beta_all))) {best_fit_beta_all <- best_fit_beta} else {
        best_fit_beta_all <- bind_rows(best_fit_beta_all,best_fit_beta) }
    list_effic_betafit[[interv_type]][[k_dis]] <- c(shape1=best_fit_beta$alphaval,shape2=best_fit_beta$betaval,
                                                      shift_fit=best_fit_beta$shift_fit,
                                                      scale_fit=best_fit_beta$scale_fit)
    }
  }
  list(list_effic_betafit,best_fit_beta_all)
}

### calc rates of exp decay
fcn_exp_waning_rate <- function(effic_list,n_row){
  exp_waning_param<-data.frame()
  for (k_interv in 1:length(effic_list)) {
    alpha_par<-log(2)/effic_list[[k_interv]]$half_life; n_dur=effic_list[[k_interv]]$duration
    for (k_dis in 1:3){
      if (names(effic_list[[k_interv]])[k_dis] %in% c("sympt_disease","hospit","severe")){
        mean_VE_trial<-as.numeric(effic_list[[k_interv]][[k_dis]]["mean"])
        c_const<-mean_VE_trial*alpha_par*n_dur/(1-exp(-n_dur*alpha_par))
        mean_VE_calc<-mean(c_const*exp(-alpha_par*((1:n_row)-1))[1:n_dur])
        if (c_const>1) {
          mean_VE_calc <- mean(c_const*exp(-alpha_par*((1:n_row)-1))[1:n_dur])
          while (abs(mean_VE_calc-mean_VE_trial)>0.001) {
            c_const<-1 
            # alpha_par*5/(1-exp(-5*alpha_par)) # # exp_scaling = -(1/(dur_prot_infant-1))*log(2)*((1:n_row)-1)
            exp_decay<-c_const*exp(-alpha_par*((1:n_row)-1)); mean_VE_calc <- mean(exp_decay[1:n_dur])
            alpha_par<-alpha_par*as.numeric(ifelse(mean_VE_trial>mean_VE_calc,0.99,1.01)) } # while
        }  else {# if c_const>1
          if (abs(mean_VE_calc-mean_VE_trial)>0.1/100) {c_const=c_const*(mean_VE_trial/mean_VE_calc);
          mean_VE_calc<-mean(c_const*exp(-alpha_par*((1:n_row)-1))[1:n_dur])} }
        # print(c(names(effic_list)[k_interv],names(effic_list[[k_interv]])[k_dis]))
        exp_w_output<-data.frame(interv=names(effic_list)[k_interv],dis_type=names(effic_list[[k_interv]])[k_dis],
                                 "c_const"=c_const,"exp_decay_rate"=alpha_par,"half_life_month"=log(2)/alpha_par,"mean_VE_calcul"=mean_VE_calc)
        if (k_interv==1&k_dis==1) {exp_waning_param<-exp_w_output} else { exp_waning_param <- rbind(exp_waning_param,exp_w_output)  }
      } # only for disease types
    } # k_dis
  }
  # list_exp_decay
  list_exp_decay <- list(mat_vacc=list(sympt_disease=exp_waning_param %>% 
                                         filter(interv=="mat_vacc"&dis_type=="sympt_disease") %>% 
          dplyr::select(c_const,exp_decay_rate),
        hospit=exp_waning_param %>% filter(interv=="mat_vacc" & dis_type=="hospit") %>% 
          dplyr::select(c_const,exp_decay_rate),
        severe=exp_waning_param %>% filter(interv=="mat_vacc" & dis_type=="severe") %>% 
          dplyr::select(c_const,exp_decay_rate)),
      monocl_ab=list(sympt_disease=exp_waning_param %>% filter(interv=="monocl_ab" & dis_type=="sympt_disease") %>% 
          dplyr::select(c_const,exp_decay_rate),
        hospit=exp_waning_param %>% filter(interv=="monocl_ab" & dis_type=="hospit") %>% 
                                          dplyr::select(c_const,exp_decay_rate)))
  
  df_exp_waning_param <- data.frame(t(data.frame(list_exp_decay))) %>% rownames_to_column %>% 
    separate(rowname,c("interv","dis_type","param"),"\\.") %>% mutate(across(where(is.numeric),round,4))
  colnames(df_exp_waning_param)[ncol(df_exp_waning_param)]="value"
  df_exp_waning_param <- df_exp_waning_param %>% pivot_wider(names_from=param)
  list(list_exp_decay,df_exp_waning_param)
}

# # MV efficacy
# ci95_vals<-c(0.216,0.976); rate_treat<-3/405; rate_plac<-5/103; vacc_eff<-1-rate_treat/rate_plac
# # mAb efficacy
# ci95_vals<-c(0.496,0.871); rate_treat<-12/994; rate_plac<-5/103; vacc_eff<-1-rate_treat/rate_plac
# #
# alpha_inputs <- 10^seq(log10(5),1.5,by=1/50) # c(0.1,0.2,10^seq(-1/2,1/2,by=1/20),4,5,8,10)
# betafit_ci95 <- data.frame(t(sapply(alpha_inputs, function(x)
#     c(quantile(rbeta(n=1e5,shape1=x,shape2=x*(1-vacc_eff)/vacc_eff),probs=c(2.5,97.5)/100),
#     alpha=x,beta=x*(1-vacc_eff)/vacc_eff,mean=mean(rbeta(n=1e6,shape1=x,shape2=x*(1-vacc_eff)/vacc_eff))) ))) %>% 
#   rename(ci95_low=`X2.5.`,ci95_up=`X97.5.`) %>% 
#   mutate(MSE=((ci95_low-ci95_vals[1])^2+(ci95_up-ci95_vals[2])^2)/2,
#          MAE=(abs(ci95_low-ci95_vals[1])+(ci95_up-ci95_vals[2]))/2,
#          MPE=(abs(ci95_low-ci95_vals[1])+(ci95_up-ci95_vals[2]))/2 )
# # plot
# ggplot(betafit_ci95 %>% select(c(alpha,MSE,MAE,MPE)) %>% pivot_longer(!alpha),
#        aes(x=alpha,y=value,color=name)) + geom_line() + geom_point() + facet_wrap(~name,scales="free_y",nrow=3) +
#   scale_x_log10(breaks=round(alpha_inputs,2)) + ylab("SSE(CI95_modeled,CI95_trial)") + 
#   theme_bw() + standard_theme + theme(axis.text.x=element_text(hjust=1,vjust=1/2))
# # extract best estim
# betadist_pars <- as.numeric((betafit_ci95 %>% filter(MSE==min(MSE)))[c("alpha","beta")])
# mean(rbeta(n=1e4,shape1=betadist_pars[1],shape2=betadist_pars[2]))
# quantile(rbeta(n=1e4,shape1=betadist_pars[1],shape2=betadist_pars[2]),probs=c(2.5,97.5)/100)
# 
# # efficacy for severe LRTI has negative value for CI95 lower bound
# vacc_eff_sev_LRTI <- 0.915; vacc_eff_sev_LRTI_ci95 <- c(-5.6,99.8)/100
# 
# mean(rsn(n=1e3, xi=0.915, omega=1, alpha=0))
# quantile(rsn(n=1e3, xi=0.915, omega=1, alpha=0),probs=c(2.5,97.5)/100)
# 
# solve_skewnorm_distr_pars <- function(x,p){
#   mean_sn <- mean(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]))
#   ci95_down <- quantile(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]),probs=c(2.5,97.5)/100)[1]
#   ci95_up <- quantile(rsn(n=1e3, xi=x[1], omega=x[2], alpha=x[3]),probs=c(2.5,97.5)/100)[2]
#   return(c(p[1]-mean_sn,p[2]-ci95_down,p[3]-ci95_up))
# }
# 
# # resol<-seq(-4,(-2.5),by=1/10)
# par_combs <- expand.grid((0:10)/10,seq(-1,1,by=1/10),(-3:2)*6)
# for (k_p in 1:nrow(par_combs)) {
#   sol_beta_fit <- nleqslv(x=as.numeric(par_combs[k_p,]),
#                           fn=solve_skewnorm_distr_pars,p=c(vacc_eff_sev_LRTI,vacc_eff_sev_LRTI_ci95))
#   rsn_out <- rsn(n=1e4, xi=sol_beta_fit$par[1], omega=sol_beta_fit$par[2], alpha=sol_beta_fit$par[3]) # nleq_output
#   mean_sn <- mean(rsn_out); ci95_down <- quantile(rsn_out,probs=c(2.5,97.5)/100)[1]
#   ci95_up <- quantile(rsn_out,probs=c(2.5,97.5)/100)[2]
#   nleq_output <- c(sol_beta_fit$par,mean_sn,ci95_down,ci95_up,sol_beta_fit$termcd,sum(sol_beta_fit$fvec^2))
#   if (k_p==1) { sol_df <- data.frame(t(nleq_output)) } else {
#     sol_df <- bind_rows(sol_df,data.frame(t(nleq_output)))
#   }
#   if (k_p %% 100 == 0) {print(k_p)}
# }
# colnames(sol_df) <- c(paste0("x",1:3),"mean_sn","ci95_down","ci95_up","exit_flag","SSE")
# sol_df <- sol_df %>% rownames_to_column(var = "n")
# 
# sol_df %>% filter(SSE<0.5) %>% select(c(x1,x2,x3)) %>% pivot_longer(cols = c(x1,x2,x3)) %>%
#   ggplot(aes(value)) + geom_histogram() + facet_wrap(~name,scales = "free") + theme_bw()
#
# 
# # VE=1-IRR, IRR=illness_vacc_group/illness_placebo_group
# n=405+103; y=3+5
# rate_treat<-3/405; rate_plac<-5/103
# vacc_eff<-1-rate_treat/rate_plac
# theta_ve<-(1-vacc_eff)/(2-vacc_eff)
# # theta is the relative rate of infx in the treatment group 
# theta_ve <- rate_treat/(rate_treat+rate_plac)
# # prior: assuming 30% VE and uninformative beta=1
# efficacy_prior<-3/10; beta_prior<-1
# # theta = (1-VE)/(1-2*VE)
# theta_prior<-(1-efficacy_prior)/(2-efficacy_prior)
# # theta_prior is the mean of beta(alpha_prior,beta_prior=1) -> theta_prior=alpha_prior/(alpha_prior+1)
# alpha_prior<-theta_prior/(1-theta_prior)
# # prior: beta(alpha_prior,beta_prior)
# (theta_prior^(alpha_prior-1))*(1-theta_prior)^(beta_prior-1)
# theta_prior^(-0.3) # *k proportionality factor
# # LIKELIHOOD: p(y|theta) ~ binom_coeff*(theta^y)*(1-theta)^(n-y) (binomial sampling distrib)
# # choose(n,y)*(theta^y)*(1-theta)^(n-y)
# # posterior = prior*likelihood= Beta(alpha_prior + y, beta_prior + n-y)
# # (theta^(alpha_prior-1+y))*(1-theta)^(n-y+beta_prior-1)
# theta_vals<-seq(0.01,0.99,by=1/1e2) # rm(theta)
# theta_posterior_pdf<-(theta_vals^(alpha_prior-1+y))*(1-theta_vals)^(n-y+beta_prior-1)
# # CI95 for theta
# # quantile((1-theta_posterior_pdf)/(1-theta_posterior_pdf),probs=c(2.5,97.5)/100)
# theta_ci95<-qbeta(c(0.025, 0.975), shape1=alpha_prior+y, shape2 = beta_prior + n-y)
# ve_ci95<-(1-2*theta_ci95[2:1])/(1-theta_ci95[2:1])
# 
# # quantiles
# quantile(theta_posterior_pdf,probs=c(5,95)/100)
# # integral of beta distribution:
# (x^alpha)*dhyper(alpha,1-betaval,alpha+1,x)/(alpha*beta(alpha,betaval))
# obj function for integral
# solve_beta_distrib_pars <- function(x,p){
#   vacc_eff_val <- p[3]
#   alphaval<-x; betaval<-alphaval*(1-vacc_eff_val)/vacc_eff_val
#   # objective fcn is the integral
#   norm_const <- beta(alphaval,betaval)
#   integral_CI95_up <- (p[2]^alphaval)*phyper(m=alphaval,n=1-betaval,k=alphaval+1,q=p[2])/(alphaval*norm_const)
#   integral_CI95_down <- (p[1]^alphaval)*phyper(m=alphaval,n=1-betaval,k=alphaval+1,q=p[1])/(alphaval*norm_const)
#   print(c(integral_CI95_up,integral_CI95_down))
#   return(integral_CI95_up - integral_CI95_down - 0.95)
# }
