# Pre-specified change point locations ------------------------------------
###
# For one setting
###

### Generate data
## Input:
# B: number iterations
# N: number subjects per group
# snr: Signal-to-noise ratio
# effect: effect size, i.e., Bar(e), for each condition (vector of two)
# hrf_basis_fun: a matrix containing the basis functions used, in our case: 3 FLOBS basis functions (available in FSL)
#               each row represents one basis function.
## Remark:
# We demonstrate for snr=2 and effect=c(1,1.5)
# the hrf_basis_fun have been defined using "Make_FLOBS" with pre-defined values in FSL, for more information see: http://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/centos6/fsl_507/doc/wiki/FLOBS.html
require(parallel)
require(doParallel)
require(foreach)
B <- 1000
env.aic <- new.env()
env.aic$hrf_basis_fun <- hrf_basis_fun
env.aic$N <- 30
env.aic$snr <- 2
env.aic$effect <- c(1,1.5)

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
clusterExport(cl,
              c("hrf_basis_fun", "N","snr", "effect"), envir = env.aic)
clusterExport(cl, c("glm_prewhitened", "cor_res", "pre_whiten_BV", "simulate_BOLD", "est_betas_var", "est_var_sp_subject",
                    "var_est_subject", "reml_var_estimation", "fmri_group_KH"))
snr_2_e_1_15 <- foreach(i= 1:B, .packages=c("neuRosim", "MASS"))%dopar%{
  sim_known_cp(effect=effect, N=N, snr1=snr, seed=i, basis_function=hrf_basis_fun, var_change=c(0.1), beta_basis =c(3.2, -6.4, 3.2),misspec = 0)
}
stopCluster(cl)

### Analysis of generated data
## Remark:
# get_TreeBH_selections: function to control the sFDR, available at: https://github.com/cbpeterson/TreeBH
# slighly updated to account for arbitrary results ()
e_1_15_snr2 <- known_cp_results(dat_considered = snr_2_e_1_15, snr1=2, eff=c(1,1.5), B=1000)


# Pre-specified change point locations: Graphics ------------------------------------
###
# Figure 3: Rejection rates
###

## Input:
# hrf_basis_fun: The (FLOBS) basis functions used for the analysis
# rejected_means_full: data frame composed of "rejected_means" data frames (i.e., known_cp_results$rejected_means) for all considered settings

require(ggplot2)
require(data.table)
## Initialize data frame as input for ggplot:
rejection_long <- melt(setDT(rejected_means_full),
                       id.vars = c("snr", "effect", "test", "mis"),
                       variable.name = "Shape_Parameter")
rejection_long$Shape_Parameter <- as.character(rejection_long$Shape_Parameter)
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="PM_r"] <- "PM"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="NA_r"] <- "NA"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="TTP_r"] <- "TTP"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="TPN_r"] <- "TPN"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="FWHM_r"] <- "FWHM"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="FWHN_r"] <- "FWHN"
rejection_long$Shape_Parameter[rejection_long$Shape_Parameter=="AUC_r"] <- "AUC"
rejection_long$SP_fac <- as.factor(rejection_long$Shape_Parameter)
rejection_long$SNR <- factor(rejection_long$snr, levels=c(1, 2),
                             ordered=T, labels=c(expression(paste(SNR, "=1")),expression(paste(SNR, "=2"))))
rejection_long$Mis <- factor(rejection_long$mis, levels=c(0, 5),
                             ordered=T, labels=c(expression(paste(abs(hat(Psi)-Psi),"=0")),expression(paste(abs(hat(Psi)-Psi),"<=5"))))
## Graphic:
ggplot(rejection_long, aes(x=effect, y=value,group=Shape_Parameter))+
  theme_bw()+
  theme(text=element_text(size=20),
        legend.key.size=unit(4, "line"))+
  geom_line(aes(linetype=Shape_Parameter, color=Shape_Parameter),linewidth=1)+
  scale_color_manual(breaks= c("PM", "NA", "TTP", "TPN", "FWHM", "FWHN", "AUC"),
                     values=c("PM"="cyan1", "NA"="cyan3", "TTP"="pink", "TPN"="pink2", "FWHM"="pink3","FWHN"="pink4", "AUC"="cyan4"))+
  scale_linetype_manual(breaks= c("PM", "NA", "TTP", "TPN", "FWHM", "FWHN","AUC"),
                        values=c("PM"=1, "NA"=2, "TTP"=5, "TPN"=4, "FWHM"=3,"FWHN"=18, "AUC"=88))+ 
  ylab("Rejection rates")+
  xlab("Mean effect size")+
  facet_grid(test~ Mis+SNR, labeller=label_parsed)+
  coord_cartesian(ylim=c(0,1))



# Unknown change point locations ------------------------------------------
###
# For one setting
###

### Generate data
## Input:
# Bn: the number of subjects \mathcal{N}= Bn*Np
# Np: number of subjects to be generated per iteration of the for-loop
# ncl: number of clusters to be used in parallel computation
## Remark:
# Since we use a time limit when computing the post-selection variances, some iterations will return NULL
# Therefore, Bn might need to be adjusted to account for that.
# Using Np<\mathcal{N} (i.e., Bn>1) and merging the newly simulated data with the already simulated data allows for saving in between each iteration of the for-loop. 
require(parallel)
require(doParallel)
require(foreach)

Np <- 30
set.seed(40)
env.aic <- new.env()
ncl <- 6
cd_snr_2_e_1_full <- vector("list",0)
for(r in 1:Bn){
  sim_BOLD <- sim_BOLD_canonical_one_condition(mean_change=1,
                                               var_change=0.1,
                                               nscan=250, nstim=60,rho_noise=0.2,
                                               min_stat=10, snr=2, N=Np)
  env.aic$sim_BOLD <- sim_BOLD
  # initialize using multiple kernels:
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  clusterExport(cl,
                c("sim_BOLD"), envir = env.aic)
  clusterExport(cl, c("glm_prewhitened", "cor_res", "pre_whiten_BV", "sim_BOLD_canonical_one_condition", 
                      "est_betas_cd", "suff_stat_posi","tryCatch_sigma","U_star_posi","suff_stat_star_posi",
                      "loglik_known_var_glm","loglik_known_var_direct","loglik_unknown_var_glm", "loglik_unknown_var_direct",
                      "selected_model", "t_obs_star_posi", "t_obs_star_posi_seed", "F_t_given_u_posi_2"))
  cd_snr_2_e_1 <- foreach(i= 1:Np, .packages=c("neuRosim", "MASS", "nnet"))%dopar%{
    func_parallel_cd_tl(BOLD_signal= sim_BOLD$BOLD[i,], sigma_w=sim_BOLD$sigma2_w[i], rho_noise=sim_BOLD$rho, 
                        onsets=sim_BOLD$onsets, nscan=sim_BOLD$nscan, nstim=sim_BOLD$nstim, min_seq=10, last.length.theta.star=50,
                        true_cp_loc=sim_BOLD$true_cp_loc[i],Var_known=T)
  }
  stopCluster(cl)
  #remove the NULL elements that are due to the time limit being reached:
  id.null <- which(sapply(cd_snr_2_e_1, is.null))
  if(length(id.null)!=0){
    cd_snr_2_e_1 <- cd_snr_2_e_1[-id.null]
  }
  # merge the newly simulated data with the already simulated data
  if(length(cd_snr_2_e_1)!=0){
    cd_snr_2_e_1_full <- c(cd_snr_2_e_1_full_new,cd_snr_2_e_1)
  }
  # If wanted: save cd_snr_2_e_1_full here
}

### Analysis of generated data
e1_snr2 <- unknown_cp_results(dat_cons=cd_snr_2_e_1_full, eff=1, snr1=2, B=1000, N=30)

# Unknown change point locations: Graphics --------------------------------
library(ggplot2)
### Rejection rates
## Input:
# summary_results_full: merged dataframes "summary" from  function "unknown_cp_results"; across the simulated scenarios 
library(data.table)
rejection_long <- melt(setDT(summary_results_full),
                       id.vars = c("snr", "effect", "procedure", "mean_est_eta", "var_est_eta", "mean_est_var_B", "var_est_var_B"),
                       variable.name = "test")
rejection_long$Test[rejection_long$test=="rejection_rate_KH"] <- "KH"
rejection_long$Test[rejection_long$test=="rejection_rate_W"] <- "Wald"
rejection_long$Test <- factor(rejection_long$Test, levels = c("KH", "Wald"))
rejection_long$Estimates <- factor(rejection_long$procedure, levels=c("naive", "posi_05", "posi_E", "posi_OLS"))
# Plot of the rejection rates:
rej_rate <- ggplot(rejection_long, aes(x=effect, y=value,group=Estimates))+
  theme_bw()+
  theme(text=element_text(size=25),legend.key.size=unit(4, "line"))+
  geom_line(aes(linetype=Estimates, color=Estimates),linewidth=1)+
  scale_color_manual(breaks= c("naive", "posi_05", "posi_E", "posi_OLS"),
                     values=c("naive"="cyan4", "posi_05"="seagreen2", "posi_E"="seagreen", "posi_OLS"="cyan2"))+ 
  scale_linetype_manual(breaks= c("naive", "posi_05", "posi_E", "posi_OLS"),
                        values=c("naive"=1, "posi_05"=3, "posi_E"=4, "posi_OLS"=2))+ #legend in linetype, do not show m_1
  scale_x_continuous(expression(eta), labels = as.character(c(0, 0.5, 1)), breaks = c(0,0.5,1))+
  ylab("Rejection rates")+
  facet_wrap(vars(Test))+
  coord_cartesian(ylim=c(0,0.5))

### Boxplot comparing the estimated etas (naive, posi_05, posi_E, posi_OLS)
## Input:
# results_full: merged dataframes "results" from  function "unknown_cp_results"; across the simulated scenarios 
dat_bp_eta <- melt(setDT(results_full),
                   id.vars=c("name", "snr", "effect"),
                   variable.name="estimates")
dat_bp_eta$Estimates <- factor(dat_bp_eta$estimates, levels=c("naive", "posi_05", "posi_E", "posi_OLS"))

id.eta <- which(dat_bp_eta$name=="est_eta")
bp_eta <- ggplot()+
  theme_bw()+
  theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5))+
  geom_boxplot(data=dat_bp_eta[id.eta,], 
               aes(x=effect,
                   y=value,
                   group=interaction(effect,Estimates),
                   fill=Estimates),
               outlier.shape = 21, outlier.fill = NULL)+
  scale_x_continuous(expression(eta), labels = as.character(c(0, 0.5, 1)), breaks = c(0,0.5,1))+
  ylab(expression(hat(eta)))+
  scale_fill_manual(values=c("cyan4", "seagreen2", "seagreen", "cyan2"))
bp_eta
ggsave(filename = "BP_eta_unknown_cp.png",height=10, width=25, units="cm", dpi=1000)
ggsave(filename = "BP_eta_unknown_cp_low.png",height=10, width=25, units="cm", dpi=300)

### Boxplot comparing the estimated within subject variances? (posi vs. naive)
## Input:
# est_var_W_dat: dataframe summarizing "est_var_W" from the function "unknown_cp_results" across all scenarios
#                the columnes of the data frame are values: est_var_W
#                                                   snr: the SNR corresponding to the est_var_W entries
#                                                   effect: eta corresponding to the est_var_W entries
#                                                   procedure: either "naive" or "posi"
est_var_W_dat$Estimates <- factor(est_var_W_dat$procedure, levels=c("naive", "posi"))
bp_var_W <- ggplot()+
  theme_bw()+
  theme(text=element_text(size=25), plot.title = element_text(hjust = 0.5))+
  geom_boxplot(data=est_var_W_dat, 
               aes(x=effect,
                   y=values,
                   group=interaction(effect,Estimates),
                   fill=Estimates),
               outlier.shape = 21, outlier.fill = NULL)+
  scale_x_continuous(expression(eta), labels = as.character(c(0, 0.5, 1)), breaks = c(0,0.5,1))+
  ylab(expression(hat(sigma)[hat(beta)]^2))+
  scale_fill_manual(values=c("cyan4", "seagreen2"))
