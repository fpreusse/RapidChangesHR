# Simulation study functions ----------------------------------------------
####
# Pre-specified change point location
####

### Simulate BOLD signal; one change point per condition; AR(1) process
## Input:
# mean_change: bar(e), i.e., the mean difference (across the sample) of the change, corresponds to eta on group level
# var_change: the variance of changes within the sample, i.e., between subject variance
# basis_func: matrix containing the HRF's basis functions, each column= one basis function
# basis_beta: value of the beta(s) have before the change? Defines the shape of the HR before the change point. Is equal for both conditions
# rho: AR(1) parameter
# snr: Signal-to-noise ratio; defines the within-subject variance
# N: number of subjects in a sample
# nscan: number of observations
# nstim: number of stimuli for both conditions; the number of stimuli per condition are equal
# min_stat: minimum number of condition onsets per stationary segment
# inter_time: time in seconds at which the basis functions are measured, it is fixed at 0.05
# intercept: baseline of the HR signal for all subjects
# nROI: number of considered ROI per subject
# rho_ROI: correlation coefficient of the ROIs, assumed is an equi-correlation structure
## Output:
# BOLD: simulated BOLD signal
# true_cp_loc: true change point location per subject
# onsets: matrix containing the condition onset time series
# nscan: nscan
# true_dif_shape: list containing true difference between the shape parameters of the HRs to condition 1 and condition 2
# no_cp: number of change points, fixed at c(1,1)
simulate_BOLD_basis_func <- function(mean_change=c(0.5,0.5),
                                   var_change=c(0.1),
                                   basis_func=NULL,
                                   basis_beta= c(3.2, -6.4, 3.2),
                                   rho=0,
                                   snr,
                                   nscan, 
                                   nstim ,
                                   min_stat=10, 
                                   number_cp=c(1,1), N,
                                   inter_time=0.05,
                                   nROI, rho_ROI=0.2,
                                  intercept){
  # compute the covariance between the changes in the ROIs
  rho_ROI <- var_change*rho_ROI
  
  # check that mean_change has the same length as number of change points
  if(sum(number_cp)!= length(mean_change))stop("check that mean change has length of the number of change points")
  
  ## Initialize Stimulus time series: 
  # Random interstimulus intervals with mean such that all stimuli fit in there
  mean_itis <- floor(nscan/nstim)
  # use while loop to ensure that max onset is not larger than nscan.
  stim_done <- 0
  while(stim_done ==0){
    itis <- sample(c(mean_itis-1,mean_itis,mean_itis+1), size=nstim, replace=T) 
    os <- numeric(nstim)
    os[1] <- 1
    for(k in 2:nstim){
      os[k] <- os[k-1]+itis[k-1]
    }
    # check that max(os) is not larger than nscan
    if(max(os)<=nscan){
      stim_done <- 1
    }
  }
  
  # Randomize the order of the stimuli
  # Ensure that the same amount of stimuli are presented for both types of stimuli
  os1 <- numeric(nstim/2)
  os2 <- numeric(nstim/2)
  count <- 1
  ind1 <- 1
  ind2 <- 1
  while(count <=nstim){
    if(ind1 <= nstim/2 && ind2<=nstim/2){
      l <- rbinom(1,1,0.5) 
      if(l==1){
        os1[ind1] <- os[count]
        ind1 <- ind1+1
      }else{
        os2[ind2] <- os[count]
        ind2 <- ind2+1
      }
    }else{
      if(ind1<=nstim/2){
        os1[ind1] <- os[count]
        ind1 <- ind1+1
      }
      if(ind2<=nstim/2){
        os2[ind2] <- os[count]
        ind2 <- ind2+1
      }
    }
    count <- count +1
  }
  stim_1 <- stim_2 <- numeric(nscan)
  stim_1[os1] <- stim_2[os2] <- 1
  
  ## Initialize change points
  # Choose randomly one change point excluding the very beginning and very end (defined by min_stat)
  
  id_cp_2 <- sort(os2) # time of stimulus presentation
  
  id_cp_1 <- sort(os1) # time of stimulus presentation
  
  cp_loc <- matrix(c(sample(id_cp_1[(min_stat+1):(length(id_cp_1)-min_stat)],N, replace=T),
                     sample(id_cp_2[(min_stat+1):(length(id_cp_2)-min_stat)],N, replace=T)),
                   nrow = N, byrow=F)
  stim_1_1 <- stim_1_2 <-stim_2_1 <- stim_2_2 <-  matrix(0,ncol=nscan, nrow=N)
  for(i in 1:N){
    stim_1_1[i,1:(cp_loc[i,1]-1)] <- stim_1[1:(cp_loc[i,1]-1)]
    stim_1_2[i,cp_loc[i,1]:nscan] <- stim_1[cp_loc[i,1]:nscan]
    stim_2_1[i,1:(cp_loc[i,2]-1)] <- stim_2[1:(cp_loc[i,2]-1)]
    stim_2_2[i,cp_loc[i,2]:nscan] <- stim_2[cp_loc[i,2]:nscan]
  }
  
  true_cp <- cp_loc-1 #true cp returns the last observation of a stationary block.
  
  ## Initialize the true shape of the two HRFs
  
  # Intercept: not important for shape of HRF
   beta_intercept <- matrix(intercept, ncol=nROI, nrow=N)
  
  ## For changing HRF:
  # initialize the HRF before the change point:
  hrf_bc <- basis_beta[1]*basis_func[,1]+basis_beta[2]*basis_func[,2]+basis_beta[3]*basis_func[,3]
  
  # multiply by the value such that beta_after[1]= effect_size+beta_before.
  # this way I can ensure that the effect size follows a normal distribution (with very low variance though!)
  sig <- matrix(rho_ROI, nrow=nROI, ncol=nROI)
  diag(sig) <- var_change
  changes1 <- mvrnorm(N, rep(mean_change[1], nROI), sig)
  changes2 <- mvrnorm(N, rep(mean_change[2], nROI), sig)
  # these are two matrices with ncol=nROI and nrow=N
  
  m_ratio1 <- (basis_beta[1]+changes1)/basis_beta[1]
  m_ratio2 <- (basis_beta[1]+changes2)/basis_beta[1]
  
 
  ## compute the true HRF
  hrf_mat <- design_mat <- lapply(rep("list", nROI), vector, length=N)
  onsets <- matrix(c(stim_1, stim_2), ncol=2, byrow=F)
  Y<- Signals_all <- lapply(rep(0, nROI), matrix, ncol=nscan, nrow=N) 
  # Y is a list of length nROI, containing a matrix with nrow=N and ncol=nscan (BOLD signal for N participants per ROI)
  basis_hrf <- basis_beta[1]*basis_func[,1]+basis_beta[2]*basis_func[,2]+basis_beta[3]*basis_func[,3]
  # number of basis function needed for that:
  nbf <- ncol(basis_func)
  for(i in 1:N){
    # initialize onset time series for each subject. This is the same across all ROIs.
    os <- matrix(0, ncol=4, nrow=nscan)
    os[,c(1,3)] <- onsets
    os[,2] <- onsets[,1]
    os[,4] <- onsets[,2]
    os[1:true_cp[i,1],2] <- os[(true_cp[i,1]+1):nscan,1] <- os[1:true_cp[i,2],4] <- os[(true_cp[i,2]+1):nscan,3] <-0
    
    # add design matrix, this is also the same across all ROIs
    hrf_matrix <- matrix(nrow= nscan, ncol= 4*nbf)
    # only take values of the HR every TR seconds:
    time_sample <- 2/inter_time
    
    for(l in 1:4){
      for(p in 1:nbf){
        hrf_matrix[,(l-1)*nbf+p] <- convolve(os[,l], rev(basis_func[seq(1,559,by=time_sample),p]), type="open")[1:nscan]
      }
    }
    for(j in 1:nROI){
      # true HR for subject i in ROI j
      hrf_helper <- matrix(0, ncol=3, nrow= nrow(basis_func))
      hrf_helper[,1] <- basis_hrf
      hrf_helper[,2] <- hrf_helper[,1]*m_ratio1[i,j]
      hrf_helper[,3] <- hrf_helper[,1]*m_ratio2[i,j]
      hrf_mat[[j]][[i]] <- hrf_helper
      
      # get the subject and ROI specific model parameters:
      beta_all <- c(beta_intercept,basis_beta, m_ratio1[i,j]*basis_beta,basis_beta, m_ratio2[i,j]*basis_beta)
      # compute the subject and ROI specific signal
      hrf_matrix <- cbind(rep(1,nscan) , hrf_matrix)
      Signals_all[[j]][i,] <- signal <- hrf_matrix%*%beta_all
      
      ## add noise 
      if(rho==0){
        noise_sd <- abs(mean(signal)/snr)
        Y[[j]][i,] <- signal+rnorm(nscan, mean=0, sd=noise_sd)
      }else{
        var_epsilon <- (mean(signal)/snr)^2
        sd_ar <- sqrt(var_epsilon*(1-rho^2))
        Y[[j]][i,] <- signal+arima.sim(model=list(ar=rho),n=nscan, sd=sd_ar)
      }
    }
  }
  ## compute the true shapes and shape parameter:
  shape_per_hrf <- matrix(ncol=10, nrow=3)
  colnames(shape_per_hrf) <- c("PM", "NA", "IUA", "TTP", "TPN", "FWHM", "FWHN", "AUC", "AC_P", "DC_P")
  dif_1 <- dif_2 <- lapply(rep(0, nROI), matrix, ncol=10, nrow=N)
  
  for(j in 1:nROI){
    colnames(dif_1[[j]]) <- colnames(dif_2[[j]]) <- c("PM", "NA", "IUA", "TTP", "TPN", "FWHM", "FWHN", "AUC", "AC_P", "DC_P")
    for(i in 1:N){
      for(q in 1:3){
        shape_per_hrf[q,] <- shape_parameters2(hrf_mat[[j]][[i]][,q], baseline=0, tres=0.1, max.ttp=12)
      }
      dif_1[[j]][i,] <- shape_per_hrf[1,]-shape_per_hrf[2,] 
      dif_2[[j]][i,] <- shape_per_hrf[1,]-shape_per_hrf[3,]
    }
  }
 
  output <- list(BOLD=Y,
                 true_cp_loc = true_cp,
                 onsets = onsets,
                 nscan=nscan,
                 true_difference_shape = list("condition_1"= dif_1,
                                              "condition_2"= dif_2),
                 no_cp = number_cp)
  return(output)
}
### Estimate regression coefficients and their variance for one subject and one ROI
## Input:
# BOLD_signal: observed BOLD signal for one subject
# true_cp_location: true location of the change point(s) per condition
# onsets: condition onset time series
# nscan: number of observations
# number_cp: number of change points per condition, should be a vector of length 2 (b.c. we consider two conditions)
# hrf_basis_fun: matrix with columns corresponding to the basis functions of the HRF
# pre_whiten: Do we need to pre-whiten the data
# misspec: larges absolute difference between true change point location and change point location used in analysis
## Output:
# full_model: output of the GLM function used to estimate the model parameters
# est_beta: estimated regression coefficients corresponding to the (non-stationary) HRs
# var_within: within subject variance of est_beta
## Remark:
# This function includes initialization of the design matrix used to analyse the observed BOLD signal
est_betas_var <- function(BOLD_signal, true_cp_location, onsets, nscan, number_cp, hrf_basis_fun,  full_covar=F, pre_whiten=F, misspec=0){
  # initialize the design matrix
  nbf <- ncol(hrf_basis_fun)
  mis <- sample(c(-misspec:misspec), size=sum(number_cp!=0) ,replace=T)
  id.onset_1 <- which(onsets[,1]==1)
  id.onset_2 <- which(onsets[,2]==1)
  if(number_cp[2]==0){
    os <- matrix(0, ncol=3, nrow=nscan)
    os[,c(1,3)] <- onsets
    os[,2] <- onsets[,1]
    id.true_cp <- which.min(abs(id.onset_1-true_cp_location))
    if(id.true_cp+mis >= 50 | id.true_cp+mis <=10){
      if(sign(mis)==1){
        id.new.cp <- id.onset_1[49]
      }else{
        id.new.cp <- id.onset_1[10]
      }
    }else{
      id.new.cp <- id.onset_1[id.true_cp+mis]
    }
    os[1:(id.new.cp),2] <- os[(id.new.cp+1):nscan,1] <- 0
    
    hrf_matrix <- matrix(nrow= nscan, ncol= 3*nbf)
    
    for(j in 1:3){
      for(p in 1:nbf){
        hrf_matrix[,(j-1)*nbf+p] <- convolve(os[,j], rev(hrf_basis_fun[seq(1,559,by=40),p]), type="open")[1:nscan]
      }
    }
    cp_loc <- c(id.new.cp, 0)
    id_nonstat <- c(1:(2*nbf))
  }else if(number_cp[1]==0){
    # only change in the second condition
    os <- matrix(0, ncol=3, nrow=nscan)
    os[,c(1,2)] <- onsets
    os[,3] <- onsets[,2]
    id.true_cp <- which.min(abs(id.onset_2-true_cp_location))
    if(id.true_cp+mis >= 50 | id.true_cp+mis <=10){
      if(sign(mis)==1){
        id.new.cp <- id.onset_2[49]
      }else{
        id.new.cp <- id.onset_2[10]
      }
    }else{
      id.new.cp <- id.onset_2[id.true_cp+mis]
    }
    os[1:(id.new.cp),3] <- os[(id.new.cp+1):nscan,2] <- 0
    
    hrf_matrix <- matrix(nrow= nscan, ncol= 3*nbf)
    
    for(j in 1:3){
      for(p in 1:nbf){
        hrf_matrix[,(j-1)*nbf+p] <- convolve(os[,j], rev(hrf_basis_fun[seq(1,559,by=40),p]), type="open")[1:nscan]
      }
    }
    cp_loc <- c(0, id.new.cp)
    id_nonstat <- c((nbf+1):(3*nbf))
  }else{
    os <- matrix(0, ncol=4, nrow=nscan)
    os[,c(1,3)] <- onsets
    os[,2] <- onsets[,1]
    os[,4] <- onsets[,2]
    
    id.true_cp_1 <- which.min(abs(id.onset_1-true_cp_location[1]))
    id.true_cp_2 <- which.min(abs(id.onset_2-true_cp_location[2]))
    id.new.cp <- numeric(2)
    if(id.true_cp_1+mis[1] >= 50 | id.true_cp_1+mis[1] <=10){
      if(sign(mis[1])==1){
        id.new.cp[1] <- id.onset_1[49]
      }else{
        id.new.cp[1] <- id.onset_1[10]
      }
    }else{
      id.new.cp[1] <- id.onset_1[id.true_cp_1+mis[1]]
    }
    if(id.true_cp_2+mis[2] >= 50 | id.true_cp_2+mis[2] <=10){
      if(sign(mis[2])==2){
        id.new.cp[2] <- id.onset_2[49]
      }else{
        id.new.cp[2] <- id.onset_2[10]
      }
    }else{
      id.new.cp[2] <- id.onset_2[id.true_cp_2+mis[2]]
    }
    os[1:id.new.cp[1],2] <- os[(id.new.cp[1]+1):nscan,1] <- os[1:id.new.cp[2],4] <- os[(id.new.cp[2]+1):nscan,3] <-0
    
    hrf_matrix <- matrix(nrow= nscan, ncol= 4*nbf)
    
    for(j in 1:4){
      for(p in 1:nbf){
        hrf_matrix[,(j-1)*nbf+p] <- convolve(os[,j], rev(hrf_basis_fun[seq(1,559,by=40),p]), type="open")[1:nscan]
      }
    }
    cp_loc <- id.new.cp
    id_nonstat <- c(1:(4*nbf))
  }
  Xdesign <- cbind(hrf_matrix, rep(1, nscan))
  
  if(pre_whiten){
    mod_1 <- glm_prewhitened(X=Xdesign, Y=BOLD_signal, whitenend_data=T)
  }else{
    mod_1 <- glm.fit(x=Xdesign, y=BOLD_signal,intercept=F)
    mod_1$x <- Xdesign
  }
  
  # See what difference it makes if I account for length of the stationary segments when estimating within subject variance.
  var_betas <- var_est_subject(glm_model=mod_1, cp_location=cp_loc, nbf=nbf, number_cond=2, full_covar =full_covar)

  
  out <- list(full_model = mod_1,
              est_beta= mod_1$coefficients[id_nonstat],
              var_within = var_betas)
  return(out)
}

### Whole set up: simulate BOLD signal, subject level analysis, group level analysis:
## Input:
# effect: mean_change in simulate_BOLD 
# N: number of subjects
# snr1: Signal-to-noise ratio on subject level, used to compute the within-subject variance
# seed: fix a seed before generating BOLD signals
# basis_function: matrix in which the columns represent the basis functions used to generate BOLD and analyse BOLD 
# var_change: between subject variance 
# beta_basis: regression coefficients corresponding to the basis_function, used to define the shape of the HRs before the change points 
# misspec: larges absolute difference between true change point location and change point location used in analysis
# nROI: number of considered ROI per subject
# rho_ROI: correlation coefficient between the considered ROIs.
## Output:
# p_values_KH: p_values corresponding to the null hypotheses of no change in the different conditions and shape parameters; based on T_KH
# p_values_Wald:  p_values corresponding to the null hypotheses of no change in the different conditions and shape parameters; based on T_WALD
# est_eta: estimated difference of the shape parameter on group level
# true_difference_shape: true difference of the shape parameters per subject and ROI
# est_difference_shape: estimated difference between the shape parameters per subject and ROI
# between_subject_variance: between subject variance
# variance_group_level: variance on group level
# within_subject_variance: within subject variance
## Remark:
# The considered shape parameters are: PM, NA, TTP, TPN, FWHM, FWHN, AUC
# For the simulation study this function is run repeatedly. 
# The set of rejected null hypotheses is defined through the application of hierarchical testing methods on the p-values 
# The FLOBS basis functions have been used.
# To make computation faster, we have set rho=0, so that pre-whitening was not necessary.
sim_known_cp <- function(effect, N, snr1, seed, basis_function, var_change=c(0.1), beta_basis =c(3.2, -6.4, 3.2),misspec = 0,rho1, nROI, rho_ROI){
  set.seed(seed)
  sim_BOLD <- simulate_BOLD_basis_func(mean_change = effect,
                                       nscan=500,
                                       nstim=120,
                                       min_stat=15,
                                       snr=snr1,
                                       rho=rho1,
                                       number_cp = c(1,1),
                                       N=N,
                                       basis_func = basis_function,
                                       var_change=var_change,
                                       basis_beta= beta_basis,
                                       nROI=nROI, rho_ROI= rho_ROI
  )
  output_est_sp_var<- lapply(rep("list", nROI), vector, length=N)
  if(rho1==0){
      pre_whit <- F
    }else{
      pre_whit <- T
    }
  for(j in 1:nROI){
    for(i in 1:N){
      betas_var<- est_betas_var(BOLD_signal=sim_BOLD$BOLD[[j]][i,], true_cp_location=sim_BOLD$true_cp_loc[i,],
                                                           onsets= sim_BOLD$onsets, nscan=sim_BOLD$nscan, number_cp = sim_BOLD$no_cp,
                                                           hrf_basis_fun=hrf_basis_fun, misspec= misspec, prewhiten=pre_whit)
      output_est_sp_var[[j]][[i]] <- est_var_sp_subject(betas= betas_var$est_beta, covar_betas=betas_var$var_within$hat_covar, hrf=hrf_basis_fun, mc_B=10000, fix_seed=seed+14)
      # this returns the estimator per shape parameter and the covariance matrix per shape parameter
  }
  }
  rm(betas_var)
  between_subject_variance <-p_vals_1 <- p_vals_2 <-est_eta<- lapply(rep(0,nROI), matrix, ncol=7, nrow=2)
  var_group_level <-est_dif_shape_para<- lapply(rep("list", nROI), vector, length=7)
 
  for(j in 1:nROI){
    colnames(between_subject_variance[[j]]) <-colnames(p_vals_1[[j]])<- colnames(p_vals_2[[j]])<-colnames(est_eta[[j]])<- names(var_group_level[[j]]) <- names(est_dif_shape_para[[j]])<- c("PM", "NA", "IUA", "TTP", "TPN", "FWHM", "FWHN", "AUC", "AC_P", "DC_P")[c(1,2,4:8)]
     indice <- 1
    for(sp in c(1,2,4:8)){
      sp_sub_1 <- sp_sub_2 <- matrix(0, ncol=2, nrow=N)
      sp_var_sub_1 <- sp_var_sub_2 <- vector("list", N)
      for(i in 1:N){
        sp_sub_1[i,1] <- output_est_sp_var[[i]]$point_estimators[[sp]][1]
        sp_sub_1[i,2] <- output_est_sp_var[[i]]$point_estimators[[sp]][2]
        sp_var_sub_1[[i]] <- output_est_sp_var[[i]]$variances[[sp]][1:2,1:2]
      
        sp_sub_2[i,1] <- output_est_sp_var[[i]]$point_estimators[[sp]][3]
        sp_sub_2[i,2] <- output_est_sp_var[[i]]$point_estimators[[sp]][4]
        sp_var_sub_2[[i]] <- output_est_sp_var[[i]]$variances[[sp]][3:4,3:4]
      }
    
      if(REML){
        group_estimation_1 <- reml_var_estimation_2(betas_sub=sp_sub_1, var_sub=sp_var_sub_1, contrast=c(-1,1),N=N,maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F,full_covar=F)
        group_estimation_2 <- reml_var_estimation_2(betas_sub=sp_sub_2, var_sub=sp_var_sub_2, contrast=c(-1,1),N=N,maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F,full_covar=F)
      }else{
        group_estimation_1 <- em_estimation(betas_sub=sp_sub_1, var_sub=sp_var_sub_1, contrast=c(-1,1),N=N,maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F,full_covar=F)
        group_estimation_2 <- em_estimation(betas_sub=sp_sub_2, var_sub=sp_var_sub_2, contrast=c(-1,1),N=N,maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F,full_covar=F)
    }
    
      est_eta[[j]][1,indice] <- group_estimation_1$beta_group
      est_eta[[j]][2, indice] <- group_estimation_2$beta_group
      between_subject_variance[[j]][1, indice] <- group_estimation_1$sigma_between
      between_subject_variance[[j]][2, indice]<- group_estimation_2$sigma_between
      var_group_level[[j]][[indice]] <- list(condition1= group_estimation_1$var_mat,
                                             condition2= group_estimation_2$var_mat)
      est_dif_shape_para[[j]][[indice]] <- list(condition1= group_estimation_1$betas_subject,
                                                condition2= group_estimation_2$betas_subject)
    
      # compute p-values per shape parameter
    
      test_KH_1 <- fmri_group_KH(est_eta=as.vector(group_estimation_1$beta_group), 
                                 var_mat_full=group_estimation_1$var_mat, 
                                 betas_c=as.vector(group_estimation_1$betas_subject), N=N)
      test_KH_2 <- fmri_group_KH(est_eta=as.vector(group_estimation_2$beta_group), 
                                 var_mat_full=group_estimation_2$var_mat, 
                                 betas_c=as.vector(group_estimation_2$betas_subject), N=N)
    
      p_vals_1[[j]][1, indice] <- test_KH_1$p_value
      p_vals_1[[j]][2, indice] <- test_KH_2$p_value
    
      p_vals_2[[j]][1, indice] <- fmri_test_own(est_eta=as.vector(group_estimation_1$beta_group),
                                           var_mat_full=group_estimation_1$var_mat,
                                           N=N, null_eta=0)
      p_vals_2[[j]][2, indice] <- fmri_test_own(est_eta=as.vector(group_estimation_2$beta_group),
                                           var_mat_full=group_estimation_2$var_mat,
                                           N=N, null_eta=0)
    
      indice <- indice+1
    }
  }
  out <- list(p_values_KH = p_vals_1,
              p_values_Wald= p_vals_2,
              est_eta = est_eta,
              true_difference_shapes= sim_BOLD$true_difference_shape,
              est_difference_shapes = est_dif_shape_para,
              between_subject_variance= between_subject_variance,
              variance_group_level = var_group_level,
              within_subject_variances= output_est_sp_var
  )
  return(out)
}

### Analysis of simulation results
## Input:
# sim_results: output of sim_known_cp
# snr: the SNR corresponding to dat_considered
# effect: the effect corresponding to dat_considered
# nROI: number of considered ROI
# alpha: level at which the FWER is supposed to be controlled
# N: number of subjects
## Output:
# summary_eta: dataframe containing summary statistics such as the empirical average and standard deviation of the group level estiamtes
# summary_rej:list containing dataframes with the empirical rejection rates for each level of the hierarchy of the hypotheses
# full_data: Full results of the estimated eta per subject and ROI, as well as the rejections at each level of the 
analyze_known_cp <- function(sim_results, effect, snr, nROI=5, alpha=0.05,N){
  ### start with the proposed procedure:
  # report the estimated eta and it's empirical standard deviation
  B <- length(sim_results)
  est_eta_1_prop <- est_eta_2_prop <-matrix(ncol=7, nrow=B*nROI)
  names_shape <- c("PM", "NA", "TTP", "TPN", "FWHM", "FWHN", "AUC")
  for(b in 1:B){
    for(j in 1:nROI){
      est_eta_1_prop[((b-1)*nROI+j),] <- sim_results[[b]]$proposed$est_eta[[j]][1,] 
      est_eta_2_prop[((b-1)*nROI+j),] <- sim_results[[b]]$proposed$est_eta[[j]][2,]
    }
  }
  
  # Compute the average and empirical standard deviation:
  av_est_eta <- data.frame(mean= c(colMeans(est_eta_1_prop, na.rm=T), colMeans(est_eta_2_prop, na.rm=T)),
                           sd= c(apply(est_eta_1_prop, 2, sd, na.rm=T), apply(est_eta_2_prop, 2, sd, na.rm=T)),
                           shape_para= rep(names_shape, 2),
                           effect = rep(effect, each=7),
                           snr= rep(snr, 14))
  est_eta_1_prop <- as.data.frame(est_eta_1_prop)
  est_eta_2_prop <- as.data.frame(est_eta_2_prop)
  
  est_eta_1_prop$effect <- effect[1]
  est_eta_2_prop$effect <- effect[2]
  est_eta_prop <- rbind(est_eta_1_prop, est_eta_2_prop)
  colnames(est_eta_prop[,1:7]) <- names_shape
  est_eta_prop$snr <- snr
  rm(est_eta_1_prop)
  rm(est_eta_2_prop)
  
  # apply inheritance procedure to get rejection rates:
  # do one or three data.frames (one for each level)
  # have the different hypotheses as rows and indicate by zero-one if the hypothesis is rejected
  rej_ROI_KH <- rej_ROI_W <- data.frame(ROI= c(1:nROI),
                                        matrix(0,ncol = B, nrow=nROI))
  rej_cond_KH <-rej_cond_W <- data.frame(ROI= rep(c(1:nROI), each=2),
                                         effect= rep(effect, nROI),
                                         matrix(0, ncol=B, nrow=nROI*2))
  rej_shape_KH <- rej_shape_W <- data.frame(ROI= rep(c(1:nROI), each=14),
                                            effect= rep(rep(effect, each=7), nROI),
                                            sp = rep(names_shape, nROI*2),
                                            matrix(0, ncol=B, nrow=14*nROI))
  
  for(b in 1:B){
    # get the p-values for that, need to have the format:
    # ROI_1_cond_1_PM, ROI_1_cond_1_NA,...,ROI_1_cond_2_PM,... ROI_5_cond_2_AUC
    p_vals_KH <-p_vals_W <- numeric(0)
    for(j in 1:nROI){
      helper1 <- as.vector(t(sim_results[[b]]$proposed$p_values_KH[[j]]))
      helper2 <- as.vector(t(sim_results[[b]]$proposed$p_values_Wald[[j]]))
      p_vals_KH <- c(p_vals_KH, helper1)
      p_vals_W <- c(p_vals_W, helper2)
    }
    rejection_KH <- inheritance_selection(family_size=c(nROI,2,7), test=c("arbitrary"), alpha=alpha, p_vals_KH)$rejections
    rejection_W <- inheritance_selection(family_size=c(nROI,2,7), test=c("arbitrary"), alpha=alpha, p_vals_W)$rejections
    # this returns a list with rejections at each level of the hypothesis. Now safe it in the data.frames
    rej_ROI_KH[rejection_KH[[1]],b+1] <- 1
    rej_ROI_W[rejection_W[[1]],b+1] <- 1
    
    rej_cond_KH[rejection_KH[[2]], b+2] <- 1
    rej_cond_W[rejection_W[[2]], b+2] <- 1
    
    rej_shape_KH[rejection_KH[[3]],b+3] <- 1
    rej_shape_W[rejection_W[[3]],b+3] <- 1
  }
  ## report the empirical rejection rate, that is, the mean of the respective rows, without the first few rows:
  av_rej_ROI <- data.frame(ROI=c(1:nROI),
                           effect1= rep(effect[1], nROI),
                           effect2= rep(effect[2], nROI),
                           rej_rate_KH= rowMeans(rej_ROI_KH[,-1]),
                           rej_rate_W= rowMeans(rej_ROI_W[,-1]))
  av_rej_cond <- data.frame(ROI=rep(c(1:nROI),each=2),
                            effect= rep(effect,nROI),
                            rej_rate_KH= rowMeans(rej_cond_KH[,-c(1:2)]),
                            rej_rate_W= rowMeans(rej_cond_W[,-c(1:2)]))
  av_rej_shape <- data.frame(ROI= rep(c(1:nROI), each=14),
                             effect= rep(rep(effect, each=7), nROI),
                             sp = rep(names_shape, nROI*2),
                             rej_rate_KH= rowMeans(rej_shape_KH[,-c(1:3)]),
                             rej_rate_W= rowMeans(rej_shape_W[,-c(1:3)]))
  
  ## report the empirical FWER (so the (epmirical) probability of making at least one rejection)
  # add to rejection dataframes a column indicating if null hypothesis is true or false
  # do the same for the average rejection rate (since this coincides then to the FWER since all rejection rates are based on the same number of trials)
  
  rej_ROI_KH$H0 <- rej_ROI_W$H0 <- av_rej_ROI$H0 <- FALSE # since in our simulations, the ROI hypothesis is always true
  
  
  if(!0%in%effect){
    # if the effect size is not equal to zero, the condition hypotheses are false as well
    rej_cond_KH$H0 <- rej_cond_W$H0 <- av_rej_cond$H0 <- FALSE 
    
    # if the effect size is not equal to zero, the PM, NA, AUC are false, all others are true
    rej_shape_KH$H0 <- rej_shape_W$H0 <- av_rej_shape$H0 <-TRUE
    rej_shape_KH$H0[rej_shape_KH$sp%in% c("PM", "NA", "AUC")] <-rej_shape_W$H0[rej_shape_W$sp%in% c("PM", "NA", "AUC")] <- av_rej_shape$H0[av_rej_shape$sp%in%c("PM", "NA", "AUC")] <-FALSE
    
  }else{
    # when effect is zero, the null hypotheses are true for the respective condition
    rej_cond_KH$H0[rej_cond_KH$effect==0] <- rej_cond_W$H0[rej_cond_W$effect==0] <-av_rej_cond$H0[av_rej_cond$effect==0]<- TRUE #if condition 1 has effect 0: the condition hypothesis is true
    rej_cond_KH$H0[rej_cond_KH$effect!=0] <- rej_cond_W$H0[rej_cond_W$effect!=0] <-av_rej_cond$H0[av_rej_cond$effect!=0] <-FALSE
    
    # then, for the shape parameters, all hypotheses are true except for PM, NA, AUC of condition 2
    rej_shape_KH$H0 <- rej_shape_W$H0 <- av_rej_shape$H0 <- TRUE
    rej_shape_KH$H0[rej_shape_KH$sp%in% c("PM", "NA", "AUC") & rej_shape_KH$effect!=0] <-rej_shape_W$H0[rej_shape_W$sp%in% c("PM", "NA", "AUC")& rej_shape_W$effect!=0]<-av_rej_shape$H0[av_rej_shape$sp%in% c("PM", "NA", "AUC")& av_rej_shape$effect!=0] <- FALSE
  }
  rej_shape_KH$snr <- rej_shape_W$snr <- rej_cond_KH$snr <- rej_cond_W$snr <-rej_ROI_KH$snr <- rej_ROI_W$snr <-snr
  av_rej_shape$snr <- av_rej_cond$snr <- av_rej_ROI$snr <- snr
  
  ## lastly: report also the test variance and its standard deviation (utilize the group level variance I believe...)
  test_var_1_prop_KH <- test_var_1_prop_W <- test_var_2_prop_KH <- test_var_2_prop_W <- matrix(ncol=7, nrow=B*nROI)
  #true_var_1_KH <- true_var_1_W <- true_var_2_KH <- true_var_2_W <- matrix(ncol=7, nrow=B*nROI)
  for(b in 1:B){
    for(j in 1:nROI){
      for(s in 1:7){
        test_var_1_prop_KH[((b-1)*nROI+j),s] <- group_sd_KH(sim_results[[b]]$proposed$variance_group_level[[j]][[s]][[1]], 
                                                            betas_c= sim_results[[b]]$proposed$est_dif_shape_param[[j]][[s]][[1]], N) 
        test_var_2_prop_KH[((b-1)*nROI+j),s] <- group_sd_KH(sim_results[[b]]$proposed$variance_group_level[[j]][[s]][[2]], 
                                                            betas_c= sim_results[[b]]$proposed$est_dif_shape_param[[j]][[s]][[2]], N) 
        
        test_var_1_prop_W[((b-1)*nROI+j),s] <- group_sd_own(sim_results[[b]]$proposed$variance_group_level[[j]][[s]][[1]],N)
        test_var_2_prop_W[((b-1)*nROI+j),s] <- group_sd_own(sim_results[[b]]$proposed$variance_group_level[[j]][[s]][[2]],N)
      }
    }
  }
  
  # Compute the average and empirical standard deviation:
  av_test_sd <- data.frame(mean= c(colMeans(test_var_1_prop_KH, na.rm=T), 
                                   colMeans(test_var_2_prop_KH, na.rm=T),
                                   colMeans(test_var_1_prop_W, na.rm=T), 
                                   colMeans(test_var_2_prop_W, na.rm=T)),
                           sd= c(apply(test_var_1_prop_KH, 2, sd, na.rm=T), 
                                 apply(test_var_2_prop_KH, 2, sd, na.rm=T),
                                 apply(test_var_1_prop_W, 2, sd, na.rm=T), 
                                 apply(test_var_2_prop_W, 2, sd, na.rm=T)),
                           shape_para= rep(names_shape, 4),
                           effect = rep(rep(effect, each=7),2),
                           snr= rep(snr, 14*2),
                           method= rep(c("KH", "Wald"), each=14))
  test_var_1_prop_KH <- as.data.frame(test_var_1_prop_KH)
  test_var_2_prop_KH <- as.data.frame(test_var_2_prop_KH)
  test_var_1_prop_W <- as.data.frame(test_var_1_prop_W)
  test_var_2_prop_W <- as.data.frame(test_var_2_prop_W)
  
  test_var_1_prop_KH$effect <-test_var_1_prop_W$effect <- effect[1]
  test_var_2_prop_KH$effect <- test_var_2_prop_W$effect <-effect[2]
  test_var_1_prop_KH$method <- test_var_2_prop_KH$method <- "KH"
  test_var_1_prop_W$method <- test_var_2_prop_W$method <- "Wald"
  
  test_sd_prop <- rbind(test_var_1_prop_KH, test_var_2_prop_KH, test_var_1_prop_W, test_var_2_prop_W)
  colnames(test_sd_prop[,1:7]) <- names_shape
  test_sd_prop$snr <- snr
  rm(test_var_1_prop_KH)
  rm(test_var_2_prop_KH)
  rm(test_var_1_prop_W)
  rm(test_var_2_prop_W)
  
  out_proposed <- list(summary_eta= av_est_eta,
                       summary_rej = list(ROI= av_rej_ROI,
                                          condition= av_rej_cond,
                                          shape= av_rej_shape),
                       full_data = list(full_eta= est_eta_prop,
                                        full_test_sd = test_sd_prop,
                                        full_rejection_KH= list(ROI= rej_ROI_KH,
                                                                condition= rej_cond_KH,
                                                                shape= rej_shape_KH),
                                        full_rejection_W=list(ROI= rej_ROI_W,
                                                              condition= rej_cond_W,
                                                              shape= rej_shape_W)))
  return(out_proposed)
}
### Functions to compute the bias and MSE of the estimated differences at the group level:
## Input:
# basis_func: basis functions of the HRF model
# basis_beta: basis beta used to simulate the BOLD signal
# effect: mean change of the simulated BOLD signal
## Output:
# true difference in the considered shape parameters at the group level
true_changes <- function(basis_func=hrf_basis_fun, basis_beta=c(3.2, -6.4, 3.2), effect){
  hrf_bc <- basis_beta[1]*basis_func[,1]+basis_beta[2]*basis_func[,2]+basis_beta[3]*basis_func[,3]
  m_ratio <- (basis_beta[1]+effect)/basis_beta[1]
  hrf_ac <- hrf_bc*m_ratio
  shape_bc <- shape_parameters2(hrf_bc, baseline=0, tres=0.1, max.ttp=12)[c(1,2,4:8)]
  shape_ac <- shape_parameters2(hrf_ac, baseline=0, tres=0.1, max.ttp=12)[c(1,2,4:8)]
  return(shape_ac-shape_bc)
}
## Input
# summary_eta_prop: summary_eta; output of the analyze_known_cp
# effect: considered effect size
# mis: number indicating if the true change point location is misspecified
analyze_eta <- function(summary_eta_prop, effect, mis){
  summary_eta_prop$cp <- TRUE
  true_effect <- c(true_changes(effect=effect[1]), true_changes(effect=effect[2]))
  summary_eta_prop$bias <- summary_eta_prop$mean-true_effect
  summary_eta_prop$mse <- (summary_eta_prop$bias)^2+summary_eta_prop$sd^2
  summary_eta_prop$procedure <- "proposed"
  return(rbind(summary_eta_prop))
}

### Compute the average rejection rates at the shape parameter level
## Input
# effect: considered effect size
# snr: considered SNR
# summary_rej: summary_rej from output of the analyze_known_cp
# mis: number indicating if the true change point location is misspecified
## Output
# dataframe containing the average rejection rates per shape parameter, test statistic and effect size
average_rej_prop <- function(effect, snr, summary_rej, mis=0){
  rej_shape <- data.frame(test=rep(c("KH", "Wald"), each=7*2),
                          effect= rep(rep(effect,each=7),2),
                          shape= rep(c("PM", "NA","TTP", "TPN", "FWHM", "FWHN", "AUC"),4),
                          rej_rate= c(mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="PM"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="NA"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="TTP"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="TPN"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="FWHM"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="FWHN"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="AUC"&summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="PM"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="NA"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="TTP"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="TPN"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="FWHM"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="FWHN"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_KH[summary_rej$shape$sp=="AUC"&summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="PM"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="NA"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="TTP"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="TPN"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="FWHM"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="FWHN"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="AUC"& summary_rej$shape$effect==effect[1]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="PM"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="NA"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="TTP"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="TPN"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="FWHM"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="FWHN"& summary_rej$shape$effect==effect[2]]),
                                      mean(summary_rej$shape$rej_rate_W[summary_rej$shape$sp=="AUC"& summary_rej$shape$effect==effect[2]])),
                          snr= rep(snr, 7*4))
  rej_ROI$mis <- rej_cond$mis <- rej_shape$mis <- mis
  rej_ROI$procedure <- rej_cond$procedure <- rej_shape$procedure <- "proposed"
  return(rej_shpae)
}

### Compute the empricial FWER
## Input
# summary_rej: summary_rej from output of the analyze_known_cp
# effect: considered effect size
# snr: considered SNR
# mis: number indicating if the true change point location is misspecified
## Output
# dataframe containing the empirical FWER per test statistic and considered effect sizes/ scenario
empirical_FWER <- function(summary_rej, effect, snr,mis=0){
  out <- data.frame(test= c("KH", "Wald"),
                         effect1 = rep(effect[1],2),
                         effect2= rep(effect[2],2),
                         snr = rep(snr,2),
                         level= rep("shape",2),
                         FWER= round(c(mean(summary_rej_prop$shape$rej_rate_KH[summary_rej_prop$shape$H0]),
                                 mean(summary_rej_prop$shape$rej_rate_W[summary_rej_prop$shape$H0])),4))
  out$procedure <- "proposed"
  out$mis <- mis
  return(out)
}


####
# Unknown change point location
####

### canoncial HRF as a function of t.
## Input:
# t: Time
# a1, a2, b1, b2, d: parameters defining the shape of the canoncial HRF
## Output:
# canonical HRF at time t.
hrf_canonical <- function(t, a1=6, a2=12, b1=0.9, b2=0.9, d=0.35){
  out <- t^(a1-1)*b1^(a1)*exp(-b1*t)/gamma(a1)-d*t^(a2-1)*b2^(a2)*exp(-b2*t)/gamma(a2)
  return(out)
}

### Simulate BOLD signal based on canonical HRF; one change point; one condition
## Input:
# mean_change: average change in the PM, corresponds to eta on group level
# var_change: variance of the change, corresponds to between subject variance
# nscan, nstim: number of scans and stimulus onsets, respectively
# rho_noise: parameter of the AR(1) process
# min_stat: minimum length of the stationary segment
# snr: signal to noise ratio
# N: number of subjects per group
## Output:
# BOLD: simulated BOLD signal for N subjects
# true_cp_loc: true location of the change points per subject
# true_changes: true effect sizes per subject
# onsets: condition onset time series per subject
# nscan, nstim: number of scans and stimulus onsets, respectively
# SNR: Signal to noise ratio
# rho: parameter of the AR(1) process
# sigma2_w: Variance of the AR(1) process, corresponds to sigma^2_{epsilon_i}
## Remarks:
# Note that the baseline of the BOLD signal per subjects corresponds to rnorm(N, 10, 1)
sim_BOLD_canonical_one_ROI <- function(mean_change=c(0.5),
                                             var_change=c(0.1),
                                             nscan=250, nstim=60,rho_noise,
                                             min_stat=10, snr, N){
  ## Initialize Stimulus time series: 
  # Random interstimulus intervals with mean defined by nscan and nstim:
  mean_itis <- floor(nscan/nstim)
  stim_done <- 0
  while(stim_done ==0){
    itis <- sample(c(mean_itis-1,mean_itis,mean_itis+1), size=nstim, replace=T) 
    os <- numeric(nstim)
    os[1] <- 1
    for(k in 2:nstim){
      os[k] <- os[k-1]+itis[k-1]
    }
    # check that max(os) is not larger than nscan
    if(max(os)<=nscan){
      stim_done <- 1
    }
  }
  
  stim_1 <-numeric(nscan)
  stim_1[os] <- 1
  
  ## Initialize change points
  # Choose randomly one change point excluding the very beginning and very end (defined by min_stat)
  id_cp_1 <- sort(os) # time of stimulus presentation
  id_cp_len_1 <- length(os)
  
  cp_loc <- sample(id_cp_1[(min_stat+1):(id_cp_len_1-min_stat)],N, replace=T)
  stim_1_1 <- stim_1_2 <- matrix(0,ncol=nscan, nrow=N)
  for(i in 1:N){
    stim_1_1[i,1:(cp_loc[i]-1)] <- stim_1[1:(cp_loc[i]-1)]
    stim_1_2[i,cp_loc[i]:nscan] <- stim_1[cp_loc[i]:nscan]
  }
  
  true_cp <- cp_loc 
  # true cp returns the first observation of a stationary block, should correspond to an onset.
  
  ## Initialize the true shape of the two HRFs
  
  # Intercept: not important for shape of HRF, assumed to not change. Add it anyways for estimation
  beta_intercept <- rnorm(N, 1, 1) 
  
  # For changing HRF:
  # use canonical HRF with changing effect size (influences PM, NA, does not influence TTP, don't know about FWHM)
  # get instantly all N participants
  # set effect size to 1 for both stimuli
  changes <- rnorm(N, mean_change,sqrt(var_change))
  
  hrf_mat <- matrix(0, nrow=nscan, ncol=2)
  for(i in 1:N){
    # onsets has to be given in seconds, currently given in TR, so multiply everything by 2:
    list_onsets <- list(first=which(stim_1_1[i,]!=0)*2,
                        second= which(stim_1_2[i,] !=0)*2)
    design_mat[[i]] <- simprepTemporal(onsets= list_onsets, duration=list(1,1), totaltime=nscan*2, TR=2, 
                                       effectsize=list(1,(1+changes[i])), hrf="double-gamma")
  }
  ## Define Y:
  # add noise (temporal noise and white noise):
  Y <- Signals_all <-  matrix(nrow=N, ncol=nscan)
  var_ar <- numeric(N)
  for(i in 1:N){
    # signal as the convolution of the condition onset time series with HR
    hrf <- hrf_canonical(seq(1,40, by=2)) # HR is less than 40 seconds long; we sample by TR=2
    hrf_mat[,1] <- convolve(stim_1_1[i,], rev(hrf), type="open")[1:nscan]
    hrf_mat[,1] <- convolve(stim_1_2[i,], rev(hrf), type="open")[1:nscan]
    beta_all <- c(1, 1+changes[i])
    Signals_all[i,] <- Signal <- hrf_mat%*%beta_all + beta_intercept[i]
    var_epsilon <- (mean(Signal)/snr)^2
    sd_ar <- var_ar[i] <- sqrt(var_epsilon*(1-rho_noise^2))
    Y[i,] <- Signal+arima.sim(model=list(ar=rho_noise),n=nscan, sd=sd_ar)
  }
  var_ar <- var_ar^2
  output <- list(BOLD=Y,
                 true_cp_loc = true_cp,
                 true_changes= changes,
                 onsets = stim_1,
                 nscan=nscan,
                 nstim= nstim, 
                 SNR=snr,
                 rho=rho_noise,
                 sigma2_w = var_ar)
  return(output)
}
### Simulate BOLD signal based on canonical HRF; one change point; one condition, multiple ROI
## Input:
# mean_change: average change in the PM, corresponds to eta on group level
# var_change: variance of the change, corresponds to between subject variance
# nscan, nstim: number of scans and stimulus onsets, respectively
# rho_noise: parameter of the AR(1) process
# min_stat: minimum length of the stationary segment
# snr: signal to noise ratio
# N: number of subjects per group
# nROI: number of considered ROI
# rho_ROI: correlation  between the ROIs
# start_number: necessary for output, used later to identify which ROI are dependent (i.e., belong to one subject)
## Output:
# BOLD: simulated BOLD signal for N subjects
# true_cp_loc: true location of the change points per subject
# true_changes: true effect sizes per subject
# onsets: condition onset time series per subject
# nscan, nstim: number of scans and stimulus onsets, respectively
# SNR: Signal to noise ratio
# rho_noise: parameter of the AR(1) process
# rho_spatial: correlation coefficient between ROIs of one subject
# sigma2_w: Variance of the AR(1) process, corresponds to sigma^2_{epsilon_i}
# counter: counts the generated ROI, can be used to later identify which ROI belong to one subject.
## Remarks:
# Note that the baseline of the BOLD signal per subjects corresponds to rnorm(N, 10, 1)
sim_BOLD_canonical_multiple_ROI <- function(mean_change=c(0.5),
                                    var_change=c(0.1),
                                    nscan=250, nstim=60,rho_noise,
                                    min_stat=10, snr, N,
                                    nROI, rho_ROI, start_number){
  rho_ROI <- var_change*rho_ROI
  ## Initialize Stimulus time series: 
  ## this is the same for all participants
  # Random interstimulus intervals with mean such that all stimuli fit in there
  mean_itis <- floor(nscan/nstim)
  # use while loop to ensure that max onset is not larger than nscan.
  stim_done <- 0
  while(stim_done ==0){
    itis <- sample(c(mean_itis-1,mean_itis,mean_itis+1), size=nstim, replace=T) 
    os <- numeric(nstim)
    os[1] <- 1
    for(k in 2:nstim){
      os[k] <- os[k-1]+itis[k-1]
    }
    # check that max(os) is not larger than nscan
    if(max(os)<=nscan){
      stim_done <- 1
    }
  }
  
  stim_1 <-numeric(nscan)
  stim_1[os] <- 1
  
  ## Initialize change points
  # Choose randomly one change point excluding the very beginning and very end (defined by min_stat)
  id_cp_1 <- sort(os) # time of stimulus presentation
  id_cp_len_1 <- length(os)
  
  cp_loc <- sample(id_cp_1[(min_stat+1):(id_cp_len_1-min_stat)],N, replace=T)
  stim_1_1 <- stim_1_2 <- matrix(0,ncol=nscan, nrow=N)
  for(i in 1:N){
    stim_1_1[i,1:(cp_loc[i]-1)] <- stim_1[1:(cp_loc[i]-1)]
    stim_1_2[i,cp_loc[i]:nscan] <- stim_1[cp_loc[i]:nscan]
  }
  
  true_cp <- rep(cp_loc, each=nROI) 
  # true cp returns the first observation of a stationary block, should correspond to an onset.

  ## Initialize the true shape of the two HRFs (now: four HRFs)
  
  # Intercept: not important for shape of HRF, assumed to not change. Add it anyways for estimation
  beta_intercept <- rnorm(N*nROI, 10, 1) 
  
  # For changing HRF:
  # use canonical HRF with changing effect size 
  # Add here the correlation between the ROI, the correlation is assumed to be in the effect and not due to noise
  sig_helper <- matrix(rho_ROI, ncol=nROI, nrow=nROI)
  diag(sig_helper) <- var_change
  
  changes1 <- as.vector(t(mvrnorm(N, mu=rep(mean_change,nROI),Sigma=sig_helper)))
  
  ## These are vectors of length N*nROI with (Sub_1_ROI_1, Sub_1_ROI_2,...,Sub_N_ROI_nROI)
  # If I wanted to have this in matrix form with ncol=nROI and nrow=N: matrix(changes, ncol=nROI, byrow=T)
  
  design_mat <- hrf_mat <- vector("list", N*nROI)
 
  for(i in 1:N){
    list_onsets <- list(first=which(stim_1_1[i,]!=0)*2,
                        second= which(stim_1_2[i,] !=0)*2)
    for(j in 1:nROI){
      design_mat[[(i-1)*nROI+j]] <- simprepTemporal(onsets= list_onsets, duration=list(1,1), totaltime=nscan*2, TR=2, 
                                         effectsize=list(1,(1+changes1[(i-1)*nROI+j])), hrf="double-gamma")
    }
  }
  ## Define Y:
  # add noise (temporal noise and white noise):
  Y <- matrix(nrow=N*nROI, ncol=nscan)
  var_ar <- numeric(N*nROI)
  for(i in 1:(N*nROI)){
    Signal <- simTSfmri(design= design_mat[[i]],
                       base= beta_intercept[i],
                       noise="none") 
    var_epsilon <- (mean(Signal)/snr)^2
    sd_ar <- var_ar[i] <- sqrt(var_epsilon*(1-rho_noise^2))
    Y[i,] <- Signal+arima.sim(model=list(ar=rho_noise),n=nscan, sd=sd_ar)
  }
 var_ar <- var_ar^2
 # BOLD signal is now a matrix with nrow with Sub_1_ROI_1, Sub_1_ROI_2,..., Sub_N_ROI_nROI
 # To ensure that at the end we can put the results in an order and remember which ROIs belong together:
 # we number the observed BOLD signal from start_number to start_number+(N*nROI-1)
 # counter <- c(start_number:(start_number+N*nROI)) (added to output)
 # the ROIs in a block of three belong together, I can always start with absurd large number if I want to make it easier.
  output <- list(BOLD=Y,
                 true_cp_loc = true_cp,
                 true_changes= changes1,
                 onsets = stim_1,
                 nscan=nscan,
                 nstim= nstim, 
                 SNR=snr,
                 rho_noise=rho_noise,
                 rho_spatial=rho_ROI,
                 sigma2_w = var_ar,
                 counter= c(start_number:(start_number+(N*nROI)-1)))
  return(output)
}

### Analysis on group level:
## Input:
# output_posi: the output of N subjects for "sim_BOLD_canonical_one_condition"
# N: number of subjects within a sample on group level
## Output:
# p_values_KH: p-values corresponding to the null hypotheses of eta=0, based on T_KH
# p_values_Wald: p-values corresponding to the null hypotheses of eta=0, based on T_Wald
# est_eta: estimated eta
# between_subject_variance: estimated between-subject variance
# variance_group_level: estimated variance on group level
group_analysis_posi <- function(output_posi, N){
  hat_dif_OLS <-hat_var_naive<-hat_dif_posi<-hat_E_posi<-hat_var_posi <-est_cp<- numeric(N)
  # initialize vectors used in further analysis:
  for(i in 1:N){
    hat_dif_OLS[i]<- output_posi[[i]]$hat_dif
    hat_var_naive[i]<- output_posi[[i]]$hat_var_beta
    hat_dif_posi[i]<- output_posi[[i]]$hat_dif_cd
    hat_E_posi[i]<- output_posi[[i]]$hat_E_cd
    hat_var_posi[i]<- output_posi[[i]]$hat_var_cd
  }
  ## Estimate eta and its variance on group level based on different subject level input:
  # Several combinations of estimates (location and scale estimates)
  group_estimation_naive <- reml_var_estimation_2_posi(b_c=hat_dif_OLS, sigmas_c= hat_var_naive, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi <- reml_var_estimation_2_posi(b_c=hat_dif_posi, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi_E <- reml_var_estimation_2_posi(b_c=hat_E_posi, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi_OLS <- reml_var_estimation_2_posi(b_c=hat_dif_OLS, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  est_eta <- c(group_estimation_naive$beta_group,group_estimation_posi$beta_group, group_estimation_posi_E$beta_group, group_estimation_posi_OLS$beta_group)
  between_subject_variance <- c(group_estimation_naive$sigma_between,group_estimation_posi$sigma_between,group_estimation_posi_E$sigma_between, group_estimation_posi_OLS$sigma_between)
  var_group_level <- list(naive= group_estimation_naive$var_mat, 
                          posi_05=group_estimation_posi$var_mat,
                          posi_E=group_estimation_posi_E$var_mat,
                          OLS_posi=group_estimation_posi_OLS$var_mat)
  ## Compute p-values corresponding to H_0: no change in PM
  p_vals_KH <- p_vals_W <- numeric(4)
  p_vals_KH[1] <- fmri_group_KH(est_eta=est_eta[1], 
                                var_mat_full=var_group_level[[1]],
                                betas_c=hat_dif_OLS, N=N)$p_value
  p_vals_KH[2] <- fmri_group_KH(est_eta=est_eta[2], 
                                var_mat_full=var_group_level[[2]],
                                betas_c=hat_dif_posi, N=N)$p_value
  p_vals_KH[3] <- fmri_group_KH(est_eta=est_eta[3], 
                                var_mat_full=var_group_level[[3]],
                                betas_c=hat_E_posi, N=N)$p_value
  p_vals_KH[4] <- fmri_group_KH(est_eta=est_eta[4], 
                                var_mat_full=var_group_level[[4]],
                                betas_c=hat_dif_OLS, N=N)$p_value
  p_vals_W[1] <- fmri_test_own(est_eta=est_eta[1], 
                               var_mat_full=var_group_level[[1]],
                               N=N)
  p_vals_W[2] <- fmri_test_own(est_eta=est_eta[2], 
                               var_mat_full=var_group_level[[2]],
                               N=N)
  p_vals_W[3] <- fmri_test_own(est_eta=est_eta[3], 
                               var_mat_full=var_group_level[[3]],
                               N=N)
  p_vals_W[4] <- fmri_test_own(est_eta=est_eta[4], 
                               var_mat_full=var_group_level[[4]],
                               N=N)
  ## Name all vectors containing results:
  names(est_eta) <-names(between_subject_variance)<-names(p_vals_KH)<- names(p_vals_W)<- c("naive", "posi_05", "posi_E", "OLS_posi")
  out <- list(p_values_KH = p_vals_KH,
              p_values_Wald= p_vals_W,
              est_eta = est_eta,
              between_subject_variance= between_subject_variance,
              variance_group_level = var_group_level)
  return(out)
}

### Full analysis of simulation results:
## Input:
# dat_cons: output of sim_BOLD_canonical_one_condition
# eff: average effect size, i.e., eta, corresponding to dat_cons
# snr1: SNR corresponding to dat_cons
# B: number of time we draw from dat_cons, i.e., number of iterations for group level analysis
# N: size of each sample for group level analysis
# sample.ids: N x B-matrix indicating which subject in dat_cons are used for group level analysis per iteration
## Output:
# samples: the sample.ids, either as given in the input or the generated ones
# summary: dataframe containing the: rejection rates based on T_KH and T_Wald,  summary_results$rejection_rate_KH <- colMeans(p_vals_KH[,4:7]<= 0.05)
                                  #  average estimated eta across all iterations 
                                  #  variance of the estimated etas
                                  #  average between subject variance
                                  #  variance of the between subject variance
#          all results are given for the four different estimator combinations: naive, posi_05, posi_E and posi_OLS.
# results: raw results; dataframe containing the estimated etas, estimated between subject variance and p-values (T_KH and T_Wald) for the B iterations
# est_var_W: matrix containing the \mathcal{N} estimates for the naive variance and posi-variance.
## Remarks:
# if no sample.ids are given, they are generated 
unknown_cp_results <- function(dat_cons, eff, snr1, B, N,sample.ids=NULL){
  set.seed(eff+1)
  p_vals_KH <- est_etas <- p_vals_W <- est_var_B <- data.frame(name= character(B),
                                                               snr=rep(snr1, B),
                                                               effect = rep(eff,B),
                                                               naive=numeric(B),
                                                               posi_05= numeric(B),
                                                               posi_E= numeric(B),
                                                               posi_OLS=numeric(B))
  cons_data <- vector("list", B)
  p_vals_KH$name <- rep("p_KH", B)
  p_vals_W$name <- rep("p_W", B)
  est_etas$name <- rep("est_eta", B)
  est_var_B$name <- rep("est_var_B", B)
  if(is.null(sample.ids)){
    sample_new <- T
    sample.ids <- matrix(0,ncol=B, nrow=N)
  }else{
    sample_new <- F
  }
  #  
  for(b in 1:B){
    ## 1) Sample from the subject level results
    if(sample_new){
      double <- 1
      while(double !=0){
        sample.ids[,b] <- sort(sample(c(1:length(dat_cons)), N)) #check here: if I do bootstrapping approach, I would sample with replacement. But if we say that our xy iterations are the population, it is more reasonable to sample without replacement?
        double <- sum(apply(sample.ids[,-b], 2, identical, y=sample.ids[,b]))
      }
    }
    ## 2) Do group level analysis
    cons_data[[b]] <- group_analysis_posi(output_posi=dat_cons[sample.ids[,b]], 
                                          N=N)
    ## 3) Initialize vectors most relevant for results
    p_vals_KH[b,4:7] <- cons_data[[b]]$p_values_KH
    p_vals_W[b,4:7] <- cons_data[[b]]$p_values_Wald
    est_etas[b,4:7] <- cons_data[[b]]$est_eta
    est_var_B[b,4:7] <- cons_data[[b]]$between_subject_variance
  }
  results <- rbind(p_vals_KH, p_vals_W, est_etas,est_var_B)
  ## Make a summary of the results (i.e., rejection rates, means)
  summary_results <- data.frame(snr=rep(snr1,4),
                                effect=rep(eff,4),
                                procedure=c("naive", "posi_05", "posi_E", "posi_OLS"),
                                rejection_rate_KH =numeric(4),
                                rejection_rate_W=numeric(4),
                                mean_est_eta=numeric(4),
                                var_est_eta= numeric(4),
                                mean_est_var_B=numeric(4),
                                var_est_var_B=numeric(4))
  
  summary_results$rejection_rate_KH <- colMeans(p_vals_KH[,4:7]<= 0.05)
  summary_results$rejection_rate_W <- colMeans(p_vals_W[,4:7]<= 0.05)
  summary_results$mean_est_eta <- colMeans(est_etas[,4:7])
  summary_results$var_est_eta <- apply(est_etas[,4:7], 2, var, na.rm=T)
  summary_results$mean_est_var_B <- colMeans(est_var_B[,4:7])
  summary_results$var_est_var_B <- apply(est_var_B[,4:7], 2, var, na.rm=T)
  
  mat_est_var <- sapply(dat_cons, est_var)
  
  out <- list(samples= sample.ids,
              summary= summary_results,
              results= results,
              est_var_W= mat_est_var)
  return(out)
}

# Subject level functions (known change points) -----------------------------------------------------------
####
# Apply pre-whitening when computing GLM
####
## Input:
# Y: observed BOLD signal
# X: design matrix
# whitened: should whitened data be returned as well?
## Output:
# coefficients: estimated regression coefficients
# residuals,
# var_mod: variance of the whitened model
# y: whitened signal(if whitened_data=T)
# x: whitened design matrix (if whitened_data=T)
# Y_original: observed signal 
# X_original: initial design matrix 
# V_1: inverse of the variance matrix V accounting for the AR(1) process
## Remark: 
# applies pre-whitening for an AR(1) process
# This function is only applied when change point locations are known;
# Output is used in var_est_subject

glm_prewhitened <- function(Y,X, whitened_data=T){
  ar_order <- 1
  # Estimate the parameter of the AR 1 process based on the residuals
  hat_beta <- solve(crossprod(X))%*%t(X)%*%Y
  res <- Y-X%*%hat_beta
  # estimate rho
  nscan <- length(Y)
  auto_cor <- arima(res, c(1,0,0))$coef
  # sigma*V as the variance of the error terms.
  V <- matrix(0, ncol=nscan, nrow=nscan)
  for(i in 1:nscan){
    for(j in 1:nscan)
      V[i,j] <- auto_cor[1]^(abs(i-j))
  }
  # compute V^-1
  V_1 <- solve(V)
  # compute updated hat_beta
  hat_beta_w <- solve(t(X)%*%V_1%*%X)%*%t(X)%*%V_1%*%Y
  # compute updated model variance
  res2 <- t(Y)%*%V_1%*%Y-t(Y)%*%V_1%*%X%*%hat_beta_w-t(hat_beta_w)%*%t(X)%*%V_1%*%Y+t(hat_beta_w)%*%t(X)%*%V_1%*%X%*%hat_beta_w
  sigma2 <- 1/(nscan-ncol(X))*res2
  if(whitened_data){
    # compute whitened data
    E <- svd(V_1)
    # then V^{-1} = E$u%*%diag(E$d)%*%t(E$u)
    # to get V^{-1/2}
    W <- E$u%*%diag(E$d)^(1/2)%*%t(E$u)
    # update both BOLD and design matrix according to Mumford et al (2006): 
    X_w <- W%*%X
    Y_w <- W%*%Y
    out <- list(coefficients= hat_beta_w,
                residuals= sqrt(res2),
                var_mod= sigma2,
                y= Y_w,
                x=X_w,
                Y_original = Y,
                X_original = X,
                V_1= V_1)
  }else{
    out <- list(coefficients= hat_beta_w,
                residuals= sqrt(res2),
                var_mod= sigma2,
                Y_original = Y,
                X_original = X,
                V_1= V_1)
  }
  
  
  return(out)
}

####
# Point and variance estimator for shape parameters
####

### Estimate variance of beta parameters
## Input:
# glm_model: estimated glm model, where y and x correspond to the observed signal and design matrix used to estimate beta coefficients
# cp number: number of change points per condition
# nbf: number of basis functions
# full_covar: return either the full covariance matrix or set manually the covariance between the segments and the different conditions to zero.
#             in the latter case, only covariance matrix for regression coefficients (stationary + nonstationary) are returned.
## Output:
# hat_var: estimated variance per subject (of model)
# hat_covar: estimated variance-covariance matrix per subject (of model parameters)
## Remarks: 
# variance has to be estimated based on whitened data, that is, 
#y and x that are actually used in the model estimation
var_est_subject <- function(glm_model, cp_number=c(0,0), nbf=3, number_cond=2, full_covar =F){
  Y <- glm_model$y
  X <- glm_model$x
  beta <- glm_model$coefficients
  Ts <- length(Y)
  p <- length(beta)
  v <- Ts-p #degrees of freedom
  var_sub <- 1/v*crossprod(x=(Y-X%*%beta))
  covar_sub <- as.vector(var_sub)*solve(crossprod(X))
  if(full_covar){
    out <- list(hat_var = as.vector(var_sub), 
                hat_covar = covar_sub)
  }else{
    nonstat_cond <- sum(cp_number!=0)
    stat_cond <- number_cond-nonstat_cond
    no_nonstat_seg <- sum(cp_number)+nonstat_cond 
    mat_covar_sub <- matrix(0, nrow=no_nonstat_seg*nbf+stat_cond*nbf, ncol=no_nonstat_seg*nbf+stat_cond*nbf)
    for(i in 1:(no_nonstat_seg+stat_cond)){
      mat_covar_sub[((i-1)*nbf+1):(nbf*i),((i-1)*nbf+1):(nbf*i)] <- covar_sub[((i-1)*nbf+1):(nbf*i),((i-1)*nbf+1):(nbf*i)]
    }
    out <- list(hat_var=as.vector(var_sub),
                hat_covar = mat_covar_sub)
  }
  return(out)
}

### Estimate location and variance of shape parameters
## Input:
# assume that beta parameters follow a normal distribution with known covariance matrix
# betas: list with betas which have been estimated
# covar_betas: corresponding covariance matrix of the betas
# hrf: basis function of the hrf
# nHR: number of HRs in model
# mc_B: on how many repetitions is the estimated variance based?
# fix_seed:fix a seed to estimate variance to ensure replicability
## Output:
# point_estimators: estimators for the shape parameters
# variances: estimated variance of the shape parameters
est_var_sp_subject <- function(betas, covar_betas,nHR=NULL, hrf=hrf_basis_fun, mc_B=10000, fix_seed=14){
  # how many HR am I interested in?
  nbetas_HR <- ncol(hrf)
  if(is.null(nHR)){
    nHR <- length(betas)/nbetas_HR
  }
  if(nHR%%1!=0)stop("Number betas and number basis function do not correspond")
  
  # initiate list in which the intermediate results are stored
  # one list element= one HR, each list element: matrix including the shape parameters
  point_est <- lapply(rep(0,10,), matrix, ncol=nHR, nrow=1)
  hr_mc <- lapply(rep(0, 10),matrix, ncol=nHR, nrow=mc_B)
  names(hr_mc) <- names(point_est) <- c("PM", "NA", "IUA", "TTP", "TPN", "FWHM", "FWHN", "AUC", "AC_P", "DC_P") #these are currently the names of all possible shape parameters computed, not all need to be investigated on group level!
  
  for(j in 1:nHR){
    # these are the betas that are relevant for one HR:
    relevant_betas <- betas[((j-1)*nbetas_HR+1):(j*nbetas_HR)]
    # combine these betas with the basis functions:
    HR <- as.vector(t(relevant_betas)%*%t(hrf))
    shape <- shape_parameters2(hrf_sampled = HR)
    for(sp in 1:10){
      point_est[[sp]][1,j] <- shape[sp]
    }
  }
  ## now compute the shape parameters
  # draw betas from multivariate normal distribution
  set.seed(fix_seed)
  betas_mc <- mvrnorm(mc_B, betas[1:(nbetas_HR*nHR)], covar_betas)
  # compute the HRs and the shape parameters
  for(i in 1:mc_B){
    # for each HR we have nbetas_HR betas
    # each column is one HR
    for(j in 1:nHR){
      # these are the betas that are relevant for one HR:
      relevant_betas <- betas_mc[i,((j-1)*nbetas_HR+1):(j*nbetas_HR)]
      # combine these betas with the basis functions:
      HR <- colSums(t(hrf)*relevant_betas)
      shape_mc <- shape_parameters2(hrf_sampled = HR)
      for(sp in 1:10){
        hr_mc[[sp]][i,j] <- shape_mc[sp]
      }
    }
  }
  # now I have a list per shape parameter with ncol=nHR and nrow=mc_B
  # now compute the covariance and variance for each HR per shape parameter
  variances <- lapply(hr_mc, var, na.rm=T)
  out <- list(point_estimators= point_est,
              variances = variances)
  return(out)
}


# Subject level functions (unknown change points) ------------------------------

### Functions to compute estimators for location and variance of the beta parameters for one person
## Input:
# BOLD_signal: observed BOLD signal
# onsets: onsets time series
# nscan, nstim: number of scans (T) and number of stimuli for one condition
# min_seg defines the minimum length of a stationary segment.
# fixed.seed: If I have a fixed seed, it has the length of iterations on which computation for one P(theta_star) is based
#             In this case, we draw Y_star directly from a multivariate normal with a fixed seed at all times.
# full_X: use complete (full) model matrix to generate samples Y_star? Must be false.
# N_iter_CD: on how many iterations for each possible theta value should P(t_obs < t_star) be based?
# N_iter_int: on how many iterations should construction of set E with considered theta-values be based?
# true_cp: should the true change point location among the considered change point locations? Only relevant for simulation study, otherwise NULL.
#          False means that the true change point location might not be considered.
# n_cons_cp: the number of change points/ models that we consider.
# min_dif_between_cp: I ensure that there are at least that many steps between the considered change point locations. 
#                     Besides around true change point location, the difference between the considered change points is a Vielfaches of min_dif_between_cp
#                     Can be set to one if I want random change point locations
#                     Only needed when true change point location is known, otherwise we use equi-distance between all possible change point locations.
# miss_spec: by how much is true change point misspecified. If true change point location is considered: set to zero.
# number: needed if several ROI are used to be able to identify ROIs belonging to one subject.
## Output:
# number: id for the considered ROI
# hat_dif: estimated beta coefficient [OLS] corresponding to the estimated change in the peak magnitude
# hat_var_beta: estimated variance of hat_dif; WITHOUT accounting for model selection 
# hat_var_cd: post-selection variance of beta
# hat_dif_cd: estimated beta coefficient [median of the post-selection confidence distribution]
# hat_E_cd: expected value of the post-selection confidence distribution
# Conf_dist: approximation of the confidence distribution
# est_cp: estimated change point location
## Remarks:
# This function is for one condition with one change point and multiple ROIs
est_betas_cd <- function(BOLD_signal, onsets, nscan, nstim,
                         min_seg,
                         Var_known=NULL,
                         last.length.theta.star=50, 
                         fixed.seed=NULL,
                         N_iter_CD=500,N_iter_int=250,
                         full_X=F, true_cp=NULL, n_cons_cp=4,
                         min_dif_between_cp=5, miss_spec=0, number){
  id_cp <- which(onsets==1)
  id_all_pos_cp <- id_cp[(min_seg+1):(nstim-min_seg+1)]
  length_all_pos_cp <- length(id_all_pos_cp)
  if(!is.null(true_cp)){
    if(true_cp%in%id_all_pos_cp){
      true_cp_1 <- true_cp}else{
        true_cp_1 <- id_all_pos_cp[min(which(id_all_pos_cp>true_cp))]
      }
    if(miss_spec!=0){
      id_true <- which(id_all_pos_cp==true_cp_1)
      id_mis <- id_true+sample(seq(-miss_spec,miss_spec,1),1)
      cp_mis <- id_all_pos_cp[id_mis]
      if(min_dif_between_cp==1){
        id_all_pos_cp <- id_all_pos_cp[-c(id_true, id_mis)]
        considered_cp <- sort(c(cp_mis, id_all_pos_cp[sample(c(1:(length_all_pos_cp-2)), size=n_cons_cp-1)]))
      }else{
        seq_cons <- seq(from=1, to=length_all_pos_cp, by=min_dif_between_cp) #sequence of considered cp with a minimum of 5 between them; crude approximation
        if(id_true%in%seq_cons| cp_mis%in%seq_cons){
          # remove the true_cp and/or cp_mis from the other considered change point locations
          considered_cp <- sort(c(cp_mis, id_all_pos_cp[sample(seq_cons[-which(seq_cons%in%c(id_true, id_mis))], size=n_cons_cp-1)]))
        }else{
          seq_helper <- c(-5, seq_cons, max(seq_cons)+10)
          # using the seq_helper ensures that at least one value is larger/ smaller than cp_mis, even if cp_mis is at an extreme position.
          id_remove <- c(max(which(seq_helper<cp_mis)), min(which(seq_helper>cp_mis)))
          # if id_remove contains the first or last entry in seq_helper: remove it
          # note: both cannot happen at the same time!
          if(1%in%id_remove){
            id_remove <- id_remove[2]
          }
          if(length(seq_helper)%in% id_remove){
            id_remove <- id_remove[1]
          }
          considered_cp <- sort(c(cp_mis, id_all_pos_cp[sample(seq_cons[-(id_remove-1)], size=n_cons_cp-1)]))
        }
      }
    }else{
      if(min_dif_between_cp==1){
        id_all_pos_cp <- id_all_pos_cp[-which(id_all_pos_cp==true_cp_1)]
        considered_cp <- sort(c(true_cp_1, id_all_pos_cp[sample(c(1:(length_all_pos_cp-1)), size=n_cons_cp-1)]))
      }else{
        id.true.cp <- which(id_all_pos_cp==true_cp_1)
        seq_cons <- seq(from=1, to=length_all_pos_cp, by=min_dif_between_cp) #sequence of considered cp with a minimum of 5 between them; crude approximation
        if(id.true.cp%in%seq_cons){
          considered_cp <- sort(c(true_cp_1, id_all_pos_cp[sample(seq_cons[-which(seq_cons==id.true.cp)], size=n_cons_cp-1)]))
        }else{
          seq_helper <- c(-5, seq_cons, max(seq_cons)+10)
          # using the seq_helper ensures that at least one value is larger/ smaller than the id.true.cp, even if it is at an extreme position.
          id_remove <- c(max(which(seq_helper<id.true.cp)), min(which(seq_helper>id.true.cp)))
          # if id_remove contains the first or last entry in seq_helper: remove it
          # note: both cannot happen at the same time!
          if(1%in%id_remove){
            id_remove <- id_remove[2]
          }
          if(length(seq_helper)%in% id_remove){
            id_remove <- id_remove[1]
          }
          considered_cp <- sort(c(true_cp_1, id_all_pos_cp[sample(seq_cons[-(id_remove-1)], size=n_cons_cp-1)]))
        }
      }
    }
    
  }else{
    considered_cp <- id_all_pos_cp[round(seq(from=1, to=length_all_pos_cp, length.out=n_cons_cp),0)]
  }
  
  # Initialize X_full, each considered cp is an additional condition with onsets starting from the cp_location
  list_onsets <- vector("list", n_cons_cp+1)
  models.id <- vector("list", n_cons_cp)
  
  # Onsets for the full condition
  list_onsets[[1]] <- id_cp*2
  for(i in 1:n_cons_cp){
    # the onsets for the dummy conditions
    dummy_cond <- onsets
    dummy_cond[1:(considered_cp[i]-1)] <- 0
    list_onsets[[i+1]] <- which(dummy_cond!=0)*2 #this doesn't work because the first one will always be an onset!
    # which dummy condition is considered in the different models?
    models.id[[i]] <- c(1,(i+1), (n_cons_cp+2)) 
    #at the last position of the design matrix we have the intercept, which is also included in each model
  }
  duration_list <- lapply(rep(1, n_cons_cp+1), mean)
  
  X_hrf <- specifydesign(onsets= list_onsets, durations=duration_list, 
                         totaltime=nscan*2, TR=2, 
                         effectsize=duration_list, conv="double-gamma")
  for(i in 1:n_cons_cp){
    X_hrf[1:(list_onsets[[i+1]][1]/2-1),i+1] <- 0
  }
  X_full <- cbind(X_hrf, rep(1, nscan))
  ### Select a model
  # Use the selected_model function to choose a model
  if(is.null(Var_known)){
    ar <- T
    mod_select <- selected_model(Y=BOLD_signal, X=X_full, id.models = models.id, AR=ar, Var=Var_known)
    ### Compute for this model the variance of the parameters, ignoring post-selective inference
    # I have the AR covariance matrix and the estimates of beta, therefore I only need the variance of the model
    X_selected <- X_full[,models.id[[mod_select[[1]]]]]
    var_mod <- 1/(nscan-ncol(X_selected))*(t(BOLD_signal)%*%mod_select$V_1%*%BOLD_signal-t(BOLD_signal)%*%mod_select$V_1%*%X_selected%*%mod_select$hat_beta-t(mod_select$hat_beta)%*%t(X_selected)%*%mod_select$V_1%*%BOLD_signal+t(mod_select$hat_beta)%*%t(X_selected)%*%mod_select$V_1%*%X_selected%*%mod_select$hat_beta)
    var_betas <- as.vector(var_mod)*solve(t(X_selected)%*%mod_select$V_1%*%X_selected) 
  }else{
    ar <- F
    E <- svd(Var_known$V)
    V.05 <- E$u%*%diag(E$d)^(1/2)%*%t(E$u)
    V_05 <- solve(V.05)
    BOLD_w <- V_05%*%BOLD_signal
    X_w <- V_05%*%X_full
    mod_select <- selected_model(Y=BOLD_w, X=X_w, id.models = models.id, AR=ar, Var=Var_known)
    ### Compute for this model the variance of the parameters, ignoring post-selective inference
    # I have the AR covariance matrix and the estimates of beta, therefore I only need the variance of the model
    X_selected <- X_w[,models.id[[mod_select[[1]]]]]
    var_mod <- 1/(nscan-ncol(X_selected))*(t(BOLD_w)%*%BOLD_w-t(BOLD_w)%*%X_selected%*%mod_select$hat_beta-t(mod_select$hat_beta)%*%t(X_selected)%*%BOLD_w+t(mod_select$hat_beta)%*%t(X_selected)%*%X_selected%*%mod_select$hat_beta)
    var_betas <- as.vector(var_mod)*solve(crossprod(X_selected))
  }
  
  ### Include POSI
  # Reorder X_full and id.models such that covariates of the selected model are the first entries in X_new 
  #   and id.models_new[[1]] corresponds to the selected model
  # since X_full is ordered according to timing, the selected model corresponds to the first column in X_full and the mod_select[[1]]+1 column in X_full
  id_not_select <- c(2:(n_cons_cp+1))[-(mod_select[[1]])]
  # the selected model is the full onset time series+ onset time series starting at cp_location + intercept
  # the intercept is always at the last place in X_full and therefore has to be at the last place in X_ord, otherwise id.mod would not be correct.
  if(is.null(Var_known)){
    X_ord <- X_full[,c(1,(mod_select[[1]]+1),id_not_select, (n_cons_cp+2))]
  }else{
    X_ord <- X_w[,c(1,(mod_select[[1]]+1),id_not_select, (n_cons_cp+2))]
    BOLD_signal <- BOLD_w
  }
  models.id.ord <- models.id[c(mod_select[[1]], id_not_select-1)]
  
  ### Compute the Confidence distribution
  # Run F_t_given_u with X_new and id.models_new
  # I will use the fixed seed to hopefully make the estimation of Y_star slightly faster.
  # Thus, I will need to update the seed every time for the optimization! F_t_given_u_posi
  # starting eta include the betas, rho and sigma
  # we are always interested in the change, therefore in the value of the second beta parameter, which is not included in eta then.
  if(is.null(Var_known)){
    # if variance is not known, eta includes estimates for rho and variance
    est_eta <- c(mod_select$hat_beta[-2], mod_select$rho, var_mod)
  }else{
    # if the variance is known, eta only includes the estimates for the regression parameters excluding the focus parameter
    est_eta <- c(mod_select$hat_beta[-2])
  }
  est_theta <- mod_select$hat_beta[2]
  #if we use the full design matrix, we need to adapt eta accordingly:
  if(full_X){
    if(is.null(Var_known)){
      # ensure here that rho and var_mod are at the end.
      # add besides that n_considered_cp-1
      n_eta <- length(est_eta)
      est_eta <- c(est_eta[1:(n_eta-2)], rep(0, n_cons_cp-1), est_eta[(n_eta-1):n_eta])
    }else{
      est_eta <- c(est_eta, rep(0, n_cons_cp-1))
    }
  }
  # if we do not use seed, we use est_eta for sample generation; need several of them.
  #if(is.null(fixed.seed)){
  #  est_eta <- matrix(rep(est_eta, N_iter_P), byrow=T, ncol=length(est_eta))
  #}
  P_star <- F_t_given_u_posi_2(X_obs=X_ord, 
                               Y_obs=BOLD_signal, 
                               focus_id=2, 
                               models_id=models.id.ord,
                               starting_eta= est_eta, 
                               est_theta= est_theta,
                               V_known=Var_known, 
                               last_length_theta_star=last.length.theta.star,
                               X_f = full_X,adapt_start=F,
                               N_iter_int=N_iter_int, N_iter_CD=N_iter_CD,
                               break_while=20000)
  ###
  ## Compute the variance based on P_star
  ###
  ## Compute the confidence density
  # first: remove any values that have less than 100 iterations (we cannot trust them...)
  id_remove <- which(P_star$iterations< 100)
  if(length(id_remove)!=0){
    F_cd <- 1-round(P_star$mean_theta_star_leq_theta[-id_remove],4)
    last.length.theta.star <- length(F_cd)
  }else{
    F_cd <- 1-round(P_star$mean_theta_star_leq_theta,4)
  }
  ## Check for unreasonable values; these can happen due to using Monte Carlo methods to compute the confidence distribution
  # Currently, we simplify life by substituting the unreasonable value with the mean of the values before and after
  # Repeat until there are no unreasonable values left (check after each update)
  id_larger <- which(F_cd[-1]< F_cd[-last.length.theta.star])
  while(length(id_larger)!=0){
    for(s in 1:length(id_larger)){
      F_cd[id_larger[s]+1] <- mean(c(F_cd[id_larger[s]], F_cd[id_larger[s]+2]))
    }
    id_larger <- which(F_cd[-1]< F_cd[-last.length.theta.star])
  }
  # what if the extreme values do not fit?
  
  ## Compute the confidence density; use 1-CD
  
  f_cd <- c(F_cd[1], F_cd[-1]-F_cd[-last.length.theta.star])
  
  ## Compute the expected value, as well as the expected value of theta^2
  if(length(id_remove)!=0){
    P_star_theta <- P_star$theta_star[-id_remove]
  }else{
    P_star_theta <- P_star$theta_star
  }
  
  E_cd <- sum(P_star_theta*f_cd, na.rm=T)
  Var_cd <- sum(P_star_theta^2*f_cd, na.rm=T)-E_cd^2
  
  ## Find the value for which F_cd=0.5
  # use weighted mean of values larger than 0.5 and lower than 0.5:
  F_05 <- F_cd-0.5
  id.05 <- c(min(which(F_05>0)), max(which(F_05<0)))
  dif_cd <- (1-abs(F_05[id.05[1]])/sum(abs(F_05[id.05])))*P_star_theta[id.05[1]]+(1-abs(F_05[id.05[2]])/sum(abs(F_05[id.05])))*P_star_theta[id.05[2]]
  
  ## Find the value at which confidence distribution has probability of 0.5 (roughly)
  # Since there will probably not be a single value, use the weighted mean?
  # or use the expected value? In many cases, the ex
  return(list(number=number,
              hat_dif = mod_select$hat_beta[2],
              hat_var_beta= var_betas[2,2],
              hat_var_cd = Var_cd,
              hat_dif_cd= dif_cd,
              hat_E_cd= E_cd,
              Conf_dist = P_star,
              est_cp= considered_cp[mod_select[[1]]]))
  
}


# Post-selection confidence distribution functions ---------------------------------------
require(nnet)

####
# Compute the sufficient statistics
####

### Compute the sufficient for the focus parameter when the focus parameter is not the variance!
## Input:
# Y: the observed response
# X: the full model, that is, X includes every covariate considered in at least one model.
# id.focus: which beta are we interested in, needs to corresponds to the regression coefficients of X
# id.selected: which covariates are included in the selected model? Only necessary if variance needs to be estimated.
# AR: should sufficient parameters be computed for Y following an AR(1)-process (TRUE) or for independence of Y (FALSE)
# Var: Does the variance need to be estimated?
## Output:
# focus: sufficient statistic for the focus parameter
# nuisance sufficient statistics for the nuisance parameters
suff_stat_posi <- function(Y, X,id.focus, id.selected=NULL, AR=T, Var=T){
  suff_focus <- as.vector(t(X[,id.focus])%*%Y)
  N <- length(Y)
  if(AR){
    zeros <- rep(0, N-1)
    diag_1 <- diag(1, nrow=N-1, ncol=N-1)
    G <- cbind(zeros, diag_1)
    Tilde_G <- cbind(diag_1, zeros)
    Bar_G <- rbind(c(1, rep(0,N-1)), c(rep(0,N-1),1))
    Tilde_G_Y <- Tilde_G%*%Y
    G_Y <- G%*%Y
    Bar_G_Y <- Bar_G%*%Y
    suff_nuisance <- c(crossprod(x=Tilde_G%*%X[,-id.focus], y=G_Y)+crossprod(x=G%*%X[,-id.focus],y=Tilde_G_Y),
                       2*crossprod(X[,-id.focus], Y)- crossprod(x=Bar_G%*%X[,-id.focus], y=Bar_G_Y),
                       crossprod(X[,-id.focus],Y),
                       crossprod(Tilde_G_Y, G_Y),
                       2*crossprod(Y)-crossprod(Bar_G_Y))
  }else{
    suff_nuisance <- c(crossprod(X[,-id.focus],Y))
  }
  if(Var){
    w <- Y-X[,id.selected]%*%solve(crossprod(X[,id.selected]))%*%t(X[,id.selected])%*%Y
    suff_nuisance <- c(suff_nuisance,crossprod(Y), w)
  }
  return(list(focus = suff_focus,
              nuisance= suff_nuisance))
}

### Optimize eta parameters such that eta_obs = eta_star.
## Input:
# eta: these are all the beta parameters IN THE SELECTED MODEL except the focus parameter, followed by rho, followed by the variance
# X_obs: same as X in suff_stat_pos, ensure that covariates belonging to the selected models are the first covariates in X
# w_i: only needed if AR=F or fix_seed=NULL, corresponds to uniformly distributed random variables
# theta: fixed value for focus parameter
# id.focus, AR, Var as in suff_stat_posi
# id.models: list of vectors which show which covariates are included in which model
# N: number of observations, i.e., length of Y that should be simulated
# fix_seed: if AR=T and fix_seed=NULL, we draw Y_star based on pre-whitenened mu and color it afterwards
#           if AR=T and fix_seed is an integer, we draw Y_star from a multivariate normal with Sigma given, always using the same seed before drawing from Sigma.
# Var: Either NULL  or the known variance
## Output: 
# Squared sum of the difference between the simulated sufficient statistics and the observed sufficient statistics.
## Remark: 
# When Var is known, we assume that X_obs is pre-whitenend!
# This function is used to find values for eta such that sufficient statistic for simulated Y are close to observed sufficient statistics.
U_star_posi <- function(eta, X_obs, w_i, u_obs, theta, id.focus,N, id.models, AR=T, Var=NULL, fix_seed=NULL){
  n_eta <- length(eta)
  X_selected <- X_obs[,id.models[[1]]]
  if(!is.null(Var)){
    AR <- F
  }
  if(AR){
    # Specify the variance-covariance matrix:
    # rho as the second to last entry in eta.
    # ensure that rho is between -1 and 1.
    if(abs(eta[(n_eta-1)])>=1){
      eta[(n_eta-1)] <- (eta[(n_eta-1)]-0.1)/abs(eta[(n_eta-1)])
    }
    V <- matrix(0, ncol=N, nrow=N)
    for(i in 1:N){
      for(j in 1:N){
        V[i,j] <- eta[(n_eta-1)]^(abs(i-j))
      }
    }
    if(is.null(fix_seed)){
      # draw pre-whitened sample and then color the whole sample
      # compute V^0.5 and V^-0.5
      E <- svd(V)
      V.05 <- E$u%*%diag(E$d)^(1/2)%*%t(E$u)
      V_05 <- solve(V.05)
      # whiten X to get whitened mean
      # compute the mean only based on the selected model.
      X_whitened <- V_05%*%X_selected
      mu <- X_whitened[,id.focus]*theta+X_whitened[,-id.focus]%*%eta[1:(n_eta-2)]
      Y_whitened <- qnorm(w_i, mean=mu, sd=sqrt(abs(eta[n_eta]))) 
      # color the sample again.
      Y_star <- V.05%*%Y_whitened
      # based on Y_star we can now compute the sufficient statistics.
      # the goal is to reduce this to be as close to zero as possible!
      var.suff <- is.null(Var)
      suff.stat <- suff_stat_posi(Y=Y_star, X=X_obs, id.focus=id.focus, id.selected=id.models[[1]], AR=AR, Var=var.suff)
      out <- sum((suff.stat$nuisance-u_obs)^2)
    }else{
      # with a fixed seed, I achieve the same result as fixing w_i and using the inverse. 
      # randomness comes from changing the seed for each iteration of w_i
      mu <- X_obs[,id.focus]*theta+X_selected[,-id.focus]%*%eta[1:(n_eta-2)]
      
      # this might return values that are not positive definite, this depends on the values of eta that are tested.
      # So I need to ensure that sigma is considered positive definite (based on the definition of abs)
      suff.stat <- tryCatch_sigma(eta=eta, theta=theta, mu=mu, fix_seed=fix_seed, 
                                  length_u=length(u_obs), X_obs=X_obs, id.focus=id.focus, 
                                  id.selected=id.selected, AR=AR, Var=Var)
      out <- sum((suff.stat$nuisance-u_obs)^2)
    }
  }else{
    mu <- X_obs[,id.focus]*theta+X_selected[,-id.focus]%*%eta
    Y_star <- qnorm(w_i, mean=mu, sd=sqrt(Var[[1]]))
    suff.stat <- suff_stat_posi(Y=Y_star, X=X_obs, id.focus=id.focus, id.selected=id.models[[1]], AR=F, Var=F)
    out <- sum((suff.stat$nuisance-u_obs)^2)
  }
  if(is.infinite(out)){out <- 10^16}
  return(out)
}

### Compute the sufficient statistic of the focus parameter for a simulated Y 
## Input:
# As in U_star_posi
## Output:
# suff_stat: sufficient statistic for the focus parameter for simulated Y
# mod_selected: selected model for the simulated Y
suff_stat_star_posi <- function(X_obs,eta,theta, w_i, id.focus,id.models, fix.seed, Var, full_X){
  # Initialize Y_star as in U_star_posi
  n_eta <- length(eta)
  if(full_X){
    X_selected <- X_obs
  }else{
    X_selected <- X_obs[,id.models[[1]]]
  }
  
  if(!is.null(Var)){
    AR <- F
  }else{AR <- T}
  if(AR){
    # Specify the variance-covariance matrix:
    # rho as the second to last entry in eta.
    # ensure that rho is in between -1 and 1.
    if(abs(eta[(n_eta-1)])>=1){
      eta[(n_eta-1)] <- (eta[(n_eta-1)]-0.1)/abs(eta[(n_eta-1)])
    }
    N <- dim(X_obs)[1]
    V <- matrix(0, ncol=N, nrow=N)
    for(i in 1:N){
      for(j in 1:N){
        V[i,j] <- eta[(n_eta-1)]^(abs(i-j))
      }
    }
    if(is.null(fix.seed)){
      # draw pre-whitened sample and then color the whole sample
      # compute V^0.5 and V^-0.5
      E <- svd(V)
      V.05 <- E$u%*%diag(E$d)^(1/2)%*%t(E$u)
      V_05 <- solve(V.05)
      # whiten X to get whitened mean
      # compute the mean only based on the selected model.
      X_whitened <- V_05%*%X_selected
      mu <- X_whitened[,id.focus]*theta+X_whitened[,-id.focus]%*%eta[1:(n_eta-2)]
      Y_whitened <- qnorm(w_i, mean=mu, sd=sqrt(abs(eta[n_eta]))) 
      # color the sample again.
      Y_star <- V.05%*%Y_whitened
    }else{
      # with a fixed seed, I achieve the same result as fixing w_i and using the inverse. 
      mu <- X_obs[,id.focus]*theta+X_selected[,-id.focus]%*%eta[1:(n_eta-2)]
      set.seed(fix.seed)
      Y_star <- mvrnorm(n=1, mu=mu, Sigma = abs(eta[n_eta])*V, tol= 10^-13)
    }
  }else{
    mu <- X_obs[,id.focus]*theta+X_selected[,-id.focus]%*%eta
    Y_star <- qnorm(w_i, mean=mu, sd=sqrt(Var[[1]]))
  }
  out <- list(suff_stat= as.vector(t(X_obs[,id.focus])%*%Y_star),
              mod_selected = selected_model(Y=Y_star, X=X_obs, id.models = id.models, AR=AR, Var=Var)[[1]])
  return(out)
}

## Helper
tryCatch_sigma <- function(eta, theta, mu, fix_seed, length_u, X_obs, id.focus, id.selected, AR, Var){
  tryCatch({
    set.seed(fix_seed)
    Y_star <- mvrnorm(n=1, mu=mu, Sigma = abs(eta[n_eta])*V, tol= 10^-13)
    suff.stat <- suff_stat_posi(Y=Y_star, X=X_obs, id.focus=id.focus, id.selected=id.selected, AR=AR, Var=Var)
    return(suff.stat)
  }, error=function(e2){
    suff.stat <- rep(10^16, length_u)
  }
  )
}

####
# Model selection
####

### Computation of the model likelihood for AR(1) processes for known variance
## Input:
# X_obs: design matrix
# Y_obs: observed signal
# Var_obs: known variance
## Output:
# loglik: Loglikelihood of the model
# hat_beta: estimated beta coefficients 
# V_1: inverse of the variance accounting for the AR(1) process
## Remark:
# The data is not pre-whitenend
loglik_known_var_direct <- function(X_obs,Y_obs,Var_obs){
  N <- length(Y_obs)
  V_1 <- solve(Var_obs$V)
  hat_beta <- solve(crossprod(X_obs))%*%t(X_obs)%*%Y_obs
  loglik <- 1/Var_obs$sigma*(t(Y_obs)%*%V_1%*%X_obs%*%hat_beta-1/2*t(Y_obs)%*%V_1%*%Y_obs-1/2*t(X_obs%*%hat_beta)%*%V_1%*%X_obs%*%hat_beta)-1/2*log(det(Var_obs$sigma*Var_obs$V))-N/2*log(2*pi)
  return(list(loglik=loglik,
              hat_beta=hat_beta,
              V_1=V_1))
}
### Computation of the model likelihood for independent data and known variance
## Input:
# X_w: whitenend design matrix
# Y_w: whitenend signal
# sigma_obs: known variance of the AR process
## Output:
# loglik: Loglikelihood of the model
# hat_beta: estimated beta coefficients 
## Remark:
# we use pre-whitened data, i.e., Y is independent
loglik_known_var_whitened <- function(X_w, Y_w, sigma_obs){
  N <- length(Y_w)
  hat_beta <- solve(crossprod(X_w))%*%t(X_w)%*%Y_w
  loglik <- 1/sigma_obs*(t(Y_w)%*%X_w%*%hat_beta-1/2*t(Y_w)%*%Y_w-1/2*t(X_w%*%hat_beta)%*%X_w%*%hat_beta)-log(sqrt(sigma_obs))-N/2*log(2*pi)
  return(list(loglik=loglik,
              hat_beta=hat_beta))
}
### Computation of the model likelihood for independent data and unknown variance
## Input:
# X_obs: design matrix
# Y_obs: observed signal
## Output:
# loglik: Loglikelihood of the model
# hat_beta: estimated beta coefficients 
# V_1: inverse of the estimated variance matrix accounting for the AR(1) process
# rho: estimated parameter of the AR(1) process
loglik_unknown_var_direct <- function(X_obs,Y_obs){
  N <- length(Y_obs)
  # Estimate beta
  hat_beta <- solve(crossprod(X_obs))%*%t(X_obs)%*%Y_obs
  res <- Y_obs-X_obs%*%hat_beta
  # Estimate rho and sigma^2
  ar.mod <- arima(res, order=c(1,0,0))
  V <- matrix(0, ncol=N, nrow=N)
  for(i in 1:N){
    for(j in 1:N){
      V[i,j] <- ar.mod$coef[1]^(abs(i-j))
    }
  }
  V_1 <- solve(V)
  loglik <- 1/ar.mod$sigma2*(t(Y_obs)%*%V_1%*%X_obs%*%hat_beta-1/2*t(Y_obs)%*%V_1%*%Y_obs-1/2*t(X_obs%*%hat_beta)%*%V_1%*%X_obs%*%hat_beta)-1/2*log(det(ar.mod$sigma2*V))-N/2*log(2*pi)
  return(list(loglik=loglik,
              hat_beta= hat_beta,
              V_1= V_1,
              rho=ar.mod$coef[1]))
}

### model selection
## Input:
# Y: signal, if AR=F: whitenend signal
# X: design matrix, including all covariates of all models, if AR=F: whitenend signal
# id.models: list of vectors which show which covariates are included in which model
# Var: if variance is known: list with element sigma= sigma^2 and V, where V[i,j]=rho^(abs(i-j))
# AR: Do we need to account for AR process. I.e., if False: X and Y are pre-whitenend
## Output:
# mod_select: indicator for selected model in id.models
# hat_beta: estimated regression coefficients for selected model
# V_1: inverse of the estimated variance matrix accounting for the AR(1) process for the selected model
# rho: estimated parameter of the AR(1) process for the selected model
selected_model <- function(Y, X,id.models, AR=T, Var=NULL){
  nm <- length(id.models)
  loglik <- rho <-numeric(nm)
  hat_beta <- V_1 <-vector("list", nm)
  for(i in 1:nm){
    if(AR){
      if(is.null(Var)){
        loglik_result <- loglik_unknown_var_direct(X_obs=X[,id.models[[i]]], Y_obs=Y)
        loglik[i] <- loglik_result$loglik
        V_1[[i]] <- loglik_result$V_1
        hat_beta[[i]] <- loglik_result$hat_beta
        rho[i] <- loglik_result$rho
        
      }else{
        loglik_result <- loglik_known_var_direct(X_obs=X[,id.models[[i]]], Y_obs=Y, Var_obs = Var)
        loglik[i] <- loglik_result$loglik
        V_1[[i]] <- loglik_result$V_1
        hat_beta[[i]] <- loglik_result$hat_beta
      }
    }else{
      loglik_result <- loglik_known_var_whitened(X_w=X[,id.models[[i]]], Y_w=Y, sigma_obs = Var$sigma)
      loglik[i] <- loglik_result$loglik
      hat_beta[[i]] <- loglik_result$hat_beta
    }
  }
  mod_select <- which.is.max(loglik)
  if(!is.null(Var)){
    rho[mod_select] <- Var$V[1,2]
    if(!AR){
      V_1[[mod_select]] <- solve(Var$V)
    }
  }
  return(list(mod_select = mod_select,
              hat_beta = hat_beta[[mod_select]],
              V_1 = V_1[[mod_select]],
              rho= rho[mod_select]))
}

####
# Compute the post-selection confidence distribution
####

### Function that returns for a fixed theta a matrix with both simulated t_star and the value of the optimization of the nuisance parameters.
## Input:
# eta_start: starting values for eta used in optimization of U_star_posi
# X_obs, Y_obs: design matrix (including all covariates in \mathcal(M) and observed signal, respectively
# suff_obs: observed sufficient statistics
# theta: fixed value of the focus parameter
# id.focus: indice for the focus parameter in X_obs
# id.models: list of vectors which show which covariates are included in which model
# N: length of generated Y
# Variance, seed.fix, AR, full_X as in U_star_posi or suff_stat_star_posi
# adapt_eta: should values of eta at which optimimum was found be returned?
## Output:
# suff.stat: the sufficient statistics of the focus parameter for the generated samples
# opt_value: Output of U_star_posi, i.e.: how well did optimization regarding sufficient statistics for nuisance parameter work?
# If adapt_eta=T: values of nuisance parametes used for sample generation
t_obs_star_posi <- function(eta_start, X_obs, Y_obs,suff_obs, theta, id.focus, id.models, 
                            N, Variance=NULL, seed.fix =NULL, AR=T, adapt_eta=F, full_X=F){
  w <- runif(N)
  if(is.null(Variance)){
    var <- T
    eta_opt <- optim(par=eta_start, fn=U_star_posi, 
                     X_obs=X_obs, w_i=w, 
                     u_obs=suff_obs$nuisance, theta=theta,
                     id.focus=id.focus, N=N,
                     id.models=id.models,AR=AR, 
                     Var=Variance, fix_seed=seed.fix, 
                     method=c("Nelder-Mead"))
  }else{
    var <- F
    eta_opt <- optim(par=eta_start, fn=U_star_posi, 
                       X_obs=X_obs, w_i=w, 
                       u_obs=suff_obs$nuisance, theta=theta,
                       id.focus=id.focus, N=N,
                       id.models=id.models,AR=AR, 
                       Var=Variance, fix_seed=seed.fix,
                       method=c("Nelder-Mead"))
  }
  # If optimization did not work:
  if(eta_opt$convergence!=0| eta_opt$value%in%c(10^16, 10^17)){
    out <- matrix(c(NA,NA), ncol=2) 
  }else{
    suff_stat<- suff_stat_star_posi(X_obs=X_obs,eta= eta_opt$par ,theta=theta, w_i=w, 
                                    id.focus=id.focus,id.models=id.models,
                                    fix.seed=seed.fix, Var=Variance, full_X=full_X)
    if(suff_stat$mod_selected!=1){
      out <- matrix(c(NA,NA), ncol=2)
    }else{
      out <- matrix(c(suff_stat$suff_stat,
                      eta_opt$value), ncol=2)
    }
  }
  colnames(out)<- c("suff.stat", "opt_value")
  if(adapt_eta){
    out <- c(out, eta_opt$par)
  }
  return(out)
}

### Approximate post-selection confidence distribution:
## Input:
# focus_id: corresponds to id.focus in t_obs_star_posi
# models_id: corresponds to id.models in t_obs_star_posi
# X_obs, Y_obs as in t_obs_star_posi
# last_length_theta_star: How many values do we want for empirical distribution?
# V_known corresponds to variance in t_obs_star_posi (is variance known beforehand or not?)
# X_f: should full design matrix be used to generate sample? (False is correct)
# starting_eta: estimated values of the nuisance parameters of the selected model given the observed data, matrix with N_iter rows
# est_theta: estimated value of the focus parameter of the selected model given the observed data
# N_iter_int: how many iterations used to determine interval
# N_iter_CD: on how many iterations should posi confidence distribution be based.
# break_while: break any while loop after xx iterations; both for definition of interval and approximation of confidence distribution
## Output:
# theta_star: considered values of the focus parameter, i.e., values in set E^[theta]
# mean_theta_star_leq_theta: average number of iterations in which sampled sufficient statistics are larger than observed sufficient statistic for focus parameter
# max_dif_suff: maximum of opt_value in t_obs_star_posi for each considered value for theta
# iterations: number of iterations used to compute mean_theta_star_leq_theta for each theta
## Remark: 
# This function includes estimating the interval of theta for which the confidence distribution is approximated
# We use iterations because we stop approximation of confidence distribution per theta after break_while iterations even if N_iter_CD iterations has not been reached.
F_t_given_u_posi_2 <- function(X_obs, Y_obs, focus_id, models_id,
                               starting_eta, est_theta,
                               last_length_theta_star=50, X_f=F,
                               V_known=NULL,
                               adapt_start = F,
                               break_while=50000, 
                               N_iter_int=250, N_iter_CD = 500){
  N <- length(Y_obs)
  starting_eta <- matrix(rep(starting_eta, N_iter_int), byrow=T, nrow=N_iter_int)
  ## Compute the sufficient statistics:
  # If Variance is known, then we compute the sufficient statistics based on the whitenend data
  if(is.null(V_known)){
    suff_obs <- suff_stat_posi(Y=Y_obs, X=X_obs, id.focus=focus_id, id.selected = models_id[[1]], AR=T,Var=T)
    ar <- T
  }else{
    suff_obs <- suff_stat_posi(Y=Y_obs, X=X_obs, id.focus=focus_id, id.selected = models_id[[1]], AR=F,Var=F)
    ar <- F
  }
  ## Adapt the interval for possible values of theta to reduce computation time
  # Start with small interval and adapt interval after each iteration
  theta_star <- seq(from=round(est_theta,2)-5, to=round(est_theta,2)+5, by=1)
  # Then compute post-selection confidence distribution:
  n_non_na <- 0
  counter_1 <- 0
  while(counter_1 < break_while){
    n_theta <- length(theta_star)
    p_theta_star <- data.frame(theta_star=theta_star,
                               mean_theta_star_leq_theta= numeric(n_theta),
                               max_dif_suff= numeric(n_theta),
                               iterations= numeric(n_theta))
    for(i in 1:n_theta){
      t_mat <- apply(starting_eta, 1, t_obs_star_posi, X_obs=X_obs,Y_obs=Y_obs, suff_obs=suff_obs, theta=theta_star[i],id.focus=focus_id,
                     id.models= models_id, N=N, iter_m=iter_opt, Variance=V_known, seed.fix=NULL, AR=ar, full_X=X_f, adapt_eta=F)
      
      p_theta_star$mean_theta_star_leq_theta[i] <- mean(t_mat[1,]<= suff_obs$focus, na.rm=T) #if this is NA, then we do not consider it further!
      p_theta_star$max_dif_suff[i] <- max(t_mat[2,], na.rm=T)
      p_theta_star$iterations[i] <- N_iter_int-sum(is.na(t_mat[1,]))
    }
    # Check: which p_theta_star$mean_theta_leq_theta are not NA?
    theta_not_na <- p_theta_star$theta_star[!is.na(p_theta_star$mean_theta_star_leq_theta)]
    n_non_na <- length(theta_not_na)
    if(n_non_na!=0)break
    # if we did not find any values: try again with updated intervals:
    theta_star <- seq(from=min(theta_star)+5, to=max(theta_star)+5, by= 0.25)
    counter_1 <- counter_1+1
  }
  if(n_non_na==0){
    stop(c("No theta values found for which we can compute the probability"))
  }else{
    # use while loop until we have last_length_theta_star entries that are non-NA?
    p_theta <- p_theta_star
    counter_2 <- 0
    while(n_non_na < last_length_theta_star){
      ## update theta_star based on theta_not_na
      # update lower half and upper half of interval separately:
      if(min(theta_not_na)==min(theta_star)){
        # the interval of reasonable theta values is larger then theta_star: interval needs to become broader
        theta_seq_1 <- seq(min(theta_star)-6,min(theta_star)-1, by=1)
      }else{
        # theta_star contains values that are not reasonable; interval is large enough.
        theta_seq_1 <- numeric(0)
      }
      if(max(theta_not_na)==max(theta_star)){
        theta_seq_2 <- seq(max(theta_star)+1,max(theta_star)+6, by=1)
      }else{
        theta_seq_2 <- numeric(0)
      }
      theta_star <- c(theta_seq_1, theta_seq_2)
      # If the interval has reasonable size (i.e. we would not update any longer): Stop the while-loop
      n_theta <- length(theta_star)
      if(n_theta==0)break
      
      ## update p_theta_star
      p_theta_star <- data.frame(theta_star=theta_star,
                                 mean_theta_star_leq_theta= numeric(n_theta),
                                 max_dif_suff= numeric(n_theta),
                                 iterations= numeric(n_theta))
      for(i in 1:n_theta){
        t_mat <- apply(starting_eta, 1, t_obs_star_posi, X_obs=X_obs,Y_obs=Y_obs, suff_obs=suff_obs, theta=theta_star[i],id.focus=focus_id,
                       id.models= models_id, N=N, iter_m=iter_opt, Variance=V_known, seed.fix=NULL, AR=ar, full_X=X_f, adapt_eta=F)
        
        p_theta_star$mean_theta_star_leq_theta[i] <- mean(t_mat[1,]<= suff_obs$focus, na.rm=T) #if this is NA, then we do not consider it further!
        p_theta_star$max_dif_suff[i] <- max(t_mat[2,], na.rm=T)
        p_theta_star$iterations[i] <- N_iter_int-sum(is.na(t_mat[1,]))
      }
      # merge the data frames
      p_theta <- rbind(p_theta,p_theta_star)
      theta_not_na <- p_theta$theta_star[!is.na(p_theta$mean_theta_star_leq_theta)]
      # update n_non_na; stop the while loop once we have a rather large number of theta_star_not_na
      n_non_na <- length(theta_not_na)
      # stop the while loop after break_while iterations:
      if(counter_2 > break_while)break
      counter_2 <- counter_2+1
    }
    ### now: decide for which values we compute the full probabilities (based on N_iter iterations which fullfill all criteria)
    # 1) remove all NA values from p_theta
    p_theta <- p_theta[p_theta$theta_star%in%theta_not_na,]
    # 2) sort p_theta according to theta_star
    p_theta <- p_theta[order(p_theta$theta_star),]
    
    # 3) Remove values of theta_star if there are at least 10 values that are the same (would be probabilities of zero and/or ones.)
    # then keep 5 of these values only to reduce the number of 1 and 0 in the following interval
    id_dif <- which(p_theta$mean_theta_star_leq_theta[-n_non_na]!= p_theta$mean_theta_star_leq_theta[-1])
    if(min(id_dif)>=10){
      id_min <- id_dif[1]-5
    }else{
      id_min <- 1
    }
    if(max(id_dif)<= n_non_na-10){
      id_max <- max(id_dif)+5
    }else{
      id_max <- n_non_na
    }
    
    # 4) Initialize new interval of pre-defined length
    theta_star_full <- seq(from=p_theta$theta_star[id_min],to= p_theta$theta_star[id_max], length.out= last_length_theta_star)
    p_theta_full <- data.frame(theta_star=theta_star_full,
                               mean_theta_star_leq_theta= numeric(last_length_theta_star),
                               max_dif_suff= numeric(last_length_theta_star),
                               iterations= numeric(last_length_theta_star))
    # All overlapping values (will be at least the first and last value): use the already computed values.
    p_theta_full[p_theta_full$theta_star%in%p_theta$theta_star,] <- p_theta[p_theta$theta_star%in%theta_star_full,]
    
    # 5) Compute the probabilities until we have iterations= N_iter for each chosen set in theta_star in p_theta
    # This last while loop will take the longest
    for(i in 1:last_length_theta_star){
      niter <- p_theta_full$iterations[i]
      counter_3 <- 0
      while(niter < N_iter_CD){
        # break the loop after 10.000 iterations to avoid endlessly running while loops. 
        # think about how to deal with values that have rather small value of 
        if( counter_3 > break_while)break 
        t_star_one <- t_obs_star_posi(starting_eta[1,], X_obs=X_obs,Y_obs=Y_obs, suff_obs=suff_obs, theta=theta_star_full[i],id.focus=focus_id,
                                      id.models= models_id, N=N, iter_m=iter_opt, Variance=V_known, seed.fix=NULL, AR=ar, full_X=X_f, adapt_eta=F)
        # no seed in t_obs_star_posi when computing random values w_i
        
        # if value is not na: add to p_theta_full[,i]
        ## Check here which value of t_star_one could be na! it is a list after all
        if(!is.na(t_star_one[1])){
          if(t_star_one[1]<= suff_obs$focus){
            p_theta_full$mean_theta_star_leq_theta[i] <-  1/(niter+1)*(p_theta_full$mean_theta_star_leq_theta[i]*niter+1)
          }else{
            p_theta_full$mean_theta_star_leq_theta[i] <-  1/(niter+1)*(p_theta_full$mean_theta_star_leq_theta[i]*niter)
          }
          niter <- niter+1
          p_theta_full$iterations[i] <- niter
          p_theta_full$max_dif_suff[i] <- max(p_theta_full$max_dif_suff[i], t_star_one[2])
        }
        counter_3 <- counter_3 +1 
      }
    }
    return(p_theta_full)
  }
}


# Group level functions ---------------------------------------------------

####
# Estimation of eta and variance on group level
####

### REML for contrasts (of shape parameters): pre-specified change point locations
## Input:
# betas_sub: betas (or shape parameters) per contrast, so that betas belonging to one subject are next to each other!
# var_sub:list of the variance matrices of the betas within each subject, corresponding only to the variance/covariance matrix for the betas of interest!
# N: number of subjects
# maxit: maximum number of itertions in reml estimation
# retol: if within 10 iterations the difference is at most retol, we consider that we have a convergence.
# convergence curve: values of beta and variance for each iteration are returned
## Output:
# beta_group: eta
# betas_subject: contrasts per subject
# sigma_between: between subject variance of contrasts
# sigma_within: within subject variances of contrasts
# trace_beta: considered values for eta (only if converge_curve=T)
# trace_sigma: considered between subject variances (only if converge_curve=T)
# var_mat: variance on group level
## Remark: 
# it might happen that REML doesn't work. In this case, we use the EM algorithm
reml_var_estimation <- function(betas_sub, var_sub, contrast= c(1,-1), N, maxit=500, retol= 1*10^-8, converge_curve=F, full_covar=F){
  b_c <- contrast[1]*betas_sub[,1]+contrast[2]*betas_sub[,2]
  b_c <- matrix(b_c, ncol=1) #simply ensure that b_c is a vector with specific shape!
  # these are the betas that we put in the model, that is, the estimated betas per subject
  # we already use the contrast here
  ## compute the variance of the contrasts (i.e., subject specific variance)
  if(full_covar){
    var_contrast <- function(var_mat, contrast){
      out <- var_mat[1,1]+var_mat[2,2]+2*contrast[1]*contrast[2]*var_mat[1,2]
      return(out)
    }
  }else{
    var_contrast <- function(var_mat, contrast){
      out <- var_mat[1,1]+var_mat[2,2]
      return(out)
    }
  }
  sigmas_c <- sapply(X=var_sub, FUN = var_contrast, contrast=contrast)
  # assume that these variances are close to the truth.
  
  # compute eta and the variance on group level, using REML estimator.
  # this is iterative as both depend on each other.
  beta_reml <- sigma_reml <- rep(NA, maxit+1)
  
  # Define initial values for beta_reml and sigma_reml
  # use glm without accounting for within-subject variance
  beta_reml[1] <- mean(b_c)
  sigma_reml[1] <- 1/(N-1)*sum((b_c-mean(b_c))^2)
  # Use REML, or, if infeasible, EM approach:
  values_reml <- reml_try_catch(b_c=b_c,sigmas_c=sigmas_c, beta_reml=beta_reml, sigma_reml=sigma_reml, maxit=maxit,N=N, retol=retol)
  
  b_reml <- values_reml$beta_reml
  s_reml <- values_reml$sigma_reml
  id.na <- which(is.na(b_reml))
  
  if(length(id.na)==0){
    W <- diag(x=1/(s_reml[maxit+1]+sigmas_c))
    gamma_between <- s_reml[maxit+1]
    betas <- b_reml[maxit+1]
  }else{
    W <- diag(x=1/(s_reml[min(id.na)-1]+sigmas_c))
    gamma_between <- s_reml[min(id.na)-1]
    betas <- b_reml[min(id.na)-1]
  }
  if(converge_curve){
    out <- list(beta_group= betas,
                betas_subject = b_c,
                sigma_between= gamma_between,
                sigma_within= sigmas_c,
                trace_beta = b_reml[-id.na],
                trace_sigma = s_reml[-id.na],
                var_mat = W)
  }else{
    out <- list(beta_group= betas,
                betas_subject = b_c,
                sigma_between= gamma_between,
                sigma_within=sigmas_c,
                var_mat = W)
  }
  return(out)
  
}

### REML for parameter: post-selection inference
## Input:
# b_c: beta coefficients
# sigmas_c: vector of within subject post selection variances
# N: number of subjects
# maxit: maximum number of itertions in reml estimation
# retol: if within 10 iterations the difference is at most retol, we consider that we have a convergence.
# convergence curve: values of beta and variance for each iteration are returned
## Output:
# same as for "reml_var_estimation"
## Remark: 
# If REML is infeasible, we automatically use the EM algorithm
reml_var_estimation_2_posi<- function(b_c, sigmas_c, N, maxit=500, retol= 1*10^-8, converge_curve=F){
  beta_reml <- sigma_reml <- rep(NA, maxit+1)
  
  # Define initial values for beta_reml, sigma_reml
  # use glm without accounting for within-subject variance
  beta_reml[1] <- mean(b_c)
  sigma_reml[1] <- 1/(N-1)*sum((b_c-mean(b_c))^2)
  # use REML and, if not possible, EM approach:
  values_reml <- reml_try_catch(b_c=b_c,sigmas_c=sigmas_c, beta_reml=beta_reml, sigma_reml=sigma_reml, maxit=maxit,N=N, retol=retol)
  
  b_reml <- values_reml$beta_reml
  s_reml <- values_reml$sigma_reml
  id.na <- which(is.na(b_reml))
  
  if(length(id.na)==0){
    W <- diag(x=1/(s_reml[maxit+1]+sigmas_c))
    gamma_between <- s_reml[maxit+1]
    betas <- b_reml[maxit+1]
  }else{
    W <- diag(x=1/(s_reml[min(id.na)-1]+sigmas_c))
    gamma_between <- s_reml[min(id.na)-1]
    betas <- b_reml[min(id.na)-1]
  }
  if(converge_curve){
    out <- list(beta_group= betas,
                betas_subject = b_c,
                sigma_between= gamma_between,
                sigma_within= sigmas_c,
                trace_beta = b_reml[-id.na],
                trace_sigma = s_reml[-id.na],
                var_mat = W)
  }else{
    out <- list(beta_group= betas,
                betas_subject = b_c,
                sigma_between= gamma_between,
                sigma_within=sigmas_c,
                var_mat = W)
  }
  return(out)
  
}

## Helper function
reml_try_catch <- function(b_c,sigmas_c, beta_reml, sigma_reml, maxit,N, retol){
  tryCatch({
    stopper <- 0
    counter <- 0
    X <- matrix(data = rep(1,N), ncol=1)
    while(stopper==0 & counter <= maxit){
      counter <- counter +1
      # compute the variance based on beta_reml[counter]
      
      sigma_reml[counter+1] <- 1/(N-1)*sum((b_c-beta_reml[counter])^2)
      # check if this leads to correct vector, need to have vector with squared differences!
      
      # compute beta based on sigma_reml[counter+1]
      V <- diag(sigmas_c)+sigma_reml[counter+1]
      # then (in accordance with both Mumford & Nichols and Fahrmeir et al:)
      XV_1 <- t(X)%*%solve(V)
      beta_reml[counter+1] <- solve(XV_1%*%X)%*%XV_1%*%b_c
      
      if(counter>=9){
        if(mean(beta_reml[(counter-8):(counter+1)])<=retol & 
           mean(sigma_reml[(counter-8):(counter+1)])<=retol){
          stopper <- 1
        }
      }
    }
    values_reml <- list(beta_reml=beta_reml,
                        sigma_reml=sigma_reml)
    return(values_reml)
  },
  error=function(e2){
    stopper <- 0
    counter <- 0
    while(stopper==0 & counter <= maxit){
      counter <- counter +1
      # E-step
      s_tilde <- (1/sigma_reml[counter]+1/sigmas_c)^-1
      B_tilde <- s_tilde*(b_c/sigma_reml[counter]+beta_reml[counter]/sigmas_c)
      
      # M-step
      beta_reml[counter+1] <- mean(B_tilde)
      sigma_reml[counter+1] <- 1/N*(sum(s_tilde)+sum((B_tilde-beta_reml[counter+1])^2))
      
      if(counter>=9){
        if(mean(beta_reml[(counter-8):(counter+1)])<=retol & 
           mean(sigma_reml[(counter-8):(counter+1)])<=retol){
          stopper <- 1
        }
      }
    }
    values_reml <- list(beta_reml=beta_reml,
                        sigma_reml=sigma_reml)
    return(values_reml)
  })
}

####
# Test statistics on group level
####

### T_KH
## Input:
# est_eta: that is beta_group output (the estimated effect on group level)
# var_mat_full: that is var_mat output of reml
#               alternatively: use diag(1/(sigma_between+sigma_within))
# betas_c: observed effects under consideration (i.e., in our case the contrasts)
#          these are betas_subject for reml and em algorithm!
# N: number of subjects
## Output:
# test_statistic: T_KH
# p_value: p-value corresponding to T_KH
# observations: betas_c,
# test_value: est_eta,
# variance_test_value: updated variance used to compute T_KH
## Remark: 
# null hypothesis: eta = 0 
fmri_group_KH <- function(est_eta, var_mat_full, betas_c, N){
  X <- matrix(rep(1, N), ncol=1)
  XWX <- t(X)%*%var_mat_full%*%X
  P <- var_mat_full-var_mat_full%*%X%*%solve(XWX)%*%t(X)%*%var_mat_full
  bPb <-t(betas_c)%*%P%*%betas_c
  var_a <- 1/(N-1)*bPb*solve(XWX)
  T_KH <- est_eta/sqrt(var_a)
  p_val <- 2*(pt(-abs(T_KH), df=N-1))
  out <- list(test_statistic= T_KH,
              p_value = p_val,
              observations = betas_c,
              test_value=est_eta,
              variance_test_value= var_a)
}

### T_WALD
## Input:
# est_eta: that is beta_group output (the estimated effect on group level)
# var_mat_full: that is var_mat output of reml, i.e., diag(1/(sigma_between+sigma_within))
# N: sample size on group level
# null_eta: eta of the null hypotheses
## Output:
# p-value corresponding to the null hypotheses: eta=null_eta 
fmri_test_own <- function(est_eta, var_mat_full, N, null_eta=0){
  X <- matrix(rep(1, N), ncol=1)
  XWX <- t(X)%*%var_mat_full%*%X
  test_statistic <- (est_eta-null_eta)/(sqrt(solve(XWX))
  p_val <- 2*(1-pt(abs(test_statistic), df=N-1))
  return(p_val)
}

###
# Group level analysis based on output of est_beta_cd
###

## Input:
# output_posi: list of size N, each list element corresponds to the output of est_beta_cd for one participant 
# N: size of the group
group_analysis_posi <- function(output_posi, N){
  # initialize vectors used in further analysis:
  hat_dif_OLS <-hat_var_naive<-hat_dif_posi<-hat_E_posi<-hat_var_posi <-est_cp<- numeric(N)
  for(i in 1:N){
    hat_dif_OLS[i]<- output_posi[[i]]$hat_dif
    hat_var_naive[i]<- output_posi[[i]]$hat_var_beta
    hat_dif_posi[i]<- output_posi[[i]]$hat_dif_cd
    hat_E_posi[i]<- output_posi[[i]]$hat_E_cd
    hat_var_posi[i]<- output_posi[[i]]$hat_var_cd
  }
  ## Estimate eta and its variance on group level based on different subject level input:
  # Several combinations of estimators (location and scale estimators)
  group_estimation_naive <- reml_var_estimation_2_posi(b_c=hat_dif_OLS, sigmas_c= hat_var_naive, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi <- reml_var_estimation_2_posi(b_c=hat_dif_posi, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi_E <- reml_var_estimation_2_posi(b_c=hat_E_posi, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  group_estimation_posi_OLS <- reml_var_estimation_2_posi(b_c=hat_dif_OLS, sigmas_c= hat_var_posi, N=N, maxit=500, retol= 1*10^-8, converge_lim=10, converge_curve=F)
  # save results in vectors for further analysis:
  est_eta <- c(group_estimation_naive$beta_group,group_estimation_posi$beta_group, group_estimation_posi_E$beta_group, group_estimation_posi_OLS$beta_group)
  between_subject_variance <- c(group_estimation_naive$sigma_between,group_estimation_posi$sigma_between,group_estimation_posi_E$sigma_between, group_estimation_posi_OLS$sigma_between)
  var_group_level <- list(naive= group_estimation_naive$var_mat, 
                          posi_05=group_estimation_posi$var_mat,
                          posi_E=group_estimation_posi_E$var_mat,
                          OLS_posi=group_estimation_posi_OLS$var_mat)
  ## Compute p-values corresponding to H_0: no change in PM
  p_vals_KH <- p_vals_W <- numeric(4)
  p_vals_KH[1] <- fmri_group_KH(est_eta=est_eta[1], 
                                var_mat_full=var_group_level[[1]],
                                betas_c=hat_dif_OLS, N=N)$p_value
  p_vals_KH[2] <- fmri_group_KH(est_eta=est_eta[2], 
                                var_mat_full=var_group_level[[2]],
                                betas_c=hat_dif_posi, N=N)$p_value
  p_vals_KH[3] <- fmri_group_KH(est_eta=est_eta[3], 
                                var_mat_full=var_group_level[[3]],
                                betas_c=hat_E_posi, N=N)$p_value
  p_vals_KH[4] <- fmri_group_KH(est_eta=est_eta[4], 
                                var_mat_full=var_group_level[[4]],
                                betas_c=hat_dif_OLS, N=N)$p_value
  p_vals_W[1] <- fmri_test_own(est_eta=est_eta[1], 
                               var_mat_full=var_group_level[[1]],
                               N=N)
  p_vals_W[2] <- fmri_test_own(est_eta=est_eta[2], 
                               var_mat_full=var_group_level[[2]],
                               N=N)
  p_vals_W[3] <- fmri_test_own(est_eta=est_eta[3], 
                               var_mat_full=var_group_level[[3]],
                               N=N)
  p_vals_W[4] <- fmri_test_own(est_eta=est_eta[4], 
                               var_mat_full=var_group_level[[4]],
                               N=N)
  ## Name all vectors containing results:
  names(est_eta) <-names(between_subject_variance)<-names(p_vals_KH)<- names(p_vals_W)<- c("naive", "posi_05", "posi_E", "OLS_posi")
  
  out <- list(p_values_KH = p_vals_KH,
              p_values_Wald= p_vals_W,
              est_eta = est_eta,
              between_subject_variance= between_subject_variance,
              variance_group_level = var_group_level)
  return(out)
}

# Helper functions --------------------------------------------------------
### Computation of the shape parameters: PM, NA, IUA, TTP, TPN, FWHM, FWHN, AUC, AC, DC
shape_parameters2 <- function(hrf_sampled, baseline=0, tres=0.05, max.ttp=12, auc_ac_dp_iua=T){
  peak_amplitude <- max(hrf_sampled)-baseline
  ttp <- which.max(hrf_sampled)
  if(length(peak_amplitude)==0 | is.na(peak_amplitude)){
    peak_amplitude <- -1
    ttp <- 1
  }
  if(length(ttp)==0)ttp <- 1
  if(ttp > max.ttp/tres){
    ind_os <- max(which(hrf_sampled[1:ttp]<=0))
    peak_amplitude <- max(hrf_sampled[1:ind_os])-baseline
    if(length(peak_amplitude)==0){
      peak_amplitude <- -1
      ttp <- 1
    }else{ttp <- min(which(hrf_sampled==max(hrf_sampled[1:ind_os])))}
  }
  if(peak_amplitude<=0){
    slope_est <- sign(hrf_sampled[-1]-hrf_sampled[-length(hrf_sampled)])
    peak_help <- slope_est[-1]*slope_est[-length(slope_est)]
    ind.peaks <- which(peak_help!=1)+1
    if(length(ind.peaks) !=0){
      ind.p <- which(ind.peaks < max.ttp/tres)
      peak_amplitude <- max(hrf_sampled[ind.peaks[ind.p]])-baseline
      ttp <- which(hrf_sampled==max(hrf_sampled[ind.peaks[ind.p]]))}
    else{
      peak_amplitude <- min(hrf_sampled[-c(1:ttp)])+baseline
      ttp <- 1
    }
  }
  if(peak_amplitude==min(hrf_sampled)){
    ttp <- 1
  }
  if(ttp != 1){
    time_to_peak <- ttp*tres
  }else{time_to_peak <- 0}
  
  nadir_amplitude <- min(hrf_sampled[-c(1:ttp)])-baseline
  time_peak_to_nadir <- (which.min(hrf_sampled[-c(1:ttp)]))*tres
  if(auc_ac_dp_iua){
    initial_undershoot_amplitude <- min(hrf_sampled[c(1:ttp)])+baseline
  }
  # Compute full width half maximum
  half_max <- peak_amplitude/2
  exact_half <- 0
  if(length(which(hrf_sampled==half_max))==0){
    ind_hm_1 <- max(which(hrf_sampled[1:ttp]<= half_max))
    ind_hm_2 <- min(which(hrf_sampled[-c(1:ttp)]<= half_max))+ttp
    ind_1 <- (ind_hm_1*2+1)/2
    ind_2 <- (ind_hm_2*2-1)/2
  }else{
    if(length(which(hrf_sampled[1:ttp]==half_max))!=0){
      exact_half <- 1
      ind_1 <- ind_hm_1 <-  max(which(hrf_sampled[1:ttp]==half_max))
      if(length(which(hrf_sampled[-c(1:ttp)]==half_max))!=0){
        exact_half <- 3
        ind_2 <- ind_hm_2 <- max(which(hrf_sampled[-c(1:ttp)]==half_max))+ttp
      }else{
        ind_2 <- ind_hm_2 <- ((min(which(hrf_sampled[-c(1:ttp)]<= half_max))+ttp)*2-1)/2}
    }else{
      exact_half <- 2
      ind_2 <- ind_hm_2 <- max(which(hrf_sampled[-c(1:ttp)]==half_max))+ttp
      ind_1 <- ind_hm_1 <-  (max(which(hrf_sampled[1:ttp]<= half_max))*2+1)/2
    }
  }
  if(peak_amplitude < 0 && ttp != 1){
    id_peak <- which(ind.peaks %in% which(hrf_sampled == peak_amplitude))
    # find out when rise starts, use this as the "height" on which we use half of it!
    rise_s <- ind.peaks[id_peak-1]
    rise_e <- ind.peaks[id_peak+1]
    # define new helper baseline
    bl_help <- max(hrf_sampled[rise_s], hrf_sampled[rise_e])
    height_peak <- peak_amplitude-bl_help
    half_max <- bl_help + height_peak/2
    ind_hm_1 <- max(which(hrf_sampled[1:ttp]<= half_max))
    ind_hm_2 <- min(which(hrf_sampled[-c(1:ttp)]<= half_max))+ttp
    ind_1 <- (ind_hm_1*2+1)/2
    ind_2 <- (ind_hm_2*2-1)/2
  }
  FWHM <- (ind_2-ind_1)*tres
  if(FWHM==Inf)FWHM <- 0
  if(FWHM==-Inf)FWHM <- 0
  
  if(auc_ac_dp_iua){
    if(FWHM > 0){
      if(exact_half ==0){
        ascent_peak <- hrf_sampled[ind_hm_1+1]-hrf_sampled[ind_hm_1]
        descent_peak <- hrf_sampled[ind_hm_2]-hrf_sampled[ind_hm_2-1]
      }else if(exact_half == 1){
        ascent_peak <- (hrf_sampled[ind_hm_1+1]-hrf_sampled[ind_hm_1-1])/2
        descent_peak <- hrf_sampled[ind_hm_2]-hrf_sampled[ind_hm_2-1]
      }else if(exact_half==2){
        ascent_peak <- hrf_sampled[ind_hm_1+1]-hrf_sampled[ind_hm_1]
        descent_peak <- (hrf_sampled[ind_hm_2+1]-hrf_sampled[ind_hm_2-1])/2
      }else{
        ascent_peak <- (hrf_sampled[ind_hm_1+1]-hrf_sampled[ind_hm_1-1])/2
        descent_peak <- (hrf_sampled[ind_hm_2+1]-hrf_sampled[ind_hm_2-1])/2
      }
    }else{
      ascent_peak <- descent_peak <- 0
    }
  }
  
  
  # Full width half nadir
  tpn <- time_peak_to_nadir/tres
  half_nadir <- nadir_amplitude/2
  hrf_after_peak <- hrf_sampled[-c(1:ttp)]
  if(length(which(hrf_after_peak==half_nadir))==0){
    ind_hm_1_n <- min(which(hrf_after_peak[1:tpn]<= half_nadir))
    ind_hm_2_n <- max(which(hrf_after_peak[-c(1:tpn)]<= half_nadir))+tpn
    ind_1_n <- (ind_hm_1_n*2-1)/2
    ind_2_n <- (ind_hm_2_n*2+1)/2
  }else{
    if(length(which(hrf_after_peak[1:tpn]==half_nadir))!=0){
      ind_1_n <- ind_hm_1_n <- min(which(hrf_after_peak[1:tpn]==half_nadir))
      if(length(which(hrf_after_peak[-c(1:tpn)]==half_nadir))!=0){
        ind_2_n <- ind_hm_2_n <-max(which(hrf_after_peak[-c(1:tpn)]==half_nadir))+tpn
      }else{ind_2_n <- ind_hm_2_n <-((max(which(hrf_sampled[-c(1:tpn)]<= half_nadir))+tpn)*2+1)/2}
    }else{
      ind_2_n <-ind_hm_2_n <- max(which(hrf_after_peak[-c(1:tpn)]==half_nadir))+tpn
      ind_1_n <- ind_hm_1_n <- (min(which(hrf_after_peak[1:tpn]<= half_nadir))*2-1)/2
    }
  }
  if(peak_amplitude < 0 && ttp != 1){
    if(!is.na(bl_help) & !is.na(nadir_amplitude)){
      if(bl_help == nadir_amplitude){
        half_nadir <- nadir_amplitude + abs((nadir_amplitude+peak_amplitude)/2)
      }
    }else{
      half_nadir <- nadir_amplitude + abs((nadir_amplitude-bl_help)/2)}
    ind_hm_1_n <- min(which(hrf_after_peak[1:tpn]<= half_nadir))
    ind_hm_2_n <- max(which(hrf_after_peak[-c(1:tpn)]<= half_nadir))+tpn
    ind_1_n <- (ind_hm_1_n*2-1)/2
    ind_2_n <- (ind_hm_2_n*2+1)/2
  }
  FWHN <- (ind_2_n-ind_1_n)*tres
  if(FWHN==Inf)FWHM <- 0
  if(FWHN==-Inf)FWHM <- 0
  
  # Area under the curve
  # find the value for which the AUC starts(i.e. values larger than the baseline after nadir)
  # compute the integral from this value forward
  if(auc_ac_dp_iua){
    AUC_s <- min(which(hrf_sampled[-c(1:(tpn+ttp))]>= baseline))
    if(AUC_s!= Inf){
      AUC_s <- AUC_s+(tpn+ttp)
      AUC <- sum(hrf_sampled[-c(1:(AUC_s-1))])
    }else{AUC <- 0}
  }
  # account for errors:
  if(auc_ac_dp_iua){
    list_help <- list(peak_amplitude,
                      nadir_amplitude,
                      initial_undershoot_amplitude,
                      time_to_peak,
                      time_peak_to_nadir,
                      FWHM,
                      FWHN,
                      AUC,
                      ascent_peak,
                      descent_peak)
    id.error <- which(lapply(list_help, length)!=1)
    if(length(id.error)==0){
      out <- c(peak_amplitude, abs(nadir_amplitude), abs(initial_undershoot_amplitude), time_to_peak, time_peak_to_nadir, FWHM, FWHN, AUC, ascent_peak, descent_peak)
    }else{
      id.correct <- c(1:10)[-id.error]
      out <- rep(0,10)
      for(i in id.correct){
        out[i] <- list_help[[i]]
      }
    }
    names(out) <- c("PM", "NA", "IUA", "TTP", "TPN", "FWHM", "FWHN", "AUC", "AC_P", "DC_P")
  }else{
    list_help <- list(peak_amplitude,
                      nadir_amplitude,
                      time_to_peak,
                      time_peak_to_nadir,
                      FWHM,
                      FWHN)
    id.error <- which(lapply(list_help, length)!=1)
    if(length(id.error)==0){
      out <- c(peak_amplitude, abs(nadir_amplitude), time_to_peak, time_peak_to_nadir, FWHM, FWHN)
    }else{
      id.correct <- c(1:6)[-id.error]
      out <- rep(0,6)
      for(i in id.correct){
        out[i] <- list_help[[i]]
      }
    }
    names(out) <- c("PM", "NA", "TTP", "TPN", "FWHM", "FWHN")
  }
  id.inf <- which(out==Inf)
  id.neg.inf <- which(out==-Inf)
  if(length(id.inf)!=0){
    out[id.inf] <- 0
  }
  if(length(id.neg.inf)!=0){
    out[id.neg.inf] <- 0
  }
  return(out)
}

### Computation of the realized false discovery proportion
fdp <- function(vec_reject, id.true.null){
  n_reject <- sum(vec_reject)
  if(n_reject!=0){
    return(sum(vec_reject[id.true.null])/sum(vec_reject))
  }else{
    return(0)
  }
}

### Computation of the post-selection variance with a time limit.
## Remark:
# The time limit ensures that functions are not stuck in never ending while loop
# The time limit is not universal for all scenarios and should be tested beforehand.
# This function returns NULL if time limit is reached.
func_parallel_cd_tl <- function(BOLD_signal, sigma_w, rho_noise, onsets, nscan, nstim, min_seq, last.length.theta.star,
                                true_cp_loc, Var_known=T, time_limit=2.2*60*60){
  if(Var_known){
    Var_ar <- matrix(ncol=nscan, nrow=nscan)
    for(k in 1:nscan){
      for(r in 1:nscan){
        Var_ar[k,r] <- rho_noise^(abs(k-r))
      }
    }
    Var_sub <- list(sigma= sigma_w,
                    V = Var_ar)
  }else{
    Var_sub <- NULL
  }
  # now add the time limit: if it takes longer than the time limit: return NULL 
  # if it takes longer than e.g. 2 hours, it typically doesn't go back to normal. 
  # This I can use in the foreach loop and hopefully I get at least some answers
  # worst case: no answers at all.
  setTimeLimit(elapsed = time_limit)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  out <- tryCatch({
    est_betas_cd(BOLD_signal=BOLD_signal, 
                 onsets= onsets, 
                 nscan=nscan, 
                 nstim= nstim,
                 min_seg= min_seq,
                 last.length.theta.star=last.length.theta.star, 
                 true_cp = true_cp_loc,
                 Var_known = Var_sub)
  },error=function(e) NULL)
  return(out)
}

### Function that returns the estimated variance of eta for post-selection estimation and naive estimator
est_var <- function(output_posi){
  return(c(output_posi$hat_var_beta, output_posi$hat_var_cd))
}
