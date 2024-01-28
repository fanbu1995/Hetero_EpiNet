# 03/09/2021

# overall estimation function for iEpi data
# (with external infection sources)

outdir = './iepi_results/'

# server stuff, set seed and example data path
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
s0 = as.numeric(slurm_arrayid)
set.seed(s0)
ss = sample(1000,1)

# source the inference util funcs
source("inference_utils_2.R")

# load data from previously saved object??
dat = ....

# set up initial values
nX = 3 #? num. of covariates to try

IV = list(beta = 0.07, Eta = 1.2, gamma = 0.3, 
          b_S = rep(0,nX), phi = 0.2)

# function for inference
# 08/16/2021: tract rejection sampler accept ratio
infer_partial_data <- function(dat, initValues,
                               tmax = 7, tmin = 0,
                               numIter=100, burn=0, maxIter=20, tol=1e-6,
                               seed=42, verbose=TRUE){
  
  # dat: list that contains data files
  #     including: events, G0, I0, X, report.times, report, stage_change, ext_labels
  # initValues: list of initial values
  # numIter: num. of iterations for the stochastic EM
  # burn: burn-in perioid for the stochastic EM
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  
  set.seed(seed)
  
  # load data and extract stuff
  events = dat$events
  G0 = dat$G0
  I0 = dat$I0
  X = dat$X
  stage_change = dat$stage_change
  ext_labels = dat$ext_labels
  report.times = dat$report.times
  report = dat$report
  
  if(verbose){cat('Data loaded...')}
  
  # get an adjmat report
  G_all = get_adjmat_report(G0, events, report.times)
  if(verbose){cat('Adjmat report created.\n')}
  
  # get intervals to impute recovery times on and those people's ids
  MR = get_miss_recov(report, report.times, events)
  
  # get list of manifested people and their times
  ## not include I0
  ## 03/09/2021: exclude "ext_labels" as well -> they don't need imputation
  manifest = events %>% filter(event %in% c(9,10)) %>%
    filter(per1 != I0 & !per1 %in% ext_labels) %>%
    select(manifested = per1, times = time) %>%
    as.list()
  
  
  # get an initial conservative imputation of recovery times 
  # (everyone recovers at the end of interval)
  recovery_times = NULL
  for(r in 1:nrow(MR$intervals)){
    recov_r = data.frame(recov = MR$recover[[r]], time = MR$intervals$ub[r])
    recovery_times = rbind(recovery_times, recov_r)
  }
  
  # start iterating...
  
  ## current values of exp_eta and gam
  exp_eta = initValues$Eta
  gam = initValues$gamma
  beta = initValues$beta
  b_S = initValues$b_S
  phi = initValues$phi
  
  ## storage for the samples
  varnames = c("beta", "eta", "exp_eta", "gamma", "phi",
               "p_s", "alpha", "omega", "b_S", "b_alpha", "b_omega")
  params = list()
  nsamps = numIter - burn
  for(v in varnames){
    if(length(grep('b_',v)) > 0){
      # if the parameter is one of b_S, b_E, b_omega, b_alpha
      # then storage should be matrix that relates to X
      params[[v]] = matrix(0, nrow=nsamps, ncol=ncol(X))
    }else if(length(grep('(omega)|(alpha)',v)) > 0){
      # if param is omega or alpha
      # then length is 6 for each param
      params[[v]] = matrix(0, nrow=nsamps, ncol=6)
    }else{
      # else: storage is a vector
      params[[v]] = numeric(nsamps)
    }
  }
  
  for(s in 1:numIter){
    if(verbose){cat('Iteration',s,'...')}
    
    # 1) impute exposure times first
    imp_expo_times = get_expo_times(manifest, G_all, tmax, tmin, 
                                    report, report.times, 
                                    events, recovery_times, X, b_S,
                                    exp_eta, beta, phi)

    
    if(verbose){cat('Exposure times imputation done! ')}
    
    # 2) impute recovery times
    ## exposure time list
    exposure_times = list(exposed=manifest$manifested, times = imp_expo_times)
    ## construct local neighborhoods
    infec_nei = get_nei_expo_all(G_all, events, exposure_times, 
                                 report, report.times)
    
    ### run a sanity check: make sure everybody gets potential infectors
    for(p in names(infec_nei)){
      if(length(unlist(infec_nei[[p]])) == 0){
        cat('Person', p, 'does not have any infectors available!!!\n')
      }
    }
    
    ## sample recovery times by interval
    # recov_times = foreach(r=1:nrow(MR$intervals), .combine = 'c') %dopar% {
    #   lb = MR$intervals$lb[r]; ub = MR$intervals$ub[r]
    #   recovs = MR$recover[[r]]
    #   propose_recov_filter(lb, ub, recovs, exposure_times, infec_nei, gam, exp_eta)
    # }
    
    # sample recovery times
    recov_times = NULL
    for(r in 1:nrow(MR$intervals)){
      lb = MR$intervals$lb[r]; ub = MR$intervals$ub[r]
      recovs = MR$recover[[r]]
      cands = propose_recov_filter(lb, ub, recovs, exposure_times, infec_nei, gam, exp_eta)
      recov_times = c(recov_times, cands)
    }
    recovery_times$time = recov_times
  
    if(verbose){cat('Recovery times imputation done!\n')}
    
    # 3) combine samples with data to get augmented events and summarize those
    events_aug = combine_data(events, exposure_times, recovery_times)
    # 03/09/2021: add ext_labels here
    summs = summarize_events2(G0, I0, events_aug, stage_change, ext_labels)
    
    if(verbose){cat('Summarizing events done!')}
    
    # 4) get MLEs
    ## 11/06/2020 update with hacked eta and phi included
    estimates = solve_MLE(summs, X, maxIter, tol, initValues$Eta)
    ## adjust "eta" names
    estimates[['exp_eta']] = estimates[['eta']]
    estimates[['eta']] = ifelse(estimates[['eta']] <= 0, -Inf, log(estimates[['eta']]))
    
    if(verbose){cat('MLE obtained!')}
    
    ## update current values of exp(eta) and gamma
    exp_eta = estimates[['exp_eta']]
    gam = estimates[['gamma']]
    ## 11/08/2020: also get beta and b_S
    beta = estimates[['beta']]
    b_S = estimates[['b_S']]
    
    ## record it after burn-in
    if(s > burn){
      for(v in varnames){
        if(length(grep('(omega)|(alpha)|(_S)|(_E)',v)) > 0){
          # if the parameter is one of b_S, b_E, b_omega, b_alpha, omega, alpha
          params[[v]][s-burn,] = estimates[[v]]
        }else{
          # else
          params[[v]][s-burn] = estimates[[v]]
        }
      }
    }
    
    if(verbose){
      # output some key parameter values as well
      cat('Key estimates:\n beta=', estimates$beta, 
          'xi=', estimates$xi,
          'phi=', estimates$phi,
          'gamma=', estimates$gamma,
          'exp_eta=', estimates$exp_eta)
      cat('\nEstimates updated and saved.\n\n')
    }
  }
  
  # finally return results
  ## return mean values of parameter draws as well
  means = list()
  for(v in varnames){
    if(length(grep('(omega)|(alpha)|(_S)|(_E)',v)) > 0){
      # if the parameter is one of b_S, b_omega, b_alpha, omega, alpha
      means[[v]] = apply(params[[v]], 2, mean, na.rm=TRUE)
    }else{
      # else
      means[[v]] = mean(params[[v]], na.rm=TRUE)
    }
  }
  
  list(params = params, means=means)
}


# run and save results
res = infer_partial_data(dat, IV, tmax=14, tmin=0, seed = ss)

saveRDS(res, paste0(outdir,slurm_arrayid,'.rds'))
