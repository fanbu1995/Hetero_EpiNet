# 06/28/2021
# estimate eta and b_S only for the "both missing case"
# while keeping the other things fixed

# first check the b_S truth for a bunch of datasets
# also keep track of the N's
# commented out for server running

# setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
# data_root = '~/Documents/Research_and_References/Hetero_EpiNet_2020/'
# source('inference_utils_3.R')
# 
# ds_num = c(0:124)
# for(ds in ds_num){
#   fpath = paste0('ex',ds)
#   dat = load_data(data_root, fpath)
#   N = nrow(dat$X)
#   b_S_tru = dat$truth$b_S
#   cat(fpath,': N=', N, 'b_S = (', b_S_tru[1], b_S_tru[2], ').\n')
# }

# N = 100, for 1,3,5,7,...29 
# N = 200, for 0,2,4,...30, 31:40, 101:124
# N = 300, for 71:100
# N = 500, for 41:70

# b_S = (0,-1) for 0:10
# b_S = (0, 0) for 11:20
# b_S = (0, 1) for 21:30
# b_S = (1, 1) for 31:124


## SERVER STUFF (start running from here!)

# data directory and outdir
data_root = 'hetero_data/'

# change outdir for the eta and b_S experiments
outdir = 'eta_bS/'


# set seed and example data path
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
s0 = as.numeric(slurm_arrayid)
set.seed(s0)
ss = sample(1000,1)

# focus on one dataset only for each job...
ds = (s0-1) %/% 2
sub_dir = paste0('ex',ds)

# source the inference util funcs
source("inference_utils_3.R")


# 1. function to only solve for eta
# (adapted from solve_MLE)
solve_eta <- function(summaries, X, truth, initEta = 1, lowerBound = 1e-6){
  # initEta: initial value of exp(eta)
  # lowerBound: the lower bound used in "optim"
  
  N = nrow(X)
  pairs = which(upper.tri(matrix(nrow=N, ncol=N)), arr.ind = T)
  XX = apply(pairs, 1, sum_covariates, X)
  XX = t(XX)
  
  # extract summary statistics
  nI = summaries$counts['nE'] # take nE = nI
  nIs = summaries$counts['nIs']
  nIa = summaries$counts['nIa']
  nR = summaries$counts['nR']
  C_all = summaries$C_all
  D_all = summaries$D_all
  epi_table = summaries$epi_table
  net_c_table = summaries$net_c_table
  net_d_table = summaries$net_d_table
  I0 = summaries$I0
  
  # set all parameters to truth
  # (disregard all the network params, since they are not used here)
  beta = truth$beta
  gamma = truth$gamma
  phi = truth$phi
  p_s = truth$p_s
  b_S = truth$b_S
  
  
  # get quantities used for estimating eta
  eta = initEta
  
  # 10/24/20 fix: only get those who had infec neighbors at infec time!
  got_infected = which((epi_table$Ia_ti + epi_table$Is_ti) > 0 & epi_table$latent_time > 0)
  got_infected = got_infected[got_infected != I0] # exclude I0 here
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  
  Expo = (epi_table$Ia_expo + eta * epi_table$Is_expo)
  has_expo = which(Expo > 0) # check I0 not here!
  # only include those who had >0 exposure to infectious sources
  Y_infec = Y_infec[has_expo]
  X_infec = X[has_expo,]
  
  indiv_effects = exp(c(X_infec %*% b_S))
  
  # estimate eta
  # (optim: minimization, so have to negate things here)
  llfunc = function(eta){
    beta * sum(indiv_effects * (epi_table$Ia_expo + eta * epi_table$Is_expo)[has_expo]) -
      sum(log((epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected]))
  }
  ## 10/18/2020: de-bugged - previously this was not correct
  llgrr = function(eta){
    beta * sum(indiv_effects * epi_table$Is_expo[has_expo]) -
      sum(epi_table$Is_ti[got_infected]/(epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected])
  }
  
  # optim
  # the "strict" version: not fixing exp(eta) to > 1
  eta = optim(initEta, llfunc, llgrr,
              method = "L-BFGS-B", 
              lower = lowerBound, upper = Inf)$par
  
  # return the solved eta
  eta
}

# 2. function to only solve for eta and b_S
# (adapted from solve_MLE and solve_eta)
solve_eta_bS <- function(summaries, X, truth, initEta = 1, 
                         maxIter=10, lowerBound = 1e-6, tol = 1e-4){
  # initEta: initial value of exp(eta)
  # lowerBound: the lower bound used in "optim"
  
  N = nrow(X)
  pairs = which(upper.tri(matrix(nrow=N, ncol=N)), arr.ind = T)
  XX = apply(pairs, 1, sum_covariates, X)
  XX = t(XX)
  
  # extract summary statistics
  nI = summaries$counts['nE'] # take nE = nI
  nIs = summaries$counts['nIs']
  nIa = summaries$counts['nIa']
  nR = summaries$counts['nR']
  C_all = summaries$C_all
  D_all = summaries$D_all
  epi_table = summaries$epi_table
  net_c_table = summaries$net_c_table
  net_d_table = summaries$net_d_table
  I0 = summaries$I0
  
  # set all parameters to truth and set init b_S to truth as well
  # (disregard all the network params, since they are not used here)
  beta = truth$beta
  gamma = truth$gamma
  phi = truth$phi
  p_s = truth$p_s
  
  b_S = truth$b_S
  
  
  # get quantities used for estimating eta
  eta = initEta
  
  # 10/24/20 fix: only get those who had infec neighbors at infec time!
  got_infected = which((epi_table$Ia_ti + epi_table$Is_ti) > 0 & epi_table$latent_time > 0)
  got_infected = got_infected[got_infected != I0] # exclude I0 here
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  
  for(it in 1:maxIter){
    ## (1) estimate b_S
    Expo = (epi_table$Ia_expo + eta * epi_table$Is_expo) * beta
    has_expo = which(Expo > 0) # check I0 not here!
    # only include those who had >0 exposure to infectious sources
    Y_infec = Y_infec[has_expo]
    X_infec = X[has_expo,]
    Expo = Expo[has_expo]
    infec_poi = glm(Y_infec~X_infec-1+offset(log(Expo)), family = poisson())
    
    b_S_old = b_S
    b_S = infec_poi$coefficients %>% as.numeric()
    
    
    ## (2) estimate eta
    indiv_effects = exp(c(X_infec %*% b_S))
    
    # (optim: minimization, so have to negate things here)
    llfunc = function(eta){
      beta * sum(indiv_effects * (epi_table$Ia_expo + eta * epi_table$Is_expo)[has_expo]) -
        sum(log((epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected]))
    }
    ## 10/18/2020: de-bugged - previously this was not correct
    llgrr = function(eta){
      beta * sum(indiv_effects * epi_table$Is_expo[has_expo]) -
        sum(epi_table$Is_ti[got_infected]/(epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected])
    }
    
    eta_old = eta
    
    # optim
    # the "strict" version: not fixing exp(eta) to > 1
    eta = optim(eta_old, llfunc, llgrr,
                method = "L-BFGS-B", 
                lower = lowerBound, upper = Inf)$par
    
    
    # (3) if difference < tol, stop
    if(max(abs(b_S_old - b_S)) < tol & abs(eta_old - eta) < tol){
      break
    }
    
  }
  
  # return the solved eta and b_S
  list(eta = eta, b_S = b_S)
}



# the partial data estimation function
# (for eta and b_S only)
infer_partial_data <- function(fpath, interval = 7, 
                               miss_recov_prop = 1, miss_expo_prop = 1,
                               tmax = 7, tmin = 0,
                               numIter=100, burn=0, maxIter=20, tol=1e-6,
                               initEta = 1.22, seed = 42, verbose = TRUE,
                               eta_only = FALSE){
  
  # fpath: folder name of data files under data_root
  # miss_**_prop: missing proportion of recovery times and exposure times
  # numIter: num. of iterations for the stochastic EM
  # burn: burn-in perioid for the stochastic EM
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  # initEta: the initial exp(eta) value to try
  # seed: seed
  # verbose: if verbose
  # eta_only: if only estimate eta
  
  set.seed(seed)
  
  # load data and extract stuff
  dat = load_data(data_root, fpath)
  
  events.orig = dat$events
  G0 = dat$G0
  X = dat$X
  I0 = dat$I0
  stage_change = dat$stage_change
  truth = dat$truth
  
  # 06/28/2021: only params exp_eta, eta, and b_S
  varnames = c('exp_eta', 'eta', 'b_S')
  
  if(verbose){cat(fpath,'data loaded...')}
  
  # create missingness
  miss_dat = miss_data(events.orig, G0, I0, interval, miss_recov_prop, miss_expo_prop)
  
  events = miss_dat$events
  if(verbose){cat('Missingness created...')}
  
  # get an adjmat report
  G_all = get_adjmat_report(G0, events, miss_dat$report.times)
  if(verbose){cat('Adjmat report created.\n')}
  
  # get intervals to impute recovery times on and those people's ids
  MR = get_miss_recov(miss_dat$report, miss_dat$report.times, events)
  
  # get list of manifested people and their times
  ## not include I0
  manifest = events %>% filter(event %in% c(9,10)) %>%
    filter(per1 != I0) %>%
    select(manifested = per1, times = time) %>%
    as.list()
  
  # # (11/06/2020) get the true exposure times for those manifested people
  # if(fix_expo_prop > 0){
  #   # put together truth
  #   true_expo_times = events.orig %>% 
  #     filter(event == 1, per1 %in% manifest$manifested) %>%
  #     select(per1, time)
  #   or = sapply(manifest$manifested, function(x) which(true_expo_times$per1==x))
  #   true_expo_times = true_expo_times[or,]
  #   
  #   # select a porportion of them to always fix
  #   n_mani = nrow(true_expo_times)
  #   selected = rbernoulli(n_mani, fix_expo_prop)
  #   fixed_expo = true_expo_times$per1[selected]
  #   fixed_expo_times = true_expo_times$time[selected]
  # }
  
  
  # get an initial conservative imputation of recovery times 
  # (everyone recovers at the end of interval)
  recovery_times = NULL
  for(r in 1:nrow(MR$intervals)){
    recov_r = data.frame(recov = MR$recover[[r]], time = MR$intervals$ub[r])
    recovery_times = rbind(recovery_times, recov_r)
  }
  
  # # 12/13/2020: fix recovery times to truth if...
  # if(fix_recovery){
  #   # get the truth
  #   true_recovery_times = events.orig %>% 
  #     filter(event == 2) %>%
  #     select(recov = per1, time)
  #   
  #   # order the truth by the order of people ids in "MR"
  #   proxy_recovery_times = NULL
  #   for(r in 1:nrow(MR$intervals)){
  #     true_r = true_recovery_times %>%
  #       filter(time < MR$intervals$ub[r] & time > MR$intervals$lb[r]) %>%
  #       arrange(recov)
  #     proxy_recovery_times = rbind(proxy_recovery_times, true_r)
  #   }
  #   
  #   recovery_times = proxy_recovery_times
  # }
  
  # start iterating...
  
  ## current values of exp_eta and gam
  ## start eta at the initial value and 
  ## others at the truth
  exp_eta = initEta
  gam = truth$gamma
  beta = truth$beta
  b_S = truth$b_S  
  phi = truth$phi
  
  ## storage for the samples
  params = list()
  nsamps = numIter - burn
  for(v in varnames){
    if(length(grep('(omega)|(alpha)|S',v)) > 0){
      # if the parameter is one of b_S, b_omega, b_alpha, omega, alpha
      params[[v]] = matrix(0, nrow=nsamps, ncol=length(truth[[v]]))
    }else{
      # else
      params[[v]] = numeric(nsamps)
    }
  }
  
  for(s in 1:numIter){
    if(verbose){cat('Iteration',s,'...')}
    
    # 1) impute exposure times first
    imp_expo_times = get_expo_times(manifest, G_all, tmax, tmin, 
                                    miss_dat$report, miss_dat$report.times, 
                                    events, recovery_times, X, b_S,
                                    exp_eta, beta, phi)
    
    # ## fix some of them to the truth if...
    # if(fix_expo_prop > 0){
    #   imp_expo_times[selected] = fixed_expo_times
    # }
    
    if(verbose){cat('Exposure times imputation done! ')}
    
    # 2) impute recovery times
    ## exposure time list
    exposure_times = list(exposed=manifest$manifested, times = imp_expo_times)
    ## construct local neighborhoods
    infec_nei = get_nei_expo_all(G_all, events, exposure_times, 
                                 miss_dat$report, miss_dat$report.times)
    
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
    
    ## sample recovery times
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
    summs = summarize_events2(G0, I0, events_aug, stage_change)
    
    if(verbose){cat('Summarizing events done!')}
    
    # 4) get MLEs (for eta and b_S)
    
    if(eta_only){
      estimate = solve_eta(summs, X, truth, initEta = initEta)
      estimates = list(exp_eta = estimate)
      estimates$eta = ifelse(estimate <= 0, -Inf, log(estimate))
      estimates$b_S = b_S
    }else{
      estimates = solve_eta_bS(summs, X, truth, initEta = initEta, 
                               maxIter=maxIter, tol = tol)
      estimates[['exp_eta']] = estimates[['eta']]
      estimates[['eta']] = ifelse(estimates[['eta']] <= 0, -Inf, log(estimates[['eta']]))
    }
    
    if(verbose){cat('MLE obtained!')}
    
    ## update current values of exp(eta) and b_S
    exp_eta = estimates[['exp_eta']]
    b_S = estimates[['b_S']]
    
    ## record it after burn-in
    if(s > burn){
      for(v in varnames){
        if(length(grep('(omega)|(alpha)|S',v)) > 0){
          # if the parameter is one of b_S, b_omega, b_alpha, omega, alpha
          params[[v]][s-burn,] = estimates[[v]]
        }else{
          # else
          params[[v]][s-burn] = estimates[[v]]
        }
      }
    }
    
    if(verbose){
      # output some key parameter values as well
      cat('Estimates:\n', 
          'exp_eta=', estimates$exp_eta, 'b_S=', estimates$b_S)
      cat('\nEstimates updated and saved.\n\n')
    }
  }
  
  # finally return results
  ## return mean values of parameter draws as well
  means = list()
  for(v in varnames){
    if(length(grep('(omega)|(alpha)|S',v)) > 0){
      # if the parameter is one of b_S, b_omega, b_alpha, omega, alpha
      means[[v]] = apply(params[[v]], 2, mean, na.rm=TRUE)
    }else{
      # else
      means[[v]] = mean(params[[v]], na.rm=TRUE)
    }
  }
  
  ## get the biases as well
  errors = list()
  for(v in varnames){
    errors[[v]] = truth[[v]] - means[[v]]
  }
  
  # 06/28/2020: also record pop size N
  list(params = params, means=means, truth = truth, errors = errors,
       N=nrow(dat$X))
}


# run things on server
# the "eta_only" indicator
EO = as.logical((s0-1) %% 2)

REP = 5
f1 = 13
f2 = 29

for(rep in 1:REP){
  this.seed = ss + rep*f1 %% f2
  part_res = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                                tmin = 0,
                                numIter = 80, burn = 50,
                                maxIter = 30, seed = this.seed,
                                eta_only = EO)
  
  ## save results
  EO_ind = ifelse(EO, 'eta_only', 'both')
  savepath = paste0(outdir,'N=',part_res$N, sub_dir, '_', EO_ind,'_',rep,'.rds')
  saveRDS(part_res, file=savepath)
}
