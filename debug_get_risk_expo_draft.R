
# 12/13/2020
# try to debug expo time imputation functions

setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
source('inference_utils_1.R')

data_root = './'

## old code...
## deep debugging of the get_expo_risk_i function

# manifest_times = manifest
# report = miss_dat$report
# times = miss_dat$report.times
# 
# manifested = manifest$manifested
# 
# for(i in manifested){
#   cat('Person',i,': ')
#   t_i = manifest_times$times[manifested == i]
#   tE_i = get_expo_risk_i(i, t_i, G_all, tmax, tmin, report, times, 
#                   events, recovery_times, exp_eta, details = FALSE)
#   cat('Imputed exposure time:',tE_i,'\n')
# }


# a slightly modified version of the original inference function (for partial data)
infer_partial_data <- function(fpath, interval = 7, 
                               miss_recov_prop = 1, miss_expo_prop = 1,
                               tmax = 7, tmin = 0,
                               numIter=100, burn=0, maxIter=20, tol=1e-6,
                               initEta = 1.22, initGam = 0.1, 
                               initBeta = 0.2, initPhi = 0.2,
                               initB_S = c(0,0),
                               seed=42, verbose=TRUE,
                               hack0 = FALSE, hack_eta_phi = FALSE,
                               hackb_S = FALSE,
                               fix_expo_prop = 0,
                               fix_recovery = FALSE){
  
  # fpath: folder name of data files under data_root
  # miss_**_prop: missing proportion of recovery times and exposure times
  # numIter: num. of iterations for the stochastic EM
  # burn: burn-in perioid for the stochastic EM
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  # initEta: the initial exp(eta) value to try
  # initGam: the initial gamma value to use
  # initBeta: the initial beta value to use
  # initB_S: the initial b_S vector to use
  # hack0: whether or not to hack b_* coefficients to zero (default FALSE)
  # hack_eta_phi: whether or not to hack eta and phi to truth (default FALSE)
  # hack_b_S: whether or not to hack b_S to truth (default FALSE)
  # fix_expo_prop: proportion of exposure times to fix to truth (default 0)
  # fix_recovery: whether or not to hack recovery times to truth (default FALSE)
  
  set.seed(seed)
  
  # load data and extract stuff
  dat = load_data(data_root, fpath)
  
  events.orig = dat$events
  G0 = dat$G0
  X = dat$X
  I0 = dat$I0
  stage_change = dat$stage_change
  truth = dat$truth
  varnames = names(truth)
  
  if(verbose){cat('Data loaded...')}
  
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
  
  # (11/06/2020) get the true exposure times for those manifested people
  # put together truth
  true_expo_times = events.orig %>% 
    filter(event == 1, per1 %in% manifest$manifested) %>%
    select(per1, time)
  # get the order correct
  or = sapply(manifest$manifested, function(x) which(true_expo_times$per1==x))
  true_expo_times = true_expo_times[or,]
  if(fix_expo_prop > 0){
    # select a porportion of them to always fix
    n_mani = nrow(true_expo_times)
    selected = rbernoulli(n_mani, fix_expo_prop)
    fixed_expo = true_expo_times$per1[selected]
    fixed_expo_times = true_expo_times$time[selected]
  }
  
  
  # get an initial conservative imputation of recovery times 
  # (everyone recovers at the end of interval)
  recovery_times = NULL
  for(r in 1:nrow(MR$intervals)){
    recov_r = data.frame(recov = MR$recover[[r]], time = MR$intervals$ub[r])
    recovery_times = rbind(recovery_times, recov_r)
  }
  
  # 12/13/2020: fix recovery times to truth if...
  if(fix_recovery){
    # get the truth
    true_recovery_times = events.orig %>% 
      filter(event == 2) %>%
      select(recov = per1, time)
    
    # order the truth by the order of people ids in "MR"
    proxy_recovery_times = NULL
    for(r in 1:nrow(MR$intervals)){
      true_r = true_recovery_times %>%
        filter(time < MR$intervals$ub[r] & time > MR$intervals$lb[r]) %>%
        arrange(recov)
      proxy_recovery_times = rbind(proxy_recovery_times, true_r)
    }
    
    recovery_times = proxy_recovery_times
  }
  
  # start iterating...
  
  ## current values of exp_eta and gam
  exp_eta = initEta
  gam = initGam
  beta = initBeta
  b_S = initB_S
  phi = initPhi
  
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
    # 12/13/2020: use the debug version
    # imp_expo_times = get_expo_times(manifest, G_all, tmax, tmin, 
    #                                 miss_dat$report, miss_dat$report.times, 
    #                                 events, recovery_times, X, b_S,
    #                                 exp_eta, beta, phi)
    
    # 12/13/2020: use the debug version
    imp_expo_times = get_expo_times_debug(manifest, G_all, tmax, tmin, 
                                          miss_dat$report, miss_dat$report.times, 
                                          events, recovery_times, X, b_S,
                                          exp_eta, beta, phi, true_expo_times)
    
    ## fix some of them to the truth if...
    if(fix_expo_prop > 0){
      imp_expo_times[selected] = fixed_expo_times
    }
    
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
    
    # 12/13/2020...
    # allow fixing recovery times to truth
    if(fix_recovery){
      recovery_times = proxy_recovery_times
    }else{
      recov_times = NULL
      for(r in 1:nrow(MR$intervals)){
        lb = MR$intervals$lb[r]; ub = MR$intervals$ub[r]
        recovs = MR$recover[[r]]
        cands = propose_recov_filter(lb, ub, recovs, exposure_times, infec_nei, gam, exp_eta)
        recov_times = c(recov_times, cands)
      }
      recovery_times$time = recov_times
    }
    
    if(verbose){cat('Recovery times imputation done!\n')}
    
    # 3) combine samples with data to get augmented events and summarize those
    events_aug = combine_data(events, exposure_times, recovery_times)
    summs = summarize_events2(G0, I0, events_aug, stage_change)
    
    if(verbose){cat('Summarizing events done!')}
    
    # 4) get MLEs
    ## 11/06/2020 update with hacked eta and phi included
    estimates = solve_MLE(summs, X, maxIter, tol, initEta, 
                          hack0, hack_eta_phi, hackb_S,
                          truth$exp_eta, truth$phi, truth$b_S)
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
      cat('Key estimates:\n beta=', estimates$beta, 'phi= ', estimates$phi,
          'gamma=', estimates$gamma,
          'exp_eta=', estimates$exp_eta, 'p_s=', estimates$p_s)
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
  
  list(params = params, means=means, truth = truth)
}



# try it on "ex30":
s0 = 31
set.seed(s0)
ss = sample(1000,1)

sub_dir = paste0('ex',s0-1)

i_beta = 0.15

this.seed = ss + 1*1453 %% 41
part_res = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                              tmin = 0,
                              numIter = 150, burn = 0,
                              maxIter = 50, seed = this.seed,
                              initBeta = i_beta)


#####
## go over some "problematic" persons

dat = load_data(data_root, sub_dir)
events.orig = dat$events
G0 = dat$G0
X = dat$X
I0 = dat$I0
stage_change = dat$stage_change
truth = dat$truth
varnames = names(truth)

# create missingness
miss_dat = miss_data(events.orig, G0, I0, 7, 1, 1)

events = miss_dat$events

# get an adjmat report
G_all = get_adjmat_report(G0, events, miss_dat$report.times)

# get intervals to impute recovery times on and those people's ids
MR = get_miss_recov(miss_dat$report, miss_dat$report.times, events)

# get list of manifested people and their times
## not include I0
manifest = events %>% filter(event %in% c(9,10)) %>%
  filter(per1 != I0) %>%
  select(manifested = per1, times = time) %>%
  as.list()

# true recovery times
true_recovery_times = events.orig %>% 
  filter(event == 2) %>%
  select(recov = per1, time)

# order the truth by the order of people ids in "MR"
proxy_recovery_times = NULL
for(r in 1:nrow(MR$intervals)){
  true_r = true_recovery_times %>%
    filter(time < MR$intervals$ub[r] & time > MR$intervals$lb[r]) %>%
    arrange(recov)
  proxy_recovery_times = rbind(proxy_recovery_times, true_r)
}

recovery_times = proxy_recovery_times

# get some other version of recovery times
# randomly sample in each interval
recovery_times = NULL
for(r in 1:nrow(MR$intervals)){
  recov_r = data.frame(recov = MR$recover[[r]], 
                       time = runif(length(MR$recover[[r]]),
                                    min =MR$intervals$lb[r], 
                                    max = MR$intervals$ub[r]))
  recovery_times = rbind(recovery_times, recov_r)
}


# look at persons 102, 44, 109
i = 69
t_i = manifest$times[manifest$manifested==i]
x_i = X[i,]
tE_i = get_expo_risk_i(i, t_i, G_all, tmax=50, tmin=0, 
                       miss_dat$report, miss_dat$report.times,
                       events, recovery_times, x_i, b_S = c(-0.3,-0.2),
                       exp_eta=0.7, beta = 0.2, phi = 0.2, details=TRUE)

# check exposure local neighborhood
get_nei_expo_i(i, tE_i$samp, G0[i,], events, 
               miss_dat$report, miss_dat$report.times)
