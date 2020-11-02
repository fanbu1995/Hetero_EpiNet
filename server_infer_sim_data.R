# 10/24/2020
# run simulation data inference code on server

# 10/31/2020
# add a "hack0" mode:
# if TRUE, manually fix b_* = 0
# to focus on the main parameters

# 11/02/2020:
# try to force longer latency by setting tmin = 5...

library(ggplot2)

source("inference_utils_1.R")

# data directory and outdir
data_root = 'hetero_data/'
outdir = 'hetero_results5/'


# server stuff, set seed and example data path
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
s0 = as.numeric(slurm_arrayid)
set.seed(s0)
ss = sample(1000,1)

sub_dir = paste0('ex',s0-1)

# source the inference util funcs
source("inference_utils_1.R")

# defined functions
## complete data
infer_complete_data <- function(fpath, maxIter=10, tol=1e-4, initEta = 1,
                                hack0 = FALSE){
  # fpath: folder name of data files under data_root
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  # initEta: the initial exp(eta) value to try
  
  # load data and extract stuff
  dat = load_data(data_root, fpath)
  
  events = dat$events
  G0 = dat$G0
  X = dat$X
  I0 = dat$I0
  stage_change = dat$stage_change
  truth = dat$truth
  
  # summarize data
  summs = summarize_events2(G0, I0, events, stage_change)
  
  # get MLEs
  estimates = solve_MLE(summs, X, maxIter, tol, initEta, hack0)
  ## adjust "eta" names
  estimates[['exp_eta']] = estimates[['eta']]
  estimates[['eta']] = ifelse(estimates[['eta']] <= 0, Inf, log(estimates[['eta']]))
  
  # calculate errors
  varnames = names(truth)
  errors= list()
  for(v in varnames){
    errors[[v]] = sqrt(mean((estimates[[v]] - truth[[v]])^2))
  }
  
  # return list of estimates and errors
  list(estimates = estimates, truth = truth, errors = errors)
}

## partial data
infer_partial_data <- function(fpath, interval = 7, miss_recov_prop = 1, miss_expo_prop = 1,
                               tmax = 7, tmin = 0,
                               numIter=100, burn=0, maxIter=20, tol=1e-6,
                               initEta = 1, initGam = 0.2, seed=42, verbose=TRUE,
                               hack0 = FALSE){
  
  # fpath: folder name of data files under data_root
  # miss_**_prop: missing proportion of recovery times and exposure times
  # numIter: num. of iterations for the stochastic EM
  # burn: burn-in perioid for the stochastic EM
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  # initEta: the initial exp(eta) value to try
  # initGam: the initial gamma value to use
  
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
  
  # get an initial conservative imputation of recovery times 
  # (everyone recovers at the end of interval)
  recovery_times = NULL
  for(r in 1:nrow(MR$intervals)){
    recov_r = data.frame(recov = MR$recover[[r]], time = MR$intervals$ub[r])
    recovery_times = rbind(recovery_times, recov_r)
  }
  
  # start iterating...
  
  ## current values of exp_eta and gam
  exp_eta = initEta
  gam = initGam
  
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
                                    events, recovery_times, exp_eta=exp_eta)
    
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
    
    # 4) get MLEs
    estimates = solve_MLE(summs, X, maxIter, tol, initEta, hack0)
    ## adjust "eta" names
    estimates[['exp_eta']] = estimates[['eta']]
    estimates[['eta']] = ifelse(estimates[['eta']] <= 0, Inf, log(estimates[['eta']]))
    
    if(verbose){cat('MLE obtained!')}
    
    ## update current values of exp(eta) and gamma
    exp_eta = estimates[['exp_eta']]
    gam = estimates[['gamma']]
    
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


## plot inference results
## 11/01/2020: add traceplot as well
plot_estimates <- function(res, trace=TRUE){
  # res: a combined list of both MLEs (complete data) and stoch. EM
  
  varnames = names(comp_res$estimates)
  
  MLEs = res$MLE
  params = res$params
  truth = res$truth
  
  for(v in varnames){
    if(length(grep('(omega)|(alpha)|S|s',v)) == 0){
      # if the parameter is NOT one of b_S, b_omega, b_alpha, omega, alpha
      # also: don't plot p_s either, since this is always the same as MLE estimate
      # then plot everything else
      pdat = tibble(val = params[[v]])
      ptruth = truth[[v]]
      print(
        ggplot(pdat, aes(x=val)) + 
          geom_histogram(bins = 10, fill='skyblue') +
          geom_vline(xintercept = truth[[v]], size=1.5, color='red') +
          geom_vline(xintercept = MLEs[[v]], size=1.5, color='black') +
          labs(x=v) +
          theme_bw()
      )
      
      # add traceplot if...
      if(trace){
        pdat = pdat %>% mutate(samp = 1:n())
        print(
          ggplot(pdat, aes(y=val,x=samp)) + 
            geom_hline(yintercept = truth[[v]], size=1.5, color='red') +
            geom_hline(yintercept = MLEs[[v]], size=1.5, color='darkgray') +
            geom_line(size=0.3)+
            labs(y=v, x='iteration') +
            theme_bw()
        )
      }
      
    }
    # if the parameter has length > 1: don't make plots for now
  }
}




# inference and plot

## 10/31/2020 more change: hack0 = TRUE to focus on main parameters
hh = TRUE

## complete data results (MLE)
comp_res = infer_complete_data(sub_dir, maxIter = 50, tol = 1e-6, 
                               initEta = 1.5, hack0 = hh)

## partial data

## 10/31/2020 change: 
# 1. burn-in 100 iters and only take last 100 samples...
# 2. increase "tmax" to 50 (so expo time imputation window will be all the history!)
# 3. try 5 different seeds (thus expo time imputation will start at different spots)

## 11/01/2020 more attempt:
# set tmin = 1, to force longer latency...

REP = 5

## plot things together
pdfpath = paste0(outdir,"res_",sub_dir,".pdf")
pdf(pdfpath, width = 8, height = 6)

for(rep in 1:REP){
  this.seed = ss + rep*103 %% 11
  part_res = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                                tmin = 1,
                                numIter = 150, burn = 50,
                                maxIter = 50, seed = this.seed,
                                hack0 = hh)
  
  ## save results
  part_res$MLE = comp_res$estimates
  part_res$MLE_errors = comp_res$errors
  savepath = paste0(outdir,"_",sub_dir,"_",rep,".rds")
  saveRDS(part_res, file=savepath)
  
  ## plot this repetition result
  plot_estimates(part_res)
}




