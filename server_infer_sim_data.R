# 10/24/2020
# run simulation data inference code on server

# 08/16/2021
# track the acceptance ratio of the rejection sampler

# 10/31/2020
# add a "hack0" mode:
# if TRUE, manually fix b_* = 0
# to focus on the main parameters

# 11/02/2020:
# try to force longer latency by setting tmin = 5...

## 11/06/2020 more attempt:
# 1. fix eta and phi and only focus on beta estimation
# 2. set some fraction of expo times to truth and see what happens

## 11/08/2020:
# try fixing the expo time imputation algorithm, and
# 1. still with hack0 = TRUE, just run stuff
# 2. set some fraction of expo times to truth

## 11/10/2020
# fix eta, phi, etc. and only focus on beta again
# to test out imputation algorithm

## 11/11/2020
# fixed a little bug in exposure time imputation code
# and retry running stuff only with b_S = 0 fixed

## 11/24/2020
# try running stuff without any 0 fixing...

## 11/29/2020
# try running stuff (without fixing) on ex0 & 2-10 (where b_S = (0,-1))

## 11/30/2020
# try running new stuff (without fixing) where b_S = (0,1)
# and see if beta gets over-estimated...

## 12/13/2020
# try running with eta fixed
# and see if beta and b_S can be correctly handled

# also: try fixing all recovery times to truth (but not anything else)

## 02/24/2020:
# try with X[,2] = rnorm(N,sd=0.5)
# AND b_S = (1,1)

library(ggplot2)

# data directory and outdir
data_root = 'hetero_data/'
#outdir = 'hetero_results17/'

# 01/26/2021
# change outdir to new dir
# outdir = 'fix_some4/'

# 02/24/2021
# change outdir again
outdir = 'new_setting0/'


# server stuff, set seed and example data path
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
s0 = as.numeric(slurm_arrayid)
set.seed(s0)
ss = sample(1000,1)

sub_dir = paste0('ex',s0-1)

# source the inference util funcs
source("inference_utils_1.R")

# 08/16/2021: new utils
source("inference_utils_3.R")

# defined functions
## complete data
infer_complete_data <- function(fpath, maxIter=10, tol=1e-4, initEta = 1,
                                hack0 = FALSE, hackb_S = FALSE){
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
  estimates = solve_MLE(summs, X, maxIter, tol, initEta, hack0=hack0, 
                        hackb_S = hackb_S, true_b_S = truth$b_S)
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

## 11/06/2020: 
# 1. make things easier by starting at the truth of exp(eta) and gamma
# 2. set a "fixed_prop" to fix a fraction of exposure times to truth

## 12/13/2020
# add a hack to fix recovery times at truth
infer_partial_data <- function(fpath, interval = 7, 
                               miss_recov_prop = 1, miss_expo_prop = 1,
                               tmax = 7, tmin = 0,
                               numIter=100, burn=0, maxIter=20, tol=1e-6,
                               initEta = 1.22, initGam = 0.1, 
                               initBeta = 0.2, initPhi = 0.2,
                               initB_S = c(0,0),
                               seed=42, verbose=TRUE,
                               hack0 = FALSE, hack_eta = FALSE,
                               hack_phi = FALSE,
                               hackb_S = FALSE,
                               fix_expo_prop = 0,
                               fix_recovery = FALSE,
                               track_accept = TRUE){
  
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
  # hack_eta: whether or not to hack eta to truth (default FALSE)
  # hack_phi: whether or not to hack phi to truth (default FALSE)
  # hack_b_S: whether or not to hack b_S to truth (default FALSE)
  # fix_expo_prop: proportion of exposure times to fix to truth (default 0)
  # fix_recovery: whether or not to hack recovery times to truth (default FALSE)
  # track_accept: if track the acceptance ratios for expo time rejection sampler
  
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
  if(fix_expo_prop > 0){
    # put together truth
    true_expo_times = events.orig %>% 
      filter(event == 1, per1 %in% manifest$manifested) %>%
      select(per1, time)
    or = sapply(manifest$manifested, function(x) which(true_expo_times$per1==x))
    true_expo_times = true_expo_times[or,]
    
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
  
  # 08/16/2021: track acceptance ratio for rejection sampler
  accepts = NULL
  for(s in 1:numIter){
    if(verbose){cat('Iteration',s,'...')}
    
    # 1) impute exposure times first
    imp_expo_times = get_expo_times(manifest, G_all, tmax, tmin, 
                                    miss_dat$report, miss_dat$report.times, 
                                    events, recovery_times, X, b_S,
                                    exp_eta, beta, phi, track_accept)
    if(track_accept){
      accepts_s = imp_expo_times$accepts
      if(s==1){
        accepts = accepts_s
      }else{
        accepts = accepts + accepts_s
      }
      imp_expo_times = imp_expo_times$res
    }
    
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
                          hack0, hackb_S, hack_eta, hack_phi, 
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
  
  list(params = params, means=means, truth = truth, accepts = accepts)
}


## plot inference results
## 11/01/2020: add traceplot as well

## 01/04/2021
## optional plotting of truth and MLE
plot_estimates <- function(res, trace=TRUE, truthMLE = TRUE){
  # res: a combined list of both MLEs (complete data) and stoch. EM
  # trace: whether or not to include traceplot
  # truthMLE: whether or not there exists (or to include) truth and MLE estimates
  
  varnames = names(res$means)
  
  params = res$params
  
  if(truthMLE){
    MLEs = res$MLE
    truth = res$truth
  }
  
  
  for(v in varnames){
    if(length(grep('(omega)|(alpha)|S|s',v)) == 0){
      # if the parameter is NOT one of b_S, b_omega, b_alpha, omega, alpha
      # also: don't plot p_s either, since this is always the same as MLE estimate
      # then plot everything else
      pdat = tibble(val = params[[v]])
      if(truthMLE){
        print(
          ggplot(pdat, aes(x=val)) + 
            geom_histogram(bins = 10, fill='skyblue') +
            geom_vline(xintercept = truth[[v]], size=1.5, color='red') +
            geom_vline(xintercept = MLEs[[v]], size=1.5, color='black') +
            labs(x=v) +
            theme_bw()
        )
      }else{
        print(
          ggplot(pdat, aes(x=val)) + 
            geom_histogram(bins = 10, fill='skyblue') +
            #geom_vline(xintercept = truth[[v]], size=1.5, color='red') +
            #geom_vline(xintercept = MLEs[[v]], size=1.5, color='black') +
            labs(x=v) +
            theme_bw()
        )
      }
      
      
      # add traceplot if...
      if(trace){
        pdat = pdat %>% mutate(samp = 1:n())
        if(truthMLE){
          print(
            ggplot(pdat, aes(y=val,x=samp)) + 
              geom_hline(yintercept = truth[[v]], size=1.5, color='red') +
              geom_hline(yintercept = MLEs[[v]], size=1.5, color='darkgray') +
              geom_line(size=0.3)+
              labs(y=v, x='iteration') +
              theme_bw()
            )
        }else{
          print(
            ggplot(pdat, aes(y=val,x=samp)) + 
              #geom_hline(yintercept = truth[[v]], size=1.5, color='red') +
              #geom_hline(yintercept = MLEs[[v]], size=1.5, color='darkgray') +
              geom_line(size=0.3)+
              labs(y=v, x='iteration') +
              theme_bw()
          )
        }
      }
      
    }
    # if the parameter has length > 1: don't make plots for now
  }
}




# inference and plot

## 10/31/2020 more change: hack0 = TRUE to focus on main parameters

## 11/24/2020: no longer fixing hh, and see how that works
hh = FALSE

## 11/06/2020 change: hack_eta_phi = TRUE to focus on beta only
## then set it back...

## 11/10/2020: fix everything again and focus on beta only
## to test out the new imputation algorithm
## (11/11/2020: set it back)


## 12/01/2020: fix b_S to truth and see if things work fine
## trying to get at why the regression part works weird

## 12/13/2020: only fix eta but not b_S and see what happens

hep = FALSE

hbS = FALSE

## 12/13/2020: fix recovery times
fre = FALSE


## complete data results (MLE)
comp_res = infer_complete_data(sub_dir, maxIter = 50, tol = 1e-6, 
                               initEta = 1.5, hack0 = hh, hackb_S = hbS)

## partial data

## 10/31/2020 change: 
# 1. burn-in 100 iters and only take last 100 samples...
# 2. increase "tmax" to 50 (so expo time imputation window will be all the history!)
# 3. try 5 different seeds (thus expo time imputation will start at different spots)

## 11/01/2020 more attempt:
# set tmin = 1, to force longer latency...
# (this is reversed in later experiments)

## 11/06/2020 more attempt:
# 1. fix eta and phi and only focus on beta estimation
# 2. set some fraction of expo times to truth and see what happens

## 01/26/2021:
## fix the following params:
## "fix_some_*"
# 1. phi
# 2. eta
# 3. b_S
# 4. eta + b_S

# hp = FALSE
# he = TRUE
# hbS = TRUE
# 
# f1 = 1361
# f2 = 53


## 02/24/2021
# try again with b_S = (1,1) and X[,2] = Normal(0,0.5)
hp = FALSE
he = FALSE
hbS = FALSE

f1 = 1361
f2 = 53

REP = 5

## set initial beta according to dataset
i_beta = ifelse(s0%%2 == 0, 0.2, 0.15)

## 02/24/2021: no need to have different betas
## (for ex31-40: all N=200 and beta=0.15)
i_beta = 0.15

## plot things together
pdfpath = paste0(outdir,"res_",sub_dir,".pdf")
pdf(pdfpath, width = 8, height = 6)

for(rep in 1:REP){
  this.seed = ss + rep*f1 %% f2
  #fep = rep/5
  fep = 0
  part_res = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                                tmin = 0,
                                numIter = 100, burn = 50,
                                maxIter = 30, seed = this.seed,
                                initBeta = i_beta,
                                hack0 = hh, hack_eta = he,
                                hack_phi = hp,
                                hackb_S = hbS,
                                fix_expo_prop = fep, 
                                fix_recovery = fre)
  
  ## save results
  part_res$MLE = comp_res$estimates
  part_res$MLE_errors = comp_res$errors
  savepath = paste0(outdir,"_",sub_dir,"_",rep,".rds")
  saveRDS(part_res, file=savepath)
  
  ## plot this repetition result
  plot_estimates(part_res)
}


###### not run on server!! ########
# 08/16/2021
data_root = "~/Documents/Research_and_References/Hetero_EpiNet_2020"

dats = paste0('ex', c(31:40))
sub_dir = 'ex40'

part_res40 = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                              tmin = 0,
                              numIter = 50, burn = 20,
                              maxIter = 20, seed = 73,
                              initBeta = 0.15,
                              track_accept = TRUE)

sub_dir = 'ex101'
part_res101 = infer_partial_data(sub_dir, interval = 7, tmax=50, 
                                tmin = 0,
                                numIter = 50, burn = 20,
                                maxIter = 20, seed = 41,
                                initBeta = 0.15,
                                track_accept = TRUE)


