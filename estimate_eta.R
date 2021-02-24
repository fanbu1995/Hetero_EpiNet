# 02/24/2021:
# inference using complete data
# but ONLY focus on estimating \eta
# to see how bad it is

setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')

# source the inference util funcs
source("inference_utils_1.R")

data_root = '~/Documents/Research_and_References/Hetero_EpiNet_2020/'


# 1. adapted from solve_MLE
# but only estimate eta (while fixing everything else at the truth)
solve_eta <- function(summaries, X, truth, initEta = 1){
  
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
  beta = nI/(sum(epi_table$Ia_expo) + eta * sum(epi_table$Is_expo))
  
  # 10/24/20 fix: only get those who had infec neighbors at infec time!
  got_infected = which((epi_table$Ia_ti + epi_table$Is_ti) > 0 & epi_table$latent_time > 0)
  got_infected = got_infected[got_infected != I0]
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  
  Expo = (epi_table$Ia_expo + eta * epi_table$Is_expo)
  has_expo = which(Expo > 0) # check I0 not here!
  ## set Y_infec to the original thing before re-subsetting
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  Y_infec = Y_infec[has_expo]
  X_infec = X[has_expo,]
  
  indiv_effects = exp(c(X_infec %*% b_S))
  
  # estimate eta
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
              method = "L-BFGS-B", lower = 1e-6, upper = Inf)$par
  
  # return the solved eta
  eta
}

# 2.adapted from complete data inference
# but ONLY estimate eta
infer_eta_only <- function(fpath, initEta = 1){
  # fpath: folder name of data files under data_root
  # initEta: the initial exp(eta) value to try
  
  # load data and extract stuff
  dat = load_data(data_root, fpath)
  
  events = dat$events
  G0 = dat$G0
  X = dat$X
  I0 = dat$I0
  stage_change = dat$stage_change
  truth = dat$truth
  
  N = nrow(X)
  
  # summarize data
  summs = summarize_events2(G0, I0, events, stage_change)
  
  # get eta through MLE estimation
  estimate = solve_eta(summs, X, truth, initEta)
  
  # return both eta and exp(eta) in the estimates
  eta = ifelse(estimate <= 0, -Inf, log(estimate))
  estimates = list(exp_eta = estimate, eta = eta)
  
  # return N as well just to get a sense
  estimates$N = N
  
  # FOR NOW: return only the estimates for the dataset
  estimates
  
  # # calculate errors
  # varnames = names(estimates)
  # errors= list()
  # for(v in varnames){
  #   errors[[v]] = sqrt(mean((estimates[[v]] - truth[[v]])^2))
  # }
  # 
  # # return list of estimates and errors
  # list(estimates = estimates, 
  #      truth = list(exp_eta= truth$exp_eta, eta=truth$eta), 
  #      errors = errors)
}


# then do this for ALL the datasets I have for now (41 in total)
ds_nums = c(0:41)
init_Eta = 1.22
etas = NULL
exp_etas = NULL
Ns = NULL
for(ds in ds_nums){
  fp = paste0('ex',ds)
  res = infer_eta_only(fp, init_Eta)
  etas = c(etas, res$eta)
  exp_etas = c(exp_etas, res$exp_eta)
  Ns = c(Ns, res$N) # record the "N" of each dataset as well
  cat('Dataset',ds,'N =', res$N, 'eta =',res$eta, 'exp(eta) =',res$exp_eta,'\n')
  # save as it goes (in case my session breaks again)
  save(etas, exp_etas, Ns, file='eta_and_Ns.RData')
}

overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')


