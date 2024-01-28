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

# 2.adapted from complete data inference
# but ONLY estimate eta
infer_eta_only <- function(fpath, initEta = 1, lowerBound = 1e-6, otherTruth=TRUE){
  # fpath: folder name of data files under data_root
  # initEta: the initial exp(eta) value to try
  # lowerBound: the lower bound used in "optim" (default 1e-6, nearly 0)
  # otherTruth: if fixing the other parameters at truth
  
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
  # add an option to do full-on estimation of all parameters (without any fixing)
  if(otherTruth){
    estimate = solve_eta(summs, X, truth, initEta, lowerBound)
  }else{
    estimate = solve_MLE(summs, X, maxIter = 20, initEta = initEta)
    estimate = estimate$eta
  }
  
  
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
ds_nums = c(0:40)
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
  save(etas, exp_etas, Ns, file='eta_and_Ns_100_200.RData')
}

overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
overall

# eta_avg   eta_std  Eeta_avg  Eeta_std 
# 0.1960815 0.4666778 1.3630100 0.7011361



# try a bunch of datasets with N=500
## try just one dataset (to ensure that we can handle it)
infer_eta_only('ex41',1.22)
## okay! things didn't break...now let's do it!

ds_nums = c(41:70)
init_Eta = 1.22
etas = NULL
exp_etas = NULL
#Ns = NULL
for(ds in ds_nums){
  fp = paste0('ex',ds)
  res = infer_eta_only(fp, init_Eta)
  etas = c(etas, res$eta)
  exp_etas = c(exp_etas, res$exp_eta)
  #Ns = c(Ns, res$N) # record the "N" of each dataset as well
  cat('Dataset',ds, 'eta =',res$eta, 'exp(eta) =',res$exp_eta,'\n')
  # save as it goes (in case my session breaks again)
  save(etas, exp_etas, file='eta_500.RData')
}

overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
overall

## save "overall" to RData object as well
save(etas, exp_etas, overall, file='eta_500.RData')

# ?????
# eta_avg     eta_std    Eeta_avg    Eeta_std 
# -0.11178316  0.04552449  0.89513996  0.04114567 


# try with N = 300??
ds_nums = c(71:100)
init_Eta = 1.22
etas = NULL
exp_etas = NULL
#Ns = NULL
for(ds in ds_nums){
  fp = paste0('ex',ds)
  res = infer_eta_only(fp, init_Eta)
  etas = c(etas, res$eta)
  exp_etas = c(exp_etas, res$exp_eta)
  #Ns = c(Ns, res$N) # record the "N" of each dataset as well
  cat('Dataset',ds, 'eta =',res$eta, 'exp(eta) =',res$exp_eta,'\n')
  # save as it goes (in case my session breaks again)
  save(etas, exp_etas, file='eta_300.RData')
}

#overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
# filter out "bad eggs" - unstable performance of "optim"
weird = which(etas < -10 | etas > 10)
etas = etas[-weird]
exp_etas = exp_etas[-weird]
overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
overall

# save overall too
save(etas, exp_etas, overall, file='eta_300.RData')


# do more N=200
# use 50 datasets (all have N=200)
# LATER: try lowerBound = 1
load('eta_and_Ns_100_200.RData')
sum(Ns==200) # already have 26

ds_nums = c(which(Ns==200),101:124)
init_Eta = 1.22
LB = 1 # default: 1e-6

etas = NULL
exp_etas = NULL
#Ns = NULL
for(ds in ds_nums){
  fp = paste0('ex',ds)
  res = infer_eta_only(fp, init_Eta, LB)
  etas = c(etas, res$eta)
  exp_etas = c(exp_etas, res$exp_eta)
  #Ns = c(Ns, res$N) # record the "N" of each dataset as well
  cat('Dataset',ds, 'eta =',res$eta, 'exp(eta) =',res$exp_eta,'\n')
  # save as it goes (in case my session breaks again)
  #save(etas, exp_etas, file='eta_200.RData')
  save(etas, exp_etas, file='eta_200_LB1.RData')
}

overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
overall
# LB = 1e-6 (default)
# eta_avg    eta_std   Eeta_avg   Eeta_std 
# 0.18314651 0.09976085 1.20689739 0.12180927 

## save "overall" to RData object as well
load('eta_200.RData')
overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
save(etas, exp_etas, overall, file='eta_200.RData')

# LB = 1 (a hacky way)
# (2 of 50 were hacked to 0, the lower limit)
# eta_avg    eta_std   Eeta_avg   Eeta_std 
# 0.18436144 0.09724235 1.20809268 0.11954989 


## produce a summary table for N=200, 300, 500 together
tab = NULL
for(n in c(200,300,500)){
  load(paste0('eta_',n,'.RData'))
  # generate "overall" just to be safe
  overall = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
  tab = rbind(tab, overall)
}
rownames(tab) = NULL
tab = as.data.frame(tab)
names(tab) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
tab$N = c(200,300,500)

tab %>% select(N, eta_avg:Eeta_std)


## do another round of N=200 
## summarize eta BUT estimation is done for ALL params
load('eta_and_Ns_100_200.RData')
#sum(Ns==200) # already have 26

ds_nums = c(which(Ns==200),101:124)
init_Eta = 1.22
# LB = 1 # default: 1e-6

etas = NULL
exp_etas = NULL
#Ns = NULL
for(ds in ds_nums){
  fp = paste0('ex',ds)
  res = infer_eta_only(fp, init_Eta, LB, FALSE)
  etas = c(etas, res$eta)
  exp_etas = c(exp_etas, res$exp_eta)
  #Ns = c(Ns, res$N) # record the "N" of each dataset as well
  cat('Dataset',ds, 'eta =',res$eta, 'exp(eta) =',res$exp_eta,'\n')
  # save as it goes (in case my session breaks again)
  #save(etas, exp_etas, file='eta_200.RData')
  save(etas, exp_etas, file='eta_200_full_on.RData')
}

overall_full = c(mean(etas),sd(etas),mean(exp_etas),sd(exp_etas))
names(overall_full) = c('eta_avg','eta_std','Eeta_avg','Eeta_std')
overall_full

