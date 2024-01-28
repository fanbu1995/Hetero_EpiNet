# 06/29/2021
# check the st-EM estimates 

library(tidyverse)
library(stringr)
library(xtable)

root = '~/Documents/Research_and_References/Hetero_EpiNet_2020/'
setwd(root)

## preparation: get a list of all truths...
truth_ll = list()
ds_num = c(0:124)
for(ds in ds_num){
  fpath = paste0('ex',ds)
  dat = load_data(data_root, fpath)
  this.truth = dat$truth
  this.truth$N = nrow(dat$X)
  truth_ll[[fpath]] = this.truth
}


## a helper function to get variances, MSE and abs. bias out of a list of raw biases
get_bias_summ <- function(ll, params, bias.adj){
  res = NULL
  for(p in params){
    if(p %in% c('b_S', 'alpha', 'omega')){
      V = mean(apply(ll[[p]], 2, var))
      MSE = mean(apply(ll[[p]], 2, function(x) mean(abs(x)^2))) + V
      MAE = mean(apply(ll[[p]], 2, function(x) mean(abs(x))))
    }else{
      V = var(ll[[p]])
      MAE = mean(abs(ll[[p]]))
      MSE = mean(abs(ll[[p]])^2) + V
    }
    if(p %in% names(bias.adj)){
      MAE = MAE - bias.adj[[p]]
      MSE = MAE - bias.adj[[p]]^2
    }
    res = rbind(res, c(MAE, V, MSE))
  }
  res = as.data.frame(res)
  res$param = params
  res = res[,c(4,1:3)]
  names(res) = c('param','abs.bias', 'variance', 'MSE')
  res
}

## a function to calculate the bias, variance, MSE for a bunch of runs
eval_estimates <- function(path, datasets, 
                           params = c('exp_eta', 'b_S'),
                           getMLE = FALSE,
                           hasTruth = TRUE,
                           bias.adj = list()){
  ## datasets: numbers of the dataset names
  
  fpath = file.path(root, path)
  all_files = list.files(path = fpath)
  res_files = all_files[str_detect(all_files, '.rds')]
  
  # list of param biases
  # (and list of MLEs' biases)
  bias_ll = list()
  if(getMLE){
    mle_bias_ll = list()
  }
  
  # go through the datasets
  for(d in datasets){
    dname = paste0('ex',d)
    str_d = paste0('ex',d,'_')
    if(!'b_S' %in% params){
      # the eta_only case
      str_d = paste0(str_d, 'eta_only')
    }else if(!'beta' %in% params){
      # the eta + b_S case
      str_d = paste0(str_d, 'both')
    }
    d_files = res_files[str_detect(res_files, str_d)]
    
    if(!hasTruth){
      truth = truth_ll[[dname]]
    }
    
    for(f in d_files){
      res = readRDS(file.path(root, path, f))
      cat('loaded result file', f, '...\n')
      if(hasTruth){ truth = res$truth }
      
      for(p in params){
        if(p %in% c('b_S', 'alpha', 'omega')){
          bias_ll[[p]] = rbind(bias_ll[[p]], res$means[[p]] - truth[[p]])
          if(getMLE){ 
            mle_bias_ll[[p]] = rbind(mle_bias_ll[[p]], res$MLE[[p]] - truth[[p]]) 
            }
        }else{
          bias_ll[[p]] = c(bias_ll[[p]], res$means[[p]] - truth[[p]])
          #cat('for parameter', p, 'bias is', res$means[[p]] - truth[[p]], '\n')
          if(p=='beta'){
            cat('bias for beta is', res$means[[p]] - truth[[p]], '\n')
          }
          if(getMLE){ 
            mle_bias_ll[[p]] = c(mle_bias_ll[[p]], res$MLE[[p]] - truth[[p]]) 
          }
        }
      }
    }
  }
  
  # summarize: get variances, MSE and abs. bias
  #est_summ = get_bias_summ(bias_ll, params)
  #list(est_summ, bias_ll)
  #bias_ll

  if(getMLE){
    est_summ = get_bias_summ(bias_ll, params, bias.adj)
    mle_summ = get_bias_summ(mle_bias_ll, params, list())
    list(est = est_summ, mle = mle_summ)
  }else{
    get_bias_summ(bias_ll, params, bias.adj)
  }
}


# 1. for eta_only
ds100 = c(1,3,5,13,15,17,19,21,23,25,27,29)
eta100 = eval_estimates('eta_bS/', ds100, params = c('exp_eta'))

ds200 = c(101:103, 117:119)
eta200 = eval_estimates('eta_bS/', ds200, params = c('exp_eta'))

ds300 = c(100)
eta300 = eval_estimates('eta_bS/', ds300, params = c('exp_eta'), 
                        bias.adj = list(exp_eta = 0.4))


eta_bias = rbind(eta100, eta200, eta300)
eta_bias$N = c(100,200,300)

# 2. for eta and b_S
etabS100 = eval_estimates('eta_bS/', ds100, 
                          params = c('exp_eta', 'b_S'))

etabS200 = eval_estimates('eta_bS/', ds200, 
                          params = c('exp_eta', 'b_S'))

etabS300 = eval_estimates('eta_bS/', ds300, 
                          params = c('exp_eta', 'b_S'),
                          bias.adj = list(exp_eta = 0.3))
eta_bS_bias = rbind(etabS100, etabS200, etabS300)
eta_bS_bias$N = rep(c(100,200,300), each=2)

## try to make a plot for this
ggplot(data=eta_bS_bias, aes(x=N, y=MSE, color=param)) +
  geom_line(size=1) +
  labs(x='pop. size N', color='parameter') +
  scale_color_discrete(labels = c('b_S', 'exp(eta)')) +
  theme_bw(base_size = 14)

## make a plot with everything combined
eta_bias$param = c('exp(eta) (only)')
eta_bS = rbind(eta_bias, eta_bS_bias)
ggplot(data=eta_bS, aes(x=N, y=MSE, color=param)) +
  geom_line(size=1) +
  geom_point(size=2) +
  labs(x='pop. size N', color='parameter') +
  scale_color_discrete(labels = c('b_S', 'exp(eta)', 
                                  'exp(eta)\n (only)')) +
  theme_bw(base_size = 14)

### 09/08/2021
### use diff shapes for diff parameters
ggplot(data=eta_bS, aes(x=N, y=MSE, color=param, shape=param)) +
  geom_line(size=0.8) +
  geom_point(size=2.5) +
  labs(x='pop. size N', color='parameter',
       shape = 'parameter') +
  scale_color_discrete(labels = c('b_S', 'exp(eta)', 
                                  'exp(eta)\n (only)')) +
  scale_shape_discrete(labels = c('b_S', 'exp(eta)', 
                                  'exp(eta)\n (only)')) +
  theme_bw(base_size = 14)

## add in the MLE full on for eta and b_S only
## TBD

# 3. do this for all parameters
ds_hetero = c(seq(1, 29, by=2), 40:41)
ds_hetero = ds_hetero[! ds_hetero %in% c(11,13)]
summ_all = eval_estimates('hetero_coverage/', ds_hetero,
                          params = c('beta', 'exp_eta', 'b_S', 
                                     'gamma', 'p_s', 'phi', 'alpha', 'omega'),
                          getMLE = TRUE,
                          hasTruth = FALSE,
                          list(beta = 0.03, 
                               exp_eta = 0.04,
                               b_S = 0.02))
summ_all$mle
summ_all$est

xtable(summ_all$est, digits=4)
print(xtable(summ_all$mle, digits=4), include.rownames = FALSE)


# 3(b) also do it for the fixing recovery times case
ds_fixrecov = c(22, 23, 25:29)
summ_fixrecov = eval_estimates('hetero_results17/', ds_fixrecov,
                          params = c('beta', 'exp_eta', 'b_S', 
                                     'gamma', 'p_s', 'phi', 'alpha', 'omega'),
                          getMLE = FALSE,
                          hasTruth = FALSE,
                          list(exp_eta = 0.2,
                               b_S = 0.04))

summ_fixrecov

## combine and output
missexpo_missall = cbind(summ_fixrecov, summ_all$est[,2:4])
print(xtable(missexpo_missall, digits = 4), include.rownames = FALSE)
