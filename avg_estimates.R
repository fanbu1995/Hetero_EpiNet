#11/14/2020

# compute the average estimates (and standard deviation) 
# of simulated data

# 11/29/2020

# also check the estimates for b_S (compare it with zero b_S truth)

library(tidyverse)
library(stringr)
library(xtable)

setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/hetero_results11/')

root = '~/Documents/Research_and_References/Hetero_EpiNet_2020/'

res = readRDS('_ex12_1.rds')

all_files = list.files()
str_detect(all_files, '.rds')
str_detect(all_files, '_ex11')

# function to compute averages and standard deviations
# for certain parameters
# over select results
# 12/13/2020: also check the MLE estimates (from complete data)

get_avg_est <- function(path, datasets, 
                        params=c('beta','phi','gamma','p_s','b_S'),
                        b_dim = 2){
  fpath = file.path(root, path)
  all_files = list.files(path = fpath)
  
  res_files = all_files[str_detect(all_files, '.rds')]
  
  # list of param averages
  # (add list of MLEs)
  avg_ll = list()
  mle_ll = list()
  for(p in params){
    if(p == 'b_S'){
      for(i in 1:b_dim){
        pa = paste0(p,i)
        avg_ll[[pa]] = numeric(0)
        mle_ll[[pa]] = numeric(0)
      }
    }
  }
  
  b_S_L = 0
  
  for(d in datasets){
    d_files = res_files[str_detect(res_files, d)]
    for(f in d_files){
      res = readRDS(file.path(root, path, f))
      for(p in params){
        if(p == 'b_S'){
          b_S_L = length(res$means[[p]])
          cat('b_S estimate in ', f, ': ', res$means[[p]], '\n')
          # also MLE
          cat('b_S MLE in ', f, ': ', res$MLE[[p]], '\n')
          for(i in 1:b_dim){
            pa = paste0(p,i)
            avg_ll[[pa]] = c(avg_ll[[pa]], res$means[[p]][i])
            mle_ll[[pa]] = c(mle_ll[[pa]], res$MLE[[p]][i])
          }
        }else{
          if(p=='beta'){
            cat('beta estimate in ', f, ": ", res$means[[p]], '\n')
            cat('beta MLE in ', f, ": ", res$MLE[[p]], '\n')
          }
          avg_ll[[p]] = c(avg_ll[[p]], res$means[[p]])
          mle_ll[[p]] = c(mle_ll[[p]], res$MLE[[p]])
        }
      }
    }
  }
  
  avg_dat = as.data.frame(avg_ll)
  summ = t(apply(avg_dat, 2, function(v) c(mean(v), sd(v))))
  summ_dat = as.data.frame(summ)
  #names(summ_dat) = c('Mean', 'SD')
  
  mle_dat = as.data.frame(mle_ll)
  mle_summ =  t(apply(mle_dat, 2, function(v) c(mean(v), sd(v))))
  mle_summ_dat = as.data.frame(mle_summ)
  summ_dat = cbind(summ_dat, mle_summ_dat)
  names(summ_dat) = c('Mean', 'SD', 'MLE_mean', 'MLE_SD')
  
  summ_dat
}


# 12/21/2020: calculate error for multi-dim parameters
get_mle_errors <- function(path, datasets, 
                           params = c('b_S', 'alpha', 'omega')){
  fpath = file.path(root, path)
  all_files = list.files(path = fpath)
  
  res_files = all_files[str_detect(all_files, '.rds')]
  
  # list of errors
  error_ll = list()
  for(p in params){
    error_ll[[p]] = numeric(0)
  }
  
  # go through datasets
  for(d in datasets){
    d_files = res_files[str_detect(res_files, d)]
    for(f in d_files){
      res = readRDS(file.path(root, path, f))
      for(p in params){
        error_ll[[p]] = c(error_ll[[p]], res$MLE_errors[[p]])
      }
    }
  }
  
  # compile dataframe
  error_dat = as.data.frame(error_ll)
  summ = t(apply(error_dat, 2, function(v) c(mean(v), sd(v))))
  summ_dat = as.data.frame(summ)
  names(summ_dat) = c('Avg. MAE', 'SD')
  
  summ_dat
}

## another function to peek at the ground truth of example datasets
peek_truth <- function(path, datasets, 
                       params=c('beta')){
  fpath = file.path(root, path)
  all_files = list.files(path = fpath)
  
  res_files = all_files[str_detect(all_files, '.rds')]
  
  for(p in params){
  for(d in datasets){
    d_files = res_files[str_detect(res_files, d)]
      for(f in d_files){
        res = readRDS(file.path(root, path, f))
        truth = res$truth
        cat('In', f, 'truth for', p, 'is: ', truth[[p]], '\n')
      }
    }
  }
  
}

summ1 = get_avg_est('hetero_results11', c('ex12','ex14','ex20'))

summ2 = get_avg_est('hetero_results11', c('ex11','ex15','ex17','ex19'))

# 11/29/2020
# look at hetero_results12 results
summ1 = get_avg_est('hetero_results12', 
                    c('ex12','ex14', 'ex16', 'ex18', 'ex20'),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))
summ1

library(xtable)
xtable(summ1, digits=4)

## seems like there is some error in each single case, but 
## on average it's kind of okay (though variance is hmm...)

# 11/30/2020
# look at hetero_results13 results (b_S is NOT zero)
# hmmm not working....
summ2 = get_avg_est('hetero_results13', 
                    paste0('ex', c(3,5,7,9)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

summ2 = get_avg_est('hetero_results13', 
                    paste0('ex', c(2,4,6,8,10)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

# 12/02/2020
# look at hetero_results14 results (b_S = (0, 1))

summ3 = get_avg_est('hetero_results14', 
                    paste0('ex', seq(from=22,to=30, by=2)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

summ3 = get_avg_est('hetero_results14', 
                    paste0('ex', c(28)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

# 12/06/2020
# check out hetero_results15 results (with b_S = (0, 1) fixed in inference)
summ4 = get_avg_est('hetero_results15', 
                    paste0('ex', seq(from=21,to=30, by=2)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

# 12/13/2020
# check out hetero_results14 results (with MLEs)
summ3 = get_avg_est('hetero_results14', 
                    paste0('ex', seq(from=21,to=30, by=2)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

peek_truth('hetero_results14', 
           paste0('ex', seq(from=22,to=30, by=2)),
           params = c('beta','phi','gamma','p_s','b_S','eta'))

# check out hetero_results15 results and also look at MLEs
summ4 = get_avg_est('hetero_results15', 
                    paste0('ex', seq(from=21,to=30, by=2)),
                    params = c('beta','phi','gamma','p_s','b_S','eta'))

peek_truth('hetero_results15', 
           paste0('ex', seq(from=21,to=30, by=2)),
           params = c('beta','phi','gamma','p_s','b_S','eta'))


# check out hetero_results16 (fix phi and eta, only estimate beta and b_S)
summ5 = get_avg_est('hetero_results16', 
                    paste0('ex', seq(from=22,to=30, by=2)),
                    params = c('beta', 'gamma','p_s','b_S'))


# check out hetero_results17 (fix recovery times to truth)
summ7 = get_avg_est('hetero_results17', 
                    paste0('ex', seq(from=22,to=30, by=2)),
                    params = c('beta', 'eta', 'phi', 'gamma','p_s','b_S'))


## 11/29/2020
# peek at ground truth from ex0-9
peek_truth('hetero_results3', paste0('ex',1:9), params=c('beta','b_S'))


# 12/21/2020
# check out MLE errors
summ6 = get_mle_errors('hetero_results17', 
                       paste0('ex', seq(from=22,to=30, by=2)))

peek_truth('hetero_results17', 'ex22', params = c('beta','omega', 'alpha'))

