#11/14/2020

# compute the average estimates (and standard deviation) 
# of simulated data

library(tidyverse)
library(stringr)

setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/hetero_results11/')

root = '~/Documents/Research_and_References/Hetero_EpiNet_2020/'

res = readRDS('_ex12_1.rds')

all_files = list.files()
str_detect(all_files, '.rds')
str_detect(all_files, '_ex11')

# function to compute averages and standard deviations
# for certain parameters
# over select results
get_avg_est <- function(path, datasets, params=c('beta','phi','gamma','p_s')){
  fpath = file.path(root, path)
  all_files = list.files(path = fpath)
  
  res_files = all_files[str_detect(all_files, '.rds')]
  
  # list of param averages
  avg_ll = list()
  for(p in params){
    avg_ll[[p]] = numeric(0)
  }
  
  for(d in datasets){
    d_files = res_files[str_detect(res_files, d)]
    for(f in d_files){
      res = readRDS(file.path(root, path, f))
      for(p in params){
        avg_ll[[p]] = c(avg_ll[[p]], res$means[[p]])
      }
    }
  }
  
  avg_dat = as.data.frame(avg_ll)
  summ = t(apply(avg_dat, 2, function(v) c(mean(v), sd(v))))
  summ_dat = as.data.frame(summ)
  names(summ_dat) = c('Mean', 'SD')
  summ_dat
}

summ1 = get_avg_est('hetero_results11', c('ex12','ex14','ex20'))

summ2 = get_avg_est('hetero_results11', c('ex11','ex15','ex17','ex19'))

library(xtable)
xtable(summ1, digits=4)
