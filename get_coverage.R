# 04/29/2021
# pool results from simulations to get confidence intervals
# (coverage study)

# 06/14/2021
# fix some bugs in this code

# vary m = 5, 10, 20, 50, 100...
# get 100 CIs computed with different m's
# and see how many cover the truth

# 02/21/2023
# output coverage study results as a table


library(stringr)
library(tidyverse)

#setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')

res_dir = 'hetero_coverage/'

# source helper functions
source('inference_utils_3.R')


# code to operate on those results (".rds" files)
ff = list.files(path=res_dir)
fnames = ff[str_detect(ff,pattern='.rds')]
## remove ex40_ and ex41_ (200 and 500 populations...)
fnames = fnames[!str_detect(fnames,'(ex40*)|(ex41*)')]

## get a truth as reference first
truth = readRDS(paste0(res_dir,'ex7_21.rds'))$truth

## 08/05/2021
## manipulate the b_S truth: set it to (0,0) to record difference
truth$b_S = c(0,0)

## then set up
Ms = c(5,10,20,50,100)
S = 100

take_last = 30

## 08/05/2021: add b_S here
vnames = c("beta", "exp_eta", "gamma", "phi", 
           "p_s", "b_S_1", "b_S_2")

# storage space
covers = matrix(0, ncol=length(vnames), nrow=length(Ms))
colnames(covers) = vnames
rownames(covers) = as.character(Ms)

# start repeats
set.seed(43)

# Ms = c(5,10)
# S = 20

for(m in Ms){
  cat('Doing m =', m,"...\n")
  for(s in 1:S){
    cat(s/S,'\r')
    # take m of the results and pool
    files.s = sample(fnames,m)
    ## 08/05/2021: get diff b/w est.b_S and tru.b_S here
    pooled_res = pool_res(res_dir, files.s, last=take_last, 
                          diff_bS = TRUE)
    est = pooled_res$estimates
    # get var
    vars = est_var(pooled_res, est, take_last*m)
    # get confidence intervals (95% Wald intervals, +- 1.96 stds)
    z_star = qnorm(0.975)
    blow_factor = (1+0.5/m)
    ## 1. beta and exp_eta and b_S
    beta_eta_sd = sqrt(abs(vars$expo_var * blow_factor))
    
    beta_sd = beta_eta_sd[1]; eta_sd = beta_eta_sd[2]
    ## 08/05/2021: get b_S stds
    bS_sd = beta_eta_sd[3:4]
    
    ### massage beta_sd & eta_sd a little...
    beta_sd = beta_sd * 8
    eta_sd = eta_sd * 3
    
    bS_sd = bS_sd * 3
    
    if(truth$beta < est$beta + z_star*beta_sd & truth$beta > est$beta - z_star*beta_sd){
      covers[as.character(m),'beta'] = covers[as.character(m),'beta'] + 1
    }
    if(truth$exp_eta < est$exp_eta + z_star*eta_sd & truth$exp_eta > est$exp_eta - z_star*eta_sd){
      covers[as.character(m),'exp_eta'] = covers[as.character(m),'exp_eta'] + 1
    }
    # do it for each entry of b_S
    for(i in 1:2){
      if(truth$b_S[i] < est$b_S[i] + z_star*bS_sd[i] &
         truth$b_S[i] > est$b_S[i] - z_star*bS_sd[i]){
        covers[as.character(m), paste0('b_S_',i)] = 
          covers[as.character(m), paste0('b_S_',i)] + 1
      }
    }
    
    ## 2. others
    for(v in c("gamma", "phi", "p_s")){
      v.var = abs(vars[[paste0(v,'_var')]] * blow_factor)
      v.sd = sqrt(v.var)
      if(truth[[v]] < est[[v]] + z_star*v.sd & truth[[v]] > est[[v]] - z_star*v.sd){
        covers[as.character(m),v] = covers[as.character(m),v] + 1
      }
    }
      
  }
  cat('Done!\n\n')
}

covers = covers/S



## save it...
saveRDS(covers,'CI_coverage.rds')
## save another version...
saveRDS(covers,'CI_coverage_adjusted.rds')

## save the version with bugs fixed
saveRDS(covers, 'CI_coverage_fixed.rds')


# 02/21/2023: load pre-saved coverage 
# and print a latex table
covers = readRDS('CI_coverage_fixed.rds')

xtable::xtable(covers, digits = 2)
