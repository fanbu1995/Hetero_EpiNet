# 01/04/2024: simulation experiment with external infection parameters

setwd('~/Documents/Research/Hetero_EpiNet_2020/')

data_root = '~/Documents/Research/Hetero_EpiNet_2024rev/'

source("inference_utils_4.R")


infer_complete_data <- function(fpath, maxIter=20, tol=1e-4, initEta = 1, ...){
  # fpath: folder name of data files under data_root
  # maxIter: max. iterations for the MLE optimization algorithms
  # tol: tolerance for MLE optimization algorithms
  # initEta: the initial exp(eta) value to try
  
  # load data and extract stuff
  dat = load_data2(data_root, fpath)
  
  events = dat$events
  G0 = dat$G0
  X = dat$X
  I0 = dat$I0
  stage_change = dat$stage_change
  truth = dat$truth
  
  # summarize data
  summs = summarize_events2(G0, I0, events, stage_change, ext_labels = dat$ext_label)
  
  # get MLEs
  estimates = solve_MLE(summs, X, maxIter, tol, initEta,
                        true_b_S = truth$b_S, true_b_E = truth$b_E,
                        ...)
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


# try it out
res = infer_complete_data("ex1N100", tol = 1e-6, maxIter = 30, hackb_E = TRUE)
res = infer_complete_data("ex10N200", tol = 1e-6, maxIter = 30, hackb_E = FALSE)

# go through all 100 simulated datasets and combine
allRes = list()
for (i in 1:100){
  N = ifelse(i %% 2==0, 200, 100)
  filename = sprintf("ex%sN%d", i, N)
  
  this.res = infer_complete_data(filename, tol = 1e-6, maxIter = 30)
  
  allRes$xi = c(allRes$xi, this.res$estimates$xi)
  allRes$b_E1 = c(allRes$b_E1, this.res$estimates$b_E[1])
  allRes$b_E2 = c(allRes$b_E2, this.res$estimates$b_E[2])
  allRes$xiError = c(allRes$xiError, this.res$errors$xi)
  allRes$bEerror = c(allRes$bEerror, this.res$errors$b_E)
  
  cat(sprintf("inference done for dataset # %s \n\n", i))
}

saveRDS(allRes, "externalSim100datasets.rds")

allRes = as.data.frame(allRes)

getAvgIQR <- function(x){
  avg = mean(x)
  iqr = quantile(x, c(.25, .75))
  
  res = sprintf("%.3f (%.3f, %.3f)", avg, iqr[1], iqr[2])
}

errorSummary = 
allRes %>% mutate(N = rep(c(100,200), n()/2)) %>% 
  group_by(N) %>%
  summarize(xiSumm = getAvgIQR(xi),
            bE1Summ = getAvgIQR(b_E1),
            bE2Summ = getAvgIQR(b_E2),
            xiErrorSumm = getAvgIQR(xiError),
            bEerrorSumm = getAvgIQR(bEerror),
            xiRme = ) %>%
  ungroup()
  
