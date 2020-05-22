# 05/14/2020
# Heterogeneous EpiNet
# MLE script

# NOTES: assumptions about data structure
# Data should contain
# 1) Event list:
# event marked by number code: 1 (infection), 2 (recovery), 
#                              3-5 (connection), 6-8 (disconnection),
#                              9 (symptom onset)
# 2) G0: the initial network (a sparse matrix)
# 3) I0: label(s) of the initially infected person(s)
# 4) stage_change: time point that NPI starts
#                  (assume TWO stages for now: T0 before this, and T1 after)
#
# Also: assume an SIR model; recovered -> immune


library(tidyverse)

# 0.5 a function that takes in a pair p=(i,j) and outputs its index 
# of their entry in the vectorized representation of the upper-triangular NxN matrix
get_vec_index <- function(p,N){
  # make sure that i < j
  if(p[1] > p[2]){
    p = p[2:1]
  }
  i = p[1]; j = p[2]
  
  (2*N-i) * (i-1)/2 + j - i
}

# 0.5 another function that takes in p=(i,j) and an epid vector
# and then output a one-hot vector for link type
get_link_type <- function(p, epid){
  N = length(epid)
  
  # make sure that i < j
  if(p[1] > p[2]){
    p = p[2:1]
  }
  
  res = numeric(3)
  i = p[1]; j = p[2]
  
  res[1] = (epid[i] %in% c(0,-1)) & (epid[j] %in% c(0,-1))
  res[3] = (epid[i] %in% c(1,2)) & (epid[j] %in% c(1,2))
  
  if((res[1] == 0) & (res[3] == 0)){
    res[2] = 1
  }
  
  return(res)
}

# 1. a function to parse through complete data 
# and acquire necessary summary statistics
summarize_events <- function(G0, I0, events, stage_change){
  
  # get pop size N
  # and adjmat & epid vector
  N = nrow(G0)
  adjmat = G0
  epid = rep(0,N)
  epid[I0] = 1 
  ## in the "epid" vector: 
  ## 1=infected, 2=infected&symptomatic, 0=susceptible, -1=recovered
  
  # get nI and nR
  nI = sum(events$event==1)
  nR = sum(events$event==2)
  
  # separate two stages
  events_T0 = events %>% filter(time <= stage_change)
  events_T1 = events %>% filter(time > stage_change)
  
  # get C counts and D counts
  type_phase = c('HH0','HI0','II0','HH1','HI1','II1')
  C_all = numeric(6); D_all = numeric(6)
  
  for(type in c(3:5)){
    # connections before change
    C_all[type-2] = sum(events_T0$event==type)
    # connections after change
    C_all[type+1] = sum(events_T1$event==type)
  }
  for(type in c(6:8)){
    # disconnections before change
    D_all[type-5] = sum(events_T0$event==type)
    # disconnections
    D_all[type-2] = sum(events_T1$event==type)
  }
  
  names(C_all) = type_phase
  names(D_all) = type_phase
  
  # data storage for 
  # individual epidemic info
  # and pair-wise network info
  epi_table = matrix(0, nrow=N, ncol=5)
  net_c_table = matrix(0, nrow=N*(N-1)/2, ncol=7)
  net_d_table = matrix(0, nrow=N*(N-1)/2, ncol=7)
  
  t_pre = 0
  
  # go through the events one by one
  for(r in 1:nrow(events)){
    z = events$event[r]
    t_cur = events$time[r]
    
    # obtain cumulated time of exposure for each indiv. up to now
    susceptible = which(epid==0)
    num_of_Es = apply(adjmat[susceptible,], 1, function(l) sum(l[epid==1]))
    epi_table[susceptible,1] = epi_table[susceptible,1] + 
      num_of_Es * (t_cur - t_pre)
    num_of_Is = apply(adjmat[susceptible,], 1, function(l) sum(l[epid==2]))
    epi_table[susceptible,2] = epi_table[susceptible,2] + 
      num_of_Is * (t_cur - t_pre)
    
    # obtain cumulative time of staying connected/disconnected for each pair up to now
    
    ## connected pairs first
    connected = which(upper.tri(adjmat) & adjmat == 1, arr.ind = T)
    c_types = apply(connected, 1, get_link_type, epid)
    c_inds = apply(connected, 1, get_vec_index, N)
    
    ## then disconnected pairs
    disconnected = which(upper.tri(adjmat) & adjmat == 0, arr.ind = T)
    d_types = apply(disconnected, 1, get_link_type, epid)
    d_inds = apply(disconnected, 1, get_vec_index, N)
    
    if(t_cur < stage_change){
      # operate under T0
      net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + t(c_types) * (t_cur - t_pre)
      net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + t(d_types) * (t_cur - t_pre)
    }else if(t_pre > stage_change){
      # operate under T1
      net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + t(c_types) * (t_cur - t_pre)
      net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + t(d_types) * (t_cur - t_pre)
    }else{
      # straddling T0 and T1
      # both need update
      net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + 
        t(c_types) * (stage_change - t_pre)
      net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + 
        t(d_types) * (stage_change - t_pre)
      
      net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + 
        t(c_types) * (t_cur - stage_change)
      net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + 
        t(d_types) * (t_cur - stage_change)
    }
    
    # then process changes to the system
    if (z==1){
      # infection
      p1 = events$per1[r]
      epid[p1] = 1
      
      # record num of E and I friends at time of infection
      epi_table[p1,3] = sum(adjmat[p1,epid==1])
      epi_table[p1,4] = sum(adjmat[p1,epid==2])
      
      # also take down infection time
      epi_table[p1,5] = -1 * t_cur
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      epid[p1] = -1
      
      # calculate total infection time
      epi_table[p1,5] = epi_table[p1,5] + t_cur
    }else if (z==9){
      # symptom onset
      p1 = events$per1[r]
      epid[p1] = 2
    }else{
      # some edge stuff
      p1 = events$per1[r]
      p2 = events$per2[r]
      if(z %in% c(3:5)){
        # reconnection
        adjmat[p1,p2] = 1; adjmat[p2,p1] = 1
        
        # record connection in table
        ind = get_vec_index(c(p1,p2),N)
        net_d_table[ind,1] = net_d_table[ind,1] + 1
      }else{
        # disconnection
        adjmat[p1,p2] = 0; adjmat[p2,p1] = 0
        
        # record disconnection in table
        ind = get_vec_index(c(p1,p2),N)
        net_c_table[ind,1] = net_c_table[ind,1] + 1
      }
    }
    
    # report progress
    cat("Processing...", r/nrow(events),"\r")
    t_pre = t_cur
  }
  
  # if at the end, some one is still sick...
  if(any(epid %in% c(1,2))){
    still_sick = which(epid %in% c(1,2))
    epi_table[still_sick,5] = epi_table[still_sick,5] + t_cur
  }
  
  # make the epi network tables into dataframes
  epi_table = as.data.frame(epi_table)
  names(epi_table) = c('E_expo','I_expo','E_ti','I_ti','sick_time')
  net_c_table = as.data.frame(net_c_table)
  names(net_c_table) = c('Nc',type_phase)
  net_d_table = as.data.frame(net_d_table)
  names(net_d_table) = c('Nd',type_phase)
    
  # change the name a bit...
  # net_c_table --> net_d_table (# of connections and time spent disconnected)
  # and same for net_d_table
  return(list(nI=nI, nR=nR, C_all = C_all, D_all = D_all,
              epi_table = epi_table, net_c_table = net_d_table,
              net_d_table = net_c_table))
}

### try it out
# dat = readRDS('~/Documents/Research_and_References/EpiNet_Codes_ImmediateData/ex3dat_1.rds')
# 
# summ = summarize_events(dat$G0, dat$I0, dat$events, 10)

# 1.5 a function that sums up X_i and X_j given p=(i,j)
sum_covariates <- function(p, X){
  # make sure that i < j
  if(p[1] > p[2]){
    p = p[2:1]
  }
  i = p[1]; j = p[2]
  
  X[i,]+X[j,]
}

# 2. a function that takes in the summaries and carry out MLE
# (numerically solve)
solve_MLE <- function(summaries, X, maxIter=10, tol=1e-4, initEta = 3){
  
  # X has to be a dataframe with column names!!
  
  # wrangle X into needed shape of Xij (for network part)
  # (original X: n by p design matrix)
  # (want XX: n*(n-1)/2 by p design matrix)
  N = nrow(X)
  pairs = which(upper.tri(matrix(nrow=N, ncol=N)), arr.ind = T)
  XX = apply(pairs, 1, sum_covariates, X)
  XX = t(XX)
  
  
  # extract summary statistics
  nI = summaries$nI
  nR = summaries$nR
  C_all = summaries$C_all
  D_all = summaries$D_all
  epi_table = summaries$epi_table
  net_c_table = summaries$net_c_table
  net_d_table = summaries$net_d_table
  
  # 1. get recovery rate gamma
  gamma = nR/sum(epi_table$sick_time)
  
  # 2. get parameters for the infection side: beta, eta, b_S
  eta = initEta
  beta = nI/(sum(epi_table$E_expo) + eta * sum(epi_table$I_expo))
  
  got_infected = which(epi_table$sick_time > 0)
  Y_infec = (epi_table$sick_time > 0) %>% as.numeric()
  
  for(it in 1:maxIter){
    # 2.1 estimate b_S
    Expo = (epi_table$E_expo + eta * epi_table$I_expo) * beta
    has_expo = which(Expo > 0)
    Y_infec = Y_infec[has_expo]
    X_infec = X[has_expo,]
    Expo = Expo[has_expo]
    infec_poi = glm(Y_infec~X_infec-1+offset(log(Expo)), family = poisson())
    b_S = infec_poi$coefficients %>% as.numeric()
    
    # 2.2 estimate beta
    beta_old = beta
    indiv_effects = exp(c(X_infec %*% b_S))
    beta = nI/sum(indiv_effects * Expo/beta_old)
    
    # 2.3 estimate eta
    llfunc = function(eta){
      beta * sum(indiv_effects * (epi_table$E_expo + eta * epi_table$I_expo)[has_expo]) -
        sum(log((epi_table$E_ti + epi_table$I_ti * eta)[got_infected]))
    }
    llgrr = function(beka){
      eta * beta * sum(indiv_effects * epi_table$I_expo[has_expo]) -
        eta * sum(epi_table$I_ti[got_infected]/(epi_table$E_ti + epi_table$I_ti * eta)[got_infected])
    }
    eta_old = eta
    eta = optim(eta_old, llfunc, llgrr, 
                method = "L-BFGS-B", lower = 1e-6, upper = Inf)$par
    
    # 2.4 if difference < tol, stop
    if(abs(beta_old - beta) < tol | abs(eta_old - eta) < tol){
      break
    }
    
  }
  
  # 3. get parameters for link activation: alphas, b_alpha
  alpha = numeric(6)
  names(alpha) = c('HH0','HI0','II0','HH1','HI1','II1')
  
  for(type in names(alpha)){
    total_dur = sum(net_c_table[[type]])
    if(total_dur > 0){
      alpha[type] = C_all[type]/total_dur
    }else{
      alpha[type] = 0
    }
  }
  
  for(it in 1:maxIter){
    # 3.1 estimate b_alpha
    all_weights = c(as.matrix(net_c_table[2:7]) %*% alpha)
    has_weight = which(all_weights > 0)
    XX_alpha = XX[has_weight,]
    all_weights = all_weights[has_weight]
    alpha_poi = glm(net_c_table$Nc[has_weight]~XX_alpha-1+offset(log(all_weights)), 
                    family = poisson())
    b_alpha = alpha_poi$coefficients %>% as.numeric()
    
    # 3.2 estimate alpha
    alpha_old = alpha
    pair_effects = exp(c(XX_alpha * b_alpha))
    for(type in names(alpha)){
      total_dur = sum(net_c_table[[type]][has_weight] * pair_effects)
      if(total_dur > 0){
        alpha[type] = C_all[type]/total_dur
      }else{
        alpha[type] = 0
      }
    }
    
    # 3.3 if difference < tol, stop
    if(all(abs(alpha_old - alpha)) < tol){
      break
    }
  }

  
  # 4. get parameters for link termination: omegas, b_omega
  omega = numeric(6)
  names(omega) = c('HH0','HI0','II0','HH1','HI1','II1')
  
  for(type in names(omega)){
    total_dur = sum(net_d_table[[type]])
    if(total_dur > 0){
      omega[type] = D_all[type]/total_dur
    }else{
      omega[type] = 0
    }
  }
  
  for(it in 1:maxIter){
    # 3.1 estimate b_omega
    all_weights = c(as.matrix(net_d_table[2:7]) %*% omega)
    has_weight = which(all_weights > 0)
    XX_omega = XX[has_weight,]
    all_weights = all_weights[has_weight]
    omega_poi = glm(net_d_table$Nd[has_weight]~XX_omega-1+offset(log(all_weights)), 
                    family = poisson())
    b_omega = omega_poi$coefficients %>% as.numeric()
    
    # 3.2 estimate omega
    omega_old = omega
    pair_effects = exp(c(XX_omega * b_omega))
    for(type in names(omega)){
      total_dur = sum(net_d_table[[type]][has_weight] * pair_effects)
      if(total_dur > 0){
        omega[type] = D_all[type]/total_dur
      }else{
        omega[type] = 0
      }
    }
    
    # 3.3 if difference < tol, stop
    if(all(abs(omega_old - omega)) < tol){
      break
    }
  }
  
  
  return(list(beta=beta, eta=eta, b_S=b_S, gamma=gamma, 
              alpha=alpha, b_alpha=b_alpha, omega=omega, b_omega=b_omega))
}



