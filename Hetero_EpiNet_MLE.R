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
# 4) stage_change: a vector of time points during which there is NPI
#                  (stage_change[0]: starting point, stage_change[1]: end point)
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

# CHANGE (09.2020)
# epid vector: 0=S, -1=E, 1=Ia, 2=Is, -2=R
# S, E, R are H
# Ia, Is are I
get_link_type <- function(p, epid){
  N = length(epid)
  
  # make sure that i < j
  if(p[1] > p[2]){
    p = p[2:1]
  }
  
  res = numeric(3)
  i = p[1]; j = p[2]
  
  res[1] = (epid[i] %in% c(0,-1,-2)) & (epid[j] %in% c(0,-1,-2))
  res[3] = (epid[i] %in% c(1,2)) & (epid[j] %in% c(1,2))
  
  if((res[1] == 0) & (res[3] == 0)){
    res[2] = 1
  }
  
  return(res)
}

# 1. a function to parse through complete data 
# and acquire necessary summary statistics
# CHANGE: 09/21/2020
# epid vector: 0=S, -1=E, 1=Ia, 2=Is, -2=R
summarize_events <- function(G0, I0, events, stage_change){
  
  # get pop size N
  # and adjmat & epid vector
  N = nrow(G0)
  adjmat = G0
  epid = rep(0,N)
  epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  ## CHANGED!! see above for new coding scheme
  ## in the "epid" vector: 
  ## 1=infected, 2=infected&symptomatic, 0=susceptible, -1=recovered
  
  # get event counts
  #nI = sum(events$event==1)
  #nR = sum(events$event==2)
  nE = sum(events$event==1)
  nR = sum(events$event==2)
  nIa = sum(events$event==9)
  nIs = sum(events$event==10)
  #nI = nIa + nIs
  
  # CHANGED
  # get st and en of NPI period
  st = stage_change[1]; en = stage_change[2]
  
  # separate two stages
  events_T0 = events %>% filter(time <= st | time >= en)
  events_T1 = events %>% filter(time > st & time < en)
  
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
  ## CHANGED: add "latent_time" column in epi_table
  epi_table = matrix(0, nrow=N, ncol=6)
  net_c_table = matrix(0, nrow=N*(N-1)/2, ncol=7)
  net_d_table = matrix(0, nrow=N*(N-1)/2, ncol=7)
  
  t_pre = 0
  
  # go through the events one by one
  for(r in 1:nrow(events)){
    z = events$event[r]
    t_cur = events$time[r]
    
    # obtain cumulated time of exposure for each indiv. up to now
    susceptible = which(epid==0)
    ## Change: 09/21/2020
    if(length(susceptible) > 0){
      num_of_Ias = apply(adjmat[susceptible,,drop=F], 1, function(l) sum(l[epid==1]))
      epi_table[susceptible,1] = epi_table[susceptible,1] + 
        num_of_Ias * (t_cur - t_pre)
      num_of_Iss = apply(adjmat[susceptible,,drop=F], 1, function(l) sum(l[epid==2]))
      epi_table[susceptible,2] = epi_table[susceptible,2] + 
        num_of_Iss * (t_cur - t_pre)
    }
    # obtain cumulative time of staying connected/disconnected for each pair up to now
    
    ## connected pairs first
    connected = which(upper.tri(adjmat) & adjmat == 1, arr.ind = T)
    c_types = apply(connected, 1, get_link_type, epid)
    c_inds = apply(connected, 1, get_vec_index, N)
    
    ## then disconnected pairs
    disconnected = which(upper.tri(adjmat) & adjmat == 0, arr.ind = T)
    d_types = apply(disconnected, 1, get_link_type, epid)
    d_inds = apply(disconnected, 1, get_vec_index, N)
    
    if(t_cur < st | t_pre > en){
      # operate under T0
      net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + t(c_types) * (t_cur - t_pre)
      net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + t(d_types) * (t_cur - t_pre)
    }else if(t_pre > st & t_cur < en){
      # operate under T1
      net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + t(c_types) * (t_cur - t_pre)
      net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + t(d_types) * (t_cur - t_pre)
    }else if(t_pre < st & t_cur > en){
      # covers the T1 period: some T0 + entire T1 + some T0
      net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + t(c_types) * (st - t_pre + t_cur - en)
      net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + t(d_types) * (st - t_pre + t_cur - en)
      
      net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + t(c_types) * (en - st)
      net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + t(d_types) * (en - st)
    }else{
      # straddling T0 and T1
      # both need update
      
      ## get change point first:
      cp = ifelse(t_cur < en, st, en)
      
      if(t_cur < en){
        # T0 first, then T1
        net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + 
          t(c_types) * (cp - t_pre)
        net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + 
          t(d_types) * (cp - t_pre)
        
        net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + 
          t(c_types) * (t_cur - cp)
        net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + 
          t(d_types) * (t_cur - cp)
      }else{
        # T1 first, then T0
        net_c_table[c_inds,5:7] = net_c_table[c_inds,5:7] + 
          t(c_types) * (cp - t_pre)
        net_d_table[d_inds,5:7] = net_d_table[d_inds,5:7] + 
          t(d_types) * (cp - t_pre)
        
        net_c_table[c_inds,2:4] = net_c_table[c_inds,2:4] + 
          t(c_types) * (t_cur - cp)
        net_d_table[d_inds,2:4] = net_d_table[d_inds,2:4] + 
          t(d_types) * (t_cur - cp)
      }
    }
    
    # then process changes to the system
    if (z==1){
      ## CHANGED 09/21/2020
      # exposure
      p1 = events$per1[r]
      epid[p1] = -1
      
      # record num of Ia and Is friends at time of infection
      epi_table[p1,3] = sum(adjmat[p1,epid==1])
      epi_table[p1,4] = sum(adjmat[p1,epid==2])
      
      # output some info
      cat('when ',p1,' got infected, had ', sum(adjmat[p1,epid==1]),
          'Ia contacts and ', sum(adjmat[p1,epid==2]), 'Is contacts.\n')
      
      # also take down exposure time
      epi_table[p1,6] = -1 * t_cur
    }else if (z %in% c(9,10)){
      # manifestation: becoming Ia or Is
      p1 = events$per1[r]
      epid[p1] = ifelse(z==9, 1, 2)
      
      # calculate total latent time
      epi_table[p1,6] = epi_table[p1,6] + t_cur
      
      # also take down infection time
      epi_table[p1,5] = -1 * t_cur
    }else if (z==2){
      # recovery
      p1 = events$per1[r]
      epid[p1] = -2 # changed coding -2=R
      
      # calculate total infection time
      epi_table[p1,5] = epi_table[p1,5] + t_cur
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
  
  # if at the end, someone is still sick...
  if(any(epid %in% c(1,2))){
    still_sick = which(epid %in% c(1,2))
    epi_table[still_sick,5] = epi_table[still_sick,5] + t_cur
  }
  
  # AND if at the end, someone is still latent...
  if(any(epid == -1)){
    still_latent = which(epid == -1)
    epi_table[still_latent,6] = epi_table[still_latent,6] + t_cur
  }
  
  # make the epi network tables into dataframes
  epi_table = as.data.frame(epi_table)
  names(epi_table) = c('Ia_expo','Is_expo','Ia_ti','Is_ti',
                       'sick_time','latent_time')
  net_c_table = as.data.frame(net_c_table)
  names(net_c_table) = c('Nd',type_phase)
  net_d_table = as.data.frame(net_d_table)
  names(net_d_table) = c('Nc',type_phase)
    
  # change the name a bit...
  # net_c_table --> net_d_table (# of connections and time spent disconnected)
  # and same for net_d_table
  return(list(counts = c(nE=nE, nIa=nIa, nIs=nIs, nR=nR),
              C_all = C_all, D_all = D_all,
              epi_table = epi_table, net_c_table = net_d_table,
              net_d_table = net_c_table, I0=I0))
}

### try it out
# dat = readRDS('~/Documents/Research_and_References/EpiNet_Codes_ImmediateData/ex3dat_1.rds')
# 
# summ = summarize_events(dat$G0, dat$I0, dat$events, 10)

## 09/21/2020
# try it out
events = read_csv('hetero_ex_dat.csv')
names(events) = c('time','event','per1','per2')
I0 = 39 # directly pulled from Python console
## also need to re-label people with 1-indexing (rather than 0-indexing)
events$per1 = events$per1 + 1
events$per2 = events$per2 + 1

G0 = as.matrix(read_delim('hetero_ex_G0.txt', delim=" ", col_names = FALSE))
colnames(G0) = NULL
X = as.matrix(read_delim('hetero_ex_X.txt', delim=" ", col_names = FALSE))
colnames(X) = NULL

# the slow way...
summ = summarize_events(G0, I0, events, c(5,30))

# it's really, really, really slow!!!

# 09/27/2020
# try to profile the code
# library(profvis)
# profvis({summ2 = summarize_events(G0, I0, events, c(5,30))})

# a lot of time spent in "apply"??? Summing over adjmat takes a lot of time?? 
# but it's somewhat inevitable


# attempt at parallel...

# a helper function to get the "link type"
# used to index the c_ij or d_ij vector for i,j network events
# given sub-vector of epid[i,j] and "stage" (binary: 0 or 1)
get_sub_type <- function(epid){
  res = numeric(3)
  
  total = sum(epid)
  if(total == 0){
    # both H
    res[1] = 1
  }else if(total == 1){
    # H and I
    res[2] = 1
  }else{
    # both I
    res[3] = 1
  }
  
  # fin = numeric(7)
  # if(stage==1){
  #   fin[2:4] = res
  # }else{
  #   fin[5:7] = res
  # }
  # 
  # fin
  
  res
}

# function to process info for i,j pairs
summarize_ij <- function(i,j, G0_ij, I0, events, stage_change){
  # tmax
  tmax = max(events$time)
  
  # select all events related to i,j
  events = events %>% filter((per1 == i & per2==j) | (per1 == j & per2 == i) | 
                               (per1==i & is.na(per2)) | (per1==j & is.na(per2)))
  
  # get st and en of NPI period
  st = stage_change[1]; en = stage_change[2]
  
  # if there is no event related to i or j...
  # i,j have stayed H
  # and G_ij == G0_ij throughout
  if(nrow(events)==0){
    if(G0_ij==1){
      # they've stayed connected
      c_ij = rep(0,7)
      d_ij = rep(0,7)
      d_ij[2] = (st-0) + (tmax-en)
      d_ij[5] = en - st
    }else{
      # they've stayed disconnected
      d_ij = rep(0,7)
      c_ij = rep(0,7)
      c_ij[2] = (st-0) + (tmax-en)
      c_ij[5] = en - st
    }
  }else{
    # else: something happened for them
    epid = numeric(2)
    G_ij = G0_ij
    
    # # separate two stages
    # events_T0 = events %>% filter(time <= st | time >= en)
    # events_T1 = events %>% filter(time > st & time < en)
    
    c_ij = rep(0,7)
    d_ij = rep(0,7)
    
    t_pre = 0
    
    # !
    # attach a final fake row 
    # to account for the dwell time till the end
    fake_row  = c(tmax, 666, -1, -1)
    events = rbind(events, fake_row)
    
    # go through events one by one
    for(r in 1:nrow(events)){
      
      # output something
      
      #print(c_ij)
      #print(d_ij)
      
      z = events$event[r]
      t_cur = events$time[r]
      
      link_type = get_sub_type(epid)
      # accumulate risks (dwell times)
      if(t_cur < st | t_pre > en){
        # operate under T0
        if(G_ij == 0){
          # spent disconnected
          c_ij[2:4] = c_ij[2:4] + link_type * (t_cur - t_pre)
        }else{
          # spent connected
          d_ij[2:4] = d_ij[2:4] + link_type * (t_cur - t_pre)
        }
      }else if(t_pre > st & t_cur < en){
        # operate under T1
        if(G_ij == 0){
          # spent disconnected
          c_ij[5:7] = c_ij[5:7] + link_type * (t_cur - t_pre)
        }else{
          # spent connected
          d_ij[5:7] = d_ij[5:7] + link_type * (t_cur - t_pre)
        }
      }else if(t_pre < st & t_cur > en){
        # covers the T1 period: some T0 + entire T1 + some T0
        if(G_ij == 0){
          # spent disconnected
          c_ij[2:4] = c_ij[2:4] + link_type * (st - t_pre + t_cur - en)
          c_ij[5:7] = c_ij[5:7] + link_type * (en - st)
        }else{
          # spent connected
          d_ij[2:4] = d_ij[2:4] + link_type * (st - t_pre + t_cur - en)
          d_ij[5:7] = d_ij[5:7] + link_type * (en - st)
        }
      }else{
        # straddling T0 and T1
        # both T0 and T1 need update

        ## get change point first:
        cp = ifelse(t_cur < en, st, en)
        
        #cat('straddle case: change point =', cp, 't_pre = ', t_pre, 't_cur = ', t_cur, '\n')
        
        if(G_ij == 0){
          # spent disconnected
          if(cp==st){
            c_ij[2:4] = c_ij[2:4] + link_type * (cp - t_pre)
            c_ij[5:7] = c_ij[5:7] + link_type * (t_cur - cp)
          }else{
            c_ij[5:7] = c_ij[5:7] + link_type * (cp - t_pre)
            c_ij[2:4] = c_ij[2:4] + link_type * (t_cur - cp)
          }
        }else{
          # spent connected
          if(cp==st){
            d_ij[2:4] = d_ij[2:4] + link_type * (cp - t_pre)
            d_ij[5:7] = d_ij[5:7] + link_type * (t_cur - cp)
          }else{
            d_ij[5:7] = d_ij[5:7] + link_type * (cp - t_pre)
            d_ij[2:4] = d_ij[2:4] + link_type * (t_cur - cp)
          }
        }
      }
      
      if(z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        per1 = events$per1[r]
        if(per1==i){
          epid[1] = 1
        }else if(per1==j){
          epid[2] = 1
        }
      }else if(z == 2){
        # recovery: becoming H again
        per1 = events$per1[r]
        if(per1==i){
          epid[1] = 0
        }else if(per1==j){
          epid[2] = 0
        }
      }else if(z %in% c(3:5)){
        # connection
        G_ij = 1
        c_ij[1] = c_ij[1] + 1
      }else if(z %in% c(6:8)){
        # disconnection
        G_ij = 0
        d_ij[1] = d_ij[1] + 1
      }
      
      t_pre = t_cur
    }
    
    
    # # finally, deal with the status through the end
    # if(t_cur < st){
    #   # last event before start of T1: some T0 and then T1 (st to en) and then T0
    # }else if(t_cur > en){
    #   # last event after end of T1: some T0 (t_cur)
    # }else{
    #   # last event during T1: some T1 (t_cur to en) and then T0
    # }
    
  }
  
  # return as a longgg vector
  c(c_ij, d_ij)
}


# try it out
# i=39,j=43
i=37;j=39
G0_ij = G0[i,j]
summ_ij = summarize_ij(i, j, G0_ij, I0, events, c(5,30))
summ_ij

# check it against "ground truth"
summ$net_c_table[get_vec_index(c(i,j),N),]
summ$net_d_table[get_vec_index(c(i,j),N),]

# they are the same!! YEAH!


# combine this with "foreach"

## a helper function to generate i, j index sequences
get_ij_seq <- function(N){
  num.E = N * (N-1)/2
  I = numeric(num.E)
  J = numeric(num.E)
  
  counter = 0
  for(i in c(1:(N-1))){
    # append: (i, i+1), (i,i+2), ..., (i, N)
    re = N - i
    I[(counter+1):(counter+re)] = rep(i, re)
    J[(counter+1):(counter+re)] = c((i+1):N)
    
    counter = counter + re
  }
  
  list(I=I, J=J)
}


# try with foreach
library(foreach)
Edges = get_ij_seq(N)
I = Edges$I; J = Edges$J

net_tables = foreach(i=I, j=J, .combine = 'rbind') %do% {
  G0_ij = G0[i,j]
  summarize_ij(i, j, G0_ij, I0, events, c(5,30))
}


# function to summarize epi data for person i
summarize_epi <- function(i, G0_i, I0, events){
  
  # init epid vector to record disease status
  N = length(G0_i)
  epid = rep(0,N)
  epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G0_i
  # make sure (i,i) entry is 0
  G0_i[i] = 0

  # data storage for 
  # individual epidemic info
  epi_tab = numeric(6)
  names(epi_tab) = c("Ia_expo", "Is_expo", "Ia_ti", "Is_ti",
                     "sick_time", "latent_time")
  
  # get tmax
  tmax = max(events$time)
  
  # epidemic events for i (E, Ia/Is, R times)
  events_i = events %>% filter(per1 == i, event %in% c(1,9,10,2))
  # calculate sick time and latent time if i was ever infected
  if(nrow(events_i)>0){
    #expo_time = events_i$time[1]
    
    # get latent time
    ill_ind = which(events_i$event %in% c(9,10))
    if(length(ill_ind) > 0){
      epi_tab["latent_time"] = events_i$time[ill_ind] - events_i$time[1]
      # then get ill time
      recov_ind = which(events_i$event==2)
      if(length(recov_ind) > 0){
        epi_tab['sick_time'] = events_i$time[recov_ind] - events_i$time[ill_ind]
      }else{
        epi_tab['sick_time'] = tmax - events_i$time[ill_ind]
      }
    }else{
      epi_tab["latent_time"] = tmax - events_i$time[1]
    }
    
    # set "tmax" to exposed time if ever infected...
    tmax = events_i$time[1]
  }
  
  # get all events that are related to i (both epi and net)
  # AND all the manisfestation and recovery events (event %in% c(2,9,10))
  # ONLY need events that are before Exposure time (if no exposure then tmax)
  events = events %>% filter(time <= tmax) %>% 
    filter(per1 == i | per2 == i | event %in% c(2,9,10))
  
  # if there is nothing, then NO exposure for i (nobody ever manisfested)
  if(nrow(events) > 0){
    # go through the events one by one
    t_pre = 0
    for(r in 1:nrow(events)){
      z = events$event[r]
      t_cur = events$time[r]
      
      # obtain cumulated time of exposure up to now
      epi_tab['Ia_expo'] = epi_tab['Ia_expo'] + sum(nei * (epid==1)) * (t_cur - t_pre)
      epi_tab['Is_expo'] = epi_tab['Is_expo'] + sum(nei * (epid==2)) * (t_cur - t_pre)
      
      # then process changes to the system
      if (z==1){
        # exposure
        p1 = events$per1[r]
        epid[p1] = -1
        
        # get local neighborhood at time of exposure
        if(p1 == i){
          epi_tab['Ia_ti'] = sum(nei * (epid==1))
          epi_tab['Is_ti'] = sum(nei * (epid==2))
          
          # some info
          cat('when ',p1,' got infected, had ', epi_tab['Ia_ti'],
              'Ia contacts and ', epi_tab['Is_ti'], 'Is contacts.\n')
        }
      }else if (z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        p1 = events$per1[r]
        epid[p1] = ifelse(z==9, 1, 2)
      }else if (z==2){
        # recovery
        p1 = events$per1[r]
        epid[p1] = -2 # changed coding -2=R
      }else{
        # some edge stuff
        p1 = events$per1[r]
        p2 = events$per2[r]
        
        if(z %in% c(3:5)){
          # reconnection
          if(p1 == i){
            nei[p2] = 1
          }else if(p2 == i){
            nei[p1] = 1
          }
        }else if(z %in% c(6:8)){
          # disconnection
          if(p1 == i){
            nei[p2] = 0
          }else if(p2 == i){
            nei[p1] = 0
          }
        }
      }
      # set t_pre to t_cur
      t_pre = t_cur
    }
  }

  epi_tab
}


# try it out
i = 23
G0_i = G0[i,]
(summ_i = summarize_epi(i, G0_i, I0, events))

# compare with previously got results
summ$epi_table[i,]
# looks the same, YEAH!

# try it with foreach
# and do it together with the network events

# 09/27/2020
# do it on a bigger dataset (N=200 instead of N=50)
## 09/21/2020
# try it out
events = read_csv('hetero_ex2_dat.csv')
names(events) = c('time','event','per1','per2')
I0 = 102 # directly pulled from Python console
## also need to re-label people with 1-indexing (rather than 0-indexing)
events$per1 = events$per1 + 1
events$per2 = events$per2 + 1

G0 = as.matrix(read_delim('hetero_ex2_G0.txt', delim=" ", col_names = FALSE))
colnames(G0) = NULL
X = as.matrix(read_delim('hetero_ex2_X.txt', delim=" ", col_names = FALSE))
colnames(X) = NULL

# let's go!!
# now it takes tens of seconds to parse 4654 events for N=200
# (the pause is noticeable but not very bad...)

N = nrow(G0)
IJ = get_ij_seq(N)
I = IJ$I; J = IJ$J

epi_tables = foreach(i=1:N, .combine = 'rbind') %do% {
  G0_i = G0[i,]
  summarize_epi(i, G0_i, I0, events)
}
net_tables = foreach(i=I, j=J, .combine = 'rbind') %do% {
  G0_ij = G0[i,j]
  summarize_ij(i, j, G0_ij, I0, events, c(5,30))
}

# compare with the previously slow way....
# it's indeed slow as a nightmare!!
summ2 = summarize_events(G0, I0, events, c(5,30))


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
solve_MLE <- function(summaries, X, maxIter=10, tol=1e-4, initEta = 1){
  
  # X has to be a dataframe with column names!!
  
  # wrangle X into needed shape of Xij (for network part)
  # (original X: n by p design matrix)
  # (want XX: n*(n-1)/2 by p design matrix)
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
  
  # 0. get manifestation rate phi and symp. prob. p_s
  phi = (nIs + nIa)/sum(epi_table$latent_time)
  p_s = nIs/(nIs + nIa)
  
  # 1. get recovery rate gamma
  gamma = nR/sum(epi_table$sick_time)
  
  # 2. get parameters for the infection side: beta, eta, b_S
  eta = initEta
  beta = nI/(sum(epi_table$Ia_expo) + eta * sum(epi_table$Is_expo))
  
  # get rid of I0: no exposure on I0 at all
  got_infected = which(epi_table$latent_time > 0)
  got_infected = got_infected[got_infected != I0]
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  
  for(it in 1:maxIter){
    # 2.1 estimate b_S
    Expo = (epi_table$Ia_expo + eta * epi_table$Is_expo) * beta
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
      beta * sum(indiv_effects * (epi_table$Ia_expo + eta * epi_table$Is_expo)[has_expo]) -
        sum(log((epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected]))
    }
    llgrr = function(beka){
      eta * beta * sum(indiv_effects * epi_table$Is_expo[has_expo]) -
        eta * sum(epi_table$Is_ti[got_infected]/(epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected])
    }
    eta_old = eta
    
    ## output some info
    #cat('beta=',beta, 'eta=', eta, '\nfn(eta)=', llfunc(eta),
    #    'b_S=',b_S)
    
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
    if(all(abs(alpha_old - alpha) < tol)){
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
    if(all(abs(omega_old - omega) < tol)){
      break
    }
  }
  
  
  return(list(beta=beta, eta=eta, b_S=b_S, gamma=gamma, 
              phi=phi, p_s = p_s,
              alpha=alpha, b_alpha=b_alpha, omega=omega, b_omega=b_omega))
}

# try it out
# the infection part is somewhat working...
# the network part might be working? but precision is quite low
estimates = solve_MLE(summ, X, initEta = 1.5, maxIter = 100)



