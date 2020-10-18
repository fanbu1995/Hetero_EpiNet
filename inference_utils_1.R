# 10/17/2020

# utility functions for Hetero EpiNet

library(dplyr)
library(foreach)
library(Matrix)
registerDoParallel()

# 0. to prepare stuff
# (1) a function to generate samples from Exp(rate) truncated at (a,b)
sam_texp <- function(n, rate, a=0, b=Inf){
  quans = runif(n)
  sams = -(log(exp(-rate*a)-quans*(exp(-rate*a)-exp(-rate*b))))/rate
  return(sams)
}

# (2.1) function to get adjmat up to time point t
get_adjmat <- function(G0, events, t){
  adjmat = G0
  net_events = events %>% 
    filter(event %in% c(3:8), time < t)
  
  if(nrow(net_events) > 0){
    for(r in 1:nrow(net_events)){
      z = net_events$event[r]
      p1 = net_events$per1[r]
      p2 = net_events$per2[r]
      if(z %in% c(3:5)){
        adjmat[p1,p2] = 1
        adjmat[p2,p1] = 1
      }else{
        adjmat[p1,p2] = 0
        adjmat[p2,p1] = 0
      }
    }
  }
  
  adjmat
}

# (2) a function to save adjmat structure at some report time points
get_adjmat_report <- function(G0, events, times){
  
  adjmat = G0
  net_events = events %>% filter(event %in% c(3:8))
  
  res = list()
  
  if(times[1]==0){
    lb = times[1]
    res[[as.character(lb)]] = adjmat
  }
  
  if(nrow(net_events)==0){
    for(i in c(1:length(times))){
      lb = times[i]
      res[[as.character(lb)]] = adjmat
    }
  }else{
    t = 0; r = 1
    for(i in c(1:length(times))){
      lb = times[i]
      while((t<lb) & (r <= nrow(net_events))){
        t = net_events$time[r]
        z = net_events$event[r]
        p1 = net_events$per1[r]
        p2 = net_events$per2[r]
        
        if(t >= lb){
          res[[as.character(lb)]] = adjmat
        }
        
        if(z %in% c(3:5)){
          adjmat[p1,p2] = 1
          adjmat[p2,p1] = 1
        }else{
          adjmat[p1,p2] = 0
          adjmat[p2,p1] = 0
        }
        r = r+1
      }
      
      if(t < lb){
        res[[as.character(lb)]] = adjmat
      }
    }
  }
  
  res
}

# try this out
times = seq(0,56,by=7)
G_report = get_adjmat_report(G0, events, times)
# test out to make sure things are correct!
for(t in times){
  adj_t = get_adjmat(G0, events, t)
  cat(all(adj_t == G_report[[as.character(t)]]),"\n")
  cat(which(adj_t != G_report[[as.character(t)]], arr.ind = TRUE),'\n')
}

# events %>% filter(event %in% c(3:8), 
#                   (per1==161 & per2==145)|(per1==145 & per2==161))


# (3) a function to generate missing exposure &/ recovery times
#     and produce status reports

miss_data <- function(events, G0, I0, interval = 7, 
                      miss_recov_prop = 1, miss_expo_prop = 1){
  # `interval`: the regular period of status report
  # `miss_recov_prop`: the proportion/probability of missing recovery times
  # `miss_expo_prop`: the proportion/probability of missing exposure times
  
  N = nrow(G0)
  epid = rep(0,N); epid[I0] = 1
  tmax = max(events$time)
  report.times = seq(from = 0, to = (tmax %/% interval + 1) * interval, by = interval)
  report = epid
  
  # select all epi events
  # events.epi = events[events$event %in% c(1,2),]
  
  recov.ind = NULL
  expo.ind = NULL
  for(ix in c(2:length(report.times))){
    lb = report.times[ix-1]; ub = report.times[ix]
    st = min(which(events$time > lb)); en = max(which(events$time <= ub))
    if(st == Inf){ break }
    for(r in st:en){
      z = events$event[r]
      p = events$per1[r]
      if(z==1){
        # exposure
        epid[p] = -1
        expo.ind = c(expo.ind, r)
      }else if(z==2){
        # recovery
        epid[p] = -2
        recov.ind = c(recov.ind, r)
      }else if(z %in% c(9,10)){
        # manifestation
        epid[p] = ifelse(z==9, 1, 2)
      }
    }
    # status report
    report = rbind(report, epid)
  }
  
  # records to remove in recovery times and exposure times
  remove_recov = as.logical(rbinom(length(recov.ind), 1, miss_recov_prop))
  remove_expo = as.logical(rbinom(length(expo.ind), 1, miss_expo_prop))
  to_remove = c(recov.ind[remove_recov], expo.ind[remove_expo])
  events = events[-to_remove,]
  
  # return new dataset and status report
  row.names(report) = as.character(report.times)
  
  # if("truth" %in% names(dat)){
  #   return(list(G0=G0, I0=I0, report.times=report.times, 
  #               events=events, report = report, 
  #               truth = dat$truth))
  # }else{
  #   cat("True parameters unavailable for this dataset!\n")
  #   return(list(G0=G0, I0=I0, report.times=report.times, 
  #               events=events, report = report))
  # }
  
  list(G0=G0, I0=I0, report.times=report.times, 
       events=events, report = report)
}

# try this
miss_dat = miss_data(events, G0, I0)


# I. for MLEs

# (1) summarize epi table for each person i
# function to summarize epi data for person i
summarize_epi <- function(i, G0_i, I0, events){
  
  # init epid vector to record disease status
  N = length(G0_i)
  epid = rep(0,N)
  epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G0_i
  # make sure (i,i) entry is 0
  nei[i] = 0
  
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


# (2) summarize connection/disconnection history for (i,j) pair

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

# # try it out
# # i=39,j=43
# i=37;j=39
# G0_ij = G0[i,j]
# summ_ij = summarize_ij(i, j, G0_ij, I0, events, c(5,30))
# summ_ij
# 
# # check it against "ground truth"
# summ$net_c_table[get_vec_index(c(i,j),N),]
# summ$net_d_table[get_vec_index(c(i,j),N),]

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


# II. for missing data

# (1) missing recovery time

# i) given exposure times etc., get the local neighborhood of person i
#    at time t^E:
#    a vector of Ia people, and a vector of Is people

# this function: need to make sure i has manifested!! (had event 9 or 10)
get_nei_expo_i <- function(i, expo_time, G_i, events, report, times){
  # expo_time: the exposure time of i
  # all adjmats at beginning points of the intervals 
  # report: the matrix of everyone's status at each report point
  # times: the report time points
  
  # the left endpoint of interval that contains expo_time
  lb = max(times[times <= expo_time])
  ix = max(which(times <= expo_time))
  
  cat('lower bound for expo time', expo_time, 'is', lb,'\n')

  # init epid vector to record disease status
  # only need to start from the LB of the interval for new manifestation events
  #G_i = G_all[[as.character(lb)]][i,]
  N = length(G_i)
  epid = report[ix,]
  #epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G_i
  # make sure (i,i) entry is 0
  nei[i] = 0
  
  # storage vectors for Ia and Is people
  # Ia = NULL
  # Is = NULL
  
  # get events up to expo_time:
  # get all events that are related to i (both epi and net)
  # AND all the manisfestation events (event %in% c(9,10))
  # ONLY deal with events after lb (start point of the interval)
  events = events %>% 
    filter(time <= expo_time & time >= lb) %>%
    filter(per1 == i | per2 == i | event %in% c(9,10))
  
  # if there is nothing, then NO new 
  # activity related to i from lb to expo_time
  if(nrow(events) > 0){
    # go through the events one by one
    t_pre = 0
    for(r in 1:nrow(events)){
      z = events$event[r]
      t_cur = events$time[r]
      
      # then process changes to the system
      if (z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        # only deal with new infections after lb
        p1 = events$per1[r]
        epid[p1] = ifelse(z==9, 1, 2)
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
  
  Ia = which(nei * (epid==1) == 1)
  Is = which(nei * (epid==2) == 1)
  
  list(Ia = Ia, Is = Is)
}

# ii) get the neiborhood for all people who were ever infectious
get_nei_expo_all <- function(G_all, events, expo_times, report, times){
  # expo_times: a list of
  ## exposed: ids of all people who need exposure time imputation
  ## times: currently imputed exposure times
  exposed = expo_times$exposed
  M = ifelse(length(exposed) > 100, length(exposed), 100)
  res = foreach(i=exposed, .combine = 'list', .maxcombine = M) %dopar% {
    
    expo_time = expo_times$times[exposed == i]
    
    # the left endpoint of interval that contains expo_time
    lb = max(times[times <= expo_time])
    # get i's local adjmat at lb
    G_i = G_all[[as.character(lb)]][i,]

    get_nei_expo_i(i, expo_time, G_i, events, report, times)
  }
  names(res) = sapply(exposed, as.character)
  res
}


# 10/17/2020: try it out
# try it out
setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
events = read_csv('hetero_ex2_dat.csv')
names(events) = c('time','event','per1','per2')
I0 = 103 # directly pulled from Python console
## also need to re-label people with 1-indexing (rather than 0-indexing)
events$per1 = events$per1 + 1
events$per2 = events$per2 + 1

G0 = as.matrix(read_delim('hetero_ex2_G0.txt', delim=" ", col_names = FALSE))
colnames(G0) = NULL
X = as.matrix(read_delim('hetero_ex2_X.txt', delim=" ", col_names = FALSE))
colnames(X) = NULL

# get an `exposure_times` list
# exposure_times = events %>% filter(event == 1) %>% 
#   mutate(exposed = per1, times = time) %>%
#   select(exposed, times) %>% as.list()

# try the `get_nei_expo_all` function
infec_nei = get_nei_expo_all(G_report, events, exposure_times, 
                             miss_dat$report, miss_dat$report.times)

# check out a few to see if their real infector is in the set
events %>% filter(event==1, per1 %in% c(124,28,172))

# compare it with epi_tables to validate
# sanity check: potential infector set should be larger than real infector set!
exposed = exposure_times$exposed
for(e in exposed){
  Ia_ti = epi_tables[e,'Ia_ti']
  Is_ti = epi_tables[e,'Is_ti']
  size_a = length(infec_nei[as.character(e)]$Ia)
  size_s = length(infec_nei[as.character(e)]$Is)
  cat('Person',e,'Ia size:', Ia_ti>=size_a, 'Is size:', Is_ti >=size_s,'\n')
}

# benchmark it...
# st = Sys.time()
# for(j in 1:200){
#   infec_nei = get_nei_expo_all(G0, I0, events, exposure_times)
#   cat(j,'\r')
# }
# diff = Sys.time() - st
# cat('Time taken:',diff,'\n')
# 100 times take about 1 minute... (not super bad)


# iii) DARCI+: lower bound update is not uniform, adapted from before
propose_recov_filter <- function(lb, ub, recovers, exposure_times, nei_infec, 
                                 gam=0.2, exp_eta=1.5){
  
  # [lb, ub]: the time interval to impute recovery times on
  # recovers: ids of people to impute recovery times for
  # nei_infec: the list of Ia and Is neighbors at time of exposure for each person
  # exposure_times: a list of
  ## exposed: ids of all people who need exposure time imputation
  ## times: currently imputed exposure times
  # gam: currect value of gamma
  # exp_eta: current value of exp(eta)
  
  # pull up infection log in this interval 
  events.infec = as.data.frame(exposure_times) %>% 
    select(per1 = exposed, time = times) %>%
    filter(time > lb & time <= ub) %>%
    arrange(time)
  
  # if events.infec non-empty:
  # obtain adjusted feasible sampling lower bounds for each person in `recovers`
  # if about-to-recover people are the only sick neighbors of someone infected at t, 
  # then randomly select one of them to mandately recover after t
  # according to the prob. in manuscript
  bounds = rep(lb, length(recovers))
  
  if(nrow(events.infec) > 0){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      if(p1 %in% recovers){
        t = events.infec$time[r]
        bounds[recovers==p1] = max(bounds[recovers==p1],t)
      }
      nei_s = nei_infec[[as.character(p1)]]$Is; size_s = length(nei_s)
      nei_a = nei_infec[[as.character(p1)]]$Ia; size_a = length(nei_a)
      nei = c(nei_s, nei_a)
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        #cat("POIs for individual",p1,":\n",poi,"\n")
        t = events.infec$time[r]
        if(size_s == 0 | size_a == 0){
          # one of Ia/Is is empty set
          poi = nei
        }else{
          # both are non-empty
          # then choose with `exp_eta` factored in
          wset = sample(1:2, 1, prob = c(size_a, size_s * exp_eta))
          if(wset == 1){
            poi = nei_a
          }else{
            poi = nei_s
          }
        }
        if(length(poi)==1){
          p = poi
        }else{
          p = sample(poi, 1)
        }
        bounds[recovers==p] = max(bounds[recovers==p],t)
      }
    }
  }
  
  # sample recovery times under the adjusted bounds
  cands = sam_texp(length(recovers), gam, bounds-lb, ub-lb) + lb
  
  # cat("Interval: [",lb,",",ub,"]\n")
  # cat("To recover:",recovers,"\n")
  # cat("Feasible lower bounds:\n",bounds,"\n")
  # cat("Proposed recovery times:\n",cands,"\n")
  
  cands
}

## try this out
# seems to work!
# LB = 28; UB = 35
# recovs = events %>% filter(time > LB & time <= UB, event==2) %>%
#   select(per1) %>% pull()
# cands = propose_recov_filter(LB, UB, recovs, exposure_times, infec_nei)


# (2) missing exposure time

# get a list of recovers - times
recovery_times = events %>% filter(event==2) %>% 
  select(recov=per1, time)

# i) a function to combine two data frames and sort by time
# df_combine <- function(event1, event2){
#   n1 = nrow(event1)
#   n2 = nrow(event2)
#   if(n1 == 0){
#     event2
#   }else if(n2 == 0){
#     event1
#   }else{
#     n = n1 + n2
#     df = NULL
#     i = 1; j = 1
#     times1 = event1$time
#     times2 = event2$time
#     for(k in 1:n){
#       if(times1[i] < times2[j]){
#         df = rbind(df, )
#       }
#     }
#   }
# }

# ii) for each person i, get risk function change points and values
#    on interval [expo_time-t_max, expo_time-t_min]
#    given current recovery times and exp_eta

get_expo_risk_i <- function(i, t_i, G_all, tmax, tmin, report, times,
                            events, recovery_times, exp_eta=1.5){
  # t_i: manifestation time of i
  # G_all: all adjmats at report times
  # recovery_times: a sorted-by-time data frame of imputed recovery times (currently)
  
  st = max(0,t_i - tmax)
  en = ifelse(t_i - tmin <= 0, t_i, t_i - tmin)
  # the end of interval is set to t_i instead if t_i - tmin <= 0
  
  # the left endpoint of interval that contains st
  lb = max(times[times <= st])
  ix = max(which(times <= st))
  
  # init epid vector to record disease status
  # only need to start from the LB of the interval for new manifestation events
  G_i = G_all[[as.character(lb)]][i,]
  nei = G_i
  nei[i] = 0
  epid = report[ix,]
  
  #cat('length of nei:', length(nei), "length of epid", length(epid),'\n')
  
  # filter events related to i
  # AND all the manifestation events!!
  events = events %>% filter(time >= lb & time <= en) %>% 
    filter(per1 == i | per2 == i | event %in% c(9,10))
  
  recovered_before_st = recovery_times %>%
    filter(time >= lb & time <= st) %>%
    select(recov) %>% pull()
  
  # recovery_times = recovery_times %>% 
  #   filter(time > st & time <= en) %>%
  #   select(time=times, per1=recov)
  
  recovery_times = recovery_times %>% 
    filter(time > st & time <= en) %>%
    mutate(event=2, per2 = NA) %>%
    select(time, event, per1=recov, per2)
  
  # combine events with recovery_times
  events = rbind(events, recovery_times) %>% 
    arrange(time)
  
  # get the G_i and epid at the time of st
  t = lb
  r = 1
  if(st > lb){
    while((t < st) & (r <= nrow(events))){
      t = events$time[r]
      z = events$event[r]
      p1 = events$per1[r]
      p2 = events$per2[r]
      if(t >= st){
        break
      }
      if(z %in% c(3:5)){
        # reconnection
        if(p1 == i){
          nei[p2] = 1
        }else if(p2 == i){
          nei[p1] = 1
        }
      }else if (z %in% c(6:8)){
        # disconnection
        if(p1 == i){
          nei[p2] = 0
        }else if(p2 == i){
          nei[p1] = 0
        }
      }else if (z %in% c(9,10)){
        epid[p1] = ifelse(z==9, 1, 2)
      }
      r = r+1
    }
    
    # finally account for those who recovered in [lb,st]
    epid[recovered_before_st] = -2
    
  }
  
  # now go from st to en and register all changes
  change_times = st
  risk_val = NULL
  
  # if already exhausted events
  if(r==nrow(events)){
    change_times = c(change_time, en)
    risk_val = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta
  }else{
    prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta
    for(k in r:nrow(events)){
      t = events$time[k]
      z = events$event[k]
      p1 = events$per1[k]
      p2 = events$per2[k]
      
      if(z %in% c(3:5)){
        # reconnection
        poi = ifelse(p1==i, p2, p1)
        nei[poi] = 1
        # if this person is Is or Ia: register change point
        if(epid[poi] %in% c(1,2)){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z %in% c(6:8)){
        # disconnection
        poi = ifelse(p1==i, p2, p1)
        nei[poi] = 0
        # if this person is Is or Ia: register change point
        if(epid[poi] %in% c(1,2)){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z %in% c(9,10)){
        # manifestation
        epid[p1] = ifelse(z==9, 1, 2)
        # if p1 is connected to i
        if(nei[p1] == 1){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z == 2){
        # recovery
        epid[p1] = -2
        # if p1 is connected to i
        if(nei[p1] == 1){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }
    }
    
    risk_val = c(risk_val, prev_risk)
    change_times = c(change_times, en)
  }
  
  # get the lengths of constant-risk intervals
  # compute relative risk
  # and sample a time according to the probs
  lens = diff(change_times)
  nL = length(lens)
  if(nL == 1){
    res = runif(1, min=st, max=en)
  }else{
    weights = lens * risk_val
    A = sample(nL, 1, prob = weights)
    res = runif(1, min=change_times[A], max=change_times[A+1])
  }
  
  # return the sampled exposure time
  list(changes = change_times, lengths = lens, risks = risk_val, samp = res)
}

# try it out
i = 23
t_i = events %>% filter(event %in% c(9,10), per1==i) %>% 
  select(time) %>% pull()
er_i = get_expo_risk_i(i,t_i, G_report, tmax=7, tmin=1,
                       miss_dat$report, miss_dat$report.times, miss_dat$events,
                       recovery_times)

i = 77 # problem with this one!! --> probably fixed!
t_i = events %>% filter(event %in% c(9,10), per1==i) %>% 
  select(time) %>% pull()
er_i = get_expo_risk_i(i,t_i, G_report, tmax=7, tmin=0,
                       miss_dat$report, miss_dat$report.times, miss_dat$events,
                       recovery_times)
events %>% filter(event ==1 , per1==i)
events %>% filter(event %in% c(9,10), per1==65) # 17.8

