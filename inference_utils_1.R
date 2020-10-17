# 10/17/2020

# utility functions for Hetero EpiNet

library(foreach)
registerDoParallel()

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
#    at time t^E
#    a vector of Ia people, and a vector of Is people

# this function: need to make sure i has manifested!! (had event 9 or 10)
# ALSO: assume the events contain recovery events (and sorted by time)
get_nei_expo_i <- function(i, expo_time, G0_i, I0, events){
  # init epid vector to record disease status
  N = length(G0_i)
  epid = rep(0,N)
  epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G0_i
  # make sure (i,i) entry is 0
  nei[i] = 0
  
  # storage vectors for Ia and Is people
  Ia = NULL
  Is = NULL
  
  # get events up to expo_time:
  # get all events that are related to i (both epi and net)
  # AND all the manisfestation and recovery events (event %in% c(2,9,10))
  events = events %>% filter(time <= expo_time) %>%
    filter(per1 == i | per2 == i | event %in% c(2,9,10))
  
  # if there is nothing, then NO exposure for i (nobody ever manisfested)
  if(nrow(events) > 0){
    # go through the events one by one
    t_pre = 0
    for(r in 1:nrow(events)){
      z = events$event[r]
      t_cur = events$time[r]
      
      # then process changes to the system
      if (z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        p1 = events$per1[r]
        epid[p1] = ifelse(z==9, 1, 2)
      }else if(z==2){
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
    
    Ia = which(nei * (epid==1) == 1)
    Is = which(nei * (epid==2) == 1)
  }
    
    list(Ia = Ia, Is = Is)
}

# ii) get the neiborhood for all people who were ever infectious
get_nei_expo_all <- function(G0, I0, events, expo_times){
  # expo_times: a list of
  ## exposed: ids of all people who need exposure time imputation
  ## times: currently imputed exposure times
  exposed = expo_times$exposed
  M = ifelse(length(exposed) > 100, length(exposed), 100)
  res = foreach(i=exposed, .combine = 'list', .maxcombine = M) %dopar% {
    G0_i = G0[i,]
    expo_time = expo_times$times[exposed == i]
    get_nei_expo_i(i, expo_time, G0_i, I0, events)
  }
  names(res) = sapply(exposed, as.character)
  res
}


# 10/17/2020: try it out
# try it out
setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
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

# get an `exposure_times` list
exposure_times = events %>% filter(event == 1) %>% 
  mutate(exposed = per1, times = time) %>%
  select(exposed, times) %>% as.list()

# try the `get_nei_expo_all` function
infec_nei = get_nei_expo_all(G0, I0, events, exposure_times)

st = Sys.time()
for(j in 1:200){
  infec_nei = get_nei_expo_all(G0, I0, events, exposure_times)
  cat(j,'\r')
}
diff = Sys.time() - st
cat('Time taken:',diff,'\n')
