#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 12:26:50 2020

@author: fan
"""

#%%
import numpy as np
from numpy.random import exponential, choice
#from random import uniform

import networkx as nx

from copy import copy

import pandas as pd

#%%
# try out manipulating design matrix X into XX with (i,j) entry being Xi + Xj
X = np.random.randn(10,2)

#mani = np.zeros((2,10))
#mani[0,:] = np.ones(10)

mani = np.zeros((10,2,2))
for j in range(10):
    mani[j,:,:] = np.eye(2)

XX = X.dot(mani)    
XXT = np.transpose(XX,(1,0,2))

Xsum = XX + XXT

for i in range(10):
    for j in range(10):
        print(Xsum[i,j,0] == X[i,0]+X[j,0])
        
Xsum.dot([1,-1]).shape

#%%
# a function to conduct the manipulation above
def Xi2Xij(X):
    '''
    X: n by p design matrix; it has to be a 2-dim array!! Not a vector
    '''
    N,p = X.shape
    mani = np.zeros((N,p,p))
    for j in range(N):
        mani[j,:,:] = np.eye(p)
        
    XX = X.dot(mani)
    Xsum = XX + np.transpose(XX,(1,0,2))
    
    return Xsum
        
#%%
# define class of heterogeneous EpiNet
# define the simulation function
        
class HeteroEpiNet:
    def __init__(self, N, T1_bounds, X, Xdyad=None, params=None, G0=None, I0=None):
        '''
        N: population size
        X: n-by-p design matrx with indivdual covariates
        Xdyad: (n,n,p) design array for dyadic covariates (used for contact networks)
        params: dictionary with all parameters
            - "beta": infection rate
            - "phi": incubation rate (avg. length of incubation = 1/phi)
            - "p_s": prob. of getting symptomatic
            - "eta": log difference of infectiousness between Is and Ia, >0
            - "gamma": recovery rate
            - "alpha": baseline link activation rate, in the format of
                        [(alphaHH0, alphaHI0, alphaII0), (alphaHH1, alphaHI1, alphaII1)]
            - "omega": baseline link termination rate, in the format of
                        [(omegaHH0, omegaHI0, omegaII0), (omegaHH1, omegaHI1, omegaII1)]
            - "b_S": length-p coefficients on susceptibility
            - "b_alpha": length-p coefficients on link activation
            - "b_omega": length-p coefficients on link termination 
        T1_bounds: a tuple/list of two entries, specifying the starting and end point of intervention
            (st,en), where 0 <= st <= en <= np.inf
        G0: a numpy array representing the adjacency matrix of the initial network
        I0: a list/array of individuals already infected and symptomatic at the beginning
            ---> individuals already EXPOSED but not yet infectious (update 09/20/2020)
        '''
        if N!=X.shape[0]:
            raise ValueError("Population size N not equal to row number of X! Try again.")
            
        # setup
        self.N = N
        # the design matrix part
        self.X = X
        self.XX = Xi2Xij(X)
        # include dyadic covariates (for contact network dynamics)
        if Xdyad is not None:
            self.XXnet = np.append(self.XX, Xdyad, axis=2)
        else:
            self.XXnet = self.XX
        # the time interval of T1
        self.T1st, self.T1en = T1_bounds
        # parameters
        if params is not None:
            self.beta = params['beta']
            self.phi = params['phi']
            self.p_s = params['p_s']
            self.eta = params['eta']
            self.expEta = np.exp(params['eta']) # e^{eta}, for convenience
            self.gamma = params['gamma']
            self.alpha = params['alpha']
            self.omega = params['omega']
            self.b_S = params['b_S']
            self.b_alpha = params['b_alpha']
            self.b_omega = params['b_omega']
            #self.tmin = params['tmin']
            #self.tmax = params['tmax']
        else:
            self.beta = None
            self.phi = None
            self.p_s = None
            self.eta = None
            self.expEta = None # e^{eta}, for convenience
            self.gamma = None
            self.alpha = None
            self.omega = None
            self.b_S = None
            self.b_alpha = None
            self.b_omega = None
            self.tmin = None
            self.tmax = None
        # the network
        self.adjmat = G0
        # the epid vector
        self.epid = np.zeros(N)
        self.I0 = None
        if I0 is not None:
            self.I0 = I0
            #self.epid[I0] = 2 # "2" = infected and symptomatic
            self.epid[I0] = -1 # "-1" = exposed but not yet infectious
            
    
    # simulation function
    def simulate(self, T=100, p0 = None, G0 = None, I0 = None, params=None, initInfec=1, seed = 42, verbose=True):
        '''
        T: maximum simulation length
        p0: initial network density for Erdos-Renyi graph
        G0: initial network adjacency matrix, as given (will)
        (ONE of p0 and G0 must be not None!)
        I0: label(s) of initially ill people, if None then sample "initInfec" many to start with
        params: dictionary of parameters
        
        NOTE: 
            1. Code for status: (UPDATED 09/20/2020, different from before!)
                S: 0, E: -1, Ia: 1, Is: 2, R: -2
            2. Still haven't dealt with dynamic (time-varying) covariates -- too complicated
        '''
        
        np.random.seed(seed)
        
        # the initial network
        if G0 is None and self.adjmat is None:
            G0 = nx.erdos_renyi_graph(self.N, p0)
            self.adjmat = nx.to_numpy_array(G0)
        elif G0 is not None:
            self.adjmat = G0
        
        # initial infections
        # (set initial infections as EXPOSED: -1)
        if I0 is None and self.I0 is None:
            self.I0 = np.random.randint(0, self.N, initInfec)
            self.epid[self.I0] = -1
        elif I0 is not None:
            self.I0 = I0
            self.epid[self.I0] = -1
            
        # parameters (if specified here then override)
        if params is not None:
            self.beta = params['beta']
            self.phi = params['phi']
            self.p_s = params['p_s']
            self.eta = params['eta']
            self.expEta = np.exp(params['eta']) # e^{eta}, for convenience
            self.gamma = params['gamma']
            self.alpha = params['alpha']
            self.omega = params['omega']
            self.b_S = params['b_S']
            self.b_alpha = params['b_alpha']
            self.b_omega = params['b_omega']
        
        # set up fixed quantities
        ## infection (a column vector)
        infec_vec = self.beta * np.exp(self.X.dot(self.b_S)).reshape((self.N,1))
        ## link activation (a matrix)
        activation = np.exp(self.XXnet.dot(self.b_alpha))
        ## link termination (a matrix)
        termination = np.exp(self.XXnet.dot(self.b_omega))
        
        ## DELETED: 09/20/2020
        ## storage array for symptom onset times (from E to I)
        ## first column: time, second column: person id
        # onset_times = np.array([T+1,-1]).reshape((1,2))
        
        # storage dictionary for all the events
        event_log = {'time': [], 'event': [], 'p1': [], 'p2': []}
        
        # start simulation
        t_cur = 0
        while t_cur < T:
            ## CHANGED: 09/20/2020
            # 0. take care of symptom onset times
#            rows = np.where(onset_times[:,0] <= t_cur)[0]
#            #print(rows)
##            if rows.size > 0:
##                print(rows)          
#            for r in rows:
#                ts, ps = onset_times[r,:]
#                ps = int(ps)
#                self.epid[ps] = 2
#                
#                if verbose:
#                    print('At time {}, {} gets symptomatic!'.format(ts, ps))
#                
#                ## log this event
#                event_log['time'].append(ts)
#                event_log['event'].append(9)
#                event_log['p1'].append(ps)
#                event_log['p2'].append(np.nan)
#                
#            if rows.size > 0:
#                #onset_times = onset_times[~rows,:]
#                to_keep = [i for i in range(len(onset_times)) if i not in rows]
#                onset_times = onset_times[to_keep,:]
        
            
            # 1. compute rates
            # UPDATED: 09/20/2020
            ## 1.0 infection (exposure)
            I = np.zeros(self.N)
            I[self.epid==1] = 1; I[self.epid==2] = self.expEta
            S = np.zeros(self.N); S[self.epid==0] = 1
            infec_rate = infec_vec * S.reshape((self.N,1)) * self.adjmat * I
            tot_infec_rate = infec_rate.sum()
            
            
            ## 1.1 manifestation (--> Ia and Is)
            tot_manifest_rate = np.sum(self.epid==-1) * self.phi
            
            ## 1.2 recovery
            ## CHANGE 09/20/2020: recover asymptomatic AND symptomatic
            ## (assume Ia and Is recover at the same rate)
            tot_recover_rate = np.sum(self.epid>=1) * self.gamma
            
            ## deal with the T0 and T1 stages
            if self.T1st <= t_cur and t_cur < self.T1en:
                alpha = self.alpha[1]
                omega = self.omega[1]
            else:
                alpha = self.alpha[0]
                omega = self.omega[0]
                
            ## 1.3 link activation
            ## (assume that asymptomatic people behave like healthy!!)
            discon_type = 1-self.adjmat
            
            #healthy = np.where(self.epid<=1)[0]
            #ill = np.where(self.epid==2)[0]
            
            # CHANGED 09/20/2020: healthy = {S,E,R}, ill = {Ia, Is}
            healthy = np.where(self.epid<1)[0]
            ill = np.where(self.epid>=1)[0]
            
            if healthy.size > 0:
                
#                rr, cc = np.meshgrid(healthy, healthy)
#                rr = rr.reshape((rr.size,))
#                cc = cc.reshape((cc.size,))
#                discon_type[rr,cc] = discon_type[rr,cc] * alpha[0]
                
                #discon_type[self.epid<=1,self.epid<=1] = discon_type[[self.epid<=1,self.epid<=1] * alpha[0]
                
                # Modified 09/20/2020
                discon_type[np.ix_(healthy, healthy)] = discon_type[np.ix_(healthy, healthy)] * alpha[0]
                
                if ill.size > 0:
#                    rr, cc = np.meshgrid(healthy, ill)
#                    rr = rr.reshape((rr.size,))
#                    cc = cc.reshape((cc.size,))
#                    discon_type[rr,cc] = discon_type[rr,cc] * alpha[1]
#                    discon_type[cc,rr] = discon_type[cc,rr] * alpha[1]
                    #discon_type[self.epid<=1,self.epid==2] = discon_type[self.epid<=1,self.epid==2] * alpha[1]
                    #discon_type[self.epid==2,self.epid<=1] = discon_type[self.epid==2,self.epid<=1] * alpha[1]
                    
                    # Modified 09/20/2020
                    discon_type[np.ix_(healthy, ill)] = discon_type[np.ix_(healthy, ill)] * alpha[1]
                    discon_type[np.ix_(ill, healthy)] = discon_type[np.ix_(ill, healthy)] * alpha[1]
                    
            if ill.size > 0:
#                rr, cc = np.meshgrid(ill, ill)
#                rr = rr.reshape((rr.size,))
#                cc = cc.reshape((cc.size,))
#                discon_type[rr,cc] = discon_type[rr,cc] * alpha[2]
                #discon_type[self.epid==2,self.epid==2] = discon_type[self.epid==2,self.epid==2] * alpha[2]
            
                # Modified 09/20/2020
                discon_type[np.ix_(ill, ill)] = discon_type[np.ix_(ill, ill)] * alpha[2]

            
            acti_rate = activation * discon_type
            tot_acti_rate = acti_rate.sum()/2
            
            ## 1.4 link termination
            ## (assume that asymptomatic people behave like healthy!!)
            con_type = copy(self.adjmat)
            if healthy.size > 0:
#                rr, cc = np.meshgrid(healthy, healthy)
#                rr = rr.reshape((rr.size,))
#                cc = cc.reshape((cc.size,))
#                con_type[rr,cc] = con_type[rr,cc] * omega[0]
                
                # Modified 09/20/2020
                con_type[np.ix_(healthy, healthy)] = con_type[np.ix_(healthy, healthy)] * omega[0]

                if ill.size > 0:
#                    rr, cc = np.meshgrid(healthy, ill)
#                    rr = rr.reshape((rr.size,))
#                    cc = cc.reshape((cc.size,))
#                    con_type[rr,cc] = con_type[rr,cc] * omega[1]
#                    con_type[cc,rr] = con_type[cc,rr] * omega[1]
                    
                    
                    # Modified 09/20/2020
                    con_type[np.ix_(healthy, ill)] = con_type[np.ix_(healthy, ill)] * omega[1]
                    con_type[np.ix_(ill, healthy)] = con_type[np.ix_(ill, healthy)] * omega[1]
                    
            if ill.size > 0:
#                rr, cc = np.meshgrid(ill, ill)
#                rr = rr.reshape((rr.size,))
#                cc = cc.reshape((cc.size,))
#                con_type[rr,cc] = con_type[rr,cc] * omega[2]
#                
                # Modified 09/20/2020
                con_type[np.ix_(ill, ill)] = con_type[np.ix_(ill, ill)] * omega[2]
                
            termi_rate = termination * con_type
            tot_termi_rate = termi_rate.sum()/2
            
            # 2. sample next event time
#            tot_rate = tot_infec_rate + tot_recover_rate + tot_acti_rate + tot_termi_rate
#            t_delta = exponential(scale=1/tot_rate)
            
            # CHANGE: add manifestation 09/20/2020
            tot_rate = tot_infec_rate + tot_recover_rate + tot_acti_rate + tot_termi_rate + tot_manifest_rate
            t_delta = exponential(scale=1/tot_rate)
            
            t_next = t_cur + t_delta
            
            # 3. sample next event type
            # CHANGE: add manifestation 09/20/2020
            probs = np.array([tot_infec_rate, tot_recover_rate, tot_acti_rate, 
                              tot_termi_rate, tot_manifest_rate])/tot_rate
            event_type = choice(range(5), p=probs)
            
            if event_type == 0:
                # infection (EXPOSURE):
                # CHANGE 09/20/2020: label E as -1 instead
                indices = np.nonzero(infec_rate)
                pair_rates = np.array([infec_rate[i] for i in zip(indices[0],indices[1])])
                pair_probs = pair_rates/pair_rates.sum()
                pair = choice(range(len(pair_rates)), p=pair_probs, replace=False)
                p1, p2 = indices[0][pair], indices[1][pair] # p2 infects p1
                
                self.epid[p1] = -1
                
                # DELETED 09/20/2020
#                ## also: sample symptom onset time for this person
#                ts_min, ts_max = t_next+self.tmin, t_next+self.tmax
#                ts = uniform(ts_min, ts_max)
#                
#                onset_times = np.append(onset_times, [[ts, p1]], axis=0)
                
                #print('At time {}, {} should get symptomatic!'.format(ts, p1))
                #print(onset_times)
                
                if verbose:
                    print('At time {}, {} gets exposed by {}!'.format(t_next, p1, p2))
                ## log this event
                event_log['time'].append(t_next)
                event_log['event'].append(1)
                event_log['p1'].append(p1)
                event_log['p2'].append(p2)
                
            elif event_type == 4:
                # manifestation
                p1 = choice(np.nonzero(self.epid==-1)[0], replace=False)
                p2 = np.nan
                
                # choose the Ia vs Is type for this person
                I_type = choice([1,2], p = [1-self.p_s, self.p_s])
                
                self.epid[p1] = I_type
                
                if verbose:
                    if I_type == 1:
                        print('At time {}, {} becomes asymptomatic!'.format(t_next, p1))
                    else:
                        print('At time {}, {} becomes symptomatic!'.format(t_next, p1))
                        
                ## log this event
                event_label = 9 if I_type == 1 else 10
                event_log['time'].append(t_next)
                event_log['event'].append(event_label)
                event_log['p1'].append(p1)
                event_log['p2'].append(p2)
                
            elif event_type == 1:
                # recovery:
                # CHANGE 09/20/2020: recover both Is and Ia
                p1 = choice(np.nonzero(self.epid>=1)[0], replace=False)
                p2 = np.nan
                
                # CHANGE 09/20/2020: code for R is "-2" now
                self.epid[p1] = -2
                
                if verbose:
                    print('At time {}, {} recovers!'.format(t_next, p1))
                ## log this event
                event_log['time'].append(t_next)
                event_log['event'].append(2)
                event_log['p1'].append(p1)
                event_log['p2'].append(p2)
            elif event_type == 2:
                # link activation:
            
                indices = np.nonzero(acti_rate)
                pair_rates = np.array([acti_rate[i] for i in zip(indices[0],indices[1])])
                pair_probs = pair_rates/pair_rates.sum()
                pair = choice(range(len(pair_rates)), p=pair_probs, replace=False)
                p1, p2 = indices[0][pair], indices[1][pair]
                
                self.adjmat[p1,p2] = 1; self.adjmat[p2,p1] = 1

                
                if verbose:
                    print('At time {}, {} and {} activates their link!'.format(t_next, p1, p2))
                    
                ## log this event
                # CHANGE 09/20/2020: healthy with label <= 0, ill with label >= 1
                event_log['time'].append(t_next)
                if self.epid[p1]<=0 and self.epid[p2]<=0:
                    # HH link
                    event_log['event'].append(3)
                elif self.epid[p1]>=1 and self.epid[p2]>=1:
                    ## II link
                    event_log['event'].append(5)
                else:
                    ## HI/IH link
                    event_log['event'].append(4)
                event_log['p1'].append(p1)
                event_log['p2'].append(p2)
            else:
                # link termination
                indices = np.nonzero(termi_rate)
                pair_rates = np.array([termi_rate[i] for i in zip(indices[0],indices[1])])
                pair_probs = pair_rates/pair_rates.sum()
                pair = choice(range(len(pair_rates)), p=pair_probs, replace=False)
                p1, p2 = indices[0][pair], indices[1][pair]
                
                self.adjmat[p1,p2] = 0; self.adjmat[p2,p1] = 0
                
                if verbose:
                    print('At time {}, {} and {} terminates their link!'.format(t_next, p1, p2))
                
                ## log this event
                # CHANGE 09/20/2020: healthy with label <= 0, ill with label >= 1
                event_log['time'].append(t_next)
                if self.epid[p1]<=0 and self.epid[p2]<=0:
                    # HH link
                    event_log['event'].append(6)
                elif self.epid[p1]>=1 and self.epid[p2]>=1:
                    ## II link
                    event_log['event'].append(8)
                else:
                    ## HI/IH link
                    event_log['event'].append(7)
                event_log['p1'].append(p1)
                event_log['p2'].append(p2)
            
            # 4. forward time
            t_cur = t_next
            
            # 5. stop when nobody is infectious or exposed any more
            if np.sum(self.epid==-1)==0 and np.sum(self.epid>=1)==0:
                print('Simulation stops at {}. No one is infectious or exposed any more.'.format(t_cur))
                
                # transform event_log into a data frame
                event_log = pd.DataFrame(event_log)
                return event_log
                
                
        # transform event_log into a data frame
        event_log = pd.DataFrame(event_log)
        
        return event_log
    
    
        
#%%
# try it out
        
#params: dictionary with all parameters
#            - "beta": infection rate
#            - "eta": log difference of infectiousness between E and I, >0
#            - "gamma": recovery rate
#            - "alpha": baseline link activation rate, in the format of
#                        [(alphaHH0, alphaHI0, alphaII0), (alphaHH1, alphaHI1, alphaII1)]
#            - "omega": baseline link termination rate, in the format of
#                        [(omegaHH0, omegaHI0, omegaII0), (omegaHH1, omegaHI1, omegaII1)]
#            - "b_S": length-p coefficients on susceptibility
#            - "b_alpha": length-p coefficients on link activation
#            - "b_omega": length-p coefficients on link termination     

        
#pa = {'beta': 0.3, 'eta': 0.2, 'gamma': 0.1, 
#      'alpha':[np.array([0.001,0.001,0.001]), np.array([0.001,0.0002,0.001])],
#      'omega':[np.array([0.005,0.005,0.005]), np.array([0.005,0.05,0.005])],
#      'b_S': np.array([0,-1]),
#      'b_alpha': np.array([1,0]),
#      'b_omega': np.array([1,-1]),
#      'tmin': 1, 'tmax': 3}  

  
# 09/20/2020: new parameter settings
#        params: dictionary with all parameters
#            - "beta": infection rate
#            - "phi": incubation rate (avg. length of incubation = 1/phi)
#            - "p_s": prob. of getting symptomatic
#            - "eta": log difference of infectiousness between Is and Ia, >0
#            - "gamma": recovery rate
#            - "alpha": baseline link activation rate, in the format of
#                        [(alphaHH0, alphaHI0, alphaII0), (alphaHH1, alphaHI1, alphaII1)]
#            - "omega": baseline link termination rate, in the format of
#                        [(omegaHH0, omegaHI0, omegaII0), (omegaHH1, omegaHI1, omegaII1)]
#            - "b_S": length-p coefficients on susceptibility
#            - "b_alpha": length-p coefficients on link activation
#            - "b_omega": length-p coefficients on link termination 

pa = {'beta': 0.3, 'eta': 0.2, 'gamma': 0.1, 
      'phi': 0.2, 'p_s': 0.6,
      'alpha':[np.array([0.001,0.001,0.001]), np.array([0.001,0.0002,0.001])],
      'omega':[np.array([0.005,0.005,0.005]), np.array([0.005,0.05,0.005])],
      'b_S': np.array([0,-1]),
      'b_alpha': np.array([1,0]),
      'b_omega': np.array([1,-1]),
      'tmin': 1, 'tmax': 3}  


N = 200

X = np.random.randint(1, size=(N,2))
phase_bounds = [5,30]

EpiNet = HeteroEpiNet(N, phase_bounds, X)

res = EpiNet.simulate(T=100, p0=0.2, params=pa, verbose=True)          
        
            
            
        
            
        
        

