# 07/19/2020

# Simulate two SEIR extension models and compare

library(GillespieSSA)
library(tidyverse)
library(ggplot2); theme_set(theme_bw(base_size = 14))

# 1. SEIIR

sim_SEIIR <- function(beta, eta, phai, ps, gammaA, gammaS, 
                      init.states = c(S=99, E=1, Ia = 0, Is = 0, R = 0),
                      N=100, tmax=100, seed=NULL, plot=TRUE){
  
  ## define parameters
  parms = c(beta=beta, eta=eta, phai=phai, 
            ps=ps, gammaA=gammaA, gammaS=gammaS)
  tf = tmax
  simName = "SEIR with asymp. and symp. infectious"
  
  ## define initial state vector
  x0 = init.states
  
  ## define state-change matrix
  nu = matrix(c(-1, 0, 0, 0, 0,
                +1, -1, -1, 0 ,0,
                0, +1, 0, -1, 0,
                0, 0, +1, 0, -1,
                0, 0, 0, +1, +1),
              nrow = 5, byrow=TRUE)
  
  
  ## define propensity functions
  a <- c(
    "beta * S * (Ia + eta * Is)",
    "(1-ps) * phai * E",
    "ps * phai * E",
    "gammaA * Ia",
    "gammaS * Is"
  ) 
  
  # Run simulations with the Direct method
  if(!is.null(seed)){
    set.seed(seed)
  }
  out <- ssa(
    x0 = x0,
    a = a,
    nu = nu,
    parms = parms,
    tf = tf,
    method = ssa.d(),
    simName = simName,
    verbose = FALSE,
    consoleInterval = 1
  ) 
  if(plot){
    ssa.plot(out, show.title = TRUE, show.legend = TRUE)
  }
  
  return(out)
}

## try it out
res_SEIIR = sim_SEIIR(beta = 0.003, eta = 1.5, phai = 0.2, ps = 0.8, 
                      gammaA = 0.1, gammaS = 0.1, N=200, tmax = 100)



# 2. SEI->IR

sim_SEI2IR <- function(beta, eta, phai, phiS, gammaA, gammaS, 
                       init.states = c(S=99, E=1, Ia = 0, Is = 0, R = 0),
                       N=100, tmax=100, seed=NULL, plot=TRUE){
  
  ## define parameters
  parms = c(beta=beta, eta=eta, phai=phai, 
            phiS=phiS, gammaA=gammaA, gammaS=gammaS)
  tf = tmax
  simName = "SEIR with asymp. and symp. infectious"
  
  ## define initial state vector
  x0 = init.states
  
  ## define state-change matrix
  nu = matrix(c(-1, 0, 0, 0, 0,
                +1, -1, 0, 0 ,0,
                0, +1, -1, -1, 0,
                0, 0, +1, 0, -1,
                0, 0, 0, +1, +1),
              nrow = 5, byrow=TRUE)
  
  
  ## define propensity functions
  a <- c(
    "beta * S * (Ia + eta * Is)",
    "phai * E",
    "phiS * Ia",
    "gammaA * Ia",
    "gammaS * Is"
  ) 
  
  # Run simulations with the Direct method
  if(!is.null(seed)){
    set.seed(seed)
  }
  out <- ssa(
    x0 = x0,
    a = a,
    nu = nu,
    parms = parms,
    tf = tf,
    method = ssa.d(),
    simName = simName,
    verbose = FALSE,
    consoleInterval = 1
  ) 
  if(plot){
    ssa.plot(out, show.title = TRUE, show.legend = TRUE)
  }
  
  return(out)
}

## try it out
res_SEI2IR = sim_SEI2IR(beta = 0.003, eta = 1+0.5*4/3, 
                       phai = 0.2, phiS = 0.4, 
                       gammaA = 0.1, gammaS = 0.1*4/3, 
                       N=200, tmax = 200)


# 3. a function to extract "summary statistics" of the simulated epidemics
extract_summary <- function(resList){
  
  dat = resList$data
  
  ## 1. final R size (# of people who went through the disease)
  finalR = max(dat[,"R"])
  
  ## 2. total time of epidemic
  time = max(dat[,"t"])
  
  ## 3. peak epi size
  peak_Is = max(dat[,"Is"])
  peak_I = max(dat[,"Is"] + dat[,"Ia"])
  
  ## 4. average proportion of Is among all I people
  ratio_Is = mean(dat[,"Is"]/(dat[,"Is"]+dat[,"Ia"]), na.rm = T)
  
  return(c(finalR=finalR, time = time, 
           peak_I = peak_I, peak_Is=peak_Is,
           ratio_Is = ratio_Is))
}

## try it out
extract_summary(res_SEIIR)


# 4. a function that repeats the simulations and compares 
# summary statistics
# (no "trajectory comparison for now")

compare_models <- function(params_SEIIR, params_SEI2IR, N, tmax, 
                           Rep=50, plotSummary=TRUE, plotCurves=TRUE){
  # params_* should be named vectors or lists of parameter values
  
  ## 1. SEIIR
  summ_SEIIR = NULL
  if(plotCurves){
    events_SEIIR = NULL
  }
  cat("Simulating SEIIR...\n")
  for(r in 1:Rep){
    this_rep = 
      within(as.list(params_SEIIR),
             {
                        res = sim_SEIIR(beta = beta, eta = eta, 
                                       phai = phai, ps = ps, 
                                       gammaA = gammaA, gammaS = gammaS, 
                                       N=N, tmax = tmax, seed = r,
                                       plot = F)
                        summ = extract_summary(res)
                        })
    summ_SEIIR = rbind(summ_SEIIR, this_rep$summ)
    if(plotCurves & r <= 20){
      # record 20 event sequences at most!
      this_event = this_rep$res$data %>% as.data.frame()
      names(this_event) = c("t", "S", "E", "Ia", "Is", "R")
      this_event$rep = r
      this_event$model = "SEIIR"
      events_SEIIR = rbind(events_SEIIR, this_event)
    }
    cat(r,"\r")
  }
  cat("\n")
  rownames(summ_SEIIR) = NULL
  
  ## 2. SEI2IR
  summ_SEI2IR = NULL
  if(plotCurves){
    events_SEI2IR = NULL
  }
  cat("Simulating SEI2IR...\n")
  for(r in 1:Rep){
    this_rep = 
      within(as.list(params_SEI2IR),
             {
             res = sim_SEI2IR(beta = beta, eta = eta, 
                              phai = phai, phiS = phiS, 
                              gammaA = gammaA, gammaS = gammaS, 
                              N=N, tmax = tmax, seed = r,
                              plot = F)
             summ = extract_summary(res)
             })
    summ_SEI2IR = rbind(summ_SEI2IR, this_rep$summ)
    if(plotCurves & r <= 20){
      # record 20 event sequences at most!
      this_event = this_rep$res$data %>% as.data.frame()
      names(this_event) = c("t", "S", "E", "Ia", "Is", "R")
      this_event$rep = r + 20 # add a "20" to avoid repeat labels
      this_event$model = "SEI->IR"
      events_SEI2IR = rbind(events_SEI2IR, this_event)
    }
    cat(r,"\r")
  }
  cat("\n")
  rownames(summ_SEI2IR) = NULL
  
  ## 3. compare by KS test (approx. p values under ties...)
  m = ncol(summ_SEIIR)
  pvals = numeric(m)
  for(j in 1:m){
    pvals[j] = ks.test(summ_SEIIR[,j], summ_SEI2IR[,j])$p.value
  }
  names(pvals) = colnames(summ_SEIIR)
  print(pvals)
  
  ## 4. visualize summary statistics
  if(plotSummary){
    dat1 = summ_SEIIR %>% as.data.frame()
    names(dat1) = names(pvals)
    dat1$model = "SEIIR"
    dat2 = summ_SEI2IR %>% as.data.frame()
    names(dat2) = names(pvals)
    dat2$model = "SEI->IR"
    dat = rbind(dat1, dat2)
    for(n in names(pvals)){
      print(
        ggplot(data = dat, aes(x=model)) +
          geom_boxplot(aes_string(y=n)) +
          labs(title=n)+
          theme_bw(base_size=14)
      )
    }
  }
  
  ## 5. visualize trajectories
  if(plotCurves){
    all_events = rbind(events_SEIIR, events_SEI2IR)
    for(v in c('S','E','Ia','Is','R')){
      print(
        ggplot(data = all_events, aes(x=t, group=rep, color=model)) +
          geom_line(size=0.3, aes_string(y=v)) +
          labs(title=paste("Trajectory of",v)) +
          theme_bw(base_size=14)
      )
    }
  }
  
  return(list(model1 = summ_SEIIR, model2 = summ_SEI2IR,
              pvals = pvals))
}


## try it out
params1 = c(beta = 0.003, eta = 1.5, 
            phai = 0.2, ps = 0.8, 
            gammaA = 0.1, gammaS = 0.1)
params2 = c(beta = 0.003, eta = 1+0.5*4/3, 
            phai = 0.2, phiS = 0.4, 
            gammaA = 0.1, gammaS = 0.1*4/3)

pdf("Compare_SEIR_plots.pdf", height = 5, width = 9)
compRes = compare_models(params_SEIIR = params1, 
                          params_SEI2IR = params2, 
                          N = 200, tmax = 200,
                          Rep = 100)
dev.off()

saveRDS(compRes,"compare_SEIR_summmary.rds")

# compRes2 = compare_models(params_SEIIR = params1, 
#                          params_SEI2IR = params2, 
#                          N = 200, tmax = 200,
#                          Rep = 10)


# 5A. a function for discretization
# default: 100 days, daily counts
get_discrete <- function(resList, tmax=100, tau=1){
  dat = resList$data
  
  times = dat[,"t"]
  values = dat[,2:ncol(dat)]
  
  
  obs = values[1,] %>% as.numeric()
  ts = seq(from=0, to=tmax, by=tau)
  
  # if(tmax > ceiling(max(times))){
  #   ind.max = min(which(ts > max(times)))
  #   for(t in ts[2:ind.max]){
  #     ind = max(which(times <= t))
  #     y_t = values[ind,] %>% as.numeric()
  #     obs = rbind(obs, y_t)
  #   }
  # }
  
  for(t in ts){
    if(t!=0){
      ## find the event time t_i such that t_i <= t < t_{i+1}
      ind = max(which(times <= t))
      y_t = values[ind,] %>% as.numeric()
      obs = rbind(obs, y_t)
    }
  }
  
  res = cbind(ts, obs)
  rownames(res) = NULL
  res = res %>% as.data.frame()
  names(res) = colnames(dat)
  
  return(res)

}


# 5B. compare distribution of one thing against observation of another thing
compare_against <- function(x, y){
  # x: a vector
  # y: a number
  
  tail_p = mean(x < y)
  if(tail_p < 0.5){
    2 * tail_p
  }else{
    2 * (1- tail_p)
  }
}

# 5C. compare one discrete data with another
# dat1 < dat2
compare_discrete <- function(params_SEIIR, params_SEI2IR, N, tmax, 
                             Rep1=200, Rep2=1, 
                             disc_tmax=100, tau=1,
                             plotRuns = 300){
  
  # 08/01/2020: allow plotting more runs in the "background"
  
  if(Rep1 > 1){
    # use SEIIR results as reference
    cat("Using SEIIR model as reference model!\n")
    
    ## 1. sample one sequence from SEI2IR
    res_SEI2IR = with(as.list(params_SEI2IR), {
      sim_SEI2IR(beta = beta, eta = eta, 
                 phai = phai, phiS = phiS, 
                 gammaA = gammaA, gammaS = gammaS, 
                 N=N, tmax = tmax,
                 plot = F)
    })
    disc_SEI2IR = get_discrete(res_SEI2IR, disc_tmax, tau)
    cat("SEI2IR simulation done.\n")
    
    ## 2. sample a lot of sequences from SEIIR
    events_SEIIR = NULL
    
    for(r in 1:Rep1){
      ### A. simulate
      this_res = with(as.list(params_SEIIR), {
        sim_SEIIR(beta = beta, eta = eta, 
                  phai = phai, ps = ps, 
                  gammaA = gammaA, gammaS = gammaS, 
                  N=N, tmax = tmax, seed = r,
                  plot = F)
      })
      ### B. discretize
      this_disc = get_discrete(this_res, disc_tmax, tau)
      
      ### C. compare discrete data
      if(r == 1){
        compare = disc_SEI2IR < this_disc
      }else{
        compare = compare + (disc_SEI2IR < this_disc)
      }
      
      ### D. save events for plotting
      if(r <= plotRuns){
        # record 100 event sequences at most!
        events_SEIIR = rbind(events_SEIIR, this_disc)
      }
      cat(r,"\r")
    }
    cat("\n")
    
    ## post-process "compare"
    compare = compare[,-1]/Rep1
    for(c in 1:ncol(compare)){
      compare[,c] = sapply(compare[,c], function(p) ifelse(p<0.5, 2*p, 2*(1-p)))
    }
    
    Names = colnames(compare)
    compare = compare %>% as.data.frame()
    names(compare) = Names
    compare$t = disc_SEI2IR$t
    
    return(list(event1 = disc_SEI2IR, event2 = events_SEIIR, 
                compare = compare, reference = "SEIIR"))
  }else{
    # use SEI2IR results as reference
    cat("Using SEI2IR model as reference model!\n")
    
    ## 1. sample one sequence from SEIIR
    res_SEIIR = with(as.list(params_SEIIR), {
      sim_SEIIR(beta = beta, eta = eta, 
                phai = phai, ps = ps, 
                gammaA = gammaA, gammaS = gammaS, 
                N=N, tmax = tmax,
                plot = F)
    })
    disc_SEIIR = get_discrete(res_SEIIR, disc_tmax, tau)
    cat("SEIIR simulation done.\n")
    
    ## 2. sample a lot of sequences from SEI2IR
    events_SEI2IR = NULL
    
    for(r in 1:Rep2){
      ### A. simulate
      this_res = with(as.list(params_SEI2IR), {
        sim_SEI2IR(beta = beta, eta = eta, 
                   phai = phai, phiS = phiS, 
                   gammaA = gammaA, gammaS = gammaS, 
                   N=N, tmax = tmax, seed = r,
                   plot = F)
      })
      ### B. discretize
      this_disc = get_discrete(this_res, disc_tmax, tau)
      
      ### C. compare discrete data
      if(r == 1){
        compare = disc_SEIIR < this_disc
      }else{
        compare = compare + (disc_SEIIR < this_disc)
      }
      
      ### D. save events for plotting
      if(r <= plotRuns){
        # record 100 event sequences at most!
        events_SEI2IR = rbind(events_SEI2IR, this_disc)
      }
      cat(r,"\r")
    }
    cat("\n")
    
    ## post-process "compare"
    compare = compare[,-1]/Rep2
    for(c in 1:ncol(compare)){
      compare[,c] = sapply(compare[,c], function(p) ifelse(p<0.5, 2*p, 2*(1-p)))
    }
    
    Names = colnames(compare)
    compare = compare %>% as.data.frame()
    names(compare) = Names
    compare$t = disc_SEIIR$t
    
    return(list(event1 = disc_SEIIR, event2 = events_SEI2IR, 
                compare = compare, reference="SEI2IR"))
  }
}

## try it out
CD_res1 = compare_discrete(params_SEIIR = params1, 
                           params_SEI2IR = params2, 
                           N = 200, tmax = 200,
                           Rep1 = 300)

CD_res2 = compare_discrete(params_SEIIR = params1, 
                           params_SEI2IR = params2, 
                           N = 200, tmax = 200,
                           Rep1 = 1, Rep2 = 300)


# 6. function to visualize results from "compare_discrete"
visualize_CD <- function(CDList, 
                         plotLimits = NULL, smooth=FALSE){
  # plotLimits: a list for y-axis upper limits of the aggregate count plots
  # smooth: whether or not to smooth over the p-vals
  
  event1 = CDList$event1
  event2 = CDList$event2
  ref = CDList$reference
  
  compare = CDList$compare
  
  ## 1: plot aggregate counts over time
  if(is.null(plotLimits)){
    for(n in names(event1)[-1]){
      # a): plot line overlay the background
      print(
        ggplot(data=event2, aes_string(x="t", y=n)) +
          stat_density_2d(geom="raster",
                          aes(alpha = ..density..), contour = FALSE,
                          fill="blue")+
          #scale_alpha_continuous(range = c(0.3, 1)) +
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(expand=c(0,0))+
          geom_line(data=event1, aes_string(x="t", y=n), 
                    color="red", size=1)+
          labs(x='days',title=paste("Aggregate daily counts of",n),
               caption = paste("Reference model:", ref))+
          theme(legend.position = "none")
      )
      
      # b): show the background only
      print(
        ggplot(data=event2, aes_string(x="t", y=n)) +
          stat_density_2d(geom="raster",
                          aes(alpha = ..density..), contour = FALSE,
                          fill="blue")+
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(expand=c(0,0))+
          # geom_line(data=event1, aes_string(x="t", y=n), 
          #           color="red", size=1)+
          labs(x='days',title=paste("Aggregate daily counts of",n),
               caption = paste("Reference model:", ref))+
          theme(legend.position = "none")
      )
    }
  }else{
    for(n in names(event1)[-1]){
      # a): plot line overlay the background
      print(
        ggplot(data=event2, aes_string(x="t", y=n)) +
          stat_density_2d(geom="raster",
                          aes(alpha = ..density..), contour = FALSE,
                          fill="blue")+
          #scale_alpha_continuous(range = c(0.3, 1)) +
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(limits = c(0,plotLimits[[n]]), expand=c(0,0))+
          geom_line(data=event1, aes_string(x="t", y=n), 
                    color="red", size=1)+
          labs(x='days',title=paste("Aggregate daily counts of",n),
               caption = paste("Reference model:", ref))+
          theme(legend.position = "none")
      )
      
      # b): show the background only
      print(
        ggplot(data=event2, aes_string(x="t", y=n)) +
          stat_density_2d(geom="raster",
                          aes(alpha = ..density..), contour = FALSE,
                          fill="blue")+
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(expand=c(0,0))+
          # geom_line(data=event1, aes_string(x="t", y=n), 
          #           color="red", size=1)+
          labs(x='days',title=paste("Aggregate daily counts of",n),
               caption = paste("Reference model:", ref))+
          theme(legend.position = "none")
      )
    }
  }
  
  
  ## 2: plot "p-values" in "compare"
  
  ## a: plot S,E,R first
  compare_a = compare %>% 
    select(t, S, E, R) %>%
    gather(key='variable', value='value', -t)
  
  if(smooth){
    print(
      ggplot(data=compare_a, aes(x=t, y=value)) +
        geom_smooth(aes(color=variable), se=FALSE,
                    method='auto', span=0.1) +
        geom_hline(yintercept = 0.05, size=0.5) +
        labs(x='days', y='two-sided p-value',
             title="Empirical 'p-values' of aggregate counts",
             caption = paste("Reference model:", ref))
    )
  }else{
    print(
      ggplot(data=compare_a, aes(x=t, y=value)) +
        geom_line(aes(color=variable)) +
        geom_hline(yintercept = 0.05, size=0.5) +
        labs(x='days', y='two-sided p-value',
             title="Empirical 'p-values' of aggregate counts",
             caption = paste("Reference model:", ref))
    )
  }
  
  
  
  ## b: then plot Ia and Is
  compare_b = compare %>% 
    select(t, Ia, Is) %>%
    gather(key='variable', value='value', -t)
  
  if(!smooth){
    print(
      ggplot(data=compare_b, aes(x=t, y=value)) +
        geom_line(aes(color=variable)) +
        geom_hline(yintercept = 0.05, size=0.5) +
        labs(x='days', y='two-sided p-value',
             title="Empirical 'p-values' of aggregate counts",
             caption = paste("Reference model:", ref))
    )
  }else{
    print(
      ggplot(data=compare_b, aes(x=t, y=value)) +
        geom_smooth(aes(color=variable), se=FALSE,
                    method='auto', span=0.1) +
        geom_hline(yintercept = 0.05, size=0.5) +
        labs(x='days', y='two-sided p-value',
             title="Empirical 'p-values' of aggregate counts",
             caption = paste("Reference model:", ref))
    )
  }
  
  # p_compare = ggplot(data=compare, aes(x=t)) +
  #   geom_hline(yintercept = 0.05)
  # for(i in 2:length(names(event1))){
  #   p_compare = p_compare + 
  #     geom_line(aes_string(y=names(event1)[i]), color=i) +
  #     geom_point(aes_string(y=names(event1)[i]), color=i)
  # }
  # p_compare = p_compare + 
  #   labs(x='days',title="Empirical 'p-values' of daily aggregate counts")
  # print(p_compare)
}


## try it out 
pdf("Compare_SEIR_discrete_1.pdf", height = 5, width = 9)
visualize_CD(CD_res1)
dev.off()

pdf("Compare_SEIR_discrete_2.pdf", height = 5, width = 9)
visualize_CD(CD_res2)
dev.off()



# 08/01/2020
# do it with more simulation runs

PL = list(S=105, E=35, Ia = 20, Is = 40, R = 105)

CD_res1 = compare_discrete(params_SEIIR = params1, 
                           params_SEI2IR = params2, 
                           N = 200, tmax = 200,
                           Rep1 = 1000, 
                           plotRuns = 500)

pdf("Compare_SEIR_Ref_M1_500.pdf", height = 5, width = 9)
visualize_CD(CD_res1, plotLimits = PL)
dev.off()

CD_res2 = compare_discrete(params_SEIIR = params1, 
                           params_SEI2IR = params2, 
                           N = 200, tmax = 200,
                           Rep1 = 1, Rep2 = 1000, 
                           plotRuns = 500)

pdf("Compare_SEIR_Ref_M2_500.pdf", height = 5, width = 9)
visualize_CD(CD_res2, plotLimits = PL)
dev.off()

save(CD_res1, CD_res2, file="compare_SEIR_500.RData")


# 08/31/2020:
# re-make some visualizations
load("compare_SEIR_500.RData")

PL = list(S=105, E=35, Ia = 20, Is = 40, R = 105)

pdf("Compare_SEIR_Ref_M1_500_new.pdf", height = 4, width = 6)
visualize_CD(CD_res1, plotLimits = PL, smooth=TRUE)
dev.off()

pdf("Compare_SEIR_Ref_M2_500_new.pdf", height = 4, width = 6)
visualize_CD(CD_res2, plotLimits = PL, smooth=TRUE)
dev.off()


# plot the processes of talks
#library(xkcd)
#library(extrafont)

talks = data.frame(x=seq(from=0,to=30,by=5),
                   good = c(40,40,60, 60,70, 85,100),
                   bad = c(40, 40, 65, 65, 65, 45,35))

# good talk
ggplot(data=talks, aes(x=x,y=good)) + 
  geom_line(size=1.5, color='purple') +
  scale_y_continuous(limits = c(20,100), breaks = NULL) +
  annotate('text', x=3, y=50, label='intro', size=6) +
  annotate('text', x=15, y=70, label='build-up', size=6) +
  annotate('text', x=28, y=85, label='aha!', size=6) +
  labs(title='A typical talk is like...', x='minutes', y='excitement')

# bad talk
ggplot(data=talks, aes(x=x,y=bad)) + 
  geom_line(size=1.5, color='purple') +
  scale_y_continuous(limits = c(20,100), breaks = NULL) +
  annotate('text', x=3, y=48, label='intro', size=6) +
  annotate('text', x=15, y=66, label='kinda\n interesting', size=6) +
  annotate('text', x=28, y=50, label='hmm...', size=6) +
  labs(title='Today\'s talk might be like...', x='minutes', y='excitement')
