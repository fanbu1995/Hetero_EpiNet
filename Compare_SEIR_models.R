# 07/19/2020

# Simulate two SEIR extension models and compare

library(GillespieSSA)

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
  ## TO DO
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

compRes2 = compare_models(params_SEIIR = params1, 
                         params_SEI2IR = params2, 
                         N = 200, tmax = 200,
                         Rep = 10)

# To Do: add summary statistics visualization as well as trajectory visualization

