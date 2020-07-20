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
                       N=200, tmax = 100)



