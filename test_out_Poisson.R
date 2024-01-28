# 11/30/2020

# test out the Poisson regression model
# to estimate beta and b_S

library(tidyverse)

# num of indivs
N = 100
# dim of covars
P = 2

set.seed(131)
X = rbernoulli(n = N * P, p=0.5)
X = matrix(as.numeric(X), nrow = N, ncol = P)
b_S = c(0,1)

# tmax
tt = 10
# beta
bet = 0.2

# individual rates for Poisson process
onset_times = rexp(n = N, rate = bet * exp(X %*% b_S))

# dwell times
dwell = ifelse(onset_times < tt, onset_times, tt)

# binary outcome
Y = as.numeric(onset_times < tt)
nI = sum(Y)

# repeat in iterations
# use formulas (15) and (16)
S = 50
bet_hat = bet
for(s in 1:S){
  p.mod = glm(Y ~ X - 1 + offset(rep(log(bet_hat), N)+log(dwell)), family = poisson())
  b_S_hat = coef(p.mod)
  bet_hat = nI/sum(exp(X %*% b_S_hat) * dwell)
  cat('b_S: ', b_S_hat, " beta: ", bet_hat, '\n')
}

# estimate beta by intercept in Poisson regr
p.mod = glm(Y ~ X + offset(log(dwell)), family = poisson())
params = coef(p.mod)
(b_S_hat = params[2:3])
(bet_hat = exp(params[1]))

# both things are the same...




