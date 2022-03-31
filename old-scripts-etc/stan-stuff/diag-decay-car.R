require(rstan)
require(bayesplot)
require(tidyverse)

dc_fit <- readRDS('stan-fits/LCARk065_2021-02-26.rds')

pairs(dc_fit, pars = c('alpha', 'x0', 'tau', 'phi[1]', 'beta[3]'))
