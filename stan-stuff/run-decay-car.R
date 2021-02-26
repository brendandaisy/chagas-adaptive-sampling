require(rstan)
require(stringr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## get the locally saved model data
args <- commandArgs(trailingOnly = TRUE)
dat_model <- readRDS(args[1])
date <- Sys.Date()

dc_fit <- stan(
    file = 'decay-car.stan',
    data = dat_model,
    chains = 4,
    iter = 7000,
    warmup = 2000,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
)

saveRDS(dc_fit, paste0(str_remove(args[1], '.rds'), '_', date, '.rds'))
