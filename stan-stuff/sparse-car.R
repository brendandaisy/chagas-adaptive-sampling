require(magrittr)
require(rstan)
require(tidyverse)
require(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('../spatial-helpers.R')
source('../covariate-helpers.R')

village <- 'Paternito'
cutoff <- 100

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:infestation),
        as_factor
    )

dm <- dist_mat(dat_org, cutoff = cutoff)

## remove isolated houses for now
iso <- which(rowSums(dm) == 0)
sdm <- dm[-iso, -iso]
dat_sub <- dat_org[-iso, ]
am <- dist_mat(dat_sub)

cm <- make_dummy_mat(dat_sub)

### Fit sparse CAR with the predictors

dat_model <- list(
    N = nrow(dat_sub),
    K = ncol(cm),
    Y = as.numeric(dat_sub$infestation) - 1,
    X = cm,
    C = sdm,
    C_n = sum(sdm) %/% 2,
    tau_a = 2,
    tau_b = 1
)

sc_fit <- stan(
    file = 'sparse-car.stan',
    data = dat_model,
    chains = 4,
    iter = 7000,
    warmup = 2000,
    control = list(max_treedepth = 15)
)

saveRDS(sc_fit, 'scfit-feb5.rds')

### Basic posterior diagnosis

pairs(sc_fit, pars = c('alpha', 'tau', 'phi[1]'))
