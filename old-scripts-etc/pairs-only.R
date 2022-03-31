###
## Script for fitting and evaluating the 'pairs only' stan model
###

require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)

rstan_options(auto_write = TRUE)
source('covariate-helpers.R')
source('spatial-helpers.R')

### Read and prepare the data
village <- 'Cerrón'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:ownership),
        as_factor
    )

dat_oh <- make_dummy_mat(dat_org, as_tibble=TRUE)

tg <- make_cov_net(dat_oh, method='knn', k=1)

edge_dat <- tg %E>%
    as_tibble

## TODO: source infestation can't be treated as a double in regression
dat_model <- list(
    E = nrow(edge_dat),
    K = ncol(dat_oh) - 4, # don't count geo ids
    Y = as.numeric(pull(edge_dat, target_inf)) - 1,
    X = select(edge_dat, inter:last_col(1))
)

## Model fit
fit_pairs <- stan(
    file = 'pairs-only.stan',
    data = dat_model
)

## Diagnostics
rhats <- rhat(fit_pairs, pars='beta')
mcmc_rhat(rhats)

post_beta <- as.data.frame(fit_pairs, pars='beta') %>%
  as_tibble
names(post_beta) <- colnames(dat_model$X)

mcmc_dens(post_beta) +
    stat_function(fun = dcauchy,
                  args = list(location = 0, scale = 2.5),
                  color = "red")
