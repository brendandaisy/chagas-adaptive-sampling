###
## Script for fitting and evaluating the 'neighbors only' stan model
###

require(rstan)
require(tidyverse)
require(bayesplot)
require(brms)
require(ggthemes)
require(RColorBrewer)

rstan_options(auto_write = TRUE)
source('covariate-helpers.R')
source('spatial-helpers.R')

run_model <- function(expr, fname, reuse = TRUE) {
  path <- paste0(fname, ".Rds")
  if (reuse) {
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
  }
  else {
    fit <- eval(expr)
    saveRDS(fit, file = path)
  }
  fit
}

### Read and prepare the data
village <- 'Cerrón'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:ownership),
        as_factor
    )

contrasts(dat_org$bed_hygiene) <- contr.sum

## dat_oh <- make_dummy_mat(dat_org, inter=TRUE, as_tibble=TRUE)

tg <- make_cov_net(dat_oh, method='cutoff', h=100)

tga <- colnames(dat_oh)[6:ncol(dat_oh)] %>% # include neighbor infestation
    reduce(~add_neighbor_attr(.x, sym(.y), paste0('sum_', .y), sum), .init=tg)

dat_brm <- select(
    as_tibble(tga),
    !(c(name, long, lat, inter))
)
colnames(dat_brm) <- str_replace_all(colnames(dat_brm), ' ', '') # brms wants like this

## fit_nb <- stan('neighbors-only.stan', data = dat_model, iter=3000)
fnb <- as.formula(
    paste0(
        'infestation ~ ',
        paste(colnames(select(dat_brm, matches('sum'))), collapse='+')
    )
)

fit_test <- brm(infestation ~ bed_hygiene, data=dat_org,

fit_nb <- run_model(
    expr = brm(
        fnb,
        data = dat_brm,
        family = 'bernoulli',
        prior = prior(cauchy(0, 2.5)),
        chains = 5,
        iter = 10000,
        warmup = 2000,
        save_all_pars=TRUE
    ),
    fname = 'neighbors-only',
    reuse=FALSE
)

## Diagnostics
rhats <- rhat(fit_nb, pars='beta')
mcmc_rhat(rhats)

plot(conditional_effects(fit_nb))

post_gamma <- as.data.frame(fit_nb) %>%
    as_tibble

mcmc_dens(post_gamma) +
    stat_function(fun = dcauchy, args = list(location = 0, scale = 2.5), color = "red")

posterior_summary(fit_nb)
