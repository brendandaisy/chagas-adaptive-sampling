require(tidyverse)
require(RColorBrewer)
require(rstan)
require(ggthemes)
require(bayesplot)

source('covariate-helpers.R')
source('spatial-helpers.R')
rstan_options(auto_write = TRUE)

### read the data
village <- 'Cerrón'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:infestation),
        as_factor
    )

dat_oh <- make_dummy_mat(dat_org, contr='contr.ltfr')[,-1]
dm <- dist_mat(dat_org)

dat_model <- list(
    N = nrow(dat_org),
    K = ncol(dat_oh),
    Y = as.double(dat_org$infestation) - 1,
    X = dat_oh,
    D = dm
)

## Model fit
fit_pairs <- stan(
    file = 'neighbors-only.stan',
    data = dat_model
)
