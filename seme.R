require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)

rstan_options(auto_write = TRUE)

source('covariate-helpers.R')
source('spatial-helpers.R')

village = 'Guayabo'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
  mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor) %>%
  filter(village == !!village)

cov_matrix <- make_dummy_mat(dat_org)
dist_mat <- dist_mat_sncar(dat_org)

dat_model <- list(
  N = nrow(dat_org),
  K = ncol(cov_matrix),
  Y = dat_org$infestation,
  X = cov_matrix,
  A = dist_mat
)

fit_seme <- stan(
  file = 'seme.stan',
  data = dat_model,
  iter=20000,
)
