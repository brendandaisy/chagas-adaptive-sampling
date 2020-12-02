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
dist_mat <- dist_mat_sncar(dat_org, d_low=30)

dat_model <- list(
  N = nrow(dat_org),
  K = ncol(cov_matrix),
  Y = dat_org$infestation,
  X = cov_matrix,
  A = dist_mat,
  d_max = max(dist_mat) / 4,
  ep = .0001,
  alpha = .05
)

fit_sncar <- stan(
  file = 'sncar.stan',
  data = dat_model,
  iter=20000,
)



pairs(fit_s2, pars = c('tau', 'alpha', 'lp__', 'phi[1]', 'phi[20]'))

mcmc_trace(as.array(phi_prior), pars='phi[1]')

mcmc_neff(neff_ratio(phi_prior))

## check posterior densities, shrinkage from prior

post_phi <- as.data.frame(fit_s2, pars='phi')
names(post_phi) <- dat_org[1:nrow(adj_matrix),]$id

prior_phi <- as.data.frame(phi_prior, pars='phi') %>%
  as_tibble
names(prior_phi) <- names(post_phi)

prior_phi %>%
  select(c(8, 25, 100, 4, 350, 420, 700, 800)) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  stat_density() +
  facet_wrap(~name, nrow=2)

mcmc_dens(post_phi[,c(8, 25, 100, 4, 350, 420, 700, 800)]) +
  stat_density(data=pivot_longer(prior_phi, c(8, 25, 100, 4, 350, 420, 700, 800), values_to='Value'), 
               col='red', 
               fill=NA,
               n=1200) +
  ylim(c(0, .05))

post_beta <- as.data.frame(fit_s2, pars='beta') %>%
  as_tibble
names(post_beta) <- c('inter', colnames(cov_matrix)[-1])

mcmc_dens(post_beta) +
  stat_function(fun = dcauchy, args = list(location = 0, scale = 2.5), color = "red")
