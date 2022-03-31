require(rstan)
require(tidyverse)
require(caret)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
require(geosphere)

rstan_options(auto_write = TRUE)

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
  mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

adj_matrix <- get_neighbors(dist = 100)
iso <- which(rowSums(adj_matrix) == 0)
adj_matrix <- adj_matrix[-iso,-iso]

dat_org <- dat_org %>%
  arrange(desc(!(id %in% names(iso))))

cov_matrix <- make_dummy_mat(dat_org)

dat_model <- list(N = nrow(adj_matrix),
                  N_iso = length(iso),
                  K = ncol(cov_matrix),
                  Y = dat_org$infestation,
                  X = cov_matrix,
                  W = adj_matrix,
                  W_n = sum(adj_matrix) %/% 2)

fit_s2 <- stan(file = 's2-cauchy.stan',
            data = dat_model,
            iter=5000,
            warmup=2500,
            control=list(max_treedepth = 15))

phi_prior <- stan(file = 's2-phi-only.stan',
                  data = list(N = nrow(adj_matrix), 
                              Y = dat_org$infestation[1:nrow(adj_matrix)],
                              W = adj_matrix,
                              W_n = sum(adj_matrix) %/% 2))

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
