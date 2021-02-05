require(magrittr)
require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
require(corrr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('../spatial-helpers.R')
source('../covariate-helpers.R')

village <- 'Paternito'

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:infestation),
        as_factor
    )

dm <- dist_mat(dat_org)

## remove isolated houses for now
iso <- which(rowSums(dist_mat(dat_org, cutoff = 100)) == 0)
dm <- dm[-iso, -iso]
dat_org <- dat_org[-iso, ]

### Code to fit model
cov_mat <- make_dummy_mat(dat_org)

dat_model <- list(
    N = nrow(dat_org),
    K = ncol(cov_mat),
    Y = as.numeric.factor(dat_org$infestation),
    X = cov_mat,
    A = dm,
    tau_a = 2,
    tau_b = 1,
    kmin = .03,
    kmax = .1,
    x0min = 50,
    x0max = 800
)

dc_fit <- stan(
    file = 'glm-decay-car.stan',
    data = dat_model1,
    chains = 4,
    iter = 7000,
    warmup = 2000,
    control = list(adapt_delta = 0.9, max_treedepth = 15)
)

### From here, either load a fit or use the one above fit

dc_fit <- readRDS('dcfit-feb3.rds')

pairs(dc_fit, pars = c('x0', 'k', 'tau', 'phi[1]'))

## how have phi correlations changed?
phi_draws <- as.data.frame(dc_fit, pars = c('phi'))

phi_cor <- cor(phi_draws) %>%
    as_cordf %>%
    shave %>%
    stretch(na.rm = TRUE) %>%
    cor_info

phi_cor %>%
    filter(dist < 500) %>%
    ggplot(aes(x = dist, y = r, col = inf)) +
    geom_point(size = .9) +
    geom_vline(xintercept = median(rstan::extract(dc_fit)$x0), linetype = 'dashed') +
    labs(y = 'Corr(i, j)', x = 'Distance (meters)', col = 'Num. Infested') +
    theme_few()

ggsave('corr-full.pdf', width = 8, height = 6)

## how much autrocorrelation is left in beta?
## TODO: consider if it is better to test each beta_p seperately
zdraws <- as.data.frame(dc_fit, pars = c('Z_hat'))

zcor <- cor(zdraws) %>%
    as_cordf %>%
    shave %>%
    stretch(na.rm = TRUE) %>%
    cor_info

ggplot(zcor, aes(x = dist, y = r)) +
    geom_point(alpha = .6) +
    geom_smooth() +
    theme_few()

## posterior prediction: find MSE
ydraws <- as.data.frame(sc_fit, pars = c('Y_hat'))

mpe <- function(yhat) {
    mean((yhat - dat_model$Y)^2)
}

mpe_ppd <- apply(ydraws, 1, mpe)
summary(mpe_ppd)

## mcmc_scatter(glm_full, pars=c('alpha', 'beta[3]'))

## mcmc_areas(glm_full, pars=c('x0'))
## mcmc_intervals(glm_full, regex_pars='phi')

## fit_df <- as.data.frame(glm_full, pars=c('phi'))

## fit_df %>%
##     map_dbl(var) %>%
##     bind_cols(k = colSums(dat_model2$A)) %>%
##     ggplot(aes(x = k, y = ...1)) +
##     geom_point()

## fit_pairs <- fit_df %>%
##     cor %>%
##     as_tibble() %>%
##     mutate(i = 1:nrow(dat_org)) %>%
##     pivot_longer(-i, values_to='cor') %>%
##     mutate(
##         j = rep(1:nrow(dat_org), nrow(dat_org))
##         ## k_i = rowSums(dat_model2$C[i, ]),
##         ## k_j = rowSums(dat_model2$C[j, ])
##     ) %>%
##     select(i, j, everything(), -name) %>%
##     filter(i < j) # now, every pair unique

## fit_pairs %<>%
##     mutate(
##         inf_i = as.numeric.factor(dat_org$infestation[i]),
##         inf_j = as.numeric.factor(dat_org$infestation[j])
##     )

## fit_pairs %>%
##     pivot_longer(c(k_i, k_j)) %>% # all pairs should appear twice
##     ggplot(aes(x = value, y = cor, col = as.factor(inf_i + inf_j))) +
##     geom_point(alpha = .4, shape = 1)

## fit_pairs %>%
##     mutate(dist = map2_dbl(i, j, ~dm[.x, .y])) %>%
##     ggplot(aes(x = dist, y = cor, col = as.factor(inf_i + inf_j))) +
##     geom_point(alpha = .5) +
##     theme_minimal()

## visualize network, network uncertainty
net_parms <- as.data.frame(glm_full, pars=c('x0', 'alpha'))

cdist <- c(dm[upper.tri(dm)]) %>%
    map_dfr(~tibble(dist = .x, alpha = net_parms$alpha, x0 = net_parms$x0)) %>%
    mutate(aC = alpha *(1 - my_logistic(dist, x0 = x0, k = .065)))

sum_cdist <- cdist %>%
    group_by(dist) %>%
    summarize(maC = median(aC), laC = quantile(aC, .05), uaC = quantile(aC, .95))

sub_cdist <- cdist %>%
    group_by(dist) %>%
    filter(aC > quantile(aC, .95))

ggplot(sum_cdist, aes(x = dist, y = maC)) +
    geom_line() +
    geom_ribbon(aes(ymin=laC, ymax=uaC), alpha = .3) +
    geom_point(aes(y = aC), data = sub_cdist, alpha = .1)

mean_mat <- reduce2(
    net_parms$x0,
    net_parms$alpha,
    ~..3 * (1 - my_logistic(dm, x0 = ..2, k = .065)) + ..1,
    .init = 0
) / nrow(net_parms)

diag(mean_mat) <- 0

mean_mat %>%
    as_tibble %>%
    mutate(sid = dat_org$id) %>%
    pivot_longer(-sid, values_to='C', names_to = 'tid') %>%
    ggplot(aes(x=sid, y=tid, fill=C)) +
    geom_tile()

mean_mat %>%
    as_tibble %>%
    mutate(sid = dat_org$id) %>%
    pivot_longer(-sid, values_to='C', names_to = 'tid') %>%
    mutate(dist = map2_dbl(sid, tid, ~dm[.x, .y])) %>%
    ggplot(aes(x = dist, y = C)) +
    geom_point() +
    theme_minimal()
