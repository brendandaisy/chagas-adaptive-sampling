require(magrittr)
require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('../spatial-helpers.R')

village <- 'Paternito'

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    select(id:long, infestation)

dist_mat <- dist_mat(dat_org)

## check for reasonable decay func
tibble(x = c(dm3), y = 1 - my_logistic(x, x0 = 250, k=10^-2)) %>%
    ggplot(aes(x, y)) +
    geom_point()

## in particular, need to know combinations leading to valid graph
## cpd <- expand.grid(
##     alpha = seq(0, 1, .05),
##     x0 = seq(0, max(dist_mat) / 2, length.out=20),
##     k = seq(0, .5, .01)
## ) %>%
##     pmap_dfr(~{
##         C <- 1 - my_logistic(dist_mat, x0=..2, k=..3)
##         D <- diag(rowSums(C))
##         tibble_row(
##             alpha = ..1,
##             x0 = ..2,
##             k = ..3,
##             pd = is.positive.definite(D - ..1 * C)
##         )
##     })

## cpd %>%
##     filter(x0 == 0) %>%
##     ggplot(aes(x = k, y = alpha, fill=pd)) +
##     geom_tile()

## pick a few houses to inspect in detail
p <- c(59, 20, 8, 52)
dat_org %>%
    mutate(indx = 1:n(), col = indx %in% p) %>%
    ggplot(aes(long, lat, col = col, label = indx)) +
    geom_text(size=2.5)

## how does a house's degree change with k and x0?
seq(.01, .5, length.out=20) %>%
    map_dfr(~{
        C <- 1 - my_logistic(dist_mat, x0=200, k = .x)
        diag(C) <- 0
        tibble(id = dat_org$id, deg = colSums(C), k = .x)
    }) %>%
    ggplot(aes(x = k, y = deg, group=id)) +
    geom_line() +
    theme_minimal()

degx <- seq(10, 1000, length.out=40) %>%
    map_dfr(~{
        C <- 1 - my_logistic(dist_mat, x0=.x, k = .065)
        diag(C) <- 0
        tibble(id = dat_org$id, deg = colSums(C), x0 = .x)
    })

degx %>%
    filter(floor(x0) == 111, deg < 1) %>%
    pull(id) %>%
    unique

dat_model <- list(
    N = nrow(dat_org),
    Y = dat_org$infestation,
    A = dist_mat - 200,
    tau_a = 2,
    tau_b = 1,
    alpha = .95
)

priors_fit <- stan(
    file = 'phi-priors-only.stan',
    data = dat_model,
    chains = 1,
    iter = 5000,
    warmup = 2000,
    control = list(adapt_delta = 0.99)
)

pairs(priors_fit, pars=c('k', 'tau', str_c('phi[', p, ']'), 'phi[47]'))

mcmc_areas(priors_fit, pars=str_c('phi[', p, ']'))

fit_df <- as.data.frame(priors_fit, pars=c('phi', 'k'))

fit_df %>%
    map_dbl(var) %>%
    bind_cols(k = colSums(dat_model$A)) %>%
    ggplot(aes(x = k, y = ...1)) +
    geom_point()

fit_df %>%
    ggplot(aes(x=k, 

fit_pairs <- fit_df %>%
    select(contains('phi')) %>%
    cor %>%
    as_tibble() %>%
    mutate(i = 1:nrow(dat_org)) %>%
    pivot_longer(-i, values_to='cor') %>%
    mutate(
        j = rep(1:nrow(dat_org), nrow(dat_org)),
        k_i = rowSums(dat_model$A[i, ]),
        k_j = rowSums(dat_model$A[j, ])
    ) %>%
    select(i, j, everything(), -name) %>%
    filter(i < j) # now, every pair unique

fit_pairs %<>%
    mutate(inf_i = dat_org$infestation[i], inf_j = dat_org$infestation[j])

fit_pairs %>%
    mutate(dist = map2_dbl(i, j, ~dist_mat[.x, .y])) %>%
    ggplot(aes(x = dist, y = cor, col = as.factor(inf_i + inf_j))) +
    geom_point(alpha = .2) +
    geom_vline(xintercept = 200, col = 'red')
