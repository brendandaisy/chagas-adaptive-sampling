require(magrittr)
require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
require(corrr)

source('spatial-helpers.R')
source('covariate-helpers.R')

extract_cor <- function(stanfit) {
    phi_draws <- as.data.frame(stanfit, pars = c('phi'))

    cor(phi_draws) %>%
        as_cordf %>%
        shave %>%
        stretch(na.rm = TRUE) %>%
        cor_info
}

cor_info <- function(cordf) {
    ## append info about the pairwise neighbors onto a cordf in long form
    n1 <- as.numeric(str_extract(cordf$x, '\\d+'))
    n2 <- as.numeric(str_extract(cordf$y, '\\d+'))
    
    cordf %>%
        mutate(
            adj = map2_lgl(n1, n2, ~sdm[.x, .y] == 1),
            dist = map2_dbl(n1, n2, ~am[.x, .y]),
            inf = map2_chr(n1, n2, ~as.character(dat_model$Y[.x] + dat_model$Y[.y]))
        )
}

### Get data ready

village <- 'Paternito'
cutoff <- 100

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
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

### Load necessary stan fits

dc_fit1 <- readRDS('stan-stuff/dcfit-feb2.rds')
dc_fit2 <- readRDS('stan-stuff/dcfit-dec15.rds')
dc_fit3 <- readRDS('stan-stuff/dcfit-feb3.rds')
sc_fit <- readRDS('stan-stuff/scfit-feb5.rds')

## summarize phi correlations for all fits together

all_cors <- list(dc_fit1, dc_fit2, dc_fit3, sc_fit) %>%
    imap_dfr(~mutate(extract_cor(.x), fit = as.factor(.y)))

all_cors %>%
    filter(dist < 500) %>%
    group_by(
        distrange = cut(dist, breaks = 10, labels = FALSE, ordered_result = TRUE),
        fit,
        inf
    ) %>%
    summarize(mean_r = mean(r)) %>%
    ggplot(aes(x = distrange, y = mean_r, col = inf, linetype = fit)) +
    ## geom_point() +
    geom_line() +
    labs(
        y = 'Corr(i, j)',
        x = 'Distance (meters / 50)',
        col = 'Num. Infested',
        shape = 'Model'
    ) +
    theme_few()

## phi correlations in more detail for single model

extract_cor(dc_fit1) %>%
    filter(dist < 500) %>%
    ggplot(aes(x = dist, y = r, col = inf)) +
    geom_point(size = .9) +
    geom_vline(xintercept = median(rstan::extract(dc_fit1)$x0), linetype = 'dashed') +
    labs(y = 'Corr(i, j)', x = 'Distance (meters)', col = 'Num. Infested') +
    theme_few()

## ggsave('corr-full.pdf', width = 8, height = 6)

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
