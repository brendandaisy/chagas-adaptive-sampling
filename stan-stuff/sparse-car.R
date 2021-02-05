require(magrittr)
require(rstan)
require(tidyverse)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
## require(MASS, exclude=c('select'))
require(corrr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('../spatial-helpers.R')
source('../covariate-helpers.R')

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
    file = 'glm-sparse-car.stan',
    data = dat_model,
    chains = 4,
    iter = 7000,
    warmup = 2000,
    control = list(max_treedepth = 15)
)

## how have phi correlations changed?
phi_draws <- as.data.frame(sc_fit, pars = c('phi'))

phi_cor <- cor(phi_draws) %>%
    as_cordf %>%
    shave %>%
    stretch(na.rm = TRUE) %>%
    cor_info

phi_cor %>%
    filter(dist < 500) %>%
    ggplot(aes(x = dist, y = r, col = inf)) +
    geom_vline(xintercept = 100, linetype = 'dashed') +
    geom_point(size = .9) +
    labs(y = 'Corr(i, j)', x = 'Distance (meters)', col = 'Num. Infested') +
    theme_few()

ggsave('corr-sparse.pdf', width = 8, height = 6)

## how much autrocorrelation is left in beta?
## TODO: consider if it is better to test each beta_p seperately
zdraws <- as.data.frame(sc_fit, pars = c('Z_hat'))

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


