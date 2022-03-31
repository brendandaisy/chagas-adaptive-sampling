###
## all code for fitting and evaluating the "house only" model, including saving model
## for comparison to other methods
###

require(rstan)
require(tidyverse)
require(magrittr)
require(brms)
require(bayestestR)
require(bayesplot)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)

rstan_options(auto_write = TRUE)
source('covariate-helpers.R')

### Load the data. Note if currently using all villages or single village
dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

cov_matrix <- make_dummy_mat(dat_org)

## fit a GLM from base for later comparison
dat_bg <- cbind(cov_matrix[,-1], infestation=dat_org$infestation) %>%
    as_tibble
base_glm <- glm(infestation ~ ., data=dat_bg, family=binomial)

## data for model
dat_model <- list(
    N = nrow(cov_matrix),
    K = ncol(cov_matrix),
    Y = dat_org$infestation,
    X = cov_matrix
)

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

options(contrasts = c('contr.bayes', 'contr.poly'))

fs1 <- as.formula(
    paste0(
        'infestation ~ ',
        paste(colnames(select(dat_org, bed_hygiene:ownership)), collapse='+'),
        '+ (1 | village)'
    )
)

## TODO: make sure priors are unbiased

## pri_s1 <- brm(
##     fs1,
##     data = dat_org,
##     family = 'bernoulli',
##     prior = prior(cauchy(0, 2.5)),
##     sample_prior = 'only',
##     iter=200
## )

## emmeans(fit_s1, ~bed_hygiene)

## now actually draw from posterior
fit_s1 <- brm(
    fs1,
    data = dat_org,
    family = 'bernoulli',
    prior = c(
        prior(cauchy(0, 2.5), class='b'),
        prior(cauchy(0, 2.5), class='Intercept')
    ),
    chains = 5,
    iter = 10001,
    warmup = 2000,
    save_all_pars=TRUE
)

plot(fit_s1)

### Run MCMC estimation of the posterior
fit_s1 <- stan(
    file = 'house-only.stan',
    data = dat_model,
    iter=3000,
    control = list(max_treedepth=18)
)

### Model diagnostics
## check for divergences, n_eff, rhat
pairs(fit_s1, pars=c("r_village[Amatillo,Intercept]",
                     "r_village[Prensa,Intercept]",
                     "lp__"),
      las=1)

rhats <- rhat(fit_s1, pars='beta')
mcmc_rhat(rhats)

ratios <- neff_ratio(fit_s1, pars='beta')

ratios %>%
    enframe %>%
    mutate(
        name = as_factor(c('(Intercept)', colnames(cov_matrix)[-1])),
        warn = ifelse(value < .1, 'bad', ifelse(value < .5, 'med', 'good'))
    ) %>%
    ggplot(aes(x=fct_rev(name), y=value, fill=warn)) +
    geom_col() +
    coord_flip()

## traceplots / autocorrelation
full_post <- as.array(fit_s1)
mcmc_trace(full_post,
           pars='beta[1]') # looks okay, but whats up with low jumps in some chains?
mcmc_trace(full_post, pars='beta[1]', window=c(1000, 1300))
mcmc_acf(full_post, pars='beta[1]')

### Check posterior densities, shrinkage from prior, compare to base GLM
post_beta <- as.data.frame(fit_s1, pars='beta')
names(post_beta) <- c('(Intercept)', colnames(cov_matrix)[-1])

post_beta %>%
    pivot_longer(everything(), names_to='coef') %>%
    ggplot(aes(x=value)) +
    stat_density() +
    stat_function(fun = dcauchy, args = list(location = 0, scale = 2.5), color = "red") +
    geom_vline(aes(xintercept=bg), 
               tibble(bg = base_glm$coefficients, coef=names(post_beta)),
               col='blue', 
               linetype='dashed') +
    facet_wrap(~coef, nrow=6, scales='free')

## ggsave('beta-post-prior.pdf', g, width=20, height=10)

## Posterior prediction, confusion matrix
ppc_bars(dat_org$infestation, as.matrix(fit, pars='y_pred'))

pred_per_col <- stan_extract(fit, parms='y_pred') %>%
    t %>%
    as_tibble

conf_dist <- tibble(
    tp = as.integer(summarize_all(
        pred_per_col,
        ~sum(. == dat_org$infestation & . == 1))[1,]),
    tn = as.integer(summarize_all(
        pred_per_col,
        ~sum(. == dat_org$infestation & . == 0))[1,]),
    fp = as.integer(summarize_all(
        pred_per_col,
        ~sum(. != dat_org$infestation & . == 1))[1,]),
    fn = as.integer(summarize_all(
        pred_per_col,
        ~sum(. != dat_org$infestation & . == 0))[1,])
)

conf_dist <- conf_dist %>%
    mutate(
        `SENS/REC` = tp / (tp + fn),
        `SPEC` = tn / (tn + fp),
        `PREC/PPV` = tp / (tp + fp),
        `ACC` = (tp + tn) / dat_model$N,
        `bACC` = (`SENS/REC` + `SPEC`) / 2,
        `NPV` = tn / (tn + fn)
    )

conf_dist %>%
    pivot_longer(cols=`SENS/REC`:NPV, names_to='Metric') %>%
    ggplot(aes(y=value, x=Metric, fill=Metric)) +
    geom_violin() +
    theme_gdocs() +
    theme(axis.title.y=element_blank())

## Plot median predicted risks in Guayabo
risks <- stan_extract(fit, 'r')
colnames(risks) <- dat_org$id

risks <- risks %>%
    select(contains('Guayabo')) %>%
    pivot_longer(everything(), names_to = 'id') %>%
    group_by(id) %>%
    summarize(med = median(value),
              sd = sd(value))

gua_risks <- dat_org %>%
    filter(village == 'Guayabo') %>%
    select(id, lat, long, infestation) %>%
    left_join(risks, by='id')

g1 <- gua_risks %>%
    ggplot(aes(x=long, y=lat, col=med, shape=infestation == 1)) +
    geom_point(size=1.7) +
    scale_color_gradientn(
        colors=rev(brewer.pal(11, 'RdYlGn')),
        name='Risk') +
    labs(x = 'Longitude', y='Latitude', shape='Infestation Status') +
    theme_tufte() +
    theme(panel.grid = element_line(color='grey70'),
          axis.ticks = element_line(color='grey70'))

g2 <- gua_risks %>%
    ggplot(aes(x=long, y=lat, col=sd, shape=infestation == 1)) +
    geom_point(size=1.7) +
    scale_color_gradientn(
        colors=rev(brewer.pal(11, 'Spectral')),
        name='Std Dev') +
    labs(x = 'Longitude', y='Latitude', shape='Infestation Status') +
    theme_tufte() +
    theme(panel.grid = element_line(color='grey70'),
          axis.ticks = element_line(color='grey70'))

## ggsave('med-risks-gua.pdf', grid.arrange(g1, g2, nrow=1), width=13, height=6)

### Plot dist. of means of risks
risk_summ <- params_fit %>%
    select(contains('r[')) %>%
    pivot_longer(cols=everything()) %>%
    mutate(r_ind = str_extract(name, 'r\\[\\d+\\]')) %>%
    group_by(r_ind) %>%
    summarize(mean_val = mean(value), sd_val = sd(value))


ggplot(risk_summ, aes(x=reorder(r_ind, mean_val), y=mean_val)) +
    geom_col()
