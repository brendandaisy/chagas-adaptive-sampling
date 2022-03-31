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
cov_matrix <- make_dummy_mat(dat_org)

count_cov <- function(...) {
  nb <- which(c(...) == 1)
  if (length(nb) > 1) { 
    t <- as_tibble_row(colSums(cov_matrix[nb,-1]))
  }
  else if (length(nb) == 1) { 
    t <- as_tibble_row(cov_matrix[nb,-1])
  }
  else {
    r <- rep(0, dat_model$K - 1)
    names(r) <- colnames(cov_matrix)[-1]
    t <- as_tibble_row(r)
  }
  return(t)
}

cov_counts <- as_tibble(adj_matrix, .name_repair='minimal') %>%
  pmap_dfr(count_cov) %>%
  mutate(id = dat_org$id, nb = unname(rowSums(adj_matrix)), village = dat_org$village, inter = 1) %>%
  select(inter, everything())

## plot some examples for checking
cov_counts %>%
  bind_cols(long = dat_org$long, lat = dat_org$lat, org=dat_org$num_pigs) %>%
  filter(village == 'Amatillo') %>%
  ggplot(aes(x=long, y = lat, col=num_pigs, size=org)) +
  geom_point()

g2 <- dat_org %>%
  filter(village == 'Amatillo') %>%
  ggplot(aes(x=long, y = lat, size=num_pigs)) +
  geom_point()

grid.arrange(g1, g2, nrow=1)

dat_model <- list(N = nrow(adj_matrix),
                  K = ncol(cov_matrix),
                  Y = dat_org$infestation,
                  Z = select(cov_counts, -id, -nb, -village))

fit_s3 <- stan(file = 's3-neighbors.stan',
               data = dat_model)

mcmc_trace(as.array(fit_s3), pars='gamma[1]')

mcmc_neff(neff_ratio(fit_s3))

## check posterior densities, shrinkage from prior

post_gamma <- as.data.frame(fit_s3, pars='gamma') %>%
  as_tibble
names(post_gamma) <- c('inter', colnames(cov_matrix)[-1])

mcmc_dens(post_gamma) +
  stat_function(fun = dcauchy, args = list(location = 0, scale = 2.5), color = "red")

## posterior prediction, confusion matrix
ppc_bars(dat_org$infestation, as.matrix(fit_s3, pars='y_pred'))

pred_per_col <- as.matrix(fit_s3, pars='y_pred') %>%
  t %>%
  as_tibble

conf_dist <- tibble(
  tp = as.integer(summarize_all(pred_per_col, ~sum(. == dat_org$infestation & . == 1))[1,]),
  tn = as.integer(summarize_all(pred_per_col, ~sum(. == dat_org$infestation & . == 0))[1,]),
  fp = as.integer(summarize_all(pred_per_col, ~sum(. != dat_org$infestation & . == 1))[1,]),
  fn = as.integer(summarize_all(pred_per_col, ~sum(. != dat_org$infestation & . == 0))[1,])
)

conf_dist <- conf_dist %>%
  mutate(`SENS/REC` = tp / (tp + fn),
         `SPEC` = tn / (tn + fp),
         `PREC/PPV` = tp / (tp + fp),
         `ACC` = (tp + tn) / dat_model$N,
         `bACC` = (`SENS/REC` + `SPEC`) / 2,
         `NPV` = tn / (tn + fn))

conf_dist %>%
  pivot_longer(cols=`SENS/REC`:NPV, names_to='Metric') %>%
  ggplot(aes(y=value, x=Metric, fill=Metric)) +
  geom_violin() +
  theme_gdocs() +
  theme(axis.title.y=element_blank())

## wow...is there a way to test if the PPD is approximately the same as in S1?
## plot the 'risks' for each house from glm, S1, and S3

both_risks <- as.data.frame(fit_s3, pars='r')[,c(1, 100, 400)] %>%
  rename_all(~paste0(., 'S3')) %>%
  bind_cols(as.data.frame(fit_s1, pars='r')[1:4000, c(1, 100, 400)]) %>%
  as_tibble

both_risks %>%
  pivot_longer(everything()) %>%
  mutate(Model = ifelse(str_detect(name, 'S3'), 'S3', 'S1'),
         House_Num = str_extract(name, '\\d+')) %>%
  ggplot(aes(x=value, col=Model)) +
  geom_density(fill=NA) +
  facet_wrap(~House_Num)
