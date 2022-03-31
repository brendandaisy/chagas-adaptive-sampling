require(tidyverse)

source('../spatial-helpers.R')
source('../covariate-helpers.R')

village <- 'Paternito'
decay_func <- 'LCAR'
k_param <- 0.065
obj_name <- 'LCARkpt065'

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

cov_mat <- make_dummy_mat(dat_org)

dat_model <- list(
    N = nrow(dat_org),
    K = ncol(cov_mat),
    Y = as.numeric.factor(dat_org$infestation),
    X = cov_mat,
    A = dm,
    dc = if(decay_func == 'LCAR') 1 else 2,
    k = k_param,
    x0min = 50,
    x0max = 800
)

saveRDS(dat_model, paste0(obj_name, '.rds'))
