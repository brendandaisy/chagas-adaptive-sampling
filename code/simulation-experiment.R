# -------------------------------------------------------------------------
# simulation-experiment.R--------------------------------------------------
# -------------------------------------------------------------------------
# simulation study described in Sec. 2.6-----------------------------------
# run for each village using the given params and a variety of `alpha`s----
# -------------------------------------------------------------------------

library(future)
source('code/sequential-sampling.R')

NREP <- 1
PRED <- 'global'
THRESH <- 0.05
N_INIT <- 10
CONF_LVL <- 0.95

plan(multisession, workers = 4) # set to number of available cores
dat_org <- read_csv('anon-survey-data.csv')

pat_results <- run_simulation_study(
    filter(dat_org, village == 'Paternito'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = N_INIT,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH,
    conf_lvl = CONF_LVL
)

gua_results <- run_simulation_study(
    filter(dat_org, village == 'Guayabo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = N_INIT,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH,
    conf_lvl = CONF_LVL
)

pre_results <- run_simulation_study(
    filter(dat_org, village == 'Prensa'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = N_INIT,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH,
    conf_lvl = CONF_LVL
)

cer_results <- run_simulation_study(
    filter(dat_org, village == 'Cerr�n'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = N_INIT,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH,
    conf_lvl = CONF_LVL
)

ama_results <- run_simulation_study(
    filter(dat_org, village == 'Amatillo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = N_INIT,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH,
    conf_lvl = CONF_LVL
)
