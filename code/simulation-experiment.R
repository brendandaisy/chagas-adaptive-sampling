source('adaptive-sampling-fn.R')
source('as-fn-helpers.R')

NREP <- 2
PRED <- 'global'
THRESH <- 0.05
### DANGER!!!! <- 3.4 # change this whenever new round (not) compat with prev samples

plan(multisession, workers = 2)
dat_org <- read_csv('data/survey-data.csv') # TODO this 

bs_pat <- run_simulation_study(
    filter(dat_org, village == 'Paternito'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_pat, paste0('bootstrapss/bigg-', PRED, '-pat.rds'))

bs_gua <- run_simulation_study(
    filter(dat_org, village == 'Guayabo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_gua, paste0('bootstrapss/bigg-', PRED, '-gua.rds'))

bs_pre <- run_simulation_study(
    filter(dat_org, village == 'Prensa'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_pre, paste0('bootstrapss/bigg-', PRED, '-pre.rds'))

bs_cer <- run_simulation_study(
    filter(dat_org, village == 'Cerr�n'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_cer, paste0('bootstrapss/bigg-', PRED, '-cer.rds'))

bs_ama <- run_simulation_study(
    filter(dat_org, village == 'Amatillo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_ama, paste0('bootstrapss/bigg-', PRED, '-ama.rds'))
