source('adaptive-sampling-fn.R')
source('as-fn-helpers.R')

NREP <- 2
PRED <- 'global'
THRESH <- 0.05
### DANGER!! Change VERSION!!!!!!!!!
VERSION <- 3.4 # change this whenever new round (not) compat with prev samples

plan(multisession, workers = 2)
dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds')

bs_pat <- bigg_bs(
    filter(dat_org, village == 'Paternito'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_pat, paste0('bootstrapss/bigg-', PRED, '-v', VERSION, '-pat.rds'))

bs_gua <- bigg_bs(
    filter(dat_org, village == 'Guayabo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_gua, paste0('bootstrapss/bigg-', PRED, '-v', VERSION, '-gua.rds'))

bs_pre <- bigg_bs(
    filter(dat_org, village == 'Prensa'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_pre, paste0('bootstrapss/bigg-', PRED, '-v', VERSION, '-pre.rds'))

bs_cer <- bigg_bs(
    filter(dat_org, village == 'Cerrón'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_cer, paste0('bootstrapss/bigg-', PRED, '-v', VERSION, '-cer.rds'))

bs_ama <- bigg_bs(
    filter(dat_org, village == 'Amatillo'),
    alphas = c(0, .15, .3, .7, 1, 2),
    n_init = 10,
    n_rep = NREP,
    pred = PRED,
    tar_thresh = THRESH
)

saveRDS(bs_ama, paste0('bootstrapss/bigg-', PRED, '-v', VERSION, '-ama.rds'))
