source('adaptive-sampling-fn.R')
source('as-fn-helpers.R')

NREP <- 15
FUN <- 'comb_risk_var'
PRED <- 'known'
SEED <- 1233 # not used by rand methods
RDS <- 'combriskvarpt5-known'

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds')

bs_perf_pat <- bs_perf(
    filter(dat_org, village == 'Paternito'),
    FUN,
    n_rep = NREP,
    pred = PRED,
    seed = SEED,
    silent = FALSE
)

saveRDS(bs_perf_pat, paste0('bootstrapss/', RDS, '-pat.rds'))

bs_perf_gua <- bs_perf(
    filter(dat_org, village == 'Guayabo'),
    FUN,
    n_rep = NREP,
    pred = PRED,
    seed = SEED,
    silent = FALSE
)

saveRDS(bs_perf_gua, paste0('bootstrapss/', RDS, '-gua.rds'))

bs_perf_pre <- bs_perf(
    filter(dat_org, village == 'Prensa'),
    FUN,
    n_rep = NREP,
    pred = PRED,
    seed = SEED,
    silent = FALSE
)

saveRDS(bs_perf_pre, paste0('bootstrapss/', RDS, '-pre.rds'))

bs_perf_cer <- bs_perf(
    filter(dat_org, village == 'Cerrón'),
    FUN,
    n_rep = NREP,
    pred = PRED,
    seed = SEED,
    silent = FALSE
)

saveRDS(bs_perf_cer, paste0('bootstrapss/', RDS, '-cer.rds'))

bs_perf_ama <- bs_perf(
    filter(dat_org, village == 'Amatillo'),
    FUN,
    n_rep = NREP,
    pred = PRED,
    seed = SEED,
    silent = FALSE
)

saveRDS(bs_perf_ama, paste0('bootstrapss/', RDS, '-ama.rds'))
