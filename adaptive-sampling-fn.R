require(tidyverse)
require(INLA)
require(MASS, exclude = 'select')
require(spdep)
require(gstat)
require(precrec)

source('spatial-helpers.R')
source('covariate-helpers.R')
source('as-fn-helpers.R')

bs_perf <- function(df, strat, n_init = 10, n_rep = 10,
                    pred = c('global', 'known', 'interp', 'latent'),
                    seed = TRUE, silent = TRUE, ...) {
    dots <- list(...) # gather dots so map not confused
    map_dfr(1:n_rep, ~{
        if (!silent)
            print(paste0('Sample ', .x, '/', n_rep))
        is <- sampling_design(df, strat, n_init, pred, args = dots)
        mutate(is_score(is, df$truth, n_init, strat), rep = .x)
    })
}

is_score <- function(is, tru_inf, n_init, strat) {
    imap_dfr(c(.25, .5, .75), ~{ # score at 3 points in run
        if (str_detect(strat, 'rand'))
            ii <- .y
        else
            ii <- ceiling((.x * length(tru_inf) - n_init) / 3)
        ooss <- fit_score(
            is$fit[[ii]],
            is$obs_idx[[ii]],
            tru_inf,
            ## set thresh to 0 since not using
            thresh = 0,
            in_samp = FALSE
        )
        tibble_row(
            !!!ooss,
            size = .x,
            obs = list(is$obs_idx[[ii]])
        )
    })
}

find_best_thresh <- function(fit, obs_idx, tru_inf) {
    obs_inf <- factor(tru_inf[obs_idx], c(0, 1))
    obs_prob <- fit$summary.fitted.values$mean[obs_idx]
    folds <- createMultiFolds(obs_inf, k = 4, times = 10)
    cv_thresh <- map_dbl(folds, cv_best_thresh, y = obs_inf, p = obs_prob)
    mean(cv_thresh)
}

cv_best_thresh <- function(fold_idx, y, p) {
    y_fold <- y[fold_idx]
    p_fold <- p[fold_idx]
    scores <- map_dfr(seq(.05, .95, .025), ~{
        pred <- factor(ifelse(p_fold >= .x, 1, 0), c(0, 1))
        tibble_row(!!!bool_stats(table(pred, y_fold)), thresh = .x)
    })
    qacc <- quantile(scores$bacc, .9)
    return(arrange(filter(scores, bacc >= qacc), desc(sens))$thresh[1])
}

fit_score <- function(fit, obs_idx, tru_inf, thresh = NULL, in_samp = FALSE) {
    if (in_samp)
        idx <- obs_idx
    else
        idx <- -obs_idx
    if (is.null(thresh))
        thresh <- find_best_thresh(fit, obs_idx, tru_inf)
    ## make the predictions for accuracy targets
    pred <- map_dbl(
        fit$summary.fitted.values$mean[idx],
        ~ifelse(.x >= thresh, 1, 0)
    )
    ## prediction accuracy metrics:
    bool <- bool_stats(table(factor(pred, c(0, 1)), tru_inf[idx]))
    ## precision-recall AUC:
    prc <- evalmod(
        mode = 'rocprc',
        scores = fit$summary.fitted.values$mean[idx],
        labels = tru_inf[idx]
    )
    auc <- attr(prc$prcs[[1]], 'auc')
    roc <- attr(prc$rocs[[1]], 'auc')
    ## baseline class proportion
    prop <- sum(tru_inf[idx]) / length(tru_inf[idx])
    ## average prediction variance:
    pred_var <- mean(fit$summary.linear.predictor$sd[idx]^2)
    ## "discovery rate" always for in-sample
    disc <- sum(tru_inf[obs_idx]) / sum(tru_inf)
    ## * (length(tru_inf) / length(obs_idx))
    return(tibble_row(
        !!!bool, auc = auc, roc = roc,
        prop = prop, var = pred_var,
        disc = disc, thresh = thresh
    ))
}

bool_stats <- function(conf_mat) {
    sens <- conf_mat[2, 2] / sum(conf_mat[,2])
    spec <- conf_mat[1, 1] / sum(conf_mat[,1])
    list(
        sens = sens,
        spec = spec,
        bacc = (sens + spec) / 2
    )
}

##' Wrapper function for adaptive and random algorithms
##'
##' @title sampling_design
##' @param args A named list to pass to sampling strat
sampling_design <- function(df, strat = 'rand_unif', n_init = 10,
                               pred = c('global', 'known', 'interp', 'latent'),
                               silent = TRUE, args = NULL) {
    if (str_detect(strat, 'rand'))
        return(random_sampling(df, strat, n_init, pred, args))
    adaptive_sampling(df, strat, n_init, pred, args)
}

adaptive_sampling <- function(df, strat, n_init, pred, args) {
    ntot <- nrow(df)
    niter <- ceiling((.75 * ntot - n_init) / 3)
    ini <- sample(ntot, n_init)
    mesh <- inla_mesh(df)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    df$infestation[-ini] <- NA
    map_dfr(1:niter, ~{
        if (str_detect(strat, 'comb'))
            args <- c(t = (.x - 1) / (niter - 1), args)
        obs <- which(!is.na(df$infestation))
        fit_dat <- build_pred_dat(df, obs, pred)
        ## now using gp2!
        ft <- fit_gp2(df = fit_dat, mesh = mesh, spde = spde)
        if (is.null(args))
            sel <- do.call(strat, splice(which(is.na(df$infestation)), ft))
        else
            sel <- do.call(strat, splice(which(is.na(df$infestation)), ft, args))
        df$infestation[sel] <<- df$truth[sel]
        tibble_row(fit = ft, obs_idx = list(obs), sel = list(sel))
    })
}

random_sampling <- function(df, strat, n_init, pred, args) {
    mesh <- inla_mesh(df)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    map_dfr(1:3, ~{
        p <- sum(rep(.25, .x))
        tdf <- df # just to be safe
        m <- n_init + 3 * ceiling((p * nrow(df) - n_init) / 3 - 1)
        ## -1^^ last adapt iter not added
        if (is.null(args))
            obs <- do.call(strat, splice(tdf, m))
        else
            obs <- do.call(strat, splice(tdf, m, args))
        tdf$infestation[-obs] <- NA
        fd <- build_pred_dat(tdf, obs, pred)
        ft <- fit_gp2(df = fd, mesh = mesh, spde = spde)
        tibble_row(fit = ft, obs_idx = list(obs))
    })
}

build_pred_dat <- function(df, obs_idx, pred_type) { # df NOT dummy coded
    if (pred_type == 'known')
        df_ret <- df
    else if (pred_type == 'global')
        df_ret <- select(df, id:village, dist_perim:truth)
    else if (pred_type == 'latent')
        return(select(df, id:village, infestation:truth)) # don't rescale anything
    return(rescale_cont_vars(df_ret))
}

## dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds')
## dat_sub <- filter(dat_org, village == 'Paternito')

## is <- sampling_design(dat_sub, 'convex_comb', n_init = 10, pred = 'global', silent = FALSE)
## is_score(is, dat_sub$truth, 10, 'convex_comb')

## ft1 <- fit_old(dat_sub)
## ft2 <- fit_gp(dat_sub, mesh, spde)

## ft1$summary.hyperpar
## ft2$summary.hyperpar

## ft1$summary.fixed$mean
## ft2$summary.fixed$mean

## plot(ft1)

## ft1$summary.fitted.values$mean - ft2$summary.fitted.values$mean

## plot_obs_surface(ft1$summary.fitted.values, 1:nrow(dat_sub), dat_sub)
## plot_obs_surface(ft2$summary.fitted.values, 1:nrow(dat_sub), dat_sub)

## plot(ft1$marginals.hyperpar$`Range for idx`, type = 'l', xlim = c(0, 20))
## lines(ft2$marginals.hyperpar$`Range for z`, col = 'red')

## plot(ft1$marginals.fixed$`(Intercept)`, type = 'l')
## lines(ft2$marginals.fixed$`(Intercept)`, col = 'red')

## plot(ft1$marginals.fitted.values$fitted.Predictor.100, type = 'l', ylim = c(0, 2.8))
## lines(ft2$marginals.fitted.values$fitted.Predictor.100, col = 'red')

## plot(seq(0, 10000, 10), inla.pc.dprec(seq(0, 100, .1), u = .1, alpha = .05))

## train <- bind_rows(
##     slice_sample(filter(dat_sub, infestation == 1), n = 23),
##     slice_sample(filter(dat_sub, infestation == 0), n = 30)
## )

## test <- anti_join(dat_sub, train)
## test$infestation <- NA
## fit_dat <- build_pred_dat(bind_rows(train, test), 0, 'global')
## fit <- fit_gp(fit_dat, inla_mesh(fit_dat))

## pred <- fit$summary.fitted.values$mean[which(is.na(fit_dat$infestation))]

## dat_sub$infestation[-s] <- NA
## pd <- build_pred_dat(dat_sub, s, 'global')

## ft <- fit_gp(pd, inla_mesh(dat_sub))
## summary(ft1)
## fit_score(ft, s, dat_sub$truth)

## plot_obs_surface(ft$summary.fitted.values, s, dat_sub)
## plot_obs_surface(ft2$summary.random$idx, s, dat_sub)

## set.seed(111)
## is <- iterative_sampling(dat_sub, 'tar_risk', 30, pred = 'global')
## is_score(is, dat_sub$truth, 30, 'tar')
## fit_score(is$fit[[2]], is$obs_idx[[2]], dat_sub$truth, 0.3)

## plot_obs_surface(
##     is$fit[[ceiling((.25 * nrow(dat_sub) - 30) / 3)]]$summary.fitted.values,
##     is$obs_idx[[ceiling((.25 * nrow(dat_sub) - 30) / 3)]],
##     dat_sub
## )

## mean(dat_sub$truth[is$obs_idx[[ceiling((.25 * nrow(dat_sub) - 30) / 3)]]])
## mean(dat_sub$truth[-is$obs_idx[[ceiling((.25 * nrow(dat_sub) - 30) / 3)]]])

## oby <- is$obs_idx[[ceiling((.5 * nrow(dat_sub) - 10) / 3)]]
## oby <- sample(nrow(dat_sub), 50)
## dfy <- dat_sub
## dfy$infestation[-oby] <- NA
## fd <- build_fit_dat(dfy, kld$obs_idx[[ceiling((.5 * nrow(dat_sub) - 10) / 3)]])
## fd <- build_fit_dat(dfy, oby)
## ft <- fit_gp(fd)
## fit_score(ft, kld$obs_idx[[ceiling((.5 * nrow(dat_sub) - 10) / 3)]], dat_sub$truth)
## fit_score(ft, oby, dat_sub$truth)
## plot_obs_surface(ft$summary.random$idx, oby, dat_sub)
