require(tidyverse)
require(INLA)
require(furrr)
## require(MASS, exclude = 'select')
## require(spdep)
## require(gstat)
## require(precrec)

source('spatial-helpers.R')
source('covariate-helpers.R')
source('as-fn-helpers.R')

bigg_bs <- function(df, alphas, n_init = 10, n_rep = 10,
                    pred = c('global', 'known', 'interp', 'latent'),
                    tar_thresh = 0.05, conf_lvl = 0.95, silent = FALSE) {
    future_map_dfr(1:n_rep, ~{
        ## if (!silent)
        ##     print(paste0('Sample ', .x, '/', n_rep))
        init <- sample(1:nrow(df), n_init)
        args <- c(list(NULL), map(alphas, ~list(alpha = .x)))
        des <- map(
            args,
            ~sampling_design(df, init, pred, tar_thresh, conf_lvl, strat_arg = .x)
        )
        map_dfr(des, ~des_score(.x, df$truth)) %>%
            mutate(thresh = tar_thresh, init = list(init), alpha = c(NA, alphas))
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
}

bs_perf <- function(df, strat, n_init = 10, n_rep = 10,
                    pred = c('global', 'known', 'interp', 'latent'),
                    tar_thresh = 0.05, conf_lvl = 0.95,
                    seed = TRUE, silent = TRUE, ...) {
    dots <- list(...) # gather dots so map not confused
    map_dfr(1:n_rep, ~{
        if (!silent)
            print(paste0('Sample ', .x, '/', n_rep))
        des <- sampling_design(
            df, strat, n_init,
            pred, tar_thresh, conf_lvl,
            strat_arg = dots
        )
        mutate(des_score(des, df$truth), thresh = tar_thresh, conf = conf_lvl)
    })
}

des_score <- function(des, tru_inf) {
    end_des <- des[nrow(des),]$obs_idx[[1]]
    tibble_row(
        m = length(end_des),
        n = length(tru_inf),
        act_pct = sum(tru_inf[-end_des]) / length(tru_inf),
        sel = list(des$sel)
    )
}

##' Workhorse function for adaptive and random algorithms
##'
##' @title sampling_design
##' @param init Initial sample indices. If single no., number of initial samples instead
##' @param strat_arg A named list to pass to sampling strat
sampling_design <- function(df, init = 10,
                            pred = c('global', 'known', 'interp', 'latent'),
                            tar_thresh = 0.05, conf_lvl = 0.95,
                            silent = TRUE, strat_arg = NULL, seed = NULL) {
    if (!is.null(seed))
        set.seed(seed)
    ntot <- nrow(df)
    if (length(init) == 1)
        init <- sample(ntot, init)
    mesh <- inla_mesh(df)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    df$infestation[-init] <- NA
    done <- FALSE
    ret <- map_dfr(1:ntot, ~{
        if (done)
            return(tibble_row(fit = NULL, obs_idx = NULL, sel = NULL))
        ## if (str_detect(strat, 'comb'))
        obs <- which(!is.na(df$infestation))
        unobs <- which(is.na(df$infestation))
        fit_dat <- build_pred_dat(df, obs, pred)
        stack <- spde_stack(fit_dat, mesh, spde)
        ft <- fit_gp2(df = fit_dat, spde = spde, stack = stack)
        ## checking here implies didn't actually "visit" sel
        if (.x > 2 & length(unique(fit_dat$truth[obs])) > 1)
            done <<- check_complete(ft, unobs, tar_thresh * ntot, conf_lvl)
        if (is.null(strat_arg))
            ## sel <- do.call(strat, splice(unobs, ft))
            sel <- rand_unif(unobs, ft)
        else {
            strat_arg$t <- (length(obs) - length(init)) / (ntot - length(init))
            sel <- comb_risk_var(unobs, ft, strat_arg$t, strat_arg$alpha)
        }
        df$infestation[sel] <<- df$truth[sel]
        tibble_row(fit = ft, obs_idx = list(obs), sel = list(sel))
    })
    return(drop_na(ret))
}

check_complete <- function(fit, unobs_idx, tar_inc, conf_lvl) {
    samp_inc <- function(s) {
        s$latent %>%
            inla.link.invlogit %>%
            map_int(~sample(c(0L, 1L), size = 1, prob = c(1 - .x, .x))) %>%
            sum
    }
    if (sum(fit$summary.fitted.values$mean[unobs_idx]) > 2.1 * tar_inc)
        return(FALSE)
    samp <- inla.posterior.sample(
        5000, fit,
        selection = list(APredictor = unobs_idx),
        use.improved.mean = FALSE,
        skew.corr = FALSE
    )
    pred_inc <- map_int(samp, samp_inc)
    (sum(pred_inc < tar_inc) / length(pred_inc)) > conf_lvl
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

## des <- sampling_design(dat_sub, init = 70, pred = 'known', strat_arg = list(alpha = 0.15), tar_thresh = 0.1)
## des_score(des, dat_sub$truth)

## plot_obs_surface(
    ## slice_tail(des)$fit[[1]]$summary.fitted.values,
    ## slice_tail(des)$obs_idx[[1]],
##     ft$summary.linear.predictor[1:107,],
##     1:nrow(dat_sub),
##     dat_sub
## )

## pd <- build_pred_dat(dat_sub, 1:nrow(dat_sub), 'global')
## s <- sample(1:nrow(dat_sub), 80)
## pd$infestation[-s] <- NA
## mesh <- inla_mesh(pd)
## spde <- inla.spde2.pcmatern(
##     mesh = mesh,
##     alpha = 2,
##     prior.range = c(0.1, 0.05),
##     prior.sigma = c(3, 0.1)
## )

## A <- inla.spde.make.A(mesh = mesh, loc = sp_project(dat_sub, normalize = TRUE))
## stack <- spde_stack(pd, mesh, spde)
## ft <- fit_gp2(pd, spde, stack)

## hyp <- inla.hyperpar.sample(5000, ft, improve.marginals = TRUE)
## dm <- dist_mat(dat_sub, normalize = TRUE)
## m <- nrow(dm)

## h1 <- hyp[1,]
## cov <- matrix(0, m, m)
## for(i in 1:m) {
##     for (j in i:m) {
##         cov[i, j] <- h1[2]^2 * INLA:::inla.matern.cf(dm[i, j], h1[1], 1)
##         cov[j, i] <- cov[i, j]
##     }
## }

## mu <- A %*% ft$summary.random$z$mean
## Z <- mvrnorm(1, mu, cov)

## fixed <- ft$marginals.fixed
## inter <- inla.rmarginal(1, 

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
