# -------------------------------------------------------------------------
# sequential-sampling.R----------------------------------------------------
# -------------------------------------------------------------------------
# high-level functions for performing and comparing adaptive and random----
# sampling strategies------------------------------------------------------
# -------------------------------------------------------------------------

require(tidyverse)
require(INLA)
require(furrr)

source('code/other-helpers.R')
source('code/seq-sampling-helpers.R')

##' Run `n_rep` replications of the simulation study
##' A sampling design is computed for a village's data, for each value of `alpha` plus random sampling
##' Multi-threading is used if enabled with `plan`
##'
##' @title run_simulation_study
##' @param df The "raw" data we will be getting samples from
##' @param alpha Vector of values in [0, 1] determining adaptive strategy
##' @param n_init Number of initial random samples
##' @param n_rep Number of replications of experiment
run_simulation_study <- function(df, alphas, n_init = 10, n_rep = 10,
                    pred = c('global', 'known', 'latent'),
                    tar_thresh = 0.05, conf_lvl = 0.95) {
    future_map_dfr(1:n_rep, ~{
        init <- sample(1:nrow(df), n_init)
        args <- c(list(NULL), map(alphas, ~list(alpha = .x)))
        des <- map(
            args,
            ~sampling_design(df, init, pred, tar_thresh, conf_lvl, strat_arg = .x)
        )
        map_dfr(des, ~design_score(.x, df$truth)) %>%
            mutate(thresh = tar_thresh, init = list(init), alpha = c(NA, alphas))
    }, .options = furrr_options(seed = TRUE), .progress = TRUE)
}

##' Score the design based on # sampled and percentage of infectious houses left
##'
##' @title design_score
##' @param des The design returned by `sampling_design`
##' @param true_inf A vector of the true infectious labels
design_score <- function(des, tru_inf) {
    end_des <- des[nrow(des),]$obs_idx[[1]]
    tibble_row(
        m = length(end_des),
        n = length(tru_inf),
        act_pct = sum(tru_inf[-end_des]) / length(tru_inf),
        sel = list(des$sel)
    )
}

##' Obtain data from a sampling design strategy.
##' This is the workhorse function for adaptive and random algorithms
##'
##' @title sampling_design
##' @param df The "raw" data we will be getting samples from
##' @param init Initial sample indices. If single no., number of initial samples instead
##' @param pred The type of predictors considered available
##' @param tar_thresh The target reduction threshold
##' @param conf_lvl The confidence level for meeting the reduction target
##' @param strat_arg A named list to pass to sampling strat
##' @param seed Randomization seed
##' @return A dataframe containing the INLA fit, observed indices so far, and locations just chosen for each iteration
sampling_design <- function(
    df, init = 10, pred = c('global', 'known', 'interp', 'latent'),
     tar_thresh = 0.05, conf_lvl = 0.95, silent = TRUE, strat_arg = NULL, seed = NULL
) {

    # setup spde mesh and other stuff
    if (!is.null(seed))
        set.seed(seed)
    ntot <- nrow(df)
    if (length(init) == 1) # then init = # of random samples to start with
        init <- sample(ntot, init)
    mesh <- inla_mesh(df)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    df$infestation[-init] <- NA # "mark" unvisited nodes
    done <- FALSE
    ret <- map_dfr(1:ntot, ~{
        if (done)
            return(tibble_row(fit = NULL, obs_idx = NULL, sel = NULL))
        if (!silent)
            print(paste0("Begin iteration ", .x))

        # while not done, fit model to current observations
        obs <- which(!is.na(df$infestation))
        unobs <- which(is.na(df$infestation))
        fit_dat <- build_pred_dat(df, pred)
        stack <- spde_stack(fit_dat, mesh, spde)
        ft <- fit_gp_both(df = fit_dat, spde = spde, stack = stack)

        # check if control target has been met
        if (.x > 2 & length(unique(fit_dat$truth[obs])) > 1)
            done <<- check_complete(ft, unobs, tar_thresh * ntot, conf_lvl)

        # if not, select new locations to sample
        if (is.null(strat_arg)) {
            sel <- rand_unif(unobs, ft)
        }
        else {
            strat_arg$t <- (length(obs) - length(init)) / (ntot - length(init))
            sel <- comb_risk_var(unobs, ft, strat_arg$t, strat_arg$alpha)
        }
        df$infestation[sel] <<- df$truth[sel]
        tibble_row(fit = ft, obs_idx = list(obs), sel = list(sel))
    })
    return(drop_na(ret))
}

##' Check whether the design target has been met
##'
##' @title check_complete
##' @param fit An INLA object
##' @return A boolean indicating whether Pr(pct. infested > `tar_inc`) > `conf_lvl`
check_complete <- function(fit, unobs_idx, tar_inc, conf_lvl) {

    # Draw from the posterior predictions on unvisited locations and return # infested
    samp_inc <- function(s) {
        s$latent %>%
            inla.link.invlogit %>%
            map_int(~sample(c(0L, 1L), size = 1, prob = c(1 - .x, .x))) %>%
            sum
    }

    # Draw from the joint risk at unvisited locations
    samp <- inla.posterior.sample(
        5000, fit,
        selection = list(APredictor = unobs_idx),
        use.improved.mean = FALSE,
        skew.corr = FALSE
    )

    # Get predicted num. infested and return estimated CDF
    pred_inc <- map_int(samp, samp_inc)
    (sum(pred_inc < tar_inc) / length(pred_inc)) > conf_lvl
}
