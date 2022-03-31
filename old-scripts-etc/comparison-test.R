require(INLA)
require(tidyverse)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
require(MASS, exclude=c('select'))

sample_y <- function(phi, inter) {
    prob <- exp(inter + phi) / (1 + exp(inter + phi))
    map_int(prob, ~rbernoulli(1, .x))
}

sample_tcar <- function(df, rho, x0, inter, tau = 1, repl = 100) {
    iso <- which(rowSums(dist_mat(df, cutoff = x0)) == 0)
    dfs <- df[-iso,]
    C <- dist_mat(dfs, cutoff = x0)
    Q <- tau * (diag(rowSums(C)) - rho * C)
    phis <- mvrnorm(n = repl, mu = rep(0, nrow(dfs)), Sigma = solve(Q)) %>%
        as_tibble
    pmap_dfr(phis, function(...) tibble_row(!!!sample_y(c(...), inter)))
}

prep_samp <- function(c_sub_iso, tb_org) {
    tb_org <- filter(tb_org, id %in% names(c_sub_iso))
    
    enframe(c_sub_iso, 'id', 'y') %>%
        bind_cols(long = tb_org$long, lat = tb_org$lat)
}

plot_samp <- function(c_sub_iso, tb_org) {
    prep_samp(c_sub_iso, tb_org) %>%
        ggplot(aes(x = long, y = lat, col = as.factor(y))) +
        geom_point()
}

run_dcar <- function(dat, k = .065, dmax = 500) {
    dcar_model <- inla.rgeneric.define(
        dcar_inla,
        A = dist_mat(dat),
        k = k,
        dmax = dmax
    )

    inla(
        y ~ 1 + f(id, model = dcar_model, values = unique(dat$id)),
        family = 'binomial',
        data = dat,
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(dic = TRUE, waic = TRUE),
        control.inla = list(control.correct = list(enable = 10))
    )
}

run_tcar <- function(dat, x0) {
    tcar_model <- inla.rgeneric.define(
        tcar_inla,
        B = dist_mat(dat, cutoff = x0),
    )

    inla(
        y ~ 1 + f(id, model = tcar_model, values = unique(dat$id)),
        family = 'binomial',
        data = dat,
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(dic = TRUE, waic = TRUE),
        control.inla = list(control.correct = list(enable = 10))
    )
}

fit_dcar <- function(draw, k = .065, dmax = 500) {
    ft <- run_dcar(draw, k, dmax)

    pred <- ifelse(ft$summary.fitted.values$mean >= .5, 1, 0)
    mpe <- mean(abs(draw$y - pred))
    mx <- dcar_tmarg(ft, dmax)[[3]]
    z <- inla.zmarginal(mx, TRUE)
    
    return(
        tibble_row(
            mpe = mpe,
            dic = ft$dic$dic,
            waic = ft$waic$waic,
            meanx0 = z$mean,
            sdx0 = z$sd,
            modex0 = inla.mmarginal(mx)
        )
    )
}

fit_tcar <- function(draw, x0) {
    ft <- run_tcar(draw, x0)
    pred <- ifelse(ft$summary.fitted.values$mean >= .5, 1, 0)
    mpe <- mean(abs(draw$y - pred))
    
    return(tibble_row(mpe = mpe, dic = ft$dic$dic, waic = ft$waic$waic))
}

fit_all <- function(samples, dat, type = 'dcar', cutoff = NULL) {
    if (type == 'dcar') {
        ret <- map_dfr(
            1:nrow(samples),
            ~fit_dcar(prep_samp(unlist(samples[.x,]), dat), dmax = 250)
        )
    }
    if (type == 'tcar') {
        ret <- map_dfr(
            1:nrow(samples),
            ~fit_tcar(prep_samp(unlist(samples[.x,]), dat), cutoff)
        )
    }
    ret
}

dat_sub <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == 'Prensa') %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:ownership),
        as_factor
    )

samp1 <- sample_tcar(dat_sub, rho = .95, x0 = 75, inter = -.2, rep = 10, tau = 1)

plot_samp(unlist(samp1[3,]), dat_sub)

draw <- prep_samp(unlist(samp1[3,]), dat_sub)

ft1 <- run_dcar(draw, dmax = 300)
plot_dcar_hyper(ft1, 300)

ft2 <- run_tcar(draw, 75)
plot_tcar_hyper(ft2)

dm <- dist_mat(draw, cutoff = 75)

bym_ft <- inla(
    y ~ 1 + f(id, model = 'besagproper', graph = dm, values = unique(draw$id)),
    family = 'binomial',
    data = draw,
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE)
    ## control.inla = list(control.correct = list(enable = 10))
)

res_dcar <- fit_all(samp1, dat_sub)
res_tcar1 <- fit_all(samp1, dat_sub, 'tcar', 75)
res_tcar2 <- fit_all(samp1, dat_sub, 'tcar', 200)
