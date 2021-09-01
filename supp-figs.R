require(tidyverse)
require(tikzDevice)
require(MASS, exclude = c('select'))

source('adaptive-sampling-fn.R')

options(
    tikzLatexPackages = c(getOption( "tikzLatexPackages" ), '\\usepackage{mathrsfs}')
    ## tikzMetricPackages = c(
    ##     "\\usepackage[utf8]{inputenc}",
    ##     "\\usepackage[T1]{fontenc}",
    ##     "\\usetikzlibrary{calc}"
    ## )
)

full_vil_summ <- function(df, village, pred_type = 'global', fit_fun = fit_gp2) {
    dfs <- filter(df, village == !!village)
    fit_dat <- build_pred_dat(dfs, 1, pred_type)
    mesh <- inla_mesh(dfs)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    stack <- spde_stack(fit_dat, mesh, spde)
    ft <- fit_fun(fit_dat, spde, stack, control = list(mlik = TRUE, dic = TRUE))
    
    tibble_row(
        cpu = ft$cpu.used[4],
        dic = ft$dic$dic,
        mlik = ft$mlik[1, 1],
        rho = rlang::list2('{village}' := ft$marginals.hyperpar[[1]]),
        fix = list(ft$marginals.fixed),
        village = village
    )
}

marg_rho <- function(summ) imap_dfr(summ$rho, ~tibble(v = .y, !!!as_tibble(.x)))

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds') %>%
    mutate(village = str_replace(village, 'ó', "\\\\'o"))

### Compare proposed model to simpler ones--------------------------------------

summ0 <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_iid))
summ <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp))
summ2 <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp2))

summ0b <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, fit_fun = fit_iid, pred_type = 'known')
)

summb <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, pred_type = 'known', fit_fun = fit_gp)
)

summ2b <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, fit_fun = fit_gp2, pred_type = 'known')
)

res <- bind_rows(
    mutate(summ0, model = 'indep', pred = 'global'),
    mutate(summ, model = 'spatial', pred = 'global'),
    mutate(summ2, model = 'both', pred = 'global'),
    mutate(summ0b, model = 'indep', pred = 'all'),
    mutate(summb, model = 'spatial', pred = 'all'),
    mutate(summ2b, model = 'both', pred = 'all')
)

res %>%
    group_by(village, pred) %>%
    slice_max(mlik, with_ties = TRUE) %>%
    ungroup %>%
    count(model)

gg <- res %>%
    rename(`CPU(sec)` = cpu, DIC = dic, ML = mlik) %>%
    pivot_longer(`CPU(sec)`:ML) %>%
    ggplot(
        aes(
            interaction(pred, village), value,
            fill = factor(model, labels = c('$\\varepsilon(s)$ only', '$Z(s)$ only', 'Both'))
            ## fill = fct_relevel(model, 'indep', 'spatial')
        )
    ) +
    geom_col(position = position_dodge()) +
    facet_wrap(~fct_relevel(name, 'DIC', 'ML'), scales='free', nrow = 3) +
    ## scale_fill_manual(labels = c('$\\varepilon$ only', '$Z$ only', 'Both')) +
    labs(x = 'Predictor $\\times$ village', y = NULL, fill = NULL) +
    theme_bw()

tikz_plot(gg, 'sens-study', 9, 4.5)

### Effective range for both predictor sets-------------------------------------

gg <- bind_rows(
    mutate(marg_rho(summ2), pred = 'Global $(p = 3)$'),
    mutate(marg_rho(summ2b), pred = 'All $(p = 27)$')
) %>%
    mutate(pred = fct_relevel(pred, 'Global $(p = 3)$')) %>%
    ggplot(aes(x, y, col = v)) +
    geom_line() +
    facet_wrap(~pred, nrow = 1, scale = 'free') +
    xlim(0, 4) +
    labs(
        x = 'Effective range $\\rho$',
        ## y = 'Density $p(\\rho \\mid \\mathscr{S})$',
        y = NULL,
        col = NULL
    ) +
    theme_bw()

tikz_plot(gg, 'range-full', 7.75, 3.25)

### Fixed effects for the best model--------------------------------------------

inla_fits <- function(village, df, pred_type = 'global', fit_fun = fit_gp2) {
    dfs <- filter(df, village == !!village)
    fit_dat <- build_pred_dat(dfs, 1, pred_type)
    mesh <- inla_mesh(dfs)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    fit_fun(fit_dat, spde, spde_stack(fit_dat, mesh, spde))
}

summ_fxe <- function(summ_df, village) {
    sumv <- filter(summ_df, village == !!village)
    imap_dfr(sumv$fix[[1]][-1], ~{
        if (.y == 'bed_hygienemala')
            return(tibble())
        ci5 <- inla.hpdmarginal(0.5, .x)
        ci95 <- inla.hpdmarginal(0.95, .x)
        tibble(
            var = .y,
            m = inla.zmarginal(.x, TRUE)$mean,
            lo5 = ci5[,1],
            hi5 = ci5[,2],
            lo95 = ci95[,1],
            hi95 = ci95[,2],
            village = village
        )
    })
}

fx_all <- map_dfr(unique(dat_org$village), ~summ_fxe(summ2b, .x))

## which vars were most consistent/important?

fx_all %>%
    filter(hi5 < 0) %>%
    count(var, sort = TRUE) %>%
    filter(n >= 3)
        
gg <- fx_all %>%
    mutate(var = str_replace_all(var, '\\_', '\\\\_')) %>%
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col = 'darkorange2') +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col = 'darkorange2') +
    geom_point(size = 1.5, col = 'dodgerblue3') +
    facet_wrap(~village, nrow = 1) +
    labs(x = 'Log odds', y = NULL) +
    theme_bw() +
    theme(panel.spacing = unit(0, "lines"))

tikz_plot(gg, 'fixed-all', 7.6, 5)

### Preferential sampling example-----------------------------------------------

dperim <- function(pt) min(1 - pt$x, pt$x, 1 - pt$y, pt$y)

sim_gp <- function(pts) {
    n <- nrow(pts)
    cov <- matrix(0, n, n)
    dprm <- rep(0, n)
    for(i in 1:n) {
        dprm[i] <- dperim(pts[i,])
        for(j in i:n) {
            p1 <- pts[i,]
            p2 <- pts[j,]
            d <- sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2)
            cov[i, j] <- 0.25 * INLA:::inla.matern.cf(d, range = 0.3, nu = 0.5)
            cov[j, i] <- cov[i, j]
        }
    }
    resp <- mvrnorm(1, mu = -3.5 * dprm + rep(0.75, n), Sigma = cov) %>%
        inla.link.invlogit %>%
        map_int(~sample(c(0L, 1L), size = 1, prob = c(1 - .x, .x)))
    mutate(pts, dperim = dprm, resp = resp)
}

fit_simple_gp <- function(df, sel) {
    df$resp[-sel] <- NA
    pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
    hyper <- list(
        range = list(param = c(0.1, 0.05)),
        prec = list(initial = log(0.25), fixed = TRUE),
        nu = list(initial = log(1), fixed = TRUE)
    )
    loc <- as.matrix(select(df, x, y))
    inla(
        resp ~ dperim + f(idx, model = 'dmatern', locations = loc, hyper = hyper),
        data = df,
        family = 'binomial',
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE),
        control.fixed = pfixed
    )
}

samp_inc <- function(fit, sel) {
    samp <- inla.posterior.sample(
        5000, fit,
        selection = list(Predictor = -sel),
        use.improved.mean = TRUE,
        skew.corr = TRUE
    )
    map_int(samp, ~{
        .x$latent %>%
            inla.link.invlogit %>%
            map_int(~sample(c(0L, 1L), size = 1, prob = c(1 - .x, .x))) %>%
            sum
    })
}

set.seed(14)

dat <- tibble(x = runif(80), y = runif(80), idx = 1:80) %>%
    sim_gp

dat <- read_csv('sim-gp.csv')

ini_sel <- c(45, 18, 61, 71, 72, 75, 8)

theme <- theme_bw() +
    theme(
        panel.background = element_rect(fill = gray(0.95)),
        axis.title = element_blank()
    )

gg1 <- ggplot(dat, aes(x, y)) +
    geom_point(aes(shape = as.factor(resp))) +
    # geom_text(aes(label = idx)) +
    geom_point(data = dat[ini_sel,], shape = 1, size = 5.5) +
    scale_shape_manual(
        labels = c('$y = 0$', '$y = 1$'),
        values = c(`0` = 16, `1` = 4)
    ) +
    ## labs(title = 'Pr$(y = 1) = g(0.2 - 2 x_1(s) + Z(s))$') +
    theme +
    theme(
        legend.position = c(.09, .93),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, 'pt'), 
        legend.background = element_rect(fill = gray(0.95), color = 'black'),
        legend.key = element_rect(fill = gray(0.95)),
        ## legend.justification = 'center',
        legend.margin = margin(.5, 1, .5, 0)
    )

tikz_plot(gg1, 'fig1a', w = 4, h = 4, cd = 'figs')

fit_sel <- function(df, sel) {
    ft <- fit_simple_gp(df, sel)
    list(
        sel = sel,
        ft = ft,
        tru = sum(df$resp[-sel]),
        pred = samp_inc(ft, sel)
    )
}

plot_sel <- function(df, sel) {
    probs <- sel$ft$summary.fitted.values$mean
    dfr <- mutate(df, r = probs)

    gg <- ggplot(dfr, aes(x, y)) +
        geom_point(data = dat[sel$sel,], shape = 1, size = 7) +
        geom_point(aes(col = r, size = r)) +
        scale_color_gradient(low = 'blue', high = 'red') +
        scale_shape_manual(values = c(`Not infested` = 16, `Infested` = 4)) +
        theme +
        theme(legend.position = 'none')

    print(gg)
}

sel1a <- fit_sel(dat, c(ini_sel, 10, 45, 11, 33, 47))
sel2a <- fit_sel(dat, c(ini_sel, 44, 23, 36, 14, 42))

tikz_plot(plot_sel(dat, sel1a), 'fig1b', 5, 5)
tikz_plot(plot_sel(dat, sel2a), 'fig1c', 5, 5)

sel1b <- fit_sel(dat, c(sel1a$sel, 22, 41, 2, 37, 21, 8, 38, 50, 24))
sel2b <- fit_sel(dat, c(sel2a$sel, 12, 4, 9, 25, 7, 5, 46, 21, 13))

tikz_plot(plot_sel(dat, sel1b), 'fig1d', 5, 5)
tikz_plot(plot_sel(dat, sel2b), 'fig1e', 5, 5)

dpred <- bind_rows(
    tibble(pred = sel1b$pred, sel = 'S1', tru = sel1b$tru),
    tibble(pred = sel2b$pred, sel = 'S2', tru = sel2b$tru)
)

gg2 <- ggplot(dpred, aes(pred, col = sel)) +
    geom_density() +
    geom_vline(aes(xintercept = tru), linetype = 'dashed') +
    labs(x = 'Remaining infestations') +
    ## geom_text(aes(y = density, label = sel), predlab) +
    theme_bw() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        ## panlel.grid = element_blank(),
        panel.background = element_rect(fill = gray(0.95)),
        legend.position = 'none'
    )

tikz_plot(gg2, 'fig1f', w = 4, h = 3.4)

### Utility function examples---------------------------------------------------

ts <- expand.grid(prop = seq(0, 1, 0.01), alpha = c(0, .15, .3, .7, 1, 2)) %>%
    mutate(t = prop^alpha)

ggplot(ts, aes(prop, t, col = as.factor(alpha))) +
    geom_line() +
    labs(x = '$(m_i - m_1) / (n - m_1)$', y = '$t(i; \\alpha)$', col = NULL) +
    theme_bw()

x <- expand.grid(risk = seq(0, 1, 0.05), t = seq(0, 1, 0.05)) %>%
    mutate(score = t^1.5 * risk + (1 - t^1.5) * 1)

ggplot(x, aes(t, risk, z = score)) +
    geom_contour_filled()
