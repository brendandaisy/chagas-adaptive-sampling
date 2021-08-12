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
summ <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x))
summ2 <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp2))

summ0b <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, fit_fun = fit_iid, pred_type = 'known')
)

summb <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, pred_type = 'known')
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

gg <- res %>%
    pivot_longer(cpu:mlik) %>%
    ggplot(aes(interaction(pred, village), value, fill = model)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~name, scales='free', nrow = 3) +
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
        y = 'Density $p(\\rho \\mid \\mathscr{S})$',
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

## summ_fxe <- function(ft) {
##     imap_dfr(ft$marginals.fixed[-c(1, 3)], ~{
##         ci5 <- inla.hpdmarginal(0.5, .x)
##         ci95 <- inla.hpdmarginal(0.95, .x)
##         tibble(
##             var = .y,
##             m = inla.zmarginal(.x, TRUE)$mean,
##             lo5 = ci5[,1],
##             hi5 = ci5[,2],
##             lo95 = ci95[,1],
##             hi95 = ci95[,2]
##         )
##     })
## }

summ_fxe <- function(summ_df, village) {
    sumv <- filter(summ_df, village == !!village)
    imap_dfr(sumv$fix[[1]][-1], ~{
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

sim_gp <- function(pts) {
    dperim <- function(pt) {
        min(abs(pt$x - 1), abs(pt$x - 0), abs(pt$y - 1), abs(pt$y - 0))
    }
    n <- nrow(pts)
    cov <- matrix(0, n, n)
    dprm <- rep(0, n)
    for(i in 1:n) {
        dprm[i] <- dperim(pts[i,])
        for(j in i:n) {
            p1 <- pts[i,]
            p2 <- pts[j,]
            d <- sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2)
            cov[i, j] <- 0.25 * INLA:::inla.matern.cf(d, range = 0.3, nu = 1);
            cov[j, i] <- cov[i, j]
        }
    }
    resp <- mvrnorm(1, mu = -2 * dprm + rep(0.2, n), Sigma = cov) %>%
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

set.seed(123)

## dat <- tibble(x = runif(50), y = runif(50), idx = 1:50) %>%
##     sim_gp
dat <- read_csv('sim-gp.csv')

ini_sel <- which(dat$sel == 'init')
ini_sel <- c(ini_sel[ini_sel != 5], 3)
## dat <- mutate(dat, sel = ifelse(1:n() %in% ini_sel, 'init', 'no'))

theme <- theme_bw() +
    theme(
        panel.background = element_rect(fill = gray(0.95)),
        axis.title = element_blank(),
        ## axis.text = element_blank(),
        ## axis.ticks = element_blank()
        ## panel.grid = element_blank()
    )

p1 <- ggplot(dat, aes(x, y)) +
    geom_point(aes(shape = factor(resp, labels = c('Not infested', 'Infested')))) +
    geom_text(aes(label = idx)) +
    scale_shape_manual(values = c(`Not infested` = 16, `Infested` = 4)) +
    ## geom_point(data = dat[sel1,], shape = 1, size = 3) +
    theme +
    theme(
        legend.position = c(.13, .93),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, 'pt'), 
        legend.background = element_rect(fill = gray(0.95), color = 'black'),
        legend.key = element_rect(fill = gray(0.95)),
        ## legend.justification = 'center',
        legend.margin = margin(.5, .5, .5, .5)
    )

tikz_plot(p1, 'fig1a', w = 4, h = 4)

pref_sel1 <- c(17, 45, 10, 24, 26, 33, 44, 38, 32, 31)
sel1 <- c(pref_sel1, ini_sel)
ft1 <- fit_simple_gp(dat, sel1)
tru1 <- sum(dat$resp[-sel1])
pred1 <- samp_inc(ft1, sel1)

pref_sel2 <- c(8, 33, 47, 26, 48, 46, 2, 22, 41, 3, 23, 25, 43, 5, 36, 42)
sel2 <- c(pref_sel2, ini_sel)
ft2 <- fit_simple_gp(dat, sel2)
tru2 <- sum(dat$resp[-sel2])
pred2 <- samp_inc(ft2, sel2)

dpred <- bind_rows(
    tibble(pred = pred1, sel = 'S1', tru = tru1),
    tibble(pred = pred2, sel = 'S2', tru = tru2)
)

predlab <- tribble(
    ~pred, ~density, ~sel,
    8, 0.07, 'S1', 
    23.7, 0.07, 'S2'
)

p2 <- dat %>%
    mutate(
        r = ft1$summary.fitted.values$mean,
        sel = 1:n() %in% sel1
    ) %>%
    ggplot(aes(x, y)) +
    geom_point(data = dat[sel1,], shape = 1, size = 5.5) +
    geom_point(
        aes(
            col = r,
            size = r,
            shape = factor(resp, labels = c('Not infested', 'Infested'))
        )
    ) +
    scale_color_gradient(low = 'blue', high = 'red') +
    scale_shape_manual(values = c(`Not infested` = 16, `Infested` = 4)) +
    theme +
    theme(legend.position = 'none')

tikz_plot(p2, 'fig1b', w = 5, h = 5)

p3 <- dat %>%
    mutate(
        r = ft2$summary.fitted.values$mean,
        sel = 1:n() %in% sel1
    ) %>%
    ggplot(aes(x, y)) +
    geom_point(data = dat[sel2,], shape = 1, size = 5.5) +
    geom_point(
        aes(
            col = r,
            size = r,
            shape = factor(resp, labels = c('Not infested', 'Infested'))
        )
    ) +
    scale_color_gradient(low = 'blue', high = 'red') +
    scale_shape_manual(values = c(`Not infested` = 16, `Infested` = 4)) +
    theme +
    theme(legend.position = 'none')

tikz_plot(p3, 'fig1c', w = 5, h = 5)

p4 <- ggplot(dpred, aes(pred, col = sel)) +
    geom_density() +
    geom_vline(aes(xintercept = tru), linetype = 'dashed') +
    labs(x = 'Remaining infestations') +
    ## geom_text(aes(y = density, label = sel), predlab) +
    theme_bw() +
    theme(
        ## axis.title = element_blank(),
        ## axis.text = element_blank(),
        ## axis.ticks = element_blank(),
        ## panlel.grid = element_blank(),
        panel.background = element_rect(fill = gray(0.95)),
        legend.position = 'none'
    )

tikz_plot(p4, 'fig1d', w = 5, h = 5)

fxe_hpd <- function(ft) {
    imap_dfr(ft$marginals.fixed, ~{
        ci5 <- inla.hpdmarginal(0.5, .x)
        ci95 <- inla.hpdmarginal(0.95, .x)
        tibble(
            beta = .y,
            m = inla.zmarginal(.x, TRUE)$mean,
            lo5 = ci5[,1],
            hi5 = ci5[,2],
            lo95 = ci95[,1],
            hi95 = ci95[,2]
        )
    })
}

fx1 <- fxe_hpd(ft1) %>%
    mutate(
        beta = ifelse(beta == 'dperim', '$\\beta_1$', '$\\beta_0$'),
        samp = 's1'
    )

fx2 <- fxe_hpd(ft2) %>%
    mutate(
        beta = ifelse(beta == 'dperim', '$\\beta_1$', '$\\beta_0$'),
        samp = 's2'
    )

bind_rows(fx1, fx2) %>%
    ggplot(aes(beta, m, col = samp)) +
    ## geom_hline(xintercept = -2, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(ymin = lo95, ymax = hi95), size = 1, position = position_dodge2(.5)) +
    geom_linerange(aes(ymin = lo5, ymax = hi5), size = 2, position = position_dodge2(.5)) +
    geom_point(size = 1.5, position = position_dodge2(.5)) +
    labs(x = 'Log odds') +
    theme_bw() +
    theme(axis.title = element_blank())

## ggplot(fixed, aes(x, y, col = samp)) +
##     geom_vline(xintercept = -2, linetype = 'dashed') +
##     geom_line() +
##     xlim(-10, 10) +
##     theme_bw()

## range <- ft1$marginals.hyperpar[[1]] %>%
##     as_tibble %>%
##     mutate(samp = 's1') %>%
##     bind_rows(mutate(as_tibble(ft2$marginals.hyperpar[[1]]), samp = 's2'))

## ggplot(range, aes(x, y, col = samp)) +
##     geom_vline(xintercept = 0.3, linetype = 'dashed') +
##     geom_line() +
##     xlim(0, 4) +
##     theme_bw()

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
