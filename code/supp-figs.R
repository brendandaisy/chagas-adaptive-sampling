require(tidyverse)
require(tikzDevice)
require(MASS, exclude = c('select'))

source('adaptive-sampling-fn.R')

options(
    tikzLatexPackages = c(getOption( "tikzLatexPackages" ), '\\usepackage{bm}')
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

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds') |> 
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
    mutate(summ0, model = '$Z(s)$ removed', pred = 'global'),
    mutate(summ, model = '$\\varepsilon(s)$ removed', pred = 'global'),
    mutate(summ2, model = 'Full model', pred = 'global'),
    mutate(summ0b, model = '$Z(s)$ removed', pred = 'all'),
    mutate(summb, model = '$\\varepsilon(s)$ removed', pred = 'all'),
    mutate(summ2b, model = 'Full model', pred = 'all')
)

# count cases where each method did best
res %>%
    group_by(village, pred) %>%
    slice_max(mlik, with_ties = TRUE) %>%
    ungroup %>%
    count(model)

gg <- res |> 
    rename(`CPU(sec)` = cpu, DIC = dic, ML = mlik) %>%
    pivot_longer(`CPU(sec)`:ML) %>%
    ggplot(aes(village, value, fill = model)) +
    geom_col(position = position_dodge()) +
    facet_grid(
        fct_relevel(name, 'DIC', 'ML') ~ factor(pred, c('global', 'all'), labels = c('Global only', 'All')), 
        scales = 'free',
        switch = 'y' 
    ) +
    labs(x = 'Village', y = NULL, fill = NULL) +
    theme_bw() +
    theme(
        # aspect.ratio = 1,
        strip.background.y = element_blank(),
        strip.placement = "outside"
    )

tikz_plot(gg, 'sens-study', 8, 4, cd = 'figs')

### Effective range for both predictor sets-------------------------------------

rho_glo <- map_dfr(
    summ2$rho, 
    ~mutate(as_tibble(inla.hpdmarginal(0.95, .x)), mode = inla.mmarginal(.x))
)

rho_all <- map_dfr(
    summ2b$rho, 
    ~mutate(as_tibble(inla.hpdmarginal(0.95, .x)), mode = inla.mmarginal(.x))
)

diam <- dat_org |> 
    group_by(village) |> 
    nest() |>
    rowwise() |> 
    mutate(diam = max(dist_mat(data)))

rgm <- mutate(rho_glo, across(everything(), ~.x * diam$diam))
ram <- mutate(rho_all, across(everything(), ~.x * diam$diam))

rgm$mode - ram$mode

# lam <- -log(0.05) * 0.1

gg <- bind_rows(
    mutate(marg_rho(summ2), pred = 'Global only'),
    mutate(marg_rho(summ2b), pred = 'All')
) %>%
    ggplot(aes(x, y, col = v)) +
    geom_line() +
    # geom_function(aes(x), fun = ~lam * .x^(-2) * exp(-lam / .x), inherit.aes = FALSE) +
    facet_wrap(~fct_relevel(pred, 'Global only'), nrow = 1, scale = 'free') +
    xlim(0, 4) +
    labs(
        x = 'Effective range $\\rho$',
        y = 'Density $p(\\rho \\mid \\bm{y})$',
        y = NULL,
        col = NULL
    ) +
    theme_bw()

tikz_plot(gg, 'range-full', 7.75, 3.25, cd = 'figs')

### Fixed effects for the best model--------------------------------------------

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

fx_glo <- map_dfr(unique(dat_org$village), ~summ_fxe(summ2, .x))
fx_all <- map_dfr(unique(dat_org$village), ~summ_fxe(summ2b, .x))

## which vars were most consistent/important?
fx_all %>%
    filter(hi5 < 0) %>%
    count(var, sort = TRUE) %>%
    filter(n >= 3)

# translate the factor levels

fx_all$var <- fx_all$var |> 
    str_replace('buena', ' - good') |> 
    str_replace('si$', ' - yes') |>
    str_replace('afuera_casa', ' - outside') |> 
    str_replace('adj_inside|directly_outside', ' - adjacent') |>
    str_replace('ninguno$', ' - none') |>
    str_replace('no$', ' - no') |>
    str_replace('no_acumula', ' - no') |> 
    str_replace('7_or_more', ' - 7 or more') |> 
    str_replace('one_or_less', ' - 1 or less') |>
    str_replace('rented', ' - rented') |>
    str_replace('TRUE', ' - yes') |>
    str_replace('FALSE', ' - no') |>
    str_replace('buen_estado', ' - good') |>
    str_replace('bajareque', ' - bajareque') |>
    str_replace('palopique', ' - palopique') |> 
    str_replace('brick_block_oth', ' - brick or block') |> 
    str_replace('ladrillo', ' - other') |>
    str_replace('wood_ceram_metal', ' - wood or metal') |> 
    str_replace('directly_outside', ' - outside') |> 
    str_replace('adentro_casa', ' - inside') |> 
    str_replace('vegetal_clay', ' - vegetal or clay') |> 
    str_replace('num_', 'num. ') |> 
    str_replace('dist_perim', 'dist. perimeter') |> 
    str_to_sentence() |> 
    str_replace_all('_', ' ')

## plot ceofs for either predictor set
gg <- fx_all %>%
    mutate(var = str_replace_all(var, '\\_', '\\\\_')) %>%
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col = 'darkorange2') +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col = 'darkorange2') +
    geom_point(size = 1.5, col = 'dodgerblue3') +
    facet_wrap(~village, nrow = 1) +
    labs(x = 'Log odds', y = NULL) +
    # xlim(-5.8, 4) +
    theme_bw() +
    theme(panel.spacing = unit(.4, "lines"))

tikz_plot(gg, 'fixed-all', 7.6, 5, cd = 'figs')

### Plot of two measures in utility---------------------------------------------

dat_sub <- filter(dat_org, village == 'Paternito')
obs <- sample(nrow(dat_sub), 50)
fit_dat <- build_pred_dat(dat_sub, obs, 'global')
ft <- fit_old(df = fit_dat)

plt_dat <- dat_sub %>%
    mutate(
        mean = ft$summary.fitted.values$mean, 
        var = ft$summary.linear.predictor$sd^2, 
        obs = 1:n() %in% obs, 
        idx = 1:n()
    )

gg1 <- ggplot(plt_dat, aes(long, lat, col = as.factor(infestation), size = mean)) +
    geom_point(alpha = .7, data = filter(plt_dat, !obs), col = 'thistle4') +
    geom_point(data = filter(plt_dat, obs), shape = 10, stroke = 1.2) +
    scale_color_manual(guide = FALSE, values = c('steelblue', 'orangered1')) +
    labs(x = '', y = '', col = NULL) +
    guides(size = guide_legend(title = 'Expected prob. infestation', nrow = 1)) +
    theme_bw() +
    theme(legend.position='bottom')

gg2 <- ggplot(plt_dat, aes(long, lat, col = as.factor(infestation), size = var)) +
    geom_point(alpha = .7, data = filter(plt_dat, !obs), col = 'thistle4') +
    geom_point(data = filter(plt_dat, obs), shape = 10, stroke = 1.2) +
    scale_color_manual(guide = FALSE, values = c('steelblue', 'orangered1')) +
    labs(x = '', y = '', col = NULL) +
    guides(size = guide_legend(title = 'Variance log odds infestation', nrow = 1)) +
    theme_bw() +
    theme(legend.position='bottom')

tikz_plot(grid.arrange(gg1, gg2, nrow = 1), 'ex-util', 10.2, 5.25, 'figs')

## plot of U for these points, for var alpha

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
    geom_text(aes(label = idx)) +
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

sel_prm <- purrr::accumulate(
    pull(arrange(dat, dperim), idx)[1:30], ~c(.x, .y), .init = ini_sel
)

sel_prm <- purrr::accumulate(
    pull(arrange(dat, dperim), idx)[1:30], ~c(.x, .y), .init = ini_sel
)

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

gg <- ggplot(ts, aes(prop, t, col = as.factor(alpha))) +
    geom_line() +
    labs(x = '$(m_i - m_1) / (n - m_1)$', y = '$t(i; \\alpha)$', col = '$\\alpha$') +
    theme_bw() +
    theme(
        legend.position = c(.9, .26),
        legend.background = element_rect(color = 'gray70')
    )

tikz_plot(gg, 'util-weight', 5, 4.2, 'figs')

x <- expand.grid(risk = seq(0, 1, 0.05), t = seq(0, 1, 0.05)) %>%
    mutate(score = t^1.5 * risk + (1 - t^1.5) * 1)

ggplot(x, aes(t, risk, z = score)) +
    geom_contour_filled()
