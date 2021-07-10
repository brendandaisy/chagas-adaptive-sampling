require(tidyverse)
require(tikzDevice)

source('adaptive-sampling-fn.R')

options(
    tikzLatexPackages = c(getOption( "tikzLatexPackages" ), '\\usepackage{mathrsfs}')
    ## tikzMetricPackages = c(
    ##     "\\usepackage[utf8]{inputenc}",
    ##     "\\usepackage[T1]{fontenc}",
    ##     "\\usetikzlibrary{calc}"
    ## )
)

tikz_plot <- function(ggp, fname = 'tikz', w = 8.5, h = 4) {
    tikz(paste0(fname, '.tex'), standAlone=TRUE, width = w, height = h)
    print(ggp)
    dev.off()
    system(paste0('pdflatex ', fname, '.tex'))
}

full_vil_summ <- function(df, village, pred_type = 'global', fit_fun = fit_gp) {
    dfs <- filter(df, village == !!village)
    fit_dat <- build_pred_dat(dfs, 1, pred_type)
    mesh <- inla_mesh(dfs)
    spde <- inla.spde2.pcmatern(
        mesh = mesh,
        alpha = 2,
        prior.range = c(0.1, 0.05),
        prior.sigma = c(3, 0.1)
    )
    ft <- fit_fun(fit_dat, mesh, spde, control = list(mlik = TRUE, dic = TRUE))
    
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
    mutate(marg_rho(summ2), pred = 'global'),
    mutate(marg_rho(summ2b), pred = 'all')
) %>%
    ggplot(aes(x, y, col = v)) +
    geom_line() +
    facet_wrap(~pred, nrow = 1, scale = 'free') +
    xlim(0, 5) +
    labs(
        x = 'Effective range $\\rho$',
        y = 'Density $p(\\rho \\mid \\mathscr{S})$',
        col = NULL
    ) +
    theme_bw()

tikz_plot(gg, 'range-full', 7.75, 3.25)

fix <- imap_dfr(
    ft$marginals.fixed, ~{
        ci5 <- inla.hpdmarginal(0.5, .x)
        ci95 <- inla.hpdmarginal(0.95, .x)
        tibble(
            var = .y,
            m = inla.zmarginal(.x, TRUE)$mean,
            lo5 = ci5[,1],
            hi5 = ci5[,2],
            lo95 = ci95[,1],
            hi95 = ci95[,2]
        )
    })

ggplot(fix, aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col = 'darkorange2') +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col = 'darkorange2') +
    geom_point(size = 3, col = 'dodgerblue3') +
    labs(x = 'Log odds', y = NULL) +
    theme_bw()
    
