library(tidyverse)
library(tikzDevice)

source('code/other-helpers.R')
source('code/seq-sampling-helpers.R')

options(
    tikzLatexPackages = c(getOption( "tikzLatexPackages" ), '\\usepackage{bm}')
)

full_vil_summ <- function(df, village, pred_type = 'global', fit_fun = fit_gp_both) {
    dfs <- filter(df, village == !!village)
    fit_dat <- build_pred_dat(dfs, pred_type)
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

dat_org <- read_csv("survey-data.csv") |> 
    mutate(village = str_replace(village, 'ó', "\\\\'o"))

### Compare proposed model to simpler ones--------------------------------------

summ_rand_global <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp_random))
summ_spatial_global <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp_spatial))
summ_both_global <- map_dfr(unique(dat_org$village), ~full_vil_summ(dat_org, .x, fit_fun = fit_gp_both))

summ_rand_known <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, fit_fun = fit_gp_random, pred_type = 'known')
)

summ_spatial_known <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, pred_type = 'known', fit_fun = fit_gp_spatial)
)

summ_both_known <- map_dfr(
    unique(dat_org$village),
    ~full_vil_summ(dat_org, .x, fit_fun = fit_gp_both, pred_type = 'known')
)

res <- bind_rows(
    mutate(summ_rand_global, model = '$Z(s)$ removed', pred = 'global'),
    mutate(summ_spatial_global, model = '$\\varepsilon(s)$ removed', pred = 'global'),
    mutate(summ_both_global, model = 'Full model', pred = 'global'),
    mutate(summ_rand_known, model = '$Z(s)$ removed', pred = 'all'),
    mutate(summ_spatial_known, model = '$\\varepsilon(s)$ removed', pred = 'all'),
    mutate(summ_both_known, model = 'Full model', pred = 'all')
)

# count cases where each method did best
res %>%
    group_by(village, pred) %>%
    slice_max(dic, with_ties = TRUE) %>%
    ungroup %>%
    count(model)

## Figure 1: Model comparison
res |> 
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
    theme(strip.background.y = element_blank(), strip.placement = "outside")

### Effective range for both predictor sets-------------------------------------

## Table 1: Summary stats for effective range (not rescaled to meters here)
rho_glo <- map_dfr(
    summ_both_global$rho, 
    ~mutate(as_tibble(inla.hpdmarginal(0.95, .x)), mode = inla.mmarginal(.x))
)

rho_all <- map_dfr(
    summ_both_known$rho, 
    ~mutate(as_tibble(inla.hpdmarginal(0.95, .x)), mode = inla.mmarginal(.x))
)

## Figure S1: Posterior of effective range
bind_rows(
    mutate(marg_rho(summ_both_global), pred = 'Global only'),
    mutate(marg_rho(summ_both_known), pred = 'All')
) %>%
    ggplot(aes(x, y, col = v)) +
    geom_line() +
    facet_wrap(~fct_relevel(pred, 'Global only'), nrow = 1, scale = 'free') +
    xlim(0, 4) +
    labs(
        x = 'Effective range $\\rho$',
        y = 'Density $p(\\rho \\mid \\bm{y})$',
        y = NULL,
        col = NULL
    ) +
    theme_bw()

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

fx_glo <- map_dfr(unique(dat_org$village), ~summ_fxe(summ_both_global, .x))
fx_all <- map_dfr(unique(dat_org$village), ~summ_fxe(summ_both_known, .x))

# which vars were most consistent/important?
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

## Figure 2: Fixed effects from full village analysis
fx_all %>%
    mutate(var = str_replace_all(var, '\\_', '\\\\_')) %>%
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col = 'darkorange2') +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col = 'darkorange2') +
    geom_point(size = 1.5, col = 'dodgerblue3') +
    facet_wrap(~village, nrow = 1) +
    labs(x = 'Log odds', y = NULL) +
    theme_bw() +
    theme(panel.spacing = unit(.4, "lines"))

### Other figures not used in the paper-----------------------------------------

## Plot of two measures in utility (risk and variance)
dat_sub <- filter(dat_org, village == 'Paternito')
obs <- sample(nrow(dat_sub), 50)
fit_dat <- build_pred_dat(dat_sub, 'global')
ft <- fit_gp_spatial_dense(df = fit_dat)

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

grid.arrange(gg1, gg2, nrow = 1)

## Visualizing the proposed utility function 
ts <- expand.grid(prop = seq(0, 1, 0.01), alpha = c(0, .15, .3, .7, 1, 2)) %>%
    mutate(t = prop^alpha)

ggplot(ts, aes(prop, t, col = as.factor(alpha))) +
    geom_line() +
    labs(x = '$(m_i - m_1) / (n - m_1)$', y = '$t(i; \\alpha)$', col = '$\\alpha$') +
    theme_bw() +
    theme(
        legend.position = c(.9, .26),
        legend.background = element_rect(color = 'gray70')
    )

x <- expand.grid(risk = seq(0, 1, 0.05), t = seq(0, 1, 0.05)) %>%
    mutate(score = t^1.5 * risk + (1 - t^1.5) * 1)

ggplot(x, aes(t, risk, z = score)) +
    geom_contour_filled()
