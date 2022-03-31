require(tidyverse)
require(gridExtra)
require(ggthemes)
require(ggpubr)

source('adaptive-sampling-fn.R')
source('as-fn-helpers.R')

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds')

### Random sampling functions---------------------------------------------------

plot_sel <- function(df, r) {
    df %>%
        mutate(samp = 1:n() %in% r) %>%
        ggplot(aes(long, lat, col = samp)) +
        geom_point(size = 3, alpha = .6) +
        theme_bw()
}

dat_sub <- filter(dat_org, village == 'Cerrón')
set.seed(124)

r <- rand_inhib(dat_sub, m = ceiling(.25 * nrow(dat_sub)), neigh = 50)
r <- rand_icp(dat_sub, m = ceiling(.25 * nrow(dat_sub)), neigh = 100)
r <- rand_pair(dat_sub, ceiling(.75 * nrow(dat_sub)))
plot_sel(dat_sub, r)

## prob of weighted selection for all houses
dm <- dist_mat(dat_sub, cutoff = 100)
w <- rowSums(dm) + 1
plot(w, w^.5 / sum(w^.5))

## weighted sampling only
r <- rand_dens(dat_sub, m = ceiling(.5 * nrow(dat_sub)), alpha = .5)
plot_sel(dat_sub, r)

r <- rand_dens(dat_sub, m = ceiling(.25 * nrow(dat_sub)), alpha = 1)
plot_sel(dat_sub, r)

### Convex tradeoff function----------------------------------------------------

x <- expand.grid(risk = seq(0, 1, 0.05), t = seq(0, 1, 0.05)) %>%
    mutate(score = t^1.5 * risk + (1 - t^1.5) * 1)

ggplot(x, aes(t, risk, z = score)) +
    geom_contour_filled()

### Iteration rounding----------------------------------------------------------

dat_sub <- filter(dat_org, village == 'Paternito')
is1 <- iterative_sampling(dat_sub, 'tar_kld', 20, interp = FALSE, silent = FALSE)
is2 <- iterative_sampling(dat_sub, 'rand_dens', 20, interp = FALSE)
length(is1$obs_idx[[3]])
length(is2$obs_idx[[1]])

dat_sub <- filter(dat_org, village == 'Cerrón')
is1 <- iterative_sampling(dat_sub, 'tar_kld', 17, interp = FALSE, silent = FALSE)
is2 <- iterative_sampling(dat_sub, 'rand_dens', 17, interp = FALSE)

length(is1$obs_idx[[32]])
length(is2$obs_idx[[3]])

### Sanity check for kriging methods--------------------------------------------

dat_sub <- filter(dat_org, village == 'Paternito')
## s <- sample(nrow(dat_sub), 35)
i1 <- build_pred_dat(dat_sub, s, 'interp')

sp_norm <- sp_project(dat_sub, normalize = TRUE)
sp_obs <- sp_norm[s,]
nb <- dnearneigh(sp_obs, 0, 1 / 3) # cutoff @ 1/3 obs diameter
lw <- nb2listwdist(nb, sp_obs, type = 'exp', alpha = -log(0.05) * 3, zero.policy = TRUE)
jj <- joincount.mc(sp_obs@data$bed_hygiene, lw, 1000)

sp_norm <- sp_project(dummify_vars(dat_sub), normalize = TRUE)
sp_obs <- sp_norm[s,]
vg <- variogram(material_house_wall.adobe ~ 1, sp_obs, cutoff = 1 / 3)
vgf <- fit.variogram(vg, vgm('Mat', range = 1 / 6), fit.kappa = TRUE)
plot(vg, vgf)

x <- krige(
        material_house_wall.adobe ~ 1,
        sp_obs,
        model = vgf,
        newdata = sp_norm[-s,]
)

spplot(sp_obs, z = 'material_house_wall')
spplot(x, z = 'var1.pred')

ptheme <- theme_few() +
    theme(
        legend.background = element_rect(linetype = "solid", color = 'black'),
        legend.box.margin = margin(),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )

p1 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>%
    filter(Observed == TRUE) %>%
    ggplot(aes(long, lat, col = material_house_wall, shape = Observed)) +
    geom_point(size = 3, alpha = .85) +
    scale_shape_manual(values = c(`TRUE` = 10, `FALSE` = 19), guide = FALSE) +
    ptheme +
    theme(legend.position = c(.83, .12))

p2 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>%
    ggplot(aes(long, lat, col = floor_type, shape = Observed)) +
    geom_point(size = 3, alpha = .85) +
    scale_shape_manual(values = c(`TRUE` = 10, `FALSE` = 19), guide = FALSE) +
    ptheme +
    theme(legend.position = c(.83, .12))

p3 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>%
    ggplot(aes(long, lat, col = grain_storage_in_house, shape = Observed)) +
    geom_point(size = 3, alpha = .85) +
    scale_shape_manual(values = c(`TRUE` = 10, `FALSE` = 19), guide = FALSE) +
    ptheme +
    theme(legend.position = c(.86, .12))

p4 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>%
    pivot_longer(starts_with('sign_r')) %>% filter(value == 1) %>%
    ggplot(aes(long, lat, col = name, shape = Observed)) +
    geom_point(size = 3, alpha = .85) + scale_shape_manual(values =
    c(`TRUE` = 10, `FALSE` = 19), guide = FALSE) + ptheme +
    theme(legend.position = c(.87, .12))

p5 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>%
    pivot_longer(starts_with('materi')) %>% filter(value == 1) %>%
    ggplot(aes(long, lat, col = name, shape = Observed)) +
    geom_point(size = 3, alpha = .85) + scale_shape_manual(values =
    c(`TRUE` = 10, `FALSE` = 19), guide = FALSE) + ptheme +
    theme(legend.position = c(.84, .12))

p6 <- i1 %>%
    mutate(Observed = 1:n() %in% s) %>% ggplot(aes(long, lat, col =
    num_humans, shape = Observed)) + geom_point(size = 3, alpha = .85)
    + scale_shape_manual(values = c(`TRUE` = 10, `FALSE` = 19), guide
    = FALSE) + scale_color_viridis_c(guide = guide_colorbar(direction
    = 'horizontal'), end = .9) + ptheme + labs(col = 'Num. residents')
    + theme(legend.position = c(.8, .1), legend.title =
    element_text())

grid.arrange(
    p1, p2, p3, p4, p5, p6,
    ncol = 2,
    top = text_grob('Final interpolation based on 35 observations', face = 'bold')
)

### Compare predictor fitting methods-------------------------------------------

interp1 <- function(df, obs_idx, pval) {
    dp <- sp_project(df)
    dp_obs <- dp[obs_idx,]
    dp_new <- dp[-obs_idx,]
    md <- max(dist_mat(df))
    nb <- dnearneigh(dp_obs, 0, md / 3)
    lw <- nb2listwdist(nb, dp_obs, alpha = 0.5, zero.policy = TRUE)
    
    vars <- dp_obs@data %>%
        select(bed_hygienemala:num_pigs) %>%
        select_if(~length(unique(.x)) > 1) %>%
        imap_dfr(
            ~tibble_row(var = .y, p = moran.mc(.x, lw, 1000, zero.policy = TRUE)$p.value)
        ) %>%
        filter(p < pval) %>%
        pull(var)

    interp <- NULL
    if (length(vars) > 0) {
        interp <- map_dfc(
            vars,
            ~tibble(!!.x := krige_var1(dp_obs, dp_new, .x, md / 3)$var1.pred)
        )
    }
    
    return(list(interp, vars)) # for interp 4
}

interp2 <- function(df, obs_idx) {
    dp <- sp_project(df)
    dp_obs <- dp[obs_idx,]
    dp_new <- dp[-obs_idx,]
    md <- max(dist_mat(df))

    map_dfc(
        colnames(select(df, bed_hygienemala:num_pigs)),
        ~tibble(!!.x := krige_var1(dp_obs, dp_new, .x, md / 3)$var1.pred)
    )
}

krige_var1 <- function(sp, new_sp, var, cut) {   
    vg <- variogram(
        as.formula(paste0(var, '~1')),
        sp,
        cutoff = cut,
        width = cut / 20
    )
    vgf <- fit.variogram(vg, vgm('Sph'))
    x <- krige(
        as.formula(paste0(var, '~1')),
        sp,
        model = vgf,
        newdata = new_sp,
        debug.level = 0
    )
    return(x@data)
}

krige_var2 <- function(sp, new_sp, var, cut) {
    vg <- variogram(
        as.formula(paste0(var, '~1')),
        sp,
        cutoff = cut,
        width = cut / 20
    )
    vgf <- fit.variogram(vg, vgm('Sph'))
    if (attr(vgf, 'singular'))
        return(NULL)
    x <- krige(
        as.formula(paste0(var, '~1')),
        sp,
        model = vgf,
        newdata = new_sp,
        debug.level = 0
    )
    return(x@data)
}

interp3 <- function(df, obs_idx, var_names) {
    dp <- sp_project(df)
    dp_obs <- dp[obs_idx,]
    dp_new <- dp[-obs_idx,]
    md <- max(dist_mat(df))

    rp <- imap_dfc(
        select(df, bed_hygienemala:num_pigs) %>%
        select_if(~length(unique(.x)) > 1), ~{
            x <- krige_var2(dp_obs, dp_new, .x, md / 3)$var1.pred
            if (!is.null(x))
                tibble(!!.y := x)
            else
                NULL
        })

    rinterp(rp, var_names)
}

interp4 <- function(df, obs_idx, vars) {
    dp <- sp_project(df)
    dp_obs <- dp[obs_idx,]
    dp_new <- dp[-obs_idx,]
    md <- max(dist_mat(df))

    interp <- NULL
    if (length(vars) > 0) {
        interp <- map_dfc(vars, ~{
            x <- krige_var2(dp_obs, dp_new, .x, md / 3)$var1.pred
            if (!is.null(x))
                tibble(!!.x := x)
            else
                NULL
        })
    }
    return(interp)
}

interp_mse <- function(tru, interp, method) {
    ret <- imap_dfc(interp, ~mean((pull(tru, .y) - .x)^2))
    colnames(ret) <- colnames(interp)
    mutate(ret, nvar = ncol(interp), method = method)
}

interp_test_rep <- function(village, prop) {
    ds <- filter(dat_cont, village == !!village)
    s <- sample(1:nrow(ds), ceiling(prop * nrow(ds)))
    i2 <- interp2(ds, s)
    i1a <- interp1(ds, s, 0.02)
    i1b <- interp1(ds, s, 0.1)
    i3 <- interp3(ds, s)
    i4a <- interp4(ds, s, i1a[[2]])
    i4b <- interp4(ds, s, i1b[[2]])
    bind_rows(
        interp_mse(ds, i2, 'i2'),
        interp_mse(ds, i1a[[1]], 'i1a'),
        interp_mse(ds, i1b[[1]], 'i1b'),
        interp_mse(ds, i3, 'i3'),
        interp_mse(ds, i4a, 'i4a'),
        interp_mse(ds, i4b, 'i4b')
    ) %>%
        mutate(v = village, p = prop)
}

reps <- expand.grid(v = unique(dat_cont$village), p = seq(.1, .9, .1))

### Compare interpolation MSEs

mse <- map2_dfr(reps$v, reps$p, interp_test_rep)

mse %>%
    select(-c(num_humans:num_pigs)) %>%
    pivot_longer(-c(nvar:p), 'Variable', values_to = 'MSE') %>%
    drop_na %>%
    ## rowwise %>%
    ## mutate(med = median(c_across(bed_hygienemala:num_pigs), na.rm = TRUE)) %>%
    ggplot(aes(p, MSE, col = v)) +
    ## geom_line() +
    ## geom_histogram(bins = 20) +
    geom_point(size = 1.8, shape = 1) +
    facet_grid(method ~ Variable, scales = 'fixed') +
    scale_x_continuous(breaks = c(.1, .5)) +
    theme_bw()

ggplot(mse, aes(p, nvar, col = method)) +
    geom_line() +
    geom_point() +
    facet_wrap(~v)

### Compare interpolation variance (guessing same everywhere is bad for adaptation)

sd_test_rep <- function(village, prop) {
    ds <- filter(dat_cont, village == !!village)
    s <- sample(1:nrow(ds), ceiling(prop * nrow(ds)))
    i2 <- interp2(ds, s)
    i1a <- interp1(ds, s, 0.02)
    i1b <- interp1(ds, s, 0.1)
    i3 <- interp3(ds, s)
    i4a <- interp4(ds, s, i1a[[2]])
    i4b <- interp4(ds, s, i1b[[2]])
    imap_dfr(
        list(i2 = i2, i1a = i1a[[1]], i1b = i1b[[1]], i3 = i3, i4a = i4a, i4b = i4b),
        ~tibble_row(!!!map(.x, sd), v = village, p = prop, m = .y)
    )
}

sd_test <- map2_dfr(reps$v, reps$p, var_test_rep)

sds <- sd_test %>%
    rowwise(v, p, m) %>%
    summarize(
        med = median(c_across(), TRUE),
        min = min(c_across(), na.rm = TRUE),
        max = max(c_across(), na.rm = TRUE)
    ) %>%
    ungroup %>%
    drop_na # these are settings that didn't pick any predictors

ggplot(sds, aes(max, fill = m)) +
    geom_histogram(bins = 35)

### Compare regression fit using only interpolated fixed effects

reg_test_rep <- function(village, prop) {
    ds <- filter(dat_cont, village == !!village)
    s <- sample(1:nrow(ds), ceiling(prop * nrow(ds)))
    ds$infestation[-s] <- NA
    i2 <- build_interp_dat(ds, s, interp2(ds, s))
    i1a <- build_interp_dat(ds, s, interp1(ds, s, 0.02))
    i1b <- build_interp_dat(ds, s, interp1(ds, s, 0.1))
    i3 <- build_interp_dat(ds, s, interp3(ds, s))
    imap_dfr(list(i2 = i2, i1a = i1a, i1b = i1b, i3 = i3), ~{
        f <- fit_fixed(.x)
        bt <- find_best_thresh(f, s, ds$truth)
        fit_score(f, s, ds$truth, thresh = bt) %>%
            mutate(bt = bt, v = village, p = prop, m = .y)
    })
}

fits <- map2_dfr(reps$v, reps$p, reg_test_rep)

fits %>%
    ## pivot_longer(sens:alt, 'Metric', values_to = 'Score') %>%
    ## filter(Metric != 'alt') %>%
    ggplot(
        aes(
            x = p,
            y = bt,
            ## shape = Metric,
            col = m
        )
    ) +
    geom_line(linetype = 'dotted') +
    geom_point(size = 2.5) +
    facet_wrap(~v)

### Finding Best Threshold------------------------------------------------------

thresh_test <- function(fit, obs_idx, tru_inf) {
    map_dfr(
        seq(.1, .9, .05),
        ~fit_score(fit, obs_idx, tru_inf, thresh = .x, in_samp = TRUE)
    )
}

dat_sub <- filter(dat_org, village == 'Paternito')
s <- sample(nrow(dat_sub), 35)
dat_sub$infestation[-s] <- NA
fd <- build_pred_dat(dat_sub, s, 'global')
ft <- fit_gp(fd)
fit_score(ft, s, dat_sub$truth)
plot_obs_surface(ft$summary.fitted.values, s, dat_sub)
scores <- thresh_test(ft, s, dat_sub$truth)

scores %>%
    rename(Sensitivity = sens, Accuracy = bacc) %>%
    pivot_longer(c(Sensitivity, Accuracy)) %>%
    ggplot(aes(thresh, value, col = name)) +
    geom_vline(xintercept = find_best_thresh(ft, s, dat_sub$truth)) +
    geom_hline(yintercept = quantile(scores$bacc, .9)) +
    geom_line() +
    geom_point(size = 2) +
    labs(x = 'Threshold', y = '', col = '', title = 'In Sample Performance')

## test some fake scores to make sure intuitive
sim <- tibble(
    sens = c(rep(1, 5), rep(0, 12)),
    acc = c(rep(.3, 5), rep(.6, 12)),
    thresh = seq(.1, .9, .05)
)

sim %>%
    pivot_longer(c(sens, acc)) %>%
    ggplot(aes(thresh, value, col = name)) +
    geom_hline(yintercept = quantile(sim$acc, .9)) +
    geom_point(size = 2) +
    geom_line()

## Compare find_best_thresh to actual out of sample performance
oos_thresh_test <- function(fit, obs_idx, tru_inf) {
    map_dfr(
        seq(.1, .9, .05),
        ~fit_score(fit, obs_idx, tru_inf, thresh = .x, in_samp = FALSE)
    )
}

oos_scores <- oos_thresh_test(ft, s, dat_sub$truth)

oos_scores %>%
    rename(Sensitivity = sens, Accuracy = bacc) %>%
    pivot_longer(c(Sensitivity, Accuracy)) %>%
    ggplot(aes(thresh, value, col = name)) +
    geom_vline(xintercept = find_best_thresh(ft, s, dat_sub$truth)) +
    geom_line() +
    geom_point(size = 2) +
    labs(x = 'Threshold', y = '', col = '', title = paste0('OOS Perf. for B=', length(s)))


### KLD Selection---------------------------------------------------------------

sim_kld <- function(m, s, p = 5) {
    post <- map2(m, s, ~rnorm(10000, .x, .y))
    pri <- rnorm(10000, 0, sum(rep(.3^(-1), p)))
    k <- map_dbl(post, ~mean(KL.divergence(.x, pri, k = 5)))

    tibble(mean = m, sd = s, k = k, dist = str_c('D', 1:length(m))) %>%
        arrange(k)
}

kld_form <- function(m, s, p = 5) {
    sdp <- sum(rep(.3^(-1), p))
    (m^2 + s^2) / (2 * sdp^2) + log(sdp / s) - 0.5
}

## large support
sim_kld(c(.1, 1, 30), c(10, 10, 6)) %>%
    mutate(tru = kld_form(mean, sd))

## shifted wide vs. centered shrunk
sim_kld(c(-6.5, -6, -2, 0), c(8, 8, 4, 1)) %>%
    mutate(tru = kld_form(mean, sd))

## shifting mean with sim. variance
sim_kld(c(-2, -10, -15, -50), c(2, 2, 2, 2)) %>%
    mutate(tru = kld_form(mean, sd))

### Bootstrap Seeding-----------------------------------------------------------

plan(multisession, workers = 2)

a <- bs_perf(dat_sub, 'random', dat_sub$truth, n_rep = 2, seed = 123)
b <- bs_perf(dat_sub, 'random', dat_sub$truth, n_rep = 2, seed = 123)
c <- bs_perf(dat_sub, 'random', dat_sub$truth, n_rep = 2, seed = 1211)
d <- bs_perf(dat_sub, 'random', dat_sub$truth, n_rep = 2, seed = 1211)
e <- bs_perf(dat_sub, tar_kld, dat_sub$truth, n_rep = 2, seed = 123)
f <- bs_perf(dat_sub, tar_kld, dat_sub$truth, n_rep = 2, seed = 123)

init_samp_test <- function(nobs, seed, arg) {
    future_map(
        1:10, ~{
            a <- arg
            sample(nobs, 10)
        },
        .options = furrr_options(seed = seed)
    )
}

init_samp_test(20, 69, 'a')



