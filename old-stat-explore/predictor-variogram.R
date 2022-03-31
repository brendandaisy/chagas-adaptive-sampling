require(sf)
require(tidyverse)
require(magrittr)
require(ggthemes)
require(geosphere)
require(gstat)
require(rgdal)

source('../spatial-helpers.R')
source('../covariate-helpers.R')

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

dat_org$long <- jitter(dat_org$long, factor = .5)
dat_org$lat <- jitter(dat_org$lat, factor = .5)

dat_cont <- make_dummy_mat(
    dat_org,
    contr = 'contr.treatment',
    inter = FALSE,
    as_tibble = TRUE
) %>%
    mutate(truth = infestation)

colnames(dat_cont) <- str_replace(colnames(dat_cont), ' ', '_')
colnames(dat_cont)[21] <- "house_age.7_year_or_more"
colnames(dat_cont)[22] <- "house_age.less_than_1_year"

moran <- function(df, village) {
    dfs <- filter(df, village == !!village)

    md <- max(dist_mat(dfs))
    dp <- sp_project(dfs)
    nb <- dnearneigh(dp, 0, md / 3)
    lw <- nb2listwdist(nb, dp, alpha = 0.5)

    dfs %>%
        select(bed_hygienemala:num_pigs) %>%
        select_if(~length(unique(.x)) > 1) %>%
        imap_dfr(
            ~tibble_row(var = .y, p = moran.mc(.x, lw, 1000)$p.value) %>%
                mutate(test = ifelse(p < 0.1, 'signif', 'notsignif'), village = village)
        )
}

krige_predictor <- function(dproj, var, cut, nbin, p = 0.5) {
    s <- sample(1:nrow(dproj@data), ceiling(p * nrow(dproj@data)))
    vg <- variogram(
        as.formula(paste0(var, '~1')),
        dproj[s,],
        cutoff = cut,
        width = cut / nbin
    )
    vgf <- fit.variogram(vg, vgm('Sph'))
    x <- krige(
        as.formula(paste0(var, '~1')),
        dproj[s,],
        model = vgf,
        newdata = dproj[-s]
    )
    list(k = x, vg = vg, vgf = vgf, obs = s)
}


morans <- map_dfr(unique(dat_cont$village), ~moran(dat_cont, .x))

ggplot(morans, aes(var, p, fill = test, col = test)) +
    geom_col() +
    facet_wrap(~village, nrow=2) +
    scale_fill_manual(values= c(signif='red4', notsignif=gray(.6, .8))) +
    scale_color_manual(values= c(signif='red1', notsignif='white')) +
    theme_few() +
    labs(x='Variables', y='P-val') +
    coord_flip()

dat_proj <- sp_project(filter(dat_cont, village == 'Guayabo'))

vg <- variogram(sign_ratsTRUE ~ 1, dat_proj, cutoff = 2000, width = 2000 / 20)
vg_fit <- fit.variogram(vg, vgm('Sph'))
plot(vg, vg_fit)

kp <- krige_predictor(dat_proj, 'sign_ratsTRUE', 2000, 20)

spplot(kp$k["var1.pred"], main = "ordinary kriging predictions")
spplot(kp$k["var1.var"],  main = "ordinary kriging variance")
spplot(dat_proj['sign_ratsTRUE'], main = "truth")

kp$k["var1.var"]

ggplot(filter(dat_cont, village == 'Amatillo'), aes(long, lat, col = num_pigs, size = num_pigs)) +
    geom_point(alpha = .85)
