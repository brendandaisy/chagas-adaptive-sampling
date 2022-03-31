require(sf)
require(tidyverse)
require(magrittr)
require(ggthemes)
require(geosphere)
require(INLA)
require(gstat)
require(rgdal)

source('../spatial-helpers.R')

rs_exper <- function(df_all, form, village, p = 0.5, nrep = 10, cut = 3000) {
    dfs <- filter(df_all, village == !!village)
    map_dfr(1:nrep, ~{
        dfss <- slice_sample(dfs, prop = p)
        v <- resid_vario(dfss, form, cutoff = cut, width = cut / 25)
        mutate(v, rep = .x)
    })
}

resid_vario <- function(df, form, ...) {
    dat_proj <- sp_project(df)

    dat_proj2 <- SpatialPointsDataFrame(
        jitter(coordinates(dat_proj), factor = 0.1), dat_proj@data
    )

    variogram(form, dat_proj2, ...)
}   

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

form <- infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
    bedroom_clutter + dark_main_room + floor_type + education_level +
    construction_pile_type + firewood_location + grain_storage_in_house + house_age +
    house_hygiene + kitchen_location + land_for_agriculture + num_humans +
    num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
    condition_bedroom_wall + condition_house_wall + material_house_wall +
    material_roof + windows_in_bedroom + dist_perim

## vg_fit <- fit.variogram(vg, vgm('Sph'))

## gs <- gstat(
##     NULL, "infestation", base_f, dat_proj2,
##     model = vg_fit,
##     set = list(gls = 1)
## )

## vg_gls <- variogram(gs, cutoff = 3000, width = 3000 / 25)
## plot(vg_gls)

pre <- rs_exper(dat_org, form, 'Prensa', nrep = 6, cut = 2500) %>%
    filter(gamma < 2) # since presum was separable

ggplot(pre, aes(dist, gamma, col = as.factor(rep))) +
    geom_line() +
    geom_point()

ama <- rs_exper(dat_org, form, 'Amatillo', nrep = 10, cut = 2000) %>%
    filter(gamma < 2)

ggplot(ama, aes(dist, gamma, col = as.factor(rep))) +
    geom_line() +
    geom_point()

gua <- rs_exper(dat_org, form, 'Guayabo', nrep = 6, cut = 2000) %>%
    filter(gamma < 2)

ggplot(ama, aes(dist, gamma, col = as.factor(rep))) +
    geom_line() +
    geom_point()
