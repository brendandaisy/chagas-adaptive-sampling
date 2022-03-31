require(sf)
require(tidyverse)
require(magrittr)
require(ggthemes)
require(geosphere)
require(gstat)
require(rgdal)
require(randomForest)

source('../spatial-helpers.R')
source('../covariate-helpers.R')

dat_org <- prep_model_data('../../data-raw/gtm-tp-mf.csv') %>%
    mutate(infestation = as.factor(infestation))

dat_sub <- sp_project(filter(dat_org, village == 'Guayabo'))
s <- sample(nrow(dat_sub), nrow(dat_sub) %/% 2)
dtr <- dat_sub[s,]
dp <- dat_sub[-s,]

rf <- randomForest(
    infestation ~ . - id - long - lat - village - dist_perim,
    data = dtr,
    mtry = 3
)

dtr@data <- dtr@data %>%
    mutate(p = predict(rf, newdata = dtr, type = 'prob')[,2])

vg <- variogram(p ~ 1, dtr, cutoff = 1000, width = 1000 / 15)
plot(vg)
vg_fit <- fit.variogram(vg, vgm('Mat'))
plot(vg, vg_fit)

x <- krige(p ~ 1, proj_rf, model = vg_fit, newdata = sp_project(dat_sub[-s,]))

dat_sub$kr[s] <- proj_rf@data$p
dat_sub$kr[-s] <- x$var1.var

dat_sub %>%
    mutate(
        s = 1:n() %in% s
    ) %>%
    ggplot(aes(long, lat, col = infestation, size = kr, shape = s)) +
    geom_point()

dat_sub %>%
    mutate(
        p = predict(rf, newdata = dat_sub, type = 'prob')[,2]
    ) %>%
    ggplot(aes(long, lat, size = p, col = p)) +
    geom_point(alpha = .75)

### Plot GLS residual variogram for combined data
base_f <- infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
    bedroom_clutter + dark_main_room + floor_type + education_level +
    construction_pile_type + firewood_location + grain_storage_in_house + house_age +
    house_hygiene + kitchen_location + land_for_agriculture + num_humans +
    num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
    condition_bedroom_wall + condition_house_wall + material_house_wall +
    material_roof + windows_in_bedroom + dist_perim

dat_proj <- sp_project(dat_org)

vg <- variogram(base_f, dat_proj2, cutoff = 3000, width = 3000 / 25)
plot(vg)

vg_fit <- fit.variogram(vg, vgm('Sph'))

gs <- gstat(
    NULL, "infestation", base_f, dat_proj2,
    model = vg_fit,
    set = list(gls = 1)
)

vg_gls <- variogram(gs, cutoff = 3000, width = 3000 / 25)
plot(vg_gls)
