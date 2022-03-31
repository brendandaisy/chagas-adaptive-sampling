require(INLA)
require(tidyverse)
require(gridExtra)

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:ownership),
        as_factor
    )

base_f <- infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
        bedroom_clutter + dark_main_room + floor_type + # no education level for now
        construction_pile_type + firewood_location + grain_storage_in_house + house_age +
        house_hygiene + kitchen_location + land_for_agriculture + num_humans +
        num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
        condition_bedroom_wall + condition_house_wall + material_house_wall +
        material_roof + windows_in_bedroom

m0 <- inla(
    base_f,
    family = 'logistic',
    data = dat_sub,
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)
)

m0r <- inla(
    infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
        bedroom_clutter + dark_main_room + education_level + floor_type +
        construction_pile_type + firewood_location + grain_storage_in_house + house_age +
        house_hygiene + kitchen_location + land_for_agriculture + num_humans +
        num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
        condition_bedroom_wall + condition_house_wall + material_house_wall +
        material_roof + windows_in_bedroom + f(id),
    family = 'logistic',
    data = dat_org,
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)
)

m1 <- inla(
    infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
        bedroom_clutter + dark_main_room + education_level + floor_type +
        construction_pile_type + firewood_location + grain_storage_in_house + house_age +
        house_hygiene + kitchen_location + land_for_agriculture + num_humans +
        num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
        condition_bedroom_wall + condition_house_wall + material_house_wall +
        material_roof + windows_in_bedroom + f(village),
    family = 'logistic',
    data = dat_org,
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)
)

summary(m0)
summary(m0r)
summary(m1)
summary(m2)

iso <- which(rowSums(dist_mat(dat_org, cutoff = 100)) == 0)
dat_sub <- dat_org[-iso,] # the data for this script

dm <- dist_mat(dat_sub, cutoff = 100)

m2 <- inla(
    update(
        base_f,
        ~ . + f(id, model = 'besagproper', graph = dm, values = unique(dat_sub$id))
    ),
    family = 'binomial',
    data = dat_sub,
    control.predictor = list(link = 1),
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)
)

dat_sub <- filter(dat_org, village == 'Paternito')
iso <- which(rowSums(dist_mat(dat_sub, cutoff = 100)) == 0)
dat_sub <- dat_sub[-iso,]
dm <- dist_mat(dat_sub)
dmax <- 250

dcar_model <- inla.rgeneric.define(
    dcar_inla,
    A = dm,
    k = .065,
    dmax = dmax
)

dcar_fit <- inla(
    update(base_f, ~ . + f(id, model = dcar_model, values = unique(dat_sub$id))),
    family = 'binomial',
    data = dat_sub,
    control.predictor = list(link = 1, compute = TRUE),
    ## control.inla = list(strategy = 'laplace'),
    control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
    control.inla = list(control.correct = list(enable = 10))
)

summary(dcar_fit)

plot_hyper(dcar_fit)

