require(tidyverse)
require(spBayes)
require(posterior)

source('spatial-helpers.R')

## TODO general func for labelling houses in space, colored by value
plot_lat_long_id <- function(value) {}

fit_spam <- function(f, df_train, df_test) {
    fit <- fitme(
        update(base_f, ~ . + Matern(1|long+lat)),
        data = dat_train,
        family = binomial(),
        fixed = list(nu = 1)
    )
}

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

base_f <- infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
        bedroom_clutter + dark_main_room + floor_type + education_level +
        construction_pile_type + firewood_location + grain_storage_in_house + house_age +
        house_hygiene + kitchen_location + land_for_agriculture + num_humans +
        num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
        condition_bedroom_wall + condition_house_wall + material_house_wall +
        material_roof + windows_in_bedroom + dist_perim

dat_org$long <- jitter(dat_org$long)

dat_sub <- filter(dat_org, village == 'Paternito')
train <- sample(nrow(dat_sub), round(.5 * nrow(dat_sub)))
dat_train <- select(dat_sub[train,], -id, -village)
dat_test <- select(dat_sub[-train,], -id, -village)

dat_proj <- sp_project(dat_train)

fit <- spGLM(
    base_f,
    data = dat_proj@data,
    coords = coordinates(dat_proj),
    starting = list(phi = 1, sigma.sq = 1, 'beta' = rep(0, 38), w = 0),
    tuning = list('phi' = 0.5, 'sigma.sq' = 0.5, 'beta' = rep(.1, 38), 'w' = .5),
    priors = list(
        beta.normal = list(rep(0, 38), rep((5/3)^2, 38)),
        sigma.sq.ig = c(2, 1),
        'phi.Unif' = c(.1, 3)
    ),
    amcmc = list(
        n.batch = 50,
        batch.length = 100,
        accept.rate = .45
    ),
    cov.model = "exponential"
    ## verbose = TRUE,
    ## n.report = 10
)

hist(fit$p.beta.theta.samples[,3])

z_post <- fit$p.w.samples %>%
    as_draws_df %>%
    summarise_draws

z_post %>%
    mutate(

pp <- spPredict(fit, coordinates(dat_proj), 

yhat <- pred$p.y.predictive.samples %>%
    t %>%
    as_tibble

colnames(yhat) <- dat_sub[-train,]$id

bla <- yhat %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarize(med_r = median(value)) %>%
    mutate(
        infestation = dat_test$infestation,
        long = dat_test$long,
        lat = dat_test$lat
    )

ggplot(dat_train, aes(long, lat, shape = as.factor(infestation))) +
    geom_point(col = 'green', size = 2.5) +
    geom_point(aes(col = med_r), data = bla, size = 2.5)
