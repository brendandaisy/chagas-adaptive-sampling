require(sf)
require(tidyverse)
require(magrittr)
require(ggthemes)
require(geosphere)
require(INLA)
require(gstat)
require(rgdal)

dprim <- function(v) {
    df <- filter(dat_org, village == v)

    pts <- select(df, long, lat) %>%
        st_as_sf(coords = c('long', 'lat'), crs=4326) %>%
        st_transform("+proj=utm +zone=16 ellps=WGS84")

    buff <- pts %>%
        st_union %>%
        st_convex_hull %>%
        st_buffer(50) %>%
        st_boundary

    return(as.double(st_distance(pts, buff)))
}

sp_project <- function(tbl) {
    dfsp <- as.data.frame(tbl)
    sp::coordinates(dfsp) <- ~long+lat
    proj4string(dfsp) <- CRS("+proj=longlat")
    spTransform(dfsp, CRS("+proj=utm +zone=16 +datum=WGS84"))
}

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

### Add distance from perimeter as predictor

dprims <- map(unique(dat_org$village), dprim)
names(dprims) <- unique(dat_org$village)

dat_org <- enframe(dprims, 'village', 'dist_perim') %>%
    unnest(cols = c(dist_perim)) %>%
    mutate(id = dat_org$id) %>%
    left_join(dat_org) %>%
    relocate(dist_perim, infestation, .after = last_col())

dat_sub <- filter(dat_org, village == 'Paternito')

dat_proj <- sp_project(dat_sub)

dat_proj <- SpatialPointsDataFrame(
    jitter(coordinates(dat_proj), factor = 0.1), dat_proj@data
)

mesh <- inla.mesh.2d(
    loc = coordinates(dat_proj),
    cutoff = 50,
    max.edge = 200,
    offset = c(100, 150)
)

plot(mesh)
points(coordinates(dat_proj), col = 'red')

spde1 <- inla.spde2.matern(
    mesh = mesh,
    alpha = 2
)

A <- inla.spde.make.A(mesh = mesh, loc = coordinates(dat_proj))
pt_idx <- inla.spde.make.index(name = "spatial.field", n.spde = spde1$n.spde)

covar <- as.data.frame(model.matrix(base_f, dat_sub)[,-1])
colnames(covar) <- str_replace(colnames(covar), ' ', '_')
colnames(covar)[14] <- "house_age.7_year_or_more"
colnames(covar)[15] <- "house_age.less_than_1_year"

inla_stack <- inla.stack(
    data  = list(infestation = dat_proj@data$infestation),
    A = list(1, A),
    effects = list(
        covar, 
        c(pt_idx, list(Intercept = 1))
    ),
    tag = "obs.stack"
)

base_f <- infestation ~ bed_hygiene + bird_nests_inside + chicken_coop_location +
    bedroom_clutter + dark_main_room + floor_type + 
    construction_pile_type + firewood_location + grain_storage_in_house + house_age +
    house_hygiene + kitchen_location + land_for_agriculture + num_humans +
    num_cats + num_chickens + num_dogs + num_pigs + sign_animals + sign_rats +
    condition_bedroom_wall + condition_house_wall + material_house_wall +
    material_roof + windows_in_bedroom + dist_perim

m00 <- glm(base_f, data = dat_sub, family = binomial)

m0 <- inla(
    base_f,
    data = dat_sub,
    family = "binomial",
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(dic = TRUE)
)

m1$summary.fixed %>%
    mutate(name = rownames(.)) %>%
    ggplot(aes(name, mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
    geom_point(aes(y = mode), col = 'green', shape = 1) +
    geom_point(data = enframe(m00$coefficients, 'name', 'mean'), col = 'red') +
    coord_flip()

please <- as.formula(
    paste0(
        'infestation ~ -1 + ',
        str_c(colnames(covar), collapse = ' + '),
        ' + Intercept + f(spatial.field, model = spde)'
    )
)

m1 <- inla(
    ## update(base_f, ~ . + 0 + Intercept + f(spatial.field, model = spde)),
    please,
    data = inla.stack.data(inla_stack, spde = spde1),
    family = "binomial",
    control.predictor = list(link = 1, A = inla.stack.A(inla_stack), compute = TRUE),
    control.compute = list(dic = TRUE),
    control.inla = list(control.correct = list(enable = 10))
)

m2 <- inla(
    update(base_f, ~ . + f(id)),
    data = dat_proj@data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
)

risks00 <- predict(m00, type = 'response')
risks0 <- m0$summary.fitted.values$mean
risks1 <- m1$summary.fitted.values$mean

tibble(
    long = coordinates(dat_proj)[,1],
    lat = coordinates(dat_proj)[,2],
    u = as.matrix(A %*% m1$summary.random$s$mean)[,1]    
) %>%
    ggplot(aes(long, lat, col = u)) +
    geom_point()

spde_result1 <- inla.spde2.result(
    inla = m1,
    name = "spatial.field",
    spde = spde1,
    do.transf = TRUE
)

ggplot(as_tibble(spde_result1$marginals.range.nominal[[1]]), aes(x, y)) +
    geom_line()

inla.zmarginal(spde_result1$marginals.range.nominal[[1]])
