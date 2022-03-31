require(sf)
require(elsa)
require(tidyverse)
require(ggthemes)
require(geosphere)

dat_org <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor)

### Are infested houses closer together than random?
rperm <- function(df) { # df filtered for village
    dm <- df %>%
        transmute(lat, long, n_inf = sample(infestation)) %>%
        filter(n_inf == 1) %>%
        select(long, lat) %>%
        distm
    
    mean(dm[upper.tri(dm)])
}

perm_test <- function(df, village, n=10000) {
    sim <- map_dbl(1:n, ~rperm(filter(df, village == !!village)))
    obs <- filter(df, village == !!village) %>%
        filter(infestation == 1) %>%
        select(long, lat) %>%
        distm
    
    obs <- mean(obs[upper.tri(obs)])
    saveRDS(list(sim=sim, obs=obs), paste0(str_to_lower(village), '-inf-perm.rds'))
}

perm_dat <- map_dfr(str_to_lower(unique(dat_org$village)), ~{
    pt <- readRDS(paste0(.x, '-inf-perm.rds'))
    tibble(village = str_to_title(.x),
           sim = pt[[1]],
           obs = pt[[2]],
           p_val = sum(obs >= sim) / length(sim),
           facet = paste0(village, ', p=', p_val))
})

ggplot(perm_dat, aes(x=sim)) +
    geom_histogram() +
    geom_vline(aes(xintercept=obs), 
               data=tibble(facet=unique(perm_dat$facet), obs=unique(perm_dat$obs)),
               col='red') +
    facet_wrap(~facet, scales='free')

### Are infested houses closer to village periphery than random?
rprim <- function(df) { # df filtered for village
    ip <- df %>%
        transmute(dprim, n_inf = sample(infestation)) %>%
        filter(n_inf == 1)
    
    mean(ip$dprim)
}

prim_test <- function(df, village, n=10000) { # df filtered for village
    sim <- map_dbl(1:n, ~rprim(df))
    obs <- df %>%
        filter(infestation == 1) %>%
        pull(dprim) %>%
        mean
    
    saveRDS(list(sim=sim, obs=obs), paste0(str_to_lower(village), '-dprim.rds'))
}

## for each village, create 50m perimeter and calculate distance, then run perm test
prim_dat <- map_dfr(unique(dat_org$village), ~{
    df <- filter(dat_org, village == .x)

    pts <- select(df, long, lat) %>%
        st_as_sf(coords = c('long', 'lat'), crs=4326) %>%
        st_transform("+proj=utm +zone=16 ellps=WGS84")

    buff <- pts %>%
        st_union %>%
        st_convex_hull %>%
        st_buffer(50) %>%
        st_boundary

    df %>%
        mutate(dprim = as.double(st_distance(pts, buff))) %>%
        prim_test(village=.x, n=10000)

    pt <- readRDS(paste0(str_to_lower(.x), '-dprim.rds'))
    
    tibble(village = .x,
           sim = pt[[1]],
           obs = pt[[2]],
           p_val = sum(obs >= sim) / length(sim),
           facet = paste0(village, ', p=', p_val))
})

## plot results compared to observed mean distance
ggplot(prim_dat, aes(x=sim)) +
    geom_histogram() +
    geom_vline(aes(xintercept=obs), 
               data=tibble(facet=unique(prim_dat$facet), obs=unique(prim_dat$obs)),
               col='red') +
    facet_wrap(~facet, scales='free')

### Local autocorrelation of infestation and other covariates

append_elsa <- function(df, village, radius) {
  df <- filter(df, village == !!village)
  dfsp <- as.data.frame(df)
  sp::coordinates(dfsp) <- ~long+lat
  proj4string(dfsp) <- CRS("+proj=longlat")
  res <- spTransform(dfsp, CRS("+proj=utm +zone=16 +datum=WGS84"))
  
  df %>%
    bind_cols(elsa(res, 
                   d=dneigh(res, 0, radius, longlat=FALSE),
                   categorical=TRUE, 
                   zcol='infestation', 
                   drop=TRUE))
}

get_entro <- function(df, village) {
  df <- filter(df, village == !!village)
  dfsp <- as.data.frame(df)
  sp::coordinates(dfsp) <- ~long+lat
  proj4string(dfsp) <- CRS("+proj=longlat")
  res <- spTransform(dfsp, CRS("+proj=utm +zone=16 +datum=WGS84"))
  
  entrogram(res, zcol='infestation', width=25, longlat=FALSE)
}

ent_res <- unique(dat_org$village) %>%
    map_dfr(~{
        ent <- get_entro(dat_org, .x)@entrogram
        tibble(village=.x, !!!ent) # same window length as ELSA
    })
    
ggplot(ent_res, aes(x = distance, y = E)) +
    geom_point(col='blue') +
    facet_wrap(~village, scales='free')

## plot ELSA statistic and endogram for infestation in a village
append_elsa(dat_org, 'Amatillo', 400) %>%
    ggplot(aes(x=long, y=lat, col=ELSA)) +
    geom_point(size=3)

## for comparison, plot infestation in Amatillo with same format
filter(dat_org, village == 'Amatillo') %>%
    ggplot(aes(x=long, y=lat, col=as.factor(infestation))) +
    geom_point(size=3)

vg_cloud <- function(df, village) {
  df_dup <- filter(df, village == !!village) %>%
      mutate(long2 = long, lat2 = lat, inf2 = infestation,
             t1 = 1:n(), t2 = 1:n())

  df_dup %>%
      expand(nesting(long, lat, infestation, t1), nesting(long2, lat2, inf2, t2)) %>%
      filter(t1 < t2) %>%
      select(-t1, -t2) %>%
      transmute(
          dist = distCosine(matrix(c(long, lat), ncol=2),
                            matrix(c(long2, lat2), ncol=2)),
          semi_var = .5 * (infestation - inf2) ^ 2
      )
}

my_binr <- function(df_vg, bin_len) { # bin_len in meters
    bins <- seq(0, max(df_vg$dist) %/% 2, by=bin_len) # same max len as ELSA
    df_vg %>%
        filter(dist <= max(df_vg$dist) %/% 2) %>%
        arrange(dist) %>%
        mutate(bin = findInterval(dist, bins)) %>%
        group_by(bin) %>%
        summarize(
            n = n(),
            mean_dist = mean(dist),
            mean_semi_var = mean(semi_var)
        )
}

vg_res <- unique(dat_org$village) %>%
    map_dfr(~{
        vgc <- vg_cloud(dat_org, .x)
        tibble(village=.x, !!!my_binr(vgc, 25)) # same window length as ELSA
    }) %>%
    group_by(village) %>%
    mutate(norm_n = max(mean_semi_var) * n / max(n)) # normalize to fit nicely in plot

ggplot(vg_res, aes(x = mean_dist, y=mean_semi_var)) +
    geom_col(aes(y = norm_n), col='pink') +
    geom_point() +
    facet_wrap(~village, scales='free')

