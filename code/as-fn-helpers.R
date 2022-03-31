require(INLA)
require(gridExtra)
require(ggthemes)
require(tikzDevice)

### Selection Methods-----------------------------------------------------------

rand_unif <- function(unobs_idx, fit) sample(unobs_idx, 3)

##' Selection based on combination of risk and predictor variance
##'
##' @title comb_risk_var
##' @param t Adaptive sampling iteration
##' @param alpha Shape of tradeoff; alpha < 1 prioritizes risk sooner
##' @return Indices corresp. to rows in data
comb_risk_var <- function(unobs_idx, fit, t, alpha) {
  eta <- fit$summary.linear.predictor[unobs_idx,]
  risk <- fit$summary.fitted.values[unobs_idx,]$mean
  eta_var <- eta$sd^2
  risk <- (risk - mean(risk)) / sd(risk)
  eta_var <- (eta_var - mean(eta_var)) / sd(eta_var)
  score <- t^alpha * risk + (1 - t^alpha) * eta_var
  names(score) <- rownames(eta)
  as.double(
    str_extract(names(sort(score, decreasing = TRUE)), '\\d+')[1:3]
  )
}

### INLA Models-----------------------------------------------------------------

fit_gp2 <- function(df, spde, stack, control = list(mlik = FALSE, config = TRUE)) {
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  pid <- list(theta = list(prior = "loggamma", param = c(1, 0.01)))
  
  inla(
    update(
      make_formula(df),
      ~ . + f(z, model = spde) + f(id, model = 'iid', hyper = pid)
    ),
    data = inla.stack.data(stack, spde = spde, pid = pid),
    family = 'binomial',
    control.fixed = pfixed,
    control.predictor = list(A = inla.stack.A(stack), link = 1, compute = TRUE),
    control.compute = control
  )
}


fit_gp <- function(df, spde, stack, control = list(mlik = FALSE)) {
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  
  inla(
    update(
      make_formula(df),
      ~ . + f(z, model = spde)
    ),
    data = inla.stack.data(stack, spde = spde),
    family = 'binomial',
    control.fixed = pfixed,
    control.predictor = list(A = inla.stack.A(stack), link = 1, compute = TRUE),
    control.compute = control
  )
}

fit_old <- function(df, control = list(mlik = FALSE)) {
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  hyper <- list(
    range = list(param = c(0.1, 0.05)),
    prec = list(param = c(3, 0.1)),
    ## range = list(initial = log(0.2), fixed = TRUE),
    # prec = list(initial = log(0.25), fixed = TRUE),
    nu = list(initial = log(1), fixed = TRUE)
  )
  idx <- 1:nrow(df)
  loc <- coordinates(sp_project(df, TRUE))
  
  inla(
    # update(make_formula(df), ~ . + f(
    #   idx,
    #   model = 'dmatern',
    #   locations = loc,
    #   hyper = hyper)
    # ),
    infestation ~ density + dist_perim + f(
      idx,
      model = 'dmatern',
      locations = loc,
      hyper = hyper
    ),
    data = mutate(df, idx = idx),
    family = 'binomial',
    control.fixed = pfixed,
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = control
  )
}

fit_iid <- function(df, mesh, spde, control = list(mlik = FALSE)) {
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  pid <- list(theta = list(prior = "loggamma", param = c(1, 0.01)))
  
  inla(
    update(make_formula(df), ~ . + f(id, model = 'iid', hyper = pid)),
    data = mutate(df, Intercept = 1),
    family = 'binomial',
    control.fixed = pfixed,
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = control
  )
}

### Basic operations------------------------------------------------------------

my_logit <- function(x, inv = FALSE) {
  if (inv)
    return(exp(x) / (1 + exp(x)))
  log(x / (1 - x))
}

make_formula <- function(df) {
  trend <- str_c(colnames(df)[-c(1:4, ncol(df) - 1, ncol(df))], collapse = '+')
  as.formula(
    paste0(
      'infestation~-1 + Intercept',
      if (trend == '') '' else '+',
      trend
    )
  )
}

spde_stack <- function(df, mesh, spde) {
  df_proj <- sp_project(df, normalize = TRUE)
  A <- inla.spde.make.A(mesh = mesh, loc = coordinates(df_proj))
  z_idx <- inla.spde.make.index(name = "z", n.spde = spde$n.spde)
  inla.stack(
    data  = list(infestation = df$infestation, link = 'logit'),
    A = list(A, 1),
    effects = list(
      c(z_idx, list(Intercept = 1)),
      select(df, -infestation, -truth)
    ),
    tag = "zstack"
  )
}



inla_mesh <- function(df) {
  df_proj <- sp_project(df, normalize = TRUE)
  v <- df_proj$village[1]
  bound <- inla.nonconvex.hull(coordinates(df_proj))
  
  inla.mesh.2d(
    loc = coordinates(df_proj),
    boundary = bound,
    cutoff = 0.07,
    max.edge = if (v == 'Cerr?n') 0.15 else 0.17,
    min.angle = 30,
    offset = 0.09
  )
}

unobs_list <- function(boot_df, village) {
  df <- prep_model_data('../data-raw/gtm-tp-mf.rds') %>%
    filter(village == !!village)
  
  obs <- filter(boot_df, village == str_to_lower(str_sub(village, end = 3)))$obs
  map(obs, ~df$id[-.x])
}

bal_sample <- function(y, n_init = 10) {
  pos <- sample(which(y == 1), floor(n_init / 2))
  neg <- sample(which(y == 0), ceiling(n_init / 2))
  return(c(pos, neg))
}

### Plotting--------------------------------------------------------------------

plot_obs_surface <- function(fit_sum, obs_idx, df) {
  dfs <- df %>%
    mutate(
      mean = fit_sum$mean, 
      sd = fit_sum$sd, 
      obs = 1:n() %in% obs_idx, 
      idx = 1:n()
    )
  
  gg1 <- dfs %>%
    filter(!obs) %>%
    ggplot(aes(long, lat, col = as.factor(infestation), size = mean)) +
    geom_point(alpha = .7) +
    geom_point(data = filter(dfs, obs), shape = 10) +
    ## geom_text(aes(label = idx), col = 'grey30', size = 2.8, nudge_y=.00025)
    labs(x = '', y = '', col = 'Inf. Status') +
    theme_bw()
  
  gg2 <- dfs %>%
    filter(!obs) %>%
    ggplot(aes(long, lat, col = as.factor(infestation), size = sd)) +
    geom_point(alpha = .7) +
    geom_point(data = filter(dfs, obs), shape = 10) +
    ## geom_text(aes(label = idx), col = 'grey30', size = 2.8, nudge_y=.00025)
    labs(x = '', y = '', col = 'Inf. Status') +
    theme_bw()
  
  grid.arrange(gg1, gg2, nrow = 1)
}

plot_stat_timeseries <- function(fits, key, stat) {
  df <- map_dfr(fits, ~{
    if (vec_depth(.x[[key]]) > 2)
      svec <- flatten(.x[[key]])[[stat]]
    else
      svec <- .x[[key]][[stat]]
    names(svec) <- as.character(1:length(svec))
    tibble_row(!!!svec)
  })
  
  gg <- df %>%
    mutate(t = 1:n()) %>%
    pivot_longer(-t) %>%
    ggplot(aes(t, log(value), col = name, group = name)) +
    geom_line(alpha = .8) +
    ## scale_color_viridis_c() +
    theme(legend.position = 'none')
  
  print(gg)
}

plot_x_timeseries <- function(fits) {
  df <- map_dfr(fits, ~{
    if (vec_depth(.x[[key]]) > 2)
      svec <- flatten(.x[[key]])[[stat]]
    else
      svec <- .x[[key]][[stat]]
    names(svec) <- as.character(1:length(svec))
    tibble_row(!!!svec)
  })
  
  gg <- df %>%
    mutate(t = 1:n()) %>%
    pivot_longer(-t) %>%
    ggplot(aes(t, log(value), col = name, group = name)) +
    geom_line(alpha = .8) +
    ## scale_color_viridis_c() +
    theme(legend.position = 'none')
  
  print(gg)
}

plot_sel_eta <- function(is, iter = c(1, 5, 10, 20)) {
  par(mfrow = c(ceiling(length(iter) / 2), 2))
  for (i in iter) {
    ms <- is$fit[[i]]$marginals.linear.predictor[is$sel[[i]]]
    plot(ms[[1]])
    points(ms[[3]], col = 'green')
    ## points(is$fit[[i]]$marginals.linear.predictor[[10]], col = 'red')
    ## points(ms[[1]][,1], dnorm(ms[[1]][,1], 0, .3^(-1) + 1), col = 'green')
  }
}

### Covariate interpolation BS--------------------------------------------------

update_pred <- function(df_all, varlist) {
  df_up <- imodify(df_all, ~{
    if (!(.y %in% names(varlist)))
      return(.x)
    if (any(is.na(varlist[[.y]])))
      return(NA)
    if (is.factor(.x))
      return(fct_other(.x, drop = varlist[[.y]]))
  })
  
  select(df_up, where(~!all(is.na(.x))))
}

elim_low_count <- function(df_all, obs_idx) {
  df_obs <- df_all %>%
    slice(obs_idx) %>%
    select(-c(id:village, dist_perim:truth))
  
  compact(map(df_obs, ~{
    ret <- list()
    if (is.factor(.x)) {
      ret <- fct_count(.x) %>%
        filter(n <= 5) %>%
        pull(f)
    }
    ret
  }))
}

elim_single_lvl <- function(df_all, obs_idx) {
  df_obs <- df_all %>%
    slice(obs_idx) %>%
    select(-c(id:village, dist_perim:truth))
  
  compact(
    map(
      df_obs,
      ~if (any(fct_count(as.factor(.x))$n > nrow(df_obs) - 5)) NA else list()
    )
  )
}

elim_nonspatial <- function(df_all, obs_idx) {
  elist <- function(var) {
    ret <- c()
    if (is.factor(var)) {
      jc <- joincount.test(var, lw, zero.policy = TRUE, adjust.n = TRUE)
      for (i in seq_along(jc)) {
        if (jc[[i]]$p.value >= 0.05)
          ret <- c(ret, levels(var)[i])
      }
    }
    else { # is numeric
      mc <- moran.test(var, lw, zero.policy = FALSE, adjust.n = TRUE)
      if (mc$p.value >= 0.05)
        ret <- NA
    }
    ret
  }
  ## setup data for Moran's and kriging:
  sp_norm <- sp_project(df_all, normalize = TRUE)    
  sp_obs <- sp_norm[obs_idx,]
  sp_obs@data <- select(sp_obs@data, -c(id, village, dist_perim:truth))
  ## assign obs weights s.t. w = 0.05 at 1/3 vil diameter:
  nb <- dnearneigh(sp_obs, 0, 1 / 3) # cutoff @ 1/3 obs diameter
  lw <- nb2listwdist(
    nb, sp_obs,
    type = 'exp',
    alpha = -log(0.05) * 3,
    zero.policy = TRUE
  )
  
  return(compact(map(sp_obs@data, elist)))
}

interp_preds <- function(df_full, obs_idx) {
  ## setup spatial data for kriging:
  df_dum <- dummify_vars(df_full)
  sp_norm <- sp_project(df_dum, normalize = TRUE)
  sp_obs <- sp_norm[obs_idx,]
  sp_unobs <- sp_norm[-obs_idx,]
  
  interp <- tibble()
  ## TODO: cant treat a spatial reference var as Other w/o changing it in full data set
  ## somehow also
  if (length(vars) > 0) {
    interp <- map_dfc(
      colnames(sp_norm@data), ~{
        x <- pmin(1, pmax(0, krige_var(sp_obs, sp_unobs, .x, 1 / 3)$var1.pred))
        if (!is.null(x))
          tibble(!!.x := x)
        else
          NULL # ok to ret NULL since will be treated as Other
      }
    )
  }
  
  ## remove single lvl after 2nd removal and return
  return(rm_single_lvl_vars(interp))
}

krige_var <- function(sp, new_sp, var, cut) {
  vg <- variogram(
    as.formula(paste0(var, '~1')),
    sp,
    cutoff = cut,
    width = cut / min(nrow(sp@data) / 2, 15),
    debug.level = 0
  )
  vgf <- fit.variogram(vg, vgm('Mat', range = cut / 2), fit.kappa = TRUE, debug.level = 0)
  if (attr(vgf, 'singular'))
    return(NULL)
  print(plot(vg, vgf, main = var))
  x <- krige(
    as.formula(paste0(var, '~1')),
    sp,
    model = vgf,
    newdata = new_sp,
    debug.level = 0
  )
  print(spplot(x, z = 'var1.pred', main = var))
  return(x@data)
}

rm_single_lvl_vars <- function(df_dum) { # id cols must be removed
  vars <- suppressMessages(
    str_match(colnames(select(df_dum, contains('.'))), '(.+)\\.(.+)') %>% # only factors
      as_tibble(.name_repair = 'unique') %>%
      add_count(...2) %>%
      filter(n >= 2) %>%
      pull(...1)
  )
  select(df_dum, all_of(vars) | !contains('.')) # add cont back on
}

adj_krige_probs <- function(disc_interp) {
  adj_obs <- function(...) {
    ob <- list(...) # a row of interp dummy vars
    ## add an 'Other' group for any var that was not kriged
    ll <- map(groups, ~list(Other = 0))
    names(ll) <- groups
    ## for each dummy var, find its group and add to sublist:
    for(n in names(ob)) {
      f <- str_extract(n, '.+(?=\\.)')
      ll[[f]] <- append(ll[[f]], rlang::list2(!!n := ob[[n]]))
    }
    ## for each sublist, pick a level
    map_if(ll, ~!any(is.na(.x)), ~{
      lvl_probs <- unlist(.x) # now a named vector
      if (sum(lvl_probs) <= 0)
        return('Other')
      if (sum(lvl_probs) < 1) # if some left over prob
        lvl_probs$Other <- 1 - sum(lvl_probs)
      sample(names(lvl_probs), 1, prob = lvl_probs)
    })
  }
  groups <- unique(str_extract(colnames(disc_interp), '.+(?=\\.)'))
  
  disc_interp %>%
    pmap_dfr(function(...) tibble_row(!!!adj_obs(...)))
}


