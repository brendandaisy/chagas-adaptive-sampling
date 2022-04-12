# -------------------------------------------------------------------------
# seq-sampling-helpers.R---------------------------------------------------
# -------------------------------------------------------------------------

library(INLA)
library(tidyverse)
library(gridExtra)
library(ggthemes)

### Selection Methods-----------------------------------------------------------

##' Select 3 new locations uniformly at random
##'
##' @title rand_unif
##' @param unobs_idx row indicies not yet selected in design
##' @param fit an INLA object
##' @return Indices corresponding to selected rows in data
rand_unif <- function(unobs_idx, fit) sample(unobs_idx, 3)

##' Selection new locations based on combination of risk and predictor variance
##'
##' @title comb_risk_var
##' @param unobs_idx row indicies not yet selected in design
##' @param fit an INLA object
##' @param t Adaptive sampling iteration
##' @param alpha Shape of tradeoff; alpha < 1 prioritizes risk sooner
##' @return Indices corresponding to selected rows in data
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

##' Fit an INLA model with both spatial and iid random effects
##'
##' @title fit_gp_both
##' @param df A dataframe
##' @param spde An spde model such as \code{inla.spde2.pcmatern}
##' @param stack An INLA stack containing projected data
##' @param control A \code{list} of control variables based to \code{inla}
##' @return An INLA model object
fit_gp_both <- function(df, spde, stack, control = list(mlik = FALSE, config = TRUE)) {

  # define hyperparameters for fixed effects
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  # define loggramma prior and hypers for iid effects 
  pid <- list(theta = list(prior = "loggamma", param = c(1, 0.01)))
  
  # produce a formula with fixed and random effects and return inla fit
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

##' Fit an INLA model with spatial random effects but not iid REs
##'
##' @title fit_gp_spatial
##' @param df A dataframe
##' @param spde An spde model such as \code{inla.spde2.pcmatern}
##' @param stack An INLA stack containing projected data
##' @param control A \code{list} of control variables based to \code{inla}
##' @return An INLA model object
fit_gp_spatial <- function(df, spde, stack, control = list(mlik = FALSE)) {

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

##' Fit an INLA model with spatial effects, using the dense \code{dmatern} model
##' This is quite slow for $n>300$ or so
##'
##' @title fit_gp_spatial_dense
##' @param df A dataframe
##' @param control A \code{list} of control variables based to \code{inla}
##' @return An INLA model object
fit_gp_spatial_dense <- function(df, control = list(mlik = FALSE)) {
  pfixed <- list(mean.intercept = 0, prec.intercept = .3, mean = 0, prec = .3)
  hyper <- list(
    range = list(param = c(0.1, 0.05)),
    prec = list(param = c(3, 0.1)),
    nu = list(initial = log(1), fixed = TRUE)
  )
  idx <- 1:nrow(df)
  loc <- coordinates(sp_project(df, TRUE))
  
  inla(
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

##' Fit an INLA model with iid random effects
##'
##' @title fit_gp_random
##' @param df A dataframe
##' @param spde An spde model such as \code{inla.spde2.pcmatern}
##' @param stack An INLA stack containing projected data
##' @param control A \code{list} of control variables based to \code{inla}
##' @return An INLA model object
fit_gp_random <- function(df, mesh, spde, control = list(mlik = FALSE)) {

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



##' Make a formula of infestation ~ Intercept + covariates
##'
##' @title make_formula
##' @param df A dataframe
##' @param exclude Column indexes not to use as covariates
##' @return A formula
make_formula <- function(df, exclude=c(1:4, ncol(df) - 1, ncol(df))) {
  trend <- str_c(colnames(df)[-exclude], collapse = '+')
  as.formula(
    paste0(
      'infestation~-1 + Intercept',
      if (trend == '') '' else '+',
      trend
    ),
    env=parent.frame(1)
  )
}

##' Make a formula of infestation ~ Intercept + covariates
##'
##' @title spde_stack
##' @param df A dataframe
##' @param mesh An `inla.mesh` object
##' @param mesh An `inla.spde` object
##' @return An `inla.stack`
spde_stack <- function(df, mesh, spde, project=FALSE) {
  if (project) {
    df_proj <- sp_project(df, normalize = TRUE)
    loc <- coordinates(df_proj)
  }
  else {
    loc <- as.matrix(select(df, long, lat))
  }
  A <- inla.spde.make.A(mesh = mesh, loc = loc)
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

##' Make an `inla.mesh.2d` object from a nonconvex full around datapoints
##' Data points are pojected using a geocoordinate system and distances normalized
##'
##' @title inla_mesh
inla_mesh <- function(df, project=FALSE) {
  if (project) {
    df_proj <- sp_project(df, normalize = TRUE)
    loc <- coordinates(df_proj)
  }
  else {
    df_proj <- df
    loc <- as.matrix(select(df, long, lat))
  }
  v <- df_proj$village[1]
  bound <- inla.nonconvex.hull(loc)
  
  inla.mesh.2d(
    loc = loc,
    boundary = bound,
    cutoff = 0.07,
    max.edge = if (v == 'Cerr?n') 0.15 else 0.17,
    min.angle = 30,
    offset = 0.09
  )
}

### Plotting--------------------------------------------------------------------

##' Plot the mean and standard deviation of some field at each location
##' `fit_sum` must be something (i.e. from INLA model fit) containing `mean` and `sd` values
##' only works for models fit using `fit_gp_spatial_dense`
##'
##' @title plot_obs_surface
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
    labs(x = '', y = '', col = 'Inf. Status') +
    theme_bw()
  
  gg2 <- dfs %>%
    filter(!obs) %>%
    ggplot(aes(long, lat, col = as.factor(infestation), size = sd)) +
    geom_point(alpha = .7) +
    geom_point(data = filter(dfs, obs), shape = 10) +
    labs(x = '', y = '', col = 'Inf. Status') +
    theme_bw()
  
  grid.arrange(gg1, gg2, nrow = 1)
}