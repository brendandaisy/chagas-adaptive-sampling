require(nimble)
require(tidyverse)
require(magrittr)
require(ggthemes)
require(geosphere)
require(sf)
require(rgdal)
require(gridExtra)
require(corrr)

source('spatial-helpers.R')
source('covariate-helpers.R')
source('posterior-helpers.R')

gp_compile_nimble <- function(code, data, const, init, mon = NULL) {
    init$C <- expcov(as.matrix(dist(cbind(const$xpos, const$ypos))), init$rho, 1)
    init$z <- c(t(chol(init$C)) %*% rnorm(const$N))

    model <- nimbleModel(code, constants = const, data = data, inits = init)
    cModel <- compileNimble(model)
    conf <- configureMCMC(model)
    ## optional monitors
    if (!is.null(mon)) {
        for (m in 1:seq_along(mon))
            conf$addMonitors(mon[m])
    }
    conf$addMonitors('r0')
    conf$addMonitors('z0')
    conf$addMonitors('z')
    conf$removeSamplers('z')
    ## reduce the initial proposal covariance scale for better mixing
    conf$addSampler('z', 'RW_block', control = list(scale = 0.1))
    mcmc <- buildMCMC(conf)
    return(compileNimble(mcmc, project = cModel))
}

gp_compile_nimble <- function(code, data, const, init, mon = NULL) {
    ## setup initial spatially-correlated latent process values
    ## init$C <- matrix(0, length(data$y), length(data$y))
    ## init$s <- rep(0, length(data$y))
    ## for (v in 1:5) {
    ##     vidx <- (const$ind[v]+1):const$ind[v+1]
    ##     init$C[vidx, vidx] <- expcov(data$dist[vidx, vidx], init$rho[v])
    ##     init$s[vidx] <- c(t(chol(init$C[vidx, vidx])) %*% rnorm(length(vidx)))
    ## }
    init$C1 <- expcov(const$D1, init$rho[1], init$alpha[1])
    init$s1 <- c(t(chol(init$C1)) %*% rnorm(const$N1))
    init$C2 <- expcov(const$D2, init$rho[2], init$alpha[2])
    init$s2 <- c(t(chol(init$C2)) %*% rnorm(const$N2))
    init$C3 <- expcov(const$D3, init$rho[3], init$alpha[3])
    init$s3 <- c(t(chol(init$C3)) %*% rnorm(const$N3))
    init$C4 <- expcov(const$D4, init$rho[4], init$alpha[4])
    init$s4 <- c(t(chol(init$C4)) %*% rnorm(const$N4))
    init$C5 <- expcov(const$D5, init$rho[5], init$alpha[5])
    init$s5 <- c(t(chol(init$C5)) %*% rnorm(const$N5))

    model <- nimbleModel(code, constants = const, data = data, inits = init)
    cModel <- compileNimble(model)
    conf <- configureMCMC(model)
    ## optional monitors
    if (!is.null(mon)) {
        for (m in 1:seq_along(mon))
            conf$addMonitors(mon[m])
    }
    conf$addMonitors('s1')
    conf$removeSamplers('1')
    ## reduce the initial proposal covariance scale for better mixing
    conf$addSampler('s1', 'RW_block', control = list(scale = 0.1))
    conf$addMonitors('s2')
    conf$removeSamplers('s2')
    conf$addSampler('s2', 'RW_block', control = list(scale = 0.1))
    conf$addMonitors('s3')
    conf$removeSamplers('s3')
    conf$addSampler('s3', 'RW_block', control = list(scale = 0.1))
    conf$addMonitors('s4')
    conf$removeSamplers('s4')
    conf$addSampler('s4', 'RW_block', control = list(scale = 0.1))
    conf$addMonitors('s5')
    conf$removeSamplers('s5')
    conf$addSampler('s5', 'RW_block', control = list(scale = 0.1))
    mcmc <- buildMCMC(conf)
    return(compileNimble(mcmc, project = cModel))
}

village_dat <- function(df, cm, v) {
    dfs <- filter(df, village == v)
    dm <- dist_mat(dfs, TRUE)
    pred <- cm %>%
        select(-c(1:3), -infestation) %>%
        filter(village == v) %>%
        select(-village) %>%
        as.matrix
    list(N = nrow(dfs), dist = dm, pred = pred, inf = dfs$infestation)
}

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    mutate_at(vars(village:land_for_agriculture, sign_animals:ownership), as_factor) %>%
    arrange(village)

dat_org$long <- jitter(dat_org$long)
cm <- make_dummy_mat(dat_org, contr = 'contr.treatment', as_tibble = TRUE)

vinfo <- map(unique(dat_org$village), ~village_dat(dat_org, cm, .x))

expcov <- nimbleFunction(     
   run = function(dist = double(2), rho = double(0), alpha = double(0)) {
      returnType(double(2))
      n <- dim(dist)[1]
      result <- matrix(nrow = n, ncol = n, init = FALSE)
      ## sigma2 <- sigma*sigma
      for(i in 1:n) {
          for(j in 1:n)
              result[i, j] <- exp(-(dist[i,j]/rho)^alpha)
      }
      return(result)
   })
cexpcov <- compileNimble(expcov)

code <- nimbleCode({
    for (v in 1:5) {
        rho[v] ~ dunif(0, 10)
        alpha[v] ~ dunif(0, 2)
    }
    C1[1:N1, 1:N1] <- expcov(D1[1:N1, 1:N1], rho[1], alpha[1])
    s1[1:N1] ~ dmnorm(zeros[1:N1], cov = C1[1:N1, 1:N1])
    for (i in 1:N1) {
        logit(r1[i]) <- inprod(beta[1:p], x1[i, 1:p]) + s1[i]
        y1[i] ~ dbern(r1[i])
    }
    C2[1:N2, 1:N2] <- expcov(D2[1:N2, 1:N2], rho[2], alpha[2])
    s2[1:N2] ~ dmnorm(zeros[1:N2], cov = C2[1:N2, 1:N2])
    for (i in 1:N2) {
        logit(r2[i]) <- inprod(beta[1:p], x2[i, 1:p]) + s2[i]
        y2[i] ~ dbern(r2[i])
    }
    C3[1:N3, 1:N3] <- expcov(D3[1:N3, 1:N3], rho[3], alpha[3])
    s3[1:N3] ~ dmnorm(zeros[1:N3], cov = C3[1:N3, 1:N3])
    for (i in 1:N3) {
        logit(r3[i]) <- inprod(beta[1:p], x3[i, 1:p]) + s3[i]
        y3[i] ~ dbern(r3[i])
    }
    C4[1:N4, 1:N4] <- expcov(D4[1:N4, 1:N4], rho[4], alpha[4])
    s4[1:N4] ~ dmnorm(zeros[1:N4], cov = C4[1:N4, 1:N4])
    for (i in 1:N4) {
        logit(r4[i]) <- inprod(beta[1:p], x4[i, 1:p]) + s4[i]
        y4[i] ~ dbern(r4[i])
    }
    C5[1:N5, 1:N5] <- expcov(D5[1:N5, 1:N5], rho[5], alpha[5])
    s5[1:N5] ~ dmnorm(zeros[1:N5], cov = C5[1:N5, 1:N5])
    for (i in 1:N5) {
        logit(r5[i]) <- inprod(beta[1:p], x5[i, 1:p]) + s5[i]
        y5[i] ~ dbern(r5[i])
    }
    for (k in 1:p) {
        beta[k] ~ dt(0, 0.16, 1)
    }
})

## code <- nimbleCode({
##     for (v in 1:5) {
##         rho[v] ~ dunif(0, 10)
##         C[v, 1:ind[v], 1:ind[v]] <- expcov(dist[v, 1:ind[v], 1:ind[v]], rho[v])
##         s[v, 1:ind[v]] ~ dmnorm(zeros[1:ind[v]], cov = C[v, 1:ind[v], 1:ind[v]])
##         for (i in 1:ind[v]) {
##             logit(r[v, i]) <- inprod(beta[1:p], x[v, i, 1:p]) + s[v, i]
##             y[v, i] ~ dbern(r[v, i])
##         }
##     }
##     for (k in 1:p) {
##         beta[k] ~ dt(0, 0.16, 1)
##     }
## })


## code <- nimbleCode({
##     for (k in 1:p) {
##         beta[k] ~ dt(0, 0.16, 1)
##     }
##     for (v in 1:5) {
##         rho[v] ~ dunif(0, 10)
##         vidx <- (ind[v]+1):ind[v+1]
##         C[vidx, vidx] <- expcov(dist[vidx, vidx], rho[v])
##         s[vidx] ~ dmnorm(zeros[vidx], cov = C[vidx, vidx])
##         ## likelihood
##         for(i in vidx) {
##             lp[i] <- inprod(beta[1:p], x[i, 1:p])
##             logit(r[i]) <- lp[i] + s[i]
##             y[i] ~ dbern(r[i])
##             ## y_pp[i] ~ dbern(r[i])
##         }
##     }
## })


const <- list(
    N1 = vinfo[[1]]$N,
    D1 = vinfo[[1]]$dist,
    x1 = vinfo[[1]]$pred,
    N2 = vinfo[[2]]$N,
    D2 = vinfo[[2]]$dist,
    x2 = vinfo[[2]]$pred,
    N3 = vinfo[[3]]$N,
    D3 = vinfo[[3]]$dist,
    x3 = vinfo[[3]]$pred,
    N4 = vinfo[[4]]$N,
    D4 = vinfo[[4]]$dist,
    x4 = vinfo[[4]]$pred,
    N5 = vinfo[[5]]$N,
    D5 = vinfo[[5]]$dist,
    x5 = vinfo[[5]]$pred,
    p = ncol(vinfo[[1]]$pred),
    zeros = rep(0, 300)
)

data <- list(
    y1 = vinfo[[1]]$inf,
    y2 = vinfo[[2]]$inf,
    y3 = vinfo[[3]]$inf,
    y4 = vinfo[[4]]$inf,
    y5 = vinfo[[5]]$inf
)
init <- list(beta = rep(0, const$p), rho = rep(.2, 5), alpha = rep(1, 5))

gp_mc <- gp_compile_nimble(code, data, const, init)
samples <- runMCMC(gp_mc, niter = 10000, nburnin = 4000, nchains = 1)

cv <- runCrossValidate(conf, 6, MCMCcontrol = list(niter = 8000, nburnin = 2000))
cv0 <- runCrossValidate(conf, 6, MCMCcontrol = list(niter = 8000, nburnin = 2000))

summary(samples)
hist(samples[,'rho[5]'], breaks = 25)

pairs(samples[,c(1, 2, 41)])

## combine all chains into a tibble
post <- if (is.list(samples))
            map_dfr(1:length(samples), ~as_tibble(samples[[.x]]))
        else
            as_tibble(samples)

##
med_s <- post %>%
    select(contains('s3[')) %>%
    summarize_all(median) %>%
    pivot_longer(everything()) %>%
    bind_cols(select(filter(dat_org, village == 'Guayabo'), long, lat, infestation))

med_lp <- post %>%
    select(contains('lp[')) %>%
    summarize_all(median) %>%
    pivot_longer(everything()) %>%
    bind_cols(select(dat_sub, long, lat, infestation))

p1 <- ggplot(med_lp, aes(long, lat, size = value, col = as.factor(infestation))) +
    geom_point(alpha = .8) +
    ## scale_color_viridis_c() +
    theme(legend.position = 'bottom')

p2 <- ggplot(med_s, aes(long, lat, size = value, col = as.factor(infestation))) +
    geom_point(alpha = .8) +
    ## scale_color_viridis_c(option = 'inferno') +
    theme(legend.position = 'bottom')

grid.arrange(p1, p2, nrow = 1)

##
phi_cor <- post %>%
    select(contains('s[')) %>%
    phi_cor(dm, dat_sub$infestation)

phi_cor %>%
    filter(dist < .25) %>%
    ggplot(aes(x = dist, y = r, col = diff_inf)) +
    geom_point(size = .9) +
    labs(y = 'Corr(i, j)', x = 'Distance (meters)', col = 'Num. Infested') +
    theme_few()

##
prop <- post %>%
    select(contains('y_pp')) %>%
    map_dbl(~sum(.x) / nrow(post))

dat_sub %>%
    mutate(prop = prop) %>%
    ggplot(aes(long, lat, size = prop, col = as.factor(infestation))) +
    geom_point()

post %>%
    select(rho) %>%
    mutate(t = 1:n()) %>%
    ggplot(aes(t, rho)) +
    geom_line()

## examine betas
post %>%
    select(contains('beta')) %>%
    pivot_longer(everything()) %>%
    ggplot(aes(value)) +
    stat_function(fun = dcauchy, args = list(location = 0, scale = 2.5), color = "red") +
    geom_density(col = 'grey70', size = 1.5) +
    facet_wrap(~name, nrow = 5, scales = 'free')
