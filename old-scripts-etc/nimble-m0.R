code0 <- nimbleCode({
    for(k in 1:p) {
        beta[k] ~ dt(0, 0.16, 1)
    }
    ## likelihood
    for(i in 1:N) {
        logit(r[i]) <- inprod(beta[1:p], x[i, 1:p])
        y[i] ~ dbern(r[i])
    }
})

const0 <- list(
    N = nrow(dat_sub),
    p = ncol(cm),
    x = cm
)

m0 <- nimbleModel(
    code0,
    constants = const0,
    data = data,
    inits = list(beta = rep(0, ncol(cm)))
)
cm0 <- compileNimble(m0)
conf0 <- configureMCMC(m0)
MCMC0 <- buildMCMC(conf0)
cMCMC0 <- compileNimble(MCMC0, project = cm0)

s0 <- runMCMC(cMCMC0, niter = 10000, nburnin = 4000, nchains = 1)
