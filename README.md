# chagas-adaptive-sampling

This repository contains the code and data necessary for reproducing the results in the paper [Spatial epidemiology and adaptive targeted sampling to manage the Chagas disease vector Triatoma dimidiata](https://arxiv.org/abs/2111.05964).

## Spatial modeling with INLA

Here we provide a short tutorial on fitting the types of spatial models
used throughout the paper using the R-INLA package. INLA (Integrated
Nested Laplacian Approximation) is an alternative to MCMC for fitting a
wide variety of Bayesian models, and has become a popular option for
geospatial models because of its speed and relative ease of use (see
[Bakka et
al. (2018)](https://wires.onlinelibrary.wiley.com/doi/pdf/10.1002/wics.1443)
for a recent review).

### Load required packages and the data

``` r
library(INLA)
library(tidyverse)
source("code/seq-sampling-helpers.R")

pat <- read_csv("anon-survey-data.csv") |>
    filter(village == "Paternito")
```

### Patterns of infestation

We are going to fit a spatial model to the infestation data in El
Paternito. To begin exploring this data, we can produce a plot of the
binary infestation status of each home. Note that the latitude and
longitude are not actual coordinates here as they have been anonimized,
but we can still treat them as \(x\) and \(y\) coordinates,
respectively.

``` r
ggplot(pat, aes(long, lat, col = as.factor(infestation))) +
    geom_point() +
    labs(x = '', y = '', col = 'Inf. Status')
```

![](examplesunnamed-chunk-3-1.png)<!-- -->

We find that infested houses tend to be in middle and upper parts of the
village, and often can be found close to other infested houses.

### Traditional logistic regression

Before incorporating spatial effects, let’s fit a model using only fixed
effects from the covariate information.  
For simplicity, let’s only use the covariates `bed_hygiene`, `num_dogs`,
and `dist_perim`. To fit such a model in R-INLA, we use the `inla`
function:

``` r
fx_only <- inla(infestation ~ bed_hygiene + num_dogs + dist_perim, data = pat, family = 'binomial')
summary(fx_only)
```

    ## 
    ## Call:
    ##    c("inla.core(formula = formula, family = family, contrasts = contrasts, 
    ##    ", " data = data, quantiles = quantiles, E = E, offset = offset, ", " 
    ##    scale = scale, weights = weights, Ntrials = Ntrials, strata = strata, 
    ##    ", " lp.scale = lp.scale, link.covariates = link.covariates, verbose = 
    ##    verbose, ", " lincomb = lincomb, selection = selection, control.compute 
    ##    = control.compute, ", " control.predictor = control.predictor, 
    ##    control.family = control.family, ", " control.inla = control.inla, 
    ##    control.fixed = control.fixed, ", " control.mode = control.mode, 
    ##    control.expert = control.expert, ", " control.hazard = control.hazard, 
    ##    control.lincomb = control.lincomb, ", " control.update = 
    ##    control.update, control.lp.scale = control.lp.scale, ", " 
    ##    control.pardiso = control.pardiso, only.hyperparam = only.hyperparam, 
    ##    ", " inla.call = inla.call, inla.arg = inla.arg, num.threads = 
    ##    num.threads, ", " blas.num.threads = blas.num.threads, keep = keep, 
    ##    working.directory = working.directory, ", " silent = silent, inla.mode 
    ##    = inla.mode, safe = FALSE, debug = debug, ", " .parent.frame = 
    ##    .parent.frame)") 
    ## Time used:
    ##     Pre = 0.781, Running = 0.0642, Post = 0.0152, Total = 0.86 
    ## Fixed effects:
    ##                   mean    sd 0.025quant 0.5quant 0.975quant mode   kld
    ## (Intercept)     -2.356 0.721     -3.853   -2.329     -1.023   NA 0.004
    ## bed_hygienemala  1.624 0.566      0.586    1.600      2.806   NA 0.003
    ## num_dogs         0.238 0.132     -0.015    0.236      0.502   NA 0.001
    ## dist_perim       0.001 0.002     -0.002    0.001      0.004   NA 0.000
    ## 
    ## Marginal log-Likelihood:  -85.51 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

### Adding random effects in R-INLA

### Spatial smoothing and the SPDE approach

### Dense Gaussian Processes

