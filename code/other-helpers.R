# -------------------------------------------------------------------------
# other-helpers.R----------------------------------------------------------
# -------------------------------------------------------------------------
# miscellaneous helper functions used throughout the project---------------
# -------------------------------------------------------------------------

library(sf)
library(geosphere)
library(caret)
library(tidyverse)

# Covariate helpers--------------------------------------------------------

##' Get the predictor data matching pred type (e.g global vs. known in text)
##' Continuous variables will be normalized
##'
##' @title build_pred_dat
build_pred_dat <- function(df, pred_type) { # df NOT dummy coded
    if (pred_type == 'known')
        df_ret <- df
    else if (pred_type == 'global')
        df_ret <- select(df, id:village, dist_perim:truth)
    else if (pred_type == 'latent')
        return(select(df, id:village, infestation:truth)) # don't rescale anything
    return(rescale_cont_vars(df_ret))
}

rescale_cont_vars <- function(df) {
    df_cont <- select(df, contains('num_'), dist_perim, density)
    df_oth <- select(df, !c(contains('num_'), dist_perim, density))

    bind_cols(
        df_oth,
        predict(
            suppressWarnings(preProcess(df_cont, method = c('center', 'scale'))),
            newdata = df_cont
        )
    ) %>%
        ## put before inf for make_formula
        relocate(contains('num_'), dist_perim, density, .before = infestation)
}

# Distance and spatial objects helpers-------------------------------------

dist_mat <- function(df, normalize=FALSE, cutoff=NULL) {
  dm <- df %>%
    select(long, lat) %>%
    distm()
  
  colnames(dm) <- df$id
  rownames(dm) <- df$id
  
  if (normalize) dm <- dm / max(dm)
  if (!is.null(cutoff)) {
    dm[dm <= cutoff] <- 1
    dm[dm > cutoff] <- 0
    diag(dm) <- 0 # no self-loops
  }
  return(dm)
}

sp_project <- function(tbl, normalize = FALSE) {
    dfsp <- as.data.frame(tbl)
    sp::coordinates(dfsp) <- ~long+lat
    proj4string(dfsp) <- CRS("+proj=longlat")
    ret <- spTransform(dfsp, CRS("+proj=utm +zone=16 +datum=WGS84"))
    if (normalize) {
        dm <- dist_mat(tbl)
        wh <- which(dm == max(dm), arr.ind = TRUE)
        p1 <- coordinates(ret)[wh[1, 1],]
        p2 <- coordinates(ret)[wh[2, 1],]
        norm <- sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
        ret <- SpatialPointsDataFrame(coordinates(ret) / norm, ret@data)
    }
    return(ret)
}
