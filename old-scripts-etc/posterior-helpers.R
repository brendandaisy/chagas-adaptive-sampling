require(tidyverse)
require(corrr)

phi_cor_x_dist <- function(df_phi, dmat, inf) {
    cor(df_phi) %>%
        as_cordf %>%
        shave %>%
        stretch(na.rm = TRUE) %>%
        cor_info(dmat, inf)
}

pair_info <- function(cordf, dmat, inf) {
    ## append info about the pairwise neighbors onto a cordf in long form
    n1 <- as.numeric(str_extract(cordf$x, '\\d+'))
    n2 <- as.numeric(str_extract(cordf$y, '\\d+'))
    
    cordf %>%
        mutate(
            dist = map2_dbl(n1, n2, ~dmat[.x, .y]),
            diff_inf = map2_chr(n1, n2, ~as.character(inf[.x] + inf[.y]))
        )
}

beta_long <- function(post0, post1) {
    beta0 <- post0 %>%
        select(contains('beta')) %>%
        pivot_longer(everything(), 'coef', values_to = 'm0')

    beta1 <- post1 %>%
        select(contains('beta')) %>%
        pivot_longer(everything(), 'coef', values_to = 'm1')

    beta0 %>%
        mutate(m1 = beta1$m1) %>%
        pivot_longer(-coef)
}

hyper_intervals <- function(post, parms) {
    post %>%
        select(contains(parms
}
