require(magrittr)
require(tidyverse)
require(ggthemes)
require(RColorBrewer)
require(gridExtra)
require(MASS, exclude=c('select'))
require(corrr)
require(tikzDevice)

source('../spatial-helpers.R')

tikz_plot <- function(ggp, fname = 'tikz', w = 8.5, h = 4) {
    tikz(paste0(fname, '.tex'), standAlone=TRUE, width = w, height = h)
    print(ggp)
    dev.off()
    system(paste0('pdflatex ', fname, '.tex'))
}

prec_tcar <- function(alpha, x0, tau = 1) {
    C <- dist_mat(sDT, cutoff = x0)
    tau * (diag(rowSums(C)) - alpha * C)
}

prec_lcar <- function(alpha, k, x0, tau = 1) {
    C <- cmat_lcar(dm_comp, x0, k)
    tau * (diag(rowSums(C)) - alpha * C)
}

prec_ecar <- function(alpha, k, x0, tau = 1) {
    C <- cmat_ecar(dm_comp, x0, k)
    tau * (diag(rowSums(C)) - alpha * C)
}

cor_pairs <- function(prec_mat) {
    cov2cor(solve(prec_mat)) %>%
        as_cordf %>%
        shave %>%
        stretch(na.rm = TRUE)
}

prec2var <- function(prec_mat) {
    diag(solve(prec_mat)) %>%
        enframe('id', 'var')
}

### Setup the data

village <- 'Paternito'

DT <- read_csv('../../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    select(id:long, infestation)

## remove isolated houses from all models for now
iso <- which(rowSums(dist_mat(DT, cutoff = 100)) == 0)
sDT <- DT[-iso,] # the data for this script
dm_comp <- dist_mat(sDT) # dist_mat for LCAR/ECAR

### Fig 1: example decaying weight functions

p1 <- tibble(
    x = c(dm_comp),
    y1 = c(cmat_lcar(dm_comp, 100, .065)),
    y2 = c(cmat_lcar(dm_comp, 100, 1))
) %>%
    ggplot(aes(x, y1), size = 1.1) +
    geom_line(aes(y = y2), col = 'lightpink') +
    ## geom_point(aes(y = y2), col = 'lightpink') +
    geom_line(col = 'tomato2') +
    ## geom_point(col = 'tomato2') +
    xlim(0, 350) +
    xlab(NULL) + ylab(NULL) +
    theme_few()

p2 <- tibble(
    x = c(dm_comp),
    y1 = c(cmat_ecar(dm_comp, 100, .5)),
    y2 = c(cmat_ecar(dm_comp, 100, 2))
) %>%
    ggplot(aes(x, y1), size = 1.1) +
    geom_line(aes(y = y2), col = 'lightpink') +
    ## geom_point(aes(y = y2), col = 'lightpink') +
    geom_line(col = 'tomato2') +
    ## geom_point(col = 'tomato2') +
    xlim(0, 350) +
    xlab(NULL) + ylab(NULL) +
    theme_few()

grid.arrange(p1, p2, nrow = 1, left = '$\\mathit{\\mathbf{C}}$', bottom = 'Distance (m)') %>%
    tikz_plot('fig1')

### Fig 2: correlation matrix by distance

## first obtain correlation in long form for each ind. matrix
tcar_cor <- cor_pairs(prec_tcar(.95, 100)) %>%
    mutate(
        dist = map2_dbl(x, y, ~dm_comp[.x, .y]),
        g = 'Traditional CAR'
    ) %>%
    select(-c(x, y))

lcar_cor1 <- cor_pairs(prec_lcar(.95, .065, 100)) %>%
    mutate(
        dist = map2_dbl(x, y, ~dm_comp[.x, .y]),
        g = 'LCAR: k = 0.065'
    ) %>%
    select(-c(x, y))

lcar_cor2 <- cor_pairs(prec_lcar(.95, .01, 100)) %>%
    mutate(
        dist = map2_dbl(x, y, ~dm_comp[.x, .y]),
        g = 'LCAR: k = 0.01'
    ) %>%
    select(-c(x, y))

ecar_cor1 <- cor_pairs(prec_ecar(.95, .6, 100)) %>%
    mutate(
        dist = map2_dbl(x, y, ~dm_comp[.x, .y]),
        g = 'ECAR: k = 0.6'
    ) %>%
    select(-c(x, y))

ecar_cor2 <- cor_pairs(prec_ecar(.95, 1.1, 100)) %>%
    mutate(
        dist = map2_dbl(x, y, ~dm_comp[.x, .y]),
        g = 'ECAR: k = 1.1'
    ) %>%
    select(-c(x, y))

## now bind them all together for faceting
reduce(list(lcar_cor1, lcar_cor2, ecar_cor1, ecar_cor2), bind_rows, .init = tcar_cor) %>%
    filter(dist < 800) %>%
    ggplot(aes(x = dist, y = r, col = fct_rev(g))) +
    geom_point(alpha = .3, size = .2) +
    facet_wrap(~fct_rev(g), nrow=5) +
    labs(x = 'Distance (m)', y = 'Corr(i, j)', col = 'Model') +
    theme_few() +
    theme(legend.position = 'none')

ggsave('fig2.pdf', width = 3.7, height = 7.5)

### Fig 3: how variance changes with 'spatial strength' alpha

alphas <- c(0, seq(-.99, .99, length.out=20))

tcar_var <- alphas %>%
    map_dfr(~{
        mutate(
            as_tibble(prec2var(prec_tcar(.x, 100))),
            alpha = .x,
            deg = colSums(dist_mat(sDT, cutoff = 100))
        )
    }) %>%
    mutate(g = 'Traditional CAR')

lcar_var <- alphas %>%
    map_dfr(~{
        mutate(
            as_tibble(prec2var(prec_lcar(.x, .065, 100))),
            alpha = .x,
            deg = rowSums(cmat_lcar(dm_comp, 100, .065))
        )
    }) %>%
    mutate(g = 'LCAR, k=0.065')

ecar_var <- alphas %>%
    map_dfr(~{
        mutate(
            as_tibble(prec2var(prec_ecar(.x, 1.1, 100))),
            alpha = .x,
            deg = rowSums(cmat_ecar(dm_comp, 100, 1.1))
        )
    }) %>%
    mutate(g = 'ECAR, k=1.1')

## combine and plot
bind_rows(tcar_var, lcar_var, ecar_var) %>%
    filter(var < 7.5) %>% # o/w scaling looks bad
    ggplot(aes(x = alpha, y = var, col = deg, group = id)) +
    geom_line(alpha=.25, size = .35) +
    scale_color_viridis_c() +
    facet_wrap(~g, nrow = 1) +
    labs(x = 'alpha', y = 'Var(phi)', col = 'Degree')  +
    theme_minimal()

ggsave('fig3.pdf', width = 7.7, height = 5.2)

## how many points did we leave out above?
bind_rows(tcar_var, lcar_var, ecar_var) %>%
    filter(var >= 7.5) %>%
    count(g)

samp2 <- alphas %>%
    map_dfr(~{
        as_tibble(
            mvrnorm(n = 1000, mu = rep(0, nrow(sDT)), Sigma = solve(prec_tcar(.x, 175)))
        ) %>%
            mutate(alpha = .x)
    })

samp3 <- expand.grid(
    alpha = alphas,
    k = c(0, 10^(seq(-4, -1, .5)), .4),
    x0 = seq(0, max(dm_comp) / 3, length.out=20)
) %>%
    pmap_dfr(~{
        as_tibble(
            mvrnorm(
                n = 1000,
                mu = rep(0, nrow(sDT)),
                Sigma = solve(prec_lcar(alpha=..1, k=..2, x0=..3))
            )
        ) %>%
            mutate(alpha=..1, k=..2, x0=..3)
    })

samp4 <- expand.grid(
    alpha = alphas,
    k = c(0, 10^(seq(-4, -1, .5)), .4),
    x0 = seq(0, max(dm_comp) / 3, length.out=20)
) %>%
    pmap_dfr(~{
        as_tibble(
            mvrnorm(
                n = 1000,
                mu = rep(0, nrow(sDT)),
                Sigma = solve(prec_ecar(alpha=..1, k=..2, x0=..3))
            )
        ) %>%
            mutate(alpha=..1, k=..2, x0=..3)
    })

## summarize/contrast distributions across parameters

p1 <- samp1 %>%
    group_by(alpha) %>%
    summarize_at(-ncol(samp1), var) %>%
    pivot_longer(-1) %>%
    mutate(
        deg = rowSums(sdm100[name,]),
        lab = str_remove(name, paste0('_', village))
    ) %>%
    ggplot(aes(x = alpha, y = value, color=deg, group=name)) +
    geom_line(alpha=.25, size=.35) +
    scale_color_viridis_c() +
    labs(x = 'alpha', y = 'Var(phi)', col = 'Degree')  +
    theme_minimal()
    

p2 <- samp2 %>%
    group_by(alpha) %>%
    summarize_at(-ncol(samp2), var) %>%
    pivot_longer(-1) %>%
    mutate(
        deg = rowSums(dm2[name,]),
        lab = str_remove(name, paste0('_', village))
    ) %>%
    ggplot(aes(x = alpha, y = value, color=deg, group=name)) +
    geom_line(alpha=.25, size=.35) +
    scale_color_viridis_c() +
    labs(x = 'alpha', y = 'Var(phi)', col = 'Degree')  +
    theme_minimal()

p3 <- samp3 %>%
    filter(floor(x0) == 242 & k == 10^-2) %>%
    select(-x0, -k) %>%
    group_by(alpha) %>%
    summarize_at(-ncol(.), var) %>%
    pivot_longer(-1) %>%
    mutate(
        deg = rowSums(1 - my_logistic(dm3[name,], x0=242.37, k=10^-2)),
        lab = str_remove(name, paste0('_', village))
    ) %>%
    ggplot(aes(x = alpha, y = value, color=deg, group=name)) +
    geom_line(alpha=.25, size=.35) +
    scale_color_viridis_c() +
    labs(x = 'alpha', y = 'Var(phi)', col = 'Degree')  +
    theme_minimal()

p4 <- samp3 %>%
    filter(floor(x0) == 103 & k == .1) %>%
    select(-x0, -k) %>%
    group_by(alpha) %>%
    summarize_at(-ncol(.), var) %>%
    pivot_longer(-1) %>%
    mutate(
        deg = rowSums(1 - my_logistic(dm3[name,], x0=103.5, k=.1)),
        lab = str_remove(name, paste0('_', village))
    ) %>%
    ggplot(aes(x = alpha, y = value, color=deg, group=name)) +
    geom_line(alpha=.25, size=.35) +
    scale_color_viridis_c() +
    labs(x = 'alpha', y = 'Var(phi)', col = 'Degree')  +
    theme_minimal()

ggsave('priors-var.pdf', grid.arrange(p1, p2, p4, p3, nrow = 2), width=14, height=12)

## TODO: VVV this includes houses that aren't actually neighbors, but I think that is actually better comparison to decaying weight models
## i.e. how much of pattern decaying weights do we see just from higher-order neighbor effects?

cormat1 <- samp1 %>%
    filter(alpha == max(alpha) | alpha == min(alpha)) %>%
    split(.$alpha) %>%
    map(~stretch(correlate(select(.x, -alpha))))

cormat1 %<>%
    map2(c(-.9, .95), ~mutate(.x, alpha=.y)) %>%
    bind_rows %>%
    drop_na() %>%
    mutate(deg = rowSums(sdm100[x,]))

p3 <- cormat1 %>%
    ggplot(aes(x = deg, y = r, col = as.factor(alpha))) +
    geom_jitter(width=.25, height=0) +
    labs(x = 'Degree', y = 'Cor(phi)', col = 'alpha')  +
    theme_minimal()

grid.arrange(p1, p2)

### choose a few houses to inspect in detail for posterior comparison

## find the current 1 neighbor houses
p <- which(rowSums(sub_dist_mat) == 1)
dat_org %>%
    mutate(col = id %in% names(p)) %>%
    ggplot(aes(long, lat, col = col, label = str_remove(id, paste0('_', village)))) +
    geom_text(size=2.5)
