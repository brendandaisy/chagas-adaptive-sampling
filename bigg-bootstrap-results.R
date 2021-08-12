require(gridExtra)
require(ggthemes)
require(arules)
require(FactoMineR)

source('adaptive-sampling-fn.R')

read_boot <- function(tags, village = c('ama', 'cer', 'pat', 'pre', 'gua')) {
    vname <- as.character(unique(dat_org$village))
    tv <- expand.grid(
        village = village,
        tag = tags,
        stringsAsFactors = FALSE
    )
    
    pmap_dfr(tv, ~{
        readRDS(paste0('bootstrapss/bigg-', .y, '-', .x, '.rds')) %>%
            mutate(
                village = vname[which(str_detect(str_to_lower(vname), .x))],
                pred = if (str_detect(.y, 'known')) 'All' else 'Global only'
            )
    })
}

pctr <- function(x) paste0(trunc(x * 100), '\\%')

results <- read_boot(
    c(str_c('known-v3', c('', '.1', '.2', '.3')),
    str_c('global-v3', c('', '.1', '.2', '.3')))
)
    ## mutate(strat = ifelse(is.na(alpha), strat, str_c(strat, alpha)))

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds') %>%
    mutate(village = str_replace(village, 'ó', "\\\\'o"))

### Bootstrap performance-------------------------------------------------------

res_summ <- results %>%
    group_by(alpha, village, pred) %>%
    summarize(
        x = mean(m / n), y = mean(act_pct),
        xmin = quantile(m / n, 0.05), xmax = quantile(m / n, 0.95),
        ymin = quantile(act_pct, 0.05), ymax = quantile(act_pct, 0.95)
    ) %>%
    ungroup %>%
    mutate(
        alpha = ifelse(is.na(alpha), 'Random', paste0('$\\alpha=', alpha, '$')),
        pred = factor(pred, levels = c('Global only', 'All'))
    )

gg <- res_summ %>%
    ggplot(aes(x, y, col = alpha)) +
    geom_hline(aes(yintercept = 0.05), col = 'tomato', alpha = .7) +
    geom_linerange(aes(xmin = xmin, xmax = xmax), alpha = 0.5) +
    geom_linerange(aes(ymin = ymin, ymax = ymax), alpha = 0.5) +
    geom_point(size = 1.4) +
    facet_grid(pred ~ village, scales = 'free') +
    scale_x_continuous(labels = pctr) +
    scale_y_continuous(labels = pctr) +
    labs(x = 'Pct. observed', y = 'True prevalence', col = NULL) +
    ## xlim(0, 1) +
    theme_bw()

tikz_plot(gg, 'results-perf', 10.5, 4.2)

results %>%
    ggplot(aes(m / n, act_pct, col = as.factor(alpha))) +
    geom_hline(aes(yintercept = thresh), col = 'grey70', linetype = 'dashed', alpha = .7) +
    ## geom_point(size = 1.2, alpha = 0.3) +
    geom_point(aes(x, y), data = res_summ, size = 3) +
    facet_grid(pred ~ village, scales = 'free') +
    ## xlim(0, 1) +
    theme_bw()

### Summarizing bootstraps------------------------------------------------------

seq_sel <- function(des_sel, df_sub)
    imap_dfr(des_sel, ~tibble(idx = .x, iter = .y, inf = df_sub[.x,]$truth))

seq_pred <- function(des_sel, df_sub)
    imap_dfr(des_sel, ~mutate(df_sub[.x,], iter = .y))

seq_fit <- function(des_sel, df_sub, init) {
    obs <- accumulate(des_sel, ~c(.x, .y), .init = init)
    map(obs)
}

seq_mca <- function(des_sel, df_sub) {
    disc_pred <- discretizeDF(select(df_sub, -c(id:village, infestation:truth, num_pigs)))
    mca <- MCA(disc_pred, ncp = 12, graph = FALSE)
    imap_dfr(des_sel, ~mutate(as_tibble(mca$ind$coord[.x, 1:3]), iter = .y))
}

seq_risk <- function(des_sel, df_sub, init) {
    obs <- accumulate(des_sel, ~c(.x, .y), .init = init)
    imap_dfr(obs, ~{
        rem <- df_sub[-.x,]
        print(.x)
        tru <- sum(rem$truth) / nrow(rem)
        tibble(iter = .y, tru = tru)
    })
}

sel_map <- function(village, a, p, fun, ...) {
    dfs <- filter(dat_org, village == !!village)
    if (is.na(a))
        res <- filter(results, village == !!village, is.na(alpha), pred == p)
    else
        res <- filter(results, village == !!village, alpha == a, pred == p)
    
    map_dfr(res$sel, fun, dfs, ...)
}

###

prev <- sel_map('Paternito', 0, 'Global only', seq_risk)

ggplot(prev, aes(iter, tru)) +
    geom_smooth()

### Fixed effects

contr <- map_dfr(
    c(0, 0.7, 2, NA),
    ~mutate(sel_map('Paternito', .x, 'Global only', seq_mca), alpha = .x)
)

contr %>%
    mutate(
        alpha = factor(
            alpha,
            c(0, 0.7, 2, NA),
            c('$\\alpha=0$', '$\\alpha=0.7$', '$\\alpha=2$', 'Random'),
            exclude = NULL
        )
    ) %>%
    group_by(iter, alpha) %>%
    summarize(
        across(
            .fns = list(m = mean, lo = ~quantile(.x, 0.025), hi = ~quantile(.x, 0.975))
        )
    )

%>%
    pivot_longer(-c(iter, alpha))

    ggplot(aes(iter, value, col = name)) +
    geom_smooth() +
    facet_grid(alpha~name) +
    labs(x = 'Iteration', y = 'Principal component') +
    theme_bw() +
    theme(
        ## strip.background = element_blank(),
        ## panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing = unit(0, "lines"),
        legend.position = 'none'
    )

tikz_plot(gg, 'mca-pat', 6, 4.2)
    
###

sel_count <- sel_time_apply('Prensa', .3, 'Global only', des_sel) %>%
    mutate(inf = factor(inf, c(0, 1), labels = c('Not infested', 'Infested'))) %>%
    add_count(idx, iter, sort = TRUE)

inf_count <- sel_time_apply('Prensa', 2, 'Global only', des_sel) %>%
    mutate(inf = factor(inf, c(0, 1), labels = c('Not infested', 'Infested'))) %>%
    count(inf, iter, sort = TRUE)

ggplot(sel_count, aes(iter, as.factor(idx), fill = log(n))) +
    geom_tile(color = 'white') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_viridis_c() +
    facet_grid(inf~., scales = "free", space = "free") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.spacing = unit(0, "lines"),
        axis.text.y = element_blank()
    )

## facet grid this by village and alpha
ggplot(inf_count, aes(iter, n, fill = inf)) +
    geom_col() +
    theme_bw()

tmp %>%
    select(where(is.factor), iter) %>%
    count(bed_hygiene, iter) %>%
    ggplot(aes(iter, n)) +
    geom_line() +
    facet_wrap(~bed_hygiene)

tmpp <- tmp %>%
    select(where(is.factor), iter) %>%
    group_by(iter) %>%
    summarize(across(everything(), ~list(fct_count(.x))))

tmpp$condition_bedroom_wall %>%
    imap_dfr(~mutate(.x, iter = .y)) %>%
    ggplot(aes(iter, n, fill = f)) +
    geom_col()
