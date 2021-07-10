require(gridExtra)
require(ggthemes)

source('adaptive-sampling-fn.R')

read_boot <- function(tags, village = c('ama', 'cer', 'pat', 'pre', 'gua')) {
    tv <- expand.grid(
        village = village,
        tag = tags
    )
    
    pmap_dfr(tv, ~{
        readRDS(paste0('bootstrapss/', .y, '-', .x, '.rds')) %>%
            mutate(method = .y, village = .x)
    })
}

mwu_test <- function(res_df, target, ref = 'runif-global') {
    mwu <- function(s, v, m) {
        dat <- filter(res_df, size == s, village == v)
        w <- wilcox.test(
            pull(filter(dat, method == m), !!target),
            pull(filter(dat, method == ref), !!target),
            alternative = if (target == 'var') 'less' else 'greater',
            exact = TRUE
        )
        tibble_row(size = s, village = v, method = m, target = target, pval = w$p.value)
    }

    res_df$method <- res_df$method %>%
        as.factor %>%
        fct_relevel(ref)
    
    met <- levels(res_df$method)

    expand.grid(
        s = unique(res_df$size),
        v = unique(res_df$village),
        method = met[met != ref],
        stringsAsFactors = FALSE
    ) %>%
        pmap_dfr(~mwu(..1, ..2, ..3))
}

results <- read_boot(
    str_c(
        c('runif', 'convexpt75', 'var'),
        ## c('runif', 'var'),
        'global',
        sep = '-'
    )
) %>%
    mutate(nauc = auc - prop)

mwu <- bind_rows(
    mwu_test(results, 'nauc'),
    mwu_test(results, 'disc'),
    mwu_test(results, 'var')
)

mwu %>%
    filter(target == 'var') %>%
    ggplot(aes(interaction(size, village), log(1 + pval), fill = method)) +
    geom_col(position = position_dodge()) +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    geom_hline(yintercept = 0.1, linetype = 'dashed') +
    facet_wrap(~method, nrow = 3)

results %>%
    ggplot(
        aes(
            x = as.factor(size),
            y = sens,
            fill = method
            ## group = interaction(method, Metric, as.factor(size))
        )
    ) +
    geom_boxplot() +
    facet_wrap(~village) +
    ## scale_fill_manual(
    ##     values = c(random_known = 'white', entropy_known = 'grey80', kld_unknown = 'grey60', kld_rand_only = 'grey40')
    ## ) +
    labs(x = 'Pct. Observed', y = 'Score', shape = 'Method')

res_summ <- results %>%
    mutate(
        ## model = ifelse(str_detect(method, 'gp|interp'), 'Interp.', 'Lat. only'),
        method = str_extract(method, '\\w+')
    ) %>%
    group_by(method, size, village) %>%
    summarize(x = mean(var), y = mean(disc))

results %>%
    mutate(
        method = str_extract(method, '\\w+')
    ) %>%
    ggplot(aes(var, disc, col = method, shape = as.factor(size))) +
    geom_point(alpha = .3) +
    geom_point(aes(x, y), data = res_summ, size = 3) +
    ## stat_ellipse(level = 0.9, size = .8) +
    facet_wrap(~village, scales = 'free') +
    ## scale_linetype_manual(values = c('Lat. only' = 2, 'Interp.' = 1)) +
    theme_bw()
    ## xlim(0, 1) +
    ## ylim(0, 1)

unobs <- filter(results, method == 'rconvex1-global', size == .75, nauc < .1) %>%
    unobs_list('Guayabo')

count(tibble(id = unlist(unobs)), id, sort = TRUE) %>%
    mutate(inf = map_dbl(id, ~pull(filter(dat_org, id == .x), 'infestation'))) %>%
    filter(inf == 1)
