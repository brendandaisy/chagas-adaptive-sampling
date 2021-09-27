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
    readRDS(paste0('bootstrapss/bigg-', .y, '-', .x, '.rds')) |>
      mutate(
        village = vname[which(str_detect(str_to_lower(vname), .x))],
        pred = if (str_detect(.y, 'known')) 'All' else 'Global only'
      )
  })
}

pctr <- function(x) paste0(trunc(x * 100), '\\%')

dat_org <- prep_model_data('../data-raw/gtm-tp-mf.rds') |>
  mutate(village = str_replace(village, 'ó', "\\\\'o"))

results <- read_boot(
  c(str_c('known-v3', c('', '.1', '.2', '.3', '.4')),
    str_c('global-v3', c('', '.1', '.2', '.3', '.4')))
) |> 
  mutate(
    alpha = ifelse(is.na(alpha), 'Random', paste0('$\\alpha=', alpha, '$')),
    pred = factor(pred, levels = c('Global only', 'All'))
  )

### Bootstrap performance-------------------------------------------------------

res_summ <- results |>
  group_by(alpha, village, pred) |>
  summarize(
    x = mean(m / n), y = mean(act_pct),
    xmin = quantile(m / n, 0.05), xmax = quantile(m / n, 0.95),
    ymin = quantile(act_pct, 0.05), ymax = quantile(act_pct, 0.95)
  ) |>
  ungroup()

res_out <- results |>
  group_by(alpha, village, pred) |>
  filter(
    act_pct < quantile(act_pct, 0.05) | act_pct > quantile(act_pct, 0.95) |
    (m / n) < quantile(m / n , 0.05) | (m / n)  > quantile(m / n , 0.95)
  ) |> 
  ungroup()

gg <- res_summ |>
  ggplot(aes(x, y, col = alpha)) +
  geom_hline(aes(yintercept = 0.05), col = 'tomato', alpha = .7) +
  geom_point(
    aes(m / n, act_pct), data = res_out, 
    alpha = 0.3, size = .7
  ) +
  geom_linerange(aes(xmin = xmin, xmax = xmax), alpha = 0.5) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), alpha = 0.5) +
  geom_point(size = 1.8, shape = 5) +
  facet_grid(pred ~ village, scales = 'free') +
  scale_x_continuous(labels = pctr, guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = pctr) +
  labs(x = 'Pct. observed', y = 'True prevalence', col = NULL) +
  ## xlim(0, 1) +
  guides(col = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position='bottom')

tikz_plot(gg, 'results-perf', 10, 4.8)

### Comparison to random sampling-----------------------------------------------

rand_size <- results |> 
  filter(alpha == 'Random') |> 
  pull(m)

rand_comp <- results |> 
  mutate(diff = (rep(rand_size, each = 7) - m) / n) |> 
  filter(alpha != 'Random')

rand_comp_q <- rand_comp |> 
  filter(alpha != 'Random') |> 
  group_by(alpha, pred) |> 
  summarize(
    ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n(), min(diff), mean(diff), max(diff)
  ) |> 
  ungroup()

rand_comp |> 
  filter(alpha != 'Random') |> 
  group_by(alpha, pred) |> 
  mutate(
    ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n(), min(diff), mean(diff), max(diff)
  )
  
res_cdf <- function(d, a, p) {
  ds <- filter(rand_comp, diff <= d, alpha == a, pred == p)
  sum(ds$act_pct <= 0.05) / nrow(ds)
}

rand_comp <- rand_comp |> 
  rowwise() |> 
  mutate(succ_rate = res_cdf(diff, alpha, pred))

rc_summ95 <- rand_comp |> 
  filter(succ_rate >= 0.95) |>
  mutate(diff = diff + 0.002) |> 
  group_by(alpha, pred) |> 
  summarise(mx95 = max(diff), mn = min(c(diff, 0)))

rc_summ <- rand_comp |> 
  filter(succ_rate >= 0.8) |>
  mutate(diff = diff + 0.002) |> 
  group_by(alpha, pred) |> 
  summarise(mx80 = max(diff)) |> 
  left_join(rc_summ95)

gg <- ggplot(rand_comp, aes(diff, act_pct)) +
  geom_rect(
    aes(xmax = mx95), 
    xmin = -0.5, ymin = -.1, ymax = 0.2, 
    data = rc_summ, 
    fill = 'green', 
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmax = mx80, xmin = mx95), 
    ymin = -.1, ymax = 0.2, 
    data = rc_summ, 
    fill = 'yellow', 
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = mx80), 
    xmax = 0.5, ymin = -.1, ymax = 0.2, 
    data = rc_summ, 
    fill = 'steelblue3', 
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0.05, linetype = 'longdash', col = 'tomato', alpha = 0.3, size = 1.2) +
  geom_point(alpha = 0.5, size = 1) +
  geom_text(
    aes(label = paste0('$\\Delta_{95}=', round(100 * mx95, 1), '\\%$'), x = mx95 - 0.01), 
    y = 0, size = 2, data = rc_summ
  ) +
  geom_text(
    aes(label = paste0('$\\Delta_{80}=', round(100 * mx80, 1), '\\%$')), 
    y = .08, x = .01, size = 2, data = rc_summ
  ) +
  facet_grid(alpha ~ pred, scales = 'free_x') +
  scale_x_continuous(labels = pctr) +
  scale_y_continuous(labels = pctr, limits = c(-.01, 0.11)) +
  labs(y = 'True infestation rate', x = 'Difference in percent sampled from random') +
  theme_bw()

tikz_plot(gg, 'comp-random', 6.5, 6, cd = 'figs')

### Summarizing bootstraps------------------------------------------------------

accu_end_des <- function(sel_list, init) {
  # map(sel_list, unlist)
  purrr::reduce(
    sel_list[1:(length(sel_list) - 1)], 
    ~c(.x, unlist(.y)), .init = unlist(init)
  )
}

sum(accu_end_des(results$sel[1], results$init[1]) %in% 1:170)

end_des <- map2(results$sel, results$init, accu_end_des)

map_dfr(end_des, ~{
  tibble(sel_idx = .x, )
})

### MCA stuff

seq_mca <- function(des_sel, df_sub, global = FALSE) {
  if (global)
    disc_pred <- discretizeDF(select(df_sub, dist_perim, density))
  else {
    disc_pred <- discretizeDF(
      select(df_sub, -c(id:village, infestation:truth, num_pigs))
    )
  }
  mca <- MCA(disc_pred, ncp = 12, graph = FALSE)
  perc <- mca$eig[1:3, 'percentage of variance']
  ret <- imap_dfr(
    des_sel,
    ~mutate(as_tibble(mca$ind$coord[.x, 1:3]), iter = .y)
  )
  colnames(ret)[1:3] <- str_c(colnames(ret)[1:3], ' (', round(perc, 2), '\\%)')
  return(ret)
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

### Fixed effect MCA analysis

mca_res <- map_dfr(
  c("$\\alpha=0$", "$\\alpha=0.7$", "$\\alpha=2$", 'Random'),
  ~mutate(sel_map('Prensa', .x, 'All', seq_mca, global = FALSE), alpha = .x)
)

mca_res_summ <- mca_res |>
  pivot_longer(-c(iter, alpha), 'dim') |>
  group_by(iter, alpha, dim) |>
  summarize(
    across(
      .fns = list(m = mean, lo = ~quantile(.x, 0.025), hi = ~quantile(.x, 0.975))
    )
  ) |> 
  ungroup()

gg <- ggplot(mca_res_summ, aes(iter, value_m, col = dim)) +
  geom_linerange(aes(ymin = value_lo, ymax = value_hi), alpha = 0.5) +
  geom_point(size = 0.9) +
  facet_grid(alpha~dim) +
  labs(x = 'Iteration', y = 'Value') +
  # ylim(-.9, 2) +
  theme_bw() +
  theme(
    ## strip.background = element_blank(),
    ## panel.border = element_rect(colour = "black", fill = NA),
    # panel.spacing = unit(0, "lines"),
    legend.position = 'none'
  )

tikz_plot(gg, 'mca-pre-all', 6, 4.2)

###

sel_count <- sel_time_apply('Prensa', .3, 'Global only', des_sel) |>
  mutate(inf = factor(inf, c(0, 1), labels = c('Not infested', 'Infested'))) |>
  add_count(idx, iter, sort = TRUE)

inf_count <- sel_time_apply('Prensa', 2, 'Global only', des_sel) |>
  mutate(inf = factor(inf, c(0, 1), labels = c('Not infested', 'Infested'))) |>
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

tmp |>
  select(where(is.factor), iter) |>
  count(bed_hygiene, iter) |>
  ggplot(aes(iter, n)) +
  geom_line() +
  facet_wrap(~bed_hygiene)

tmpp <- tmp |>
  select(where(is.factor), iter) |>
  group_by(iter) |>
  summarize(across(everything(), ~list(fct_count(.x))))

tmpp$condition_bedroom_wall |>
  imap_dfr(~mutate(.x, iter = .y)) |>
  ggplot(aes(iter, n, fill = f)) +
  geom_col()
