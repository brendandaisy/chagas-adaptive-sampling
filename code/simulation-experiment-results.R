require(gridExtra)
require(ggthemes)

source('sequential-sampling.R')

read_boot <- function(tags, village = c('ama', 'cer', 'pat', 'pre', 'gua')) {
  vname <- as.character(unique(dat_org$village))
  tv <- expand.grid(
    village = village,
    tag = tags,
    stringsAsFactors = FALSE
  )
  
  pmap_dfr(tv, ~{
    readRDS(paste0('bootstrap-results/bigg-', .y, '-', .x, '.rds')) |>
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
  geom_linerange(aes(xmin = xmin, xmax = xmax), alpha = 0.8) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), alpha = 0.8) +
  geom_point(size = 2.8, shape = 18, alpha = 0.8) +
  geom_text(
    aes(x, y, label = ll),
    data = tibble(
      x = c(0.8, 0.8, 0.85),
      y = c(0.005, 0.008, 0.004),
      ll = c('target met', 'target unmet', 'more efficient'),
      # alpha = factor(rep('$\\alpha=0$', 3), levels = unique(results$alpha))
      village = 'Amatillo',
      pred = 'Global only'
    ),
    inherit.aes = FALSE,
    size = 0.8
  ) +
  facet_grid(fct_relevel(pred, 'Global only') ~ village, scales = 'free') +
  scale_x_continuous(labels = pctr, guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = pctr) +
  scale_color_manual(values = rev(colorblind_pal()(8)[-5])) +
  labs(x = 'Percent of houses observed', y = 'True infestation rate', col = NULL) +
  ## xlim(0, 1) +
  guides(col = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position='bottom')

tikz_plot(gg, 'results-perf', 10, 4.8, cd = 'figs')

### Comparison to random sampling-----------------------------------------------

rand_size <- results |> 
  filter(alpha == 'Random') |> 
  pull(m)

rand_comp <- results |> 
  mutate(diff = (rep(rand_size, each = 7) - m) / n)

rand_comp_q <- rand_comp |> 
  group_by(alpha, pred) |> 
  summarize(
    ac5 = sum(act_pct <= 0.05) / n(), 
    ac8 = sum(act_pct <= 0.08) / n(), 
    quantile(diff, 0.05), median(diff), quantile(diff, 0.95)
  ) |> 
  ungroup()

rand_comp |> 
  group_by(alpha, pred) |> 
  mutate(
    ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n()
  ) |> 
  ggplot(aes(ac5, diff, col = alpha)) +
  geom_boxplot() +
  facet_wrap(~pred, scales = 'free')

rand_comp |> 
  group_by(alpha, pred) |> 
  mutate(
    ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n()
  ) |> 
  ggplot(aes(diff, act_pct, col = alpha)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE, method = 'lm') +
  facet_wrap(~pred, scales = 'free')

ggplot(rand_comp_q, aes(ac5, `mean(diff)`, col = alpha, shape = pred)) +
  geom_point()
