# -------------------------------------------------------------------------
# simulation-experiment-results.R------------------------------------------
# -------------------------------------------------------------------------
# Summarize and plot results from the simulation study---------------------
# -------------------------------------------------------------------------

require(gridExtra)
require(ggthemes)

source('code/sequential-sampling.R')

pctr <- function(x) paste0(trunc(x * 100), '\\%')

dat_org <- read_csv("anon-survey-data.csv") |>
  mutate(village = str_replace(village, 'ó', "\\\\'o"))

results <- read_csv("simulation-experiment-results.csv")

# Results across villages--------------------------------------------------

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

## Figure 3: performance of each strat on the two goals
res_summ |>
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
  guides(col = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position='bottom')

# Comparison to random sampling--------------------------------------------

rand_size <- results |> 
  filter(alpha == 'Random') |> 
  pull(m)

rand_comp <- results |> 
  mutate(diff = (rep(rand_size, each = 7) - m) / n)

## Table 2: accuracy and size gains compared to random
rand_comp |> 
  group_by(alpha, pred) |> 
  summarize(
    ac5 = sum(act_pct <= 0.05) / n(), # acc. to 5% target
    ac8 = sum(act_pct <= 0.08) / n(), # 8% target
    quantile(diff, 0.05), median(diff), quantile(diff, 0.95) # median (95% CI) for % difference in m compared to random
  ) |> 
  ungroup()

# rand_comp |> 
#   group_by(alpha, pred) |> 
#   mutate(
#     ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n()
#   ) |> 
#   ggplot(aes(ac5, diff, col = alpha)) +
#   geom_boxplot() +
#   facet_wrap(~pred, scales = 'free')

# rand_comp |> 
#   group_by(alpha, pred) |> 
#   mutate(
#     ac5 = sum(act_pct <= 0.05) / n(), ac8 = sum(act_pct <= 0.08) / n()
#   ) |> 
#   ggplot(aes(diff, act_pct, col = alpha)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(se = FALSE, method = 'lm') +
#   facet_wrap(~pred, scales = 'free')
