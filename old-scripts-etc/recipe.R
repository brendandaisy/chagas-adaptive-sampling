require(tidyverse)
require(tidygraph)
require(ggraph)
require(magrittr)

add_step <- function(recipe, ingredients, instruction, wrap=20, wrap2=20) {
  bind_rows(
    recipe, tibble(
      from=str_wrap(ingredients, width=wrap), 
      to=str_wrap(instruction, width=wrap2), 
      weight=''
      )
    )
}

ing <- c(
  'Silken tofu',
  'Garlic',
  'Nutritional yeast',
  'Tahini/natural peanut butter',
  'Olive oil',
  'Garlic powder',
  'Dill',
  'Smoked paprika',
  'Water',
  'Salt',
  'Pepper',
  'Corn starch',
  'Fettuccini',
  'Daiya shreds',
  'Vegan butter'
)

w <- c(
  '1 pack (e.g. 16oz)',
  '5 cloves',
  '1/4 cup',
  '1/4 cup',
  '2 tbsp',
  '1 tbsp',
  '1 tsp',
  '1/4 tsp (optional)',
  '0-1/4 cup',
  '',
  '',
  '1 tsp (optional)',
  '1 lb',
  'small handfull (optional)',
  '1-2 tbsp (optional)'
)

edge_list <- tibble(from='Fettuccini Alfredo', to=str_wrap(ing, 20), weight=w) %>%
  add_step(ing[1:8], 'Combine in blender until smooth') %>%
  add_step('Fettuccini', 'Cook al dente, stir frequently') %>%
  add_step('Cook al dente, stir frequently', 'Strain') %>%
  add_step(c(ing[9:12], 'Combine in blender until smooth'), 'Simmer on low (~5 min), correct thickness and season') %>%
  add_step(
    c(ing[14:15], 'Simmer on low (~5 min), correct thickness and season', 'Strain'), 
    'Add pasta to sauce to desired ratio, add cheese and butter if using and stir until cheese just melty',
    wrap2=18
  )

as_tbl_graph(edge_list) %>%
  ggraph() +
  geom_edge_link(aes(label=weight), angle_calc='along') +
  geom_node_label(aes(label = name)) +
  coord_flip() +
  scale_y_reverse()
