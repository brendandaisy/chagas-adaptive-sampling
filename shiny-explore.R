###
## shiny app for visualizing a neighbor or neighborhood's covariates, and their
## interaction with infestation
###

require(shiny)
require(tidyverse)
require(RColorBrewer)
require(GGally)
require(ggraph)
require(tidygraph)
require(ggthemes)
require(gridExtra)
require(ggmosaic)
require(igraph)

source('covariate-helpers.R')
source('spatial-helpers.R')

### Read the data
dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
  mutate_at(
    vars(village:land_for_agriculture, sign_animals:infestation),
    as_factor
  )

### Create UI of plot output and controllers
ui <- fluidPage(
  title = "Neighborhood Explorer",  
  plotOutput('plot_obj'), # may have to change height dep on screen
  hr(),
  fluidRow(
    column(
      width=3,
      selectInput("vil" , h3("Village"), unique(dat_org$village)),
      selectInput("method" , h3("Graph type"), c('cutoff', 'knn', 'knn_rev', 'decay')),
      selectInput("pg" , h3("Page"), c(1, 2, 3, 4))
      # selectInput("fun" , h3("Neighbor function"), c('sum', 'mean'))
    ),
    column(
      width=3,
      sliderInput("h" , h3("Cutoff radius"),
                  min = 0, max = 1000, step = 5, value = 100),
      sliderInput("k" , h3("Nearest neighbors"),
                  min = 1, max = 6, step = 1, value = 3)
    ),
    column(
      width=4,
      selectInput(
        'cov',
        h3('Neighbor covariate'),
        colnames(dat_org)[5:ncol(dat_org)]
      )
    )
  )
)

page_one <- function(input, tg) {
  ## plot of graph in space
  layout <- create_layout(tg, layout='manual', x=long, y=lat)

  p1 <- layout %>%
    mutate(infestation = ifelse(infestation == '0', 1, 2)) %>%
    ggraph() +
    geom_edge_link(
      arrow=arrow(length=unit(1.5, 'mm'), type = 'closed'),
      alpha=.25
    ) +
    geom_node_point(
      aes(
        col=.data[[input$cov]],
        ## size=.data[[paste0(input$fun, '_', input$cov)]]
        shape=infestation
      ),
      size=2.5
    ) +
    ## geom_node_point(
    ##     col='grey40',
    ##     size=2,
    ##     data=filter(layout, is.na(.data[[paste0(input$fun, '_', input$cov)]]))
    ## ) +
    scale_shape_identity() +
    labs(x = 'Longitude', y= 'Latitude') +
    theme_minimal()

  if (is.factor(layout[,input$cov])) {
    ## p1 <- p1 + scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326'))

    p2 <- tg %>%
      count_edge_attr(input$cov) %>%
      ggplot() +
      geom_mosaic(aes(x=product(cov), fill=target_inf)) +
      scale_fill_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
      labs(x = paste0('source_', input$cov)) +
      theme_minimal() +
      theme(legend.position='none')
    
    p3 <- tg %>%
      count_neighbor_attr(input$cov) %>%
      pivot_longer(cols=3:ncol(.), names_to=input$cov) %>%
      ggplot(aes(x=.data[[input$cov]], y=value, col=infestation)) +
      geom_boxplot(outlier.shape=NA, alpha=.8) +
      geom_jitter(width=.25, height=0, alpha=.7) +
      scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
      theme_minimal() +
      theme(legend.position='none')
  }
  else {
    p1 <- p1 + scale_color_gradientn(
      colors=c('#f7e540', '#ffb326', '#ff5926', '#bb5de3')
    )

    p2 <- tg %>%
      count_edge_attr(input$cov) %>%
      ggplot(aes(x=cov, y=target_inf, col=target_inf)) +
      geom_boxplot(outlier.shape=NA, alpha=.8) +
      geom_jitter(width=0, height=.25, alpha=.7) +
      scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
      theme_minimal()
    
    p3 <- tg %>%
      count_neighbor_attr(input$cov) %>%
      ggplot(aes(x=infestation, y=nb_prof, col=infestation)) +
      geom_boxplot(outlier.shape=NA, alpha=.8) +
      geom_jitter(width=.25, height=0, alpha=.7) +
      scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
      labs(y=input$cov) +
      theme_minimal() +
      theme(legend.position='none')
  }

  ## define layout and return plot grid
  lay <- rbind(
    c(1, 1, 1, 2, 2),
    c(1, 1, 1, 3, 3)
  )

  grid.arrange(p1, p2, p3, layout_matrix=lay)
}

page_two <- function(input, tg) {
  tg %>%
    count_neighbor_attr(input$cov) %>%
    ggparcoord(
      columns=3:ncol(.),
      groupColumn=2,
      scale='globalminmax',
      alphaLines=.3,
      mapping=aes(shape=infestation)
    ) +
    scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
    theme_minimal()
}

page_three <- function(input, tg, df_vil) {
  dat_oh <- make_dummy_mat(df_vil, contr='contr.ltfr')
  dm <- dist_mat(df_vil)
  am <- as.matrix(get.adjacency(tg, attr='weight'))
  Z <- t(am) %*% dat_oh[,-1]
  colnames(Z) <- colnames(dat_oh)[-1]
  
  for (i in 1:ncol(Z)) Z[,i] <- Z[,i] / max(Z[,i])
  
  Z %<>% as_tibble(.name_repair='minimal') %>%
    mutate(infestation = df_vil$infestation, id = df_vil$id) %>%
    pivot_longer(!c(infestation, id))
  
  p1 <- ggplot(Z, aes(x=name, y=id, fill=value)) +
    geom_tile() +
    facet_wrap(~infestation, scales='free_y') +
    scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  p1
}

page_four <- function(input, tg) {
  meas <- tg %N>%
    mutate(
      in_degree = centrality_degree(mode = 'in'),
      transitivity = local_transitivity()
    ) %>%
    as_tibble
  
  p1 <- ggplot(meas, aes(x=in_degree, fill=infestation)) +
    geom_density() +
    facet_wrap(~infestation)
  
  p2 <- ggplot(meas, aes(x=transitivity, fill=infestation)) +
    geom_density() +
    facet_wrap(~infestation)
  
  p3 <- tg %N>%
    morph(to_components) %>%
    mutate(dist_from_lcen = node_distance_from(node_is_center())) %>%
    unmorph %>%
    as_tibble %>%
    ggplot(aes(x=dist_from_lcen, fill=infestation)) +
    geom_density() +
    facet_wrap(~infestation)
  
  grid.arrange(p1, p2, p3, ncol=1)
}

### Do work for current input and return plots
server <- function(input, output) {
  output$plot_obj <- renderPlot({ # gets called when inputs are changed
    dat_vil <- filter(dat_org, village == !!input$vil)
    tg <- make_cov_net(dat_vil, method=input$method, k=input$k, h=input$h)
    
    if (input$method == 'decay') return(page_three(input, tg, dat_vil))
    
    switch (
      input$pg,
      `1` = page_one(input, tg),
      `2` = page_two(input, tg),
      `3` = page_three(input, tg, dat_vil),
      `4` = page_four(input, tg)
    )
  })
}

### Run the app
shinyApp(ui = ui, server = server)
