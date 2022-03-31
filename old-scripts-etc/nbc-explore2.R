require(shiny)
require(tidyverse)
require(ggraph)
require(RColorBrewer)
require(gridExtra)
require(geosphere)
require(caret)
require(ggmosaic)

source('dummy-var.R')

### read the data
village <- 'Cerrón'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(vars(village:land_for_agriculture,
                   sign_animals:ownership),
              as_factor)

tg <- get_cov_net(dat_org, k=4)

## plot resulting planar graph
layout <- create_layout(tg, layout='auto')
layout$x <- layout$long
layout$y <- layout$lat

ggraph(layout) +
    geom_edge_link(arrow=arrow(length=unit(1.5, 'mm'), type = 'closed'),
                   end_cap = circle(.5, 'mm')) +
    geom_node_point(aes(col=inf))

## plot dist of distances
tg %>%
    as_tibble %>%
    ggplot(aes(x=weight)) +
    geom_histogram()

cov_matrix <- make_dummy_mat(dat_org, inter=FALSE)

## relationship between source levels and target infestation
tg %E>%
    as_tibble %>%
    select_if(is.factor) %>%
    pivot_longer(bed_hygiene:ownership) %>%
    ggplot() +
    facet_wrap(name ~ ., scales='free') +
    geom_mosaic(aes(x=product(target_inf, value), fill=target_inf)) +
    labs(x='Source attribute',
         y='Target status',
         fill='Target status')

count_cov <- function(...) {
  nb <- which(c(...) == 1)
  if (length(nb) > 1) { 
    tb <- as_tibble_row(colSums(cov_matrix[nb,]))
  }
  else if (length(nb) == 1) { 
    tb <- as_tibble_row(cov_matrix[nb,])
  }
  else {
    r <- rep(0, ncol(cov_matrix))
    names(r) <- colnames(cov_matrix)
    tb <- as_tibble_row(r)
  }
  return(tb)
}

### 1. Create sidebar to control plot input
ui <- fluidPage(
  title = "Neighborhood Explorer",  
  plotOutput('plot_obj'),
  hr(),
  fluidRow(
      column(3,
             sliderInput("h" , h3("Neighborhood radius"),
                         min = min(dm), max = max(dm), step = 1, value = 100)
             ## br(),
             ## checkboxInput('jitter', 'Jitter'),
             ## checkboxInput('smooth', 'Smooth')
             ),
    column(4, offset = 1,
      selectInput('z', 'Neighbor covariate', colnames(cov_matrix))
    ),
    column(4,
      ## selectInput('facet_row', 'Facet Row', c(None='.', names(dataset))),
      ## selectInput('facet_col', 'Facet Column', c(None='.', names(dataset)))
    )
  )
)

### 2. do work for current input
server <- function(input, output) {
    ## my_pal <- colorRampPalette(brewer.pal(7, "Blues")) # pick color palette

    main_dat <- reactive({ # get the data for main panel heatmap
        adj_matrix <- dm
        adj_matrix[dm <= input$h] <- 1
        adj_matrix[dm > input$h] <- 0
        diag(adj_matrix) <- 0
        cov_counts <- as_tibble(adj_matrix, .name_repair='minimal') %>%
            pmap_dfr(count_cov) %>%
            mutate(id = dat_org$id, long = dat_org$long, lat = dat_org$lat,
                   nb = unname(rowSums(adj_matrix)),
                   org = if (str_starts(input$z, 'num')) cov_matrix[,input$z] else as.factor(cov_matrix[,input$z]),
                   inf = dat_org$infestation == 1)
        return(cov_counts)
    })

    output$plot_obj <- renderPlot({ # dynamically create plot for given input values
        pd <- main_dat()
        g1 <- ggplot(pd, aes(x=long, y=lat, shape=inf, size=.data[[input$z]], col=.data[[input$z]])) +
            geom_point(alpha=.8) +
            scale_shape_manual(values=c(1, 2), name='Infestation') +
            guides(size=FALSE)

        g2 <- ggplot(pd, aes(y=.data[[input$z]], x=inf, col=inf)) +
            geom_boxplot() +
            coord_flip() +
            labs(x = 'Infestation', col='Infestation')

        g3 <- ggplot(pd, aes(x=inf, fill=as.factor(org))) +
            geom_bar() +
            coord_flip() +
            labs(x = 'Infestation', fill=input$z, y='Count')


        ## define layout and return plot grid
        lay <- rbind(c(1, 1, 1, 2, 2),
                     c(1, 1, 1, 3, 3))

        grid.arrange(g1, g2, g3, layout_matrix=lay)
    })
}

## run the app
shinyApp(ui = ui, server = server)
