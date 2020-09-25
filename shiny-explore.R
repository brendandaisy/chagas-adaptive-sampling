###
## shiny app for visualizing a neighbor or neighborhood's covariates, and their
## interaction with infestation
###

require(shiny)
require(tidyverse)
require(RColorBrewer)
require(ggraph)
require(tidygraph)
require(ggthemes)
require(gridExtra)
require(ggmosaic)

source('covariate-helpers.R')
source('spatial-helpers.R')

### read the data
village <- 'Cerrón'

dat_org <- read_csv('../data-raw/gtm-tp-mf.csv') %>%
    filter(village == !!village) %>%
    mutate_at(
        vars(village:land_for_agriculture, sign_animals:ownership),
        as_factor
    )

dat_oh <- make_dummy_mat(dat_org, inter=TRUE, as_tibble=TRUE)

### Create UI of plot output and controllers
ui <- fluidPage(
    title = "Neighborhood Explorer",  
    plotOutput('plot_obj', height='650px'), # may have to change height dep on screen
    hr(),
    fluidRow(
        column(
            width=3,
            sliderInput("h" , h3("Cutoff radius"),
                        min = min(dm), max = max(dm) %/% 3, step = 1, value = 100),
            sliderInput("k" , h3("Nearest neighbors"),
                        min = 1, max = 6, step = 1, value = 1)
        ),
        column(
            width=3,
            selectInput("method" , h3("Graph type"), c('knn', 'knn_rev', 'cutoff')),
            selectInput("fun" , h3("Neighbor function"), c('sum', 'mean'))
        ),
        column(
            width=4,
            selectInput(
                'cov',
                h3('Neighbor covariate'),
                colnames(dat_oh)[6:ncol(dat_oh)]
            )
        )
    )
)

### Do work for current input and return plots
server <- function(input, output) {
    output$plot_obj <- renderPlot({ # gets called when inputs are changed
        tg <- make_cov_net(dat_oh, method=input$method, k=input$k, h=input$h)

        ## apply fun to each neighborhood and add to graph
        tga <- add_neighbor_attr(
            tg,
            sym(input$cov),
            paste0(input$fun, '_', input$cov),
            get(input$fun)
        )

        ## plot of graph in space
        layout <- create_layout(tga, layout='manual', x=long, y=lat)
        if (!any(layout[,input$cov] > 1)) {
            layout[,input$cov] <- as.factor(layout[,input$cov])
        }

        p1 <- ggraph(layout) +
            geom_edge_link(
                arrow=arrow(length=unit(1.5, 'mm'), type = 'closed'),
                alpha=.25
            ) +
            geom_node_point(
                aes(
                    col=.data[[input$cov]],
                    size=.data[[paste0(input$fun, '_', input$cov)]]
                ),
                shape=1,
                data=drop_na(layout)
            ) +
            geom_node_point(
                col='grey40',
                size=2,
                data=filter(layout, is.na(.data[[paste0(input$fun, '_', input$cov)]]))
            ) +
            scale_size(range = c(3, 10)) +
            labs(x = 'Longitude', y= 'Latitude') +
            theme_minimal()

        if (is.factor(layout[,input$cov])) {
            p1 <- p1 + scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326'))
        }
        else {
            p1 <- p1 + scale_color_gradientn(colors=c('#f7e540', '#ffb326', '#ff5926', '#bb5de3'))
        }
            
        ## interaction between infestation and new column
        p2 <- tga %N>%
            as_tibble %>%
            drop_na %>%
            ggplot(
                aes(
                    y=.data[[paste0(input$fun, '_', input$cov)]],
                    x=as.factor(infestation),
                    col=as.factor(infestation)
                )
            ) +
            geom_boxplot(outlier.shape=NA, alpha=.8) +
            geom_jitter(width=.25, height=0, alpha=.7) +
            scale_color_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
            coord_flip() +
            labs(x = 'Infestation') +
            theme_minimal() +
            theme(legend.position='none')

        ## interaction between infestation and original cov, for each edge
        p3 <- tga %E>%
            as_tibble %>%
            transmute(Infestation = target_inf, cov = as.factor(.data[[input$cov]])) %>%
            ggplot() +
            geom_mosaic(aes(x=product(cov), fill=Infestation)) +
            scale_fill_manual(values = c(`1` = '#bb5de3', `0` = '#ffb326')) +
            labs(x = input$cov) +
            theme_minimal() +
            theme(legend.position='none')

        ## define layout and return plot grid
        lay <- rbind(
            c(1, 1, 1, 2, 2),
            c(1, 1, 1, 3, 3)
        )

        grid.arrange(p1, p2, p3, layout_matrix=lay)
    })
}

### Run the app
shinyApp(ui = ui, server = server)
