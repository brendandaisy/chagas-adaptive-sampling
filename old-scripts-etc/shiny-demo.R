require(shiny)
require(tidyverse)
require(gridExtra)

### Read the data, do any work that doesn't need to react to controllers
dat_all <- mtcars

### Create UI of plot output and controllers
ui <- fluidPage(
  title = "Tab Title", # title that appears in webpage tab
  plotOutput('plot_obj'), # show the output of the object plot_obj, as returned by server
  hr(),
  fluidRow( # let's make a row of two columns of controllers
    column(
      width=4,
      ## controllers take a variable name, title of controller, then other options (here vector of strings)
      selectInput("col" , h3("Drop down menu"), colnames(dat_all))
    ),
    column(
      width=3,
      sliderInput("temp" , h3("Slider"), min = 0, max = 100, step = 1, value = 10),
      checkboxInput("isfriday", "Checkbox", FALSE)
    )
  )
)

### Do work for current input and return updated output object
server <- function(input, output) {
  output$plot_obj <- renderPlot({ # gets called when inputs are changed
    message <- paste0(
      'Today is ',
      if (input$isfriday) '' else 'NOT ',
      'Friday and it is ',
      input$temp,
      ' \u00B0C outside'
    )
    
    p1 <- ggplot(dat_all, aes(x=.data[[input$col]])) +
      geom_histogram(bins=10) +
      coord_flip() +
      labs(title=message)
    
    grid.arrange(p1, nrow=1) # make complex layouts with layout_matrix
  })
}

### Run the app
shinyApp(ui = ui, server = server)
