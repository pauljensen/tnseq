
library(shiny)

create_shiny_app <- function(tnseq) {
  getunique <- function(name) {
    tnseq$insertions[[name]] %>% 
      unique() %>% 
      as.character()
  }
  
  server <- function(input, output) {
    output$insertionPlot <- renderPlot({
      if (is.null(input$gene) || nchar(input$gene) == 0) {
        gene <- ""
      }
      plot_fitness(tnseq, 
                   strain=input$strain, 
                   condition=input$condition, 
                   gene=input$gene, 
                   lib=input$library, 
                   margin=input$margin)
    })
    
    output$geneSelect <- renderUI({
      inserts <- tnseq$insertions %>% filter(strain==input$strain,
                                             condition==input$condition)
      genes <- na.omit(unique(inserts$gene))
      if (!is.null(input$gene) && input$gene %in% genes) {
        selected <- input$gene
      } else {
        selected <- NULL
      }
      selectInput("gene", "Gene", choices=c("", genes), selected=selected)
    })
  }

  ui <- fluidPage(
    titlePanel("Tn-seq Data Browser"),
    
    fluidRow(
      column(2,
        selectInput("strain", "Strain", 
                    choices=getunique("strain")),
        selectInput("condition", "Condition",
                    choices=getunique("condition")),
        selectInput("library", "Library",
                    choices=c("all", getunique("library"))),
        uiOutput("geneSelect"),
        numericInput("margin", "Margin [bp]", value=1000)
      ),
      column(10, plotOutput("insertionPlot"))
    )
  )
  
  shinyApp(ui=ui, server=server)
}

#create_shiny_app(tnseq) %>% print()

#' @export
show_gui <- function(file) {
  load(file)
  print(create_shiny_app(tnseq))
}

