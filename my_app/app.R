library(shiny)
library(ASTRA)
library(ggplot2)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Uploading Files"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a file ----
      fileInput("rds", "Input ASTRA result",
                multiple = FALSE),

      # Horizontal line ----
      tags$hr(),

      #Input: Select chromosome number
      textInput("genes", label = h3("Gene(s) to plot (separate by \",\")"), value = "MEG3,BAG5"),

      # Horizontal line ----
      tags$hr(),

      #Input: Select chromosome number
      textInput("title", label = h3("Plot title"), value = "ASTRA result"),

      actionButton("do", "Plot Me"),

      # Horizontal line ----
      tags$hr(),

      downloadButton("downloadData", "Download")
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      plotOutput("plot")
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  data <- reactiveValues()
  observeEvent(input$do, {
    output$plot <- renderPlot({

      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.

      req(input$rds)

      df_input<-readRDS(input$rds$datapath)

      isolate(data$p<-plot_ase(df=df_input,genes_to_plot=as.list(strsplit(input$genes, ","))[[1]],input$title))
      if(!is.null(data$p)){
        plot(data$p)
      }else{
        data$p<-ggplot() +
          annotate("text", x = 4, y = 25, size=8, label = "No SNPs in the selected gene(s)") +
          theme_void()
        plot(data$p)
      }
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".svg", sep="")
      },
      content = function(file) {
        ggsave(file, data$p, width=10, height=8)
      }
    )
  })


}

# Create Shiny app ----
shinyApp(ui, server)
