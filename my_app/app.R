library(shiny)
library(ASTRA)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Uploading Files"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a file ----
      fileInput("vcf", "Input het SNPs vcf (gzipped)",
                multiple = FALSE),
      fileInput("shapeit", "Input het SNPs vcf (gzipped) phased with SHAPEIT4",
                multiple = FALSE),
      fileInput("haptree", "Input het SNPs vcf (gzipped) phased with HapTree-X",
                multiple = FALSE),

      # Horizontal line ----
      tags$hr(),

      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),

      # Horizontal line ----
      tags$hr(),

      #Input: Select chromosome number
      textInput("chrom", label = h3("Chromosome number"), value = "1"),

      # Horizontal line ----
      tags$hr(),

      #Input: Select chromosome number
      textInput("sampleid", label = h3("Sample name"), value = "SampleID"),
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Data file ----
      tableOutput("contents")

    )

  )
)

# Define server logic to read selected file ----
server <- function(input, output) {

  output$contents <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$vcf)
    req(input$shapeit)
    req(input$haptree)
    df1 <- read_shapeit(input$vcf$datapath, input$shapeit$datapath,
                       sample_name = input$sampleid)
    df2 <- read_haptreex(input$vcf$datapath, input$haptree$datapath,
                       sample_name = input$sampleid)
    BM<-biomart_df()
    BMsel<-BM[BM$chromosome_name%in%input$chrom,]
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- conshap(df1[df1$CHROM==input$chrom,],df2[df2$CHROM==input$chrom,],BMsel,CHROM = input$chrom)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )

    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }

  })

}

# Create Shiny app ----
shinyApp(ui, server)
