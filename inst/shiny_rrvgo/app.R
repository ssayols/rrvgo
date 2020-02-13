library(shiny)
library(shinydashboard)
library(DT)
library(rrvgo)

shinyApp(
  #
  # Define UI -------
  #
  ui=dashboardPage(
    dashboardHeader(title="rrvgo: reduce + visualize GO"),
    dashboardSidebar(disable=TRUE),
    dashboardBody(
      fluidRow(
        column(width=9,
          tabBox(width=NULL,
            tabPanel("simMatrixPlot", div(style='overflow-x: scroll', plotOutput("simMatrixPlot"))),
            tabPanel("scatterPlot"  , div(style='overflow-x: scroll', plotOutput("scatterPlot"))),
            tabPanel("treemapPlot"  , div(style='overflow-x: scroll', plotOutput("treemapPlot"))),
            tabPanel("wordcloudPlot", div(style='overflow-x: scroll', plotOutput("wordcloudPlot")))
          ),
          tabBox(width=NULL,
            tabPanel("reducedTable",
              DTOutput("reducedTable"),
              downloadLink("downloadReducedTable", "Download reduced table")
            ),
            tabPanel("simMatrix",
              DTOutput("simMatrix"),
              downloadLink("downloadSimMatrix", "Download similarity matrix")
            )
          )
        ),
        column(width=3,
          box(width=NULL, title="GO terms", status="warning",
            textAreaInput("goterms", label="Paste here the GO terms", height="200px",
                          value=paste("# Columns are <space> separated.",
                                      "# First column is mandatory and must contain valid GO ids (GO:0009268).",
                                      "# Second column is optional and should contain scores (higher is better.",
                                      "# Other columns are ignored.", sep="\n")
            ),
            hr(),
            actionLink("reduce", "Reduce!")
          ),
          box(width=NULL, title="Options", status="warning",
            selectInput("organism", label="Organism", selected=1,
                        choices=c(Anopheles="org.Ag.eg.db",
                                  Arabidopsis="org.At.tair.db",
                                  Bovine="org.Bt.eg.db",
                                  Worm="org.Ce.eg.db",
                                  Canine="org.Cf.eg.db",
                                  Fly="org.Dm.eg.db",
                                  Zebrafish="org.Dr.eg.db",
                                  `E coli strain K12`="org.EcK12.eg.db",
                                  `E coli strain Sakai`="org.EcSakai.eg.db",
                                  Chicken="org.Gg.eg.db",
                                  Human="org.Hs.eg.db",
                                  Mouse="org.Mm.eg.db",
                                  Rhesus="org.Mmu.eg.db",
                                  Malaria="org.Pf.plasmo.db",
                                  Chimp="org.Pt.eg.db",
                                  Rat="org.Rn.eg.db",
                                  Yeast="org.Sc.sgd.db",
                                  Pig="org.Ss.eg.db",
                                  Xenopus="org.Xl.eg.db")),
            selectInput("ontology", label="Ontology", selected=1,
                        choices=c(`Biologiocal Process`="BP",
                                  `Molecular Function`="MF",
                                  `Cellular Compartment`="CC")),
            selectInput("stringency", label="Stringency", selected=2,
                        choices=c(`Large (allowed similarity=0.9)`="Large",
                                  `Medium (0.7)`="Medium",
                                  `Small (0.5)`="Small",
                                  `Tiny (0.4)`="Tiny")),
            selectInput("method", label="Distance measure", selected=1,
                        choices=c("Resnik", "Lin", "Rel", "Jiang", "Wang")),
            a("Click to open a webpage with detailed info about these measures",
              href="http://bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.htmll#semantic-similarity-measurement-based-on-go")
          )
        )
      )
    )
  ),
  #
  # Define server logic -----
  #
  server=function(input, output, session) {
    
    #
    # helpers -----
    #

    #
    # reactive components -----
    #
    goterms <- eventReactive(input$Reduce, {
      read.table(stringsAsFactors=FALSE, strip.white=TRUE, fill=TRUE,
                 text=gsub("[\\t| ]+", "\t", input$goterms))
    })
    
    simMatrix <- reactive({
      calculateSimMatrix(goterms(), org=input$organism, ont=input$ontology, method=input$method)
    })

    reducedTable <- reactive({
      reduceSimMatrix(simMatrix(), scores=goterms()[, 2])
    })
      
    #
    # UI -----
    #
    output$simMatrixPlot <- renderPlot({
    })
    
    output$scatterPlot <- renderPlot({
    })
    
    output$treemapPlot <- renderPlot({
    })
    
    output$wordcloudPlot <- renderPlot({
    })
    
    output$reducedTable <- DTOutput({
      datatable(reduceTable(simMatrix()))
    })
    
    output$simMatrix <- DTOutput({
      datatable(simMatrix())
    })
    
    #
    # Download buttons -----
    #
    output$downloadReducedTable <- downloadHandler(
      filename="reducedTable.csv",
      content=function(f) {
        write.csv(simMatrix(), file=f, row.names=FALSE)
      }
    )
    
    output$downloadSimMatrix <- downloadHandler(
      filename="similarityMatrix.csv",
      content=function(f) {
        write.csv(reducedTable(), file=f)
      }
    )
    
  }
)
