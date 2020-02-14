library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(heatmaply)
library(rrvgo)

orgdb <- c(Anopheles="org.Ag.eg.db",
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
           Xenopus="org.Xl.eg.db")

shinyApp(
  #
  # Define UI -------
  #
  ui=dashboardPage(
    dashboardHeader(title="rrvgo: reduce + visualize GO", titleWidth=300),
    dashboardSidebar(disable=TRUE),
    dashboardBody(
      fluidRow(
        column(width=9,
          tabBox(width=NULL,
            tabPanel("simMatrixPlot", div(style='overflow-x: scroll', plotlyOutput("simMatrixPlot"))),
            tabPanel("scatterPlot"  , div(style='overflow-x: scroll', plotlyOutput("scatterPlot"))),
            tabPanel("treemapPlot"  , div(style='overflow-x: scroll', plotOutput("treemapPlot"))),
            tabPanel("wordcloudPlot", div(style='overflow-x: scroll', plotOutput("wordcloudPlot")))
          ),
          tabBox(width=NULL,
            tabPanel("reducedTable",
              div(style='overflow-x: scroll', DTOutput("reducedTable")),
              downloadLink("downloadReducedTable", "Download reduced table")
            ),
            tabPanel("simMatrix",
              div(style='overflow-x: scroll', DTOutput("simMatrix")),
              downloadLink("downloadSimMatrix", "Download similarity matrix")
            )
          )
        ),
        column(width=3,
          box(width=NULL, title="GO terms", status="warning",
            textAreaInput("goterms", label="Paste here the GO terms", height="200px",
                          value=paste("# Columns are <space> or <tab> separated.",
                                      "# First column is mandatory and must contain valid GO ids (GO:0009268).",
                                      "# Second column is optional and should contain scores (higher is better.",
                                      "# Other columns are ignored.",
                                      "# Lines starting with '#' are ignored", sep="\n")
            ),
            hr(),
            fluidRow(column(6, actionLink("reduce", "Reduce!")),
                     column(6, actionLink("example", "example"), align="right"),
            ),
          ),
          box(width=NULL, title="Options", status="warning",
            selectInput("organism", label="Organism", selected="org.Hs.eg.db",
                        choices=orgdb[order(names(orgdb))]),
            selectInput("ontology", label="Ontology", selected="BP",
                        choices=c(`Biologiocal Process`="BP",
                                  `Molecular Function`="MF",
                                  `Cellular Compartment`="CC")),
            selectInput("stringency", label="Stringency", selected="Medium",
                        choices=c(`Large (allowed similarity=0.9)`="Large",
                                  `Medium (0.7)`="Medium",
                                  `Small (0.5)`="Small",
                                  `Tiny (0.4)`="Tiny")),
            selectInput("method", label="Distance measure", selected="Rel",
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
    goterms <- reactive({
      print(input$reduce)
      isolate({
        tryCatch({
          read.table(stringsAsFactors=FALSE, strip.white=TRUE, fill=TRUE, text=gsub("[\\t| ]+", "\t", input$goterms))
        }, error=function(e) NULL)
      })
    })
    
    simMatrix <- reactive({
      #validate(need(!is.null(goterms())))
      req(!is.null(goterms()))
      tryCatch(calculateSimMatrix(goterms()[, 1], org=input$organism, ont=input$ontology, method=input$method),
               error=function(e) NULL)
    })

    reducedTable <- reactive({
      #validate(need(!is.null(goterms())))
      req(!is.null(goterms()))
      scores <- if (ncol(goterms()) >  1)  setNames(goterms()[, 2], goterms()[, 1]) else NULL
      tryCatch(reduceSimMatrix(simMatrix(), scores=scores, threshold=input$stringency, orgdb=input$organism),
               error=function(e) NULL)
    })
      
    observeEvent(input$example, {
      x <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
      x$qvalue <- -log10(x$qvalue)
      x <- paste(apply(x[, c("ID", "qvalue")], 1, paste, collapse="\t"), collapse="\n")
      updateTextInput(session, "goterms", value=x)
    })
    
    #
    # UI -----
    #
    output$simMatrixPlot <- renderPlotly({
      req(simMatrix(), cancelOutput=TRUE)
      heatmaply::heatmaply(simMatrix(), plot_method="plotly")
    })
    
    output$scatterPlot <- renderPlotly({
      req(simMatrix(), reducedTable(), cancelOutput=TRUE)
      ggplotly(scatterPlot(simMatrix(), reducedTable()))
    })
    
    output$treemapPlot <- renderPlot({
      req(reducedTable(), cancelOutput=TRUE)
      treemapPlot(reducedTable())
    })
    
    output$wordcloudPlot <- renderPlot({
      req(reducedTable(), cancelOutput=TRUE)
      wordcloudPlot(reducedTable())
    })
    
    output$reducedTable <- renderDT({
      req(reducedTable(), cancelOutput=TRUE)
      datatable(reducedTable())
    })
    
    output$simMatrix <- renderDT({
      req(simMatrix(), cancelOutput=TRUE)
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
