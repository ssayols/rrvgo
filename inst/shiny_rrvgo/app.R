library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(magrittr)
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
            tabPanel("simMatrixPlot",
              fluidRow(
                column(width=10, div(style='overflow-x: scroll', plotlyOutput("simMatrixPlot"))),
                box(width=2, title="Options", status="warning",
                  checkboxInput("simMatrixDisplayDendro", "draw dendrogram", value=FALSE),
                  sliderInput("simMatrixFontSize", "font size", ticks=FALSE, min=6, max=12, value=9)
                )
              )
            ),
            tabPanel("scatterPlot"  , div(style='overflow-x: scroll', plotlyOutput("scatterPlot"))),
            tabPanel("treemapPlot"  , div(style='overflow-x: scroll', plotOutput("treemapPlot"))),
            tabPanel("wordcloudPlot",
              fluidRow(
                column(width=10, div(style='overflow-x: scroll', plotOutput("wordcloudPlot"))),
                box(width=2, title="Options", status="warning",
                  sliderInput("wordcloudMinFreq", "min frequency", ticks=FALSE, min=1, max=5, value=2)
                )
              )
            )
          ),
          tabBox(width=NULL,
            tabPanel("reducedTerms",
              div(style='overflow-x: scroll', DTOutput("reducedTerms")),
              downloadLink("downloadReducedTerms", "Download reduced terms")
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
            selectInput("stringency", label="Stringency", selected=0.7,
                        choices=c(`Large (allowed similarity=0.9)`=0.9,
                                  `Medium (0.7)`=0.7,
                                  `Small (0.5)`=0.5,
                                  `Tiny (0.4)`=0.4)),
            selectInput("method", label="Distance measure", selected="Rel",
                        choices=c("Resnik", "Lin", "Rel", "Jiang", "Wang")),
            a("Click to open a webpage with detailed info about these measures",
              href="https://www.bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#semantic-similarity-measurement-based-on-go",
              target="_blank")
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
      input$reduce
      isolate({
        tryCatch({
          read.table(stringsAsFactors=FALSE, strip.white=TRUE, fill=TRUE, text=gsub("[\\t| ]+", "\t", input$goterms))
        }, error=function(e) NULL)
      })
    })
    
    simMatrix <- reactive({
      req(!is.null(goterms()))
      withProgress(message="Calculating similarity matrix, this may take a while...", value=0, {
        tryCatch(calculateSimMatrix(goterms()[, 1], org=input$organism, ont=input$ontology, method=input$method),
                 error=function(e) NULL)
      })
    })

    reducedTerms <- reactive({
      req(!is.null(goterms()))
      withProgress(message="Reducing GO terms...", value=0, {
        scores <- if (ncol(goterms()) >  1)  setNames(goterms()[, 2], goterms()[, 1]) else NULL
        tryCatch(reduceSimMatrix(simMatrix(), scores=scores, threshold=input$stringency, orgdb=input$organism),
                 error=function(e) NULL)
      })
    })
      
    observeEvent(input$example, {
      x <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
      x$qvalue <- -log10(x$qvalue)
      x <- paste(apply(x[, c("ID", "qvalue")], 1, paste, collapse="\t"), collapse="\n")
      updateTextInput  (session, "goterms"   , value=x)
      updateSelectInput(session, "organism"  , selected="org.Hs.eg.db")
      updateSelectInput(session, "ontology"  , selected="BP")
      updateSelectInput(session, "stringency", selected=0.7)
      updateSelectInput(session, "method"    , selected="Rel")
    })
    
    #
    # UI -----
    #
    output$simMatrixPlot <- renderPlotly({
      req(simMatrix(), reducedTerms(), cancelOutput=TRUE)
      ann <- reducedTerms()$term[match(reducedTerms()$parent, reducedTerms()$go)]
      ann <- data.frame(ann[match(rownames(simMatrix()), reducedTerms()$go)])
      colnames(ann) <- ""
      heatmaply::heatmaply(simMatrix(), row_side_colors=ann, plot_method="plotly",
                           symm=TRUE, labRow=NULL, key.title="Similarity",
                           showticklabels=c(FALSE, TRUE),
                           width=1024, height=1024,
                           row_side_palette=rrvgo:::gg_color_hue,
                           show_dendrogram=rep(input$simMatrixDisplayDendro, 2),
                           fontsize_col=input$simMatrixFontSize) %>%
        colorbar(xanchor="left", yanchor="bottom", len=.2, tickfont=list(size=input$simMatrixFontSize), which=1) %>%
        colorbar(xanchor="left", yanchor="bottom", len=.5, tickfont=list(size=input$simMatrixFontSize), which=2)
    })
    
    output$scatterPlot <- renderPlotly({
      req(simMatrix(), reducedTerms(), cancelOutput=TRUE)
      
      x <- cmdscale(as.matrix(as.dist(1-simMatrix())), eig=TRUE, k=2)
      df <- cbind(as.data.frame(x$points),
                  reducedTerms()[match(rownames(x$points), reducedTerms()$go), c("term", "parent", "parentTerm", "size")])
      
      ggplotly(
        scatterPlot(simMatrix(), reducedTerms()) +
          geom_text(aes(label=parentTerm), data=subset(df, parent == rownames(df)), size=3)
      )
    })
    
    output$treemapPlot <- renderPlot({
      req(reducedTerms(), cancelOutput=TRUE)
      treemapPlot(reducedTerms(), palette=gg_color_hue(length(unique(reducedTerms()$parent))))
    })
    
    output$wordcloudPlot <- renderPlot({
      req(reducedTerms(), cancelOutput=TRUE)
      wordcloudPlot(reducedTerms(), min.freq=input$wordcloudMinFreq)
    })
    
    output$reducedTerms <- renderDT({
      req(reducedTerms(), cancelOutput=TRUE)
      datatable(reducedTerms(), rownames=FALSE, selection="none")
    })
    
    output$simMatrix <- renderDT({
      req(simMatrix(), cancelOutput=TRUE)
      datatable(simMatrix(), rownames=FALSE, selection="none")
    })
    
    #
    # Download buttons -----
    #
    output$downloadReducedTerms <- downloadHandler(
      filename="reducedTerms.csv",
      content=function(f) {
        if(!is.null(reducedTerms())) {
          write.csv(reducedTerms(), file=f)
        }
      }
    )
    
    output$downloadSimMatrix <- downloadHandler(
      filename="similarityMatrix.csv",
      content=function(f) {
        if(!is.null(simMatrix())) {
          write.csv(simMatrix(), file=f, row.names=FALSE)
        }
      }
    )
    
  }
)
