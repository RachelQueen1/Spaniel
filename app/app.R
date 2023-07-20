#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(shiny)
library(shinydashboard)
library(shinyjs)
source("appFuncs.R")

ui <- dashboardPage(
  dashboardHeader(title = paste0("Spaniel X.0.0") ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Main", tabName = "Main", icon = icon('home')),
      menuItem("Upload Object", tabName = "obj", icon = icon('upload')),
      menuItem("Quality control", tabName = "QC", icon = icon("dashboard")),
      menuItem("Dimensionality reduction", tabName = "Dimensionalityreduction", icon = icon("magnifying-glass-chart")),
      menuItem("Markers", tabName = "markers", icon = icon("table")),
      menuItem("Point Selection", tabName = "pselection", icon = icon("dna")),
      menuItem("Spaniel plot", tabName = "spotplot", icon = icon('dog'))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Main",   fluidRow(
        box(
          title = "Welcome to Spaniel!", status = "primary",
          "A really cool introduction and a biorxiv link.", width = 12
        ))),
      tabItem(tabName = "obj", fluidRow( 
        box(
          title = "Upload SpatialFeatureExperiment object", status = 'primary',
          fileInput('file1', 'Choose rds File',
                                                    accept=c('.rds'))),
        box(
          title = "How to upload", "Instructions go here"
        ))),
      tabItem(tabName = "QC", fluidRow(
        box(
          plotOutput("qc")),
        box(actionButton("QCbtn", "Go!"),
            shinyjs::hidden(p(id = "QCtext", "Processing...")) ) )),
      tabItem(tabName = "Dimensionalityreduction", fluidRow(
        box(plotOutput("dimred")),
        tabBox(id = 'maindimred', tabPanel('Main', actionButton("Dimredbtn", "Go!"),
            shinyjs::hidden(p(id = "Dimredtext", "Processing...")) ), 
            tabPanel('PCA', 'filler'), tabPanel('UMAP', 'filler')) ) ),
      tabItem(tabName = "spotplot",
              tabBox(id = 'imagetabbed', tabPanel('Cluster',
              plotOutput("image")),tabPanel('Gene', plotOutput('imagegene') )),
              box(selectInput("geneVisu", label="SELECT GENE",
                              choices = NULL, selected = "Full", multiple = FALSE),
                  actionButton("imagebtn", "Go!"),
                 shinyjs::hidden(p(id = "Imagetext", "Processing...")) ) ),
      tabItem(tabName = "markers",
              box(plotOutput("marker"), width = 12),
              box(actionButton("markersbtn", "Go!"),
                  shinyjs::hidden(p(id = "markertext", "Processing...")),
                  selectInput("markerVisu", label="SELECT GENE",
                              choices = NULL, selected = "Full", multiple = TRUE)) ),
      
      tabItem(tabName = "pselection", fluidRow(
        column(width = 4,
               plotOutput("plot1", height = 600, width = 600,
                          # Equivalent to: click = clickOpts(id = "plot_click")
                          click = "plot1_click",
                          brush = brushOpts(
                            id = "plot1_brush", resetOnNew = FALSE
                          )
               )
        )
      ),
      fluidRow(
        column(width = 6,
               h4("Selected points"),
               verbatimTextOutput("selectedPoints")
        )#,
        # column(width = 6,
        #        h4("Brushed points"),
        #        verbatimTextOutput("brush_info")
        # )
      )
              )
    )
  ), shinyjs::useShinyjs()
)

server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^3)
  myData <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- readRDS(inFile$datapath)
    return(data)
  })
  
  
  
  dataR <- reactiveVal()
  sfObj <- reactiveVal()
  markersR <- reactiveVal()
  sfpoints <- reactiveVal()
  
  plotReady <- reactiveValues(qc = FALSE, dimred = FALSE, image = FALSE, markers =FALSE, 
                              imagegene = FALSE)
  
  observe({
    shinyjs::disable("QCbtn")
    shinyjs::show("QCtext")
    print('click')
    dataIn <- QC(myData())
    print(dataIn)
    dataIn <- NormFindVarFeats(dataIn)
    print(dataIn)
    print('boop')
    dataR(dataIn)
    plotReady$qc <- TRUE
    shinyjs::enable("Dimredbtn")
  }) %>% bindEvent(.,input$QCbtn)
  
  observe({
    shinyjs::disable("Dimredbtn")
    shinyjs::show("Dimredtext")
    print('click')
    dataIn <- DimRed(dataR())
    print(dataIn)
    dataIn <- SpClusters(dataIn)
    print(dataIn)
    print('boop')
    dataR(dataIn)
    plotReady$dimred <- TRUE
  }) %>% bindEvent(.,input$Dimredbtn)
  
  observe({
    #shinyjs::disable("imagebtn")
    shinyjs::show("Imagetext")
    print('click')
    sfe <- dataR()
    sf <- colGeometries(sfe)[[1]]
    sf$mid <- sf::st_centroid(sf$geometry)
    sf$clusters <- sfe$clusters
    sfObj(sf)
    print('boop')
    plotReady$image <- TRUE
  }) %>% bindEvent(.,input$imagebtn)
  
  observeEvent(input$geneVisu, { 
    gene<<-input$geneVisu
  })
  
  observe({
    #shinyjs::disable("imagegenesbtn")
    #shinyjs::show("Imagegenestext")
    print('click')
    sfe <- dataR()
    sf <- colGeometries(sfe)[[1]]
    sf$mid <- sf::st_centroid(sf$geometry)
    sf$clusters <- sfe$clusters
    sf$gene <- sfe@assays@data$counts[input$geneVisu,]
    sfObj(sf)
    print('boop')
    plotReady$imagegene <- TRUE
  }) %>% bindEvent(.,input$imagebtn)


  observe({
    shinyjs::disable("markersbtn")
    shinyjs::show("markertext")
    print('click')
    sfe <- dataR()
    #toprank <- Markers(sfe, "clusters", "1")
    #markersR(toprank)
    plotReady$markers <- TRUE
  }) %>% bindEvent(.,input$markersbtn)
  
  observe({
    if(plotReady$qc)
      updateSelectInput(session, "markerVisu", 
                        choices = rownames(dataR()))
    updateSelectInput(session, "geneVisu", 
                      choices = rownames(dataR()))
  })
  
  
  output$qc <-renderPlot({
    print(plotReady$qc)
    if (plotReady$qc) {
      #shinyjs::enable("QC")
      shinyjs::hide("QCtext")
      print(plotReady$qc)
      print(dataR())
      plot_grid(plotColData(dataR(), x = 'sample_id', y = 'detected', colour_by = 'sample_id'),
                plotColData(dataR(), y = "total", x = "sample_id", colour_by = "sample_id"))
      plot_grid(plotColData(dataR(), y = "detected", x = "sample_id", colour_by = "sample_id"),
                plotColData(dataR(), y = "total", x = "sample_id", colour_by = "sample_id"),
                plotColData(dataR(), y = "subsets_mt_percent",x = "sample_id", colour_by = "sample_id"),
                plotColData(dataR(), y = "subsets_ribo_percent",x = "sample_id", colour_by = "sample_id"),
                ncol = 2)
    }
  })

  output$dimred <-renderPlot({
    print(plotReady$dimred)
    if (plotReady$dimred) {
      #shinyjs::enable("QC")
      shinyjs::hide("Dimredtext")
      print(plotReady$dimred)
      print(dataR())
      plotUMAP(dataR(), colour_by='clusters')
    }
  })

  output$image <-renderPlot({
    print(plotReady$image)
    if (plotReady$image) {
      #shinyjs::enable("QC")
      shinyjs::hide("Imagetext")
      ggplot(sfObj()) +
        geom_sf(aes(fill = clusters, geometry = geometry))
    }
  })
  
  output$image <-renderPlot({
    print(plotReady$image)
    if (plotReady$image) {
      #shinyjs::enable("QC")
      shinyjs::hide("Imagetext")
      #ggplot(sfObj()) +
      #  background_image(imgRaster(dataR())) + geom_sf(aes(fill = clusters, geometry = geometry))
      plotVisium(sfObj(), fill = clusters, highlight = in_tissue)
      }
  })
  
  output$imagegene <-renderPlot({
    print(plotReady$image)
    if (plotReady$image) {
      #shinyjs::enable("QC")
      #shinyjs::hide("Imagetext")
      # ggplot(sfObj()) +
      #   background_image(imgRaster(dataR())) + geom_sf(aes(fill = gene, geometry = geometry)) + scale_fill_viridis_c()
      #plotVisium(my(), fill = 'gene', highlight = in_tissue)
      plotSFE(sfObj(), 'cluster')
    }
  })
  
  output$marker <-renderPlot({
    print(plotReady$markers)
    if (plotReady$markers) {
      print(markersR())
      #shinyjs::enable("QC")
      print(input$markerVisu)
      marker_list <- as.vector(input$markerVisu)
      print(marker_list)
      shinyjs::hide("markertext")
      MarkerPlots(dataR(), marker_list, "clusters")
    }
  })
  
  output$plot1 <- renderPlot({
    ggplot(sfpoints()) +
      #geom_sf(aes(fill = clusters, geometry = mid)) +
      geom_point(data = selected(), colour = "red")
  })
  
  output$click_info <- renderPrint({
    nearPoints(sfePoints, input$plot1_click, addDist = TRUE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(sfePoints, input$plot1_brush)
  })
  
  selected <- reactive({
    # add clicked
    print(sfObj())
    sfePoints <- sf::st_coordinates(sfObj()$mid)
    sfpoints(sfePoints)
    selected_points <- sfePoints[0, ]
    selected_points <<- rbind(selected_points, nearPoints(sfePoints, input$plot1_click), brushedPoints(sfePoints, input$plot1_brush))
    # remove _all_ duplicates if any (toggle mode) 
    # http://stackoverflow.com/a/13763299/3817004
    selected_points <<- 
      selected_points[!(duplicated(selected_points) | 
                          duplicated(selected_points, fromLast = TRUE)), ]
    str(selected_points)
    return(selected_points)
  })
  
  output$selectedPoints <- renderPrint({
    selected()
  })
}

# Run app ----
shinyApp(ui, server)
