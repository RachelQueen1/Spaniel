#' @include addClusterCols.R
#' @include spanielPlot.R
#' @import shiny
#'
NULL

#' RunShinySpaniel
#' 
#' A function to visualise Spatial Transcriptomics. It requires a prepocessed 
#' Seurat Object or a SingleCellExperiment object as well as a rasterised image 
#' saved as an .rds object. There are 4 plots available in the app showing: 
#' a) the number of genes detected per spot, 
#' b) the number of reads detected per spot, 
#' c) clustering results,
#' d) the gene expression of a selected gene."
#' To view the clustering results the columns of the meta.data or colData 
#' containing clustering results must be prefixed with cluster_ . This can be 
#' done by using the markClusterCol() function included in Spaniel. 
#' 
#' @return Runs a Shiny App
#' 
#' @examples
#' ## mark the columns of metadata/colData that contain clustering 
#' ## information see ?markClusterCol for more details#'  
#' sObj <-  readRDS(file.path(system.file(package = "Spaniel"),
#'                     "extdata/sceData.rds"))
#' 
#' img <- readRDS(file.path(system.file(package = "Spaniel"),
#'                     "extdata/image.rds"))
#' 
#' ## run shinySpaniel (upload data.rds and image.rds in the shiny app)
#' ## Not Run:
#' # runShinySpaniel()
#' @export

runShinySpaniel <-function(){
    options(shiny.maxRequestSize=100*1024^2)
    ui <- pageWithSidebar(
        
        # App title ----
        headerPanel("Spatial Transcriptomics"),
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            fileInput("dataFile", "Upload Data File",
                        multiple = FALSE,
                        accept = c(".rds")),
            
            # Title: Upload image file ----  
            fileInput("imageFile", "Upload Image File",
                        multiple = FALSE,
                        accept = c(".rds")),
            
            # Extra options for cluster or gene plots
            uiOutput("plotType"),
            
            # Input: for type of plot ----
            uiOutput("moreControls"), 
            
            p(
                #"Side End"
            )
            
        ),
        
        
        # Main panel for displaying outputs ----
        mainPanel(
            #plotOutput("plotPressed"), 
            
            tabsetPanel(id = "inTabset", 
                        type = "tabs",
                        tabPanel("Getting started",
                                value = "panel1",
                                h3("Plotting Spatial Data"),
                                p("1. Upload the data.rds file and image.rds 
                                    files. It can take a couple of minutes for 
                                    the data to upload"),
                                p("2. Select the type of plot you want to 
                                    look at. There are 4 plots available 
                                    showing:
                                    a) the number of genes detected per spot, 
                                    b) the number of reads detected per spot, 
                                    c) clustering results,
                                    d) the gene expression of a selected 
                                    gene."),
                                p("3. For the cluster plot you must 
                                    also select the cluster resolution you 
                                    wish to plot
                                    (generally a lower resolution equates to 
                                    fewer clusters."),
                                p("4. For the gene plot you must select a gene 
                                    from the drop downlist. There is a bit of a 
                                    delay whilst the gene list is loading.
                                    You can jump to the gene in list by typing 
                                    the first few letters of the gene 
                                    of interest."), 
                                p("5. Click 'Plot' button in the side bar ")
                        ),
                        tabPanel(title = "View Plots",
                                value = "panel2",
                                plotOutput("plotPressed"))
                        
            )
        )    
    )    
    
    
    
    ############ Server #################################
    
    # Define server  ----
    server <- function(input, output, session) {
        
        output$summary <- renderPrint({
            "1. Upload the data.rds file and image.rds files. 
            It can take a couple of minutes for the data to upload"
        })
        
        
        ### S4 object
        Object <- reactive({
            req(input$dataFile)
            readRDS(input$dataFile$datapath)
        })
        
        ### Image object
        imageObj <- reactive({
            req(input$imageFile)
            readRDS(input$imageFile$datapath)
        })
        
        ## Choose plot type
        ## if image and seurat objects uploaded 
        output$plotType <- renderUI({
            req(input$dataFile) 
            req(input$imageFile)
            radioButtons("Type_Of_Plot", "Type of plot:",
                        c("Gene Number Per Spot Plot" = "NoGenes",
                            "Counts Per Spot Plot" = "CountsPerSpot",
                            "Cluster Plot" = "Cluster",
                            "Gene Plot" = "Gene"))
            
        })
        
        #### Cluster list object TO ADD!!
        clusterList <- reactive({
            req(Object())
            metadata = getMetadata(Object())
            colnames(metadata)[grep("cluster_", colnames(metadata))]
        })
        
        ### Extra options for Gene or Cluster plots
        output$moreControls <- renderUI({
            if (req(input$Type_Of_Plot) == "Cluster") {
                list(selectInput("noClusters", "Select clustering resolution:", 
                                clusterList()),
                    actionButton("doPlot", "Plot")
                )
                
            }
            else if (req(input$Type_Of_Plot) == "Gene") {
                s = Object()
                geneList = rownames(s)
                list(selectInput("gene", "Select gene to plot:", 
                                geneList),
                    actionButton("doPlot", "Plot")
                )
            }
            else {
                actionButton("doPlot", "Plot")
            }
            
        })
        
        
        output$plotPressed = renderPlot({
            ## seurat object
            req(input$doPlot)
            s = Object()
            
            ##create coordinates df
            # coordinates = s@meta.data[, c("x", "y")]
            # coordinates$spot = rownames(coordinates)
            metadata = getMetadata(s)
            coordinates = getCoordinates(metadata)
            
            
            ## image grob
            g = imageObj()
            
            ## plot type
            pType = input$Type_Of_Plot
            
            ## set features (NULL for all plots except Gene)
            f = NULL
            if (input$Type_Of_Plot == "Gene"){
                f = input$gene
            }
            
            ## set clusters (NULL for all plots except Cluster)
            cl = NULL
            if (input$Type_Of_Plot == "Cluster"){
                cl = input$noClusters
            }
            
            ### create plot
            spanielPlot(object = s, 
                    grob = g, 
                    plotType =  pType, 
                    gene = f, 
                    clusterRes = cl,
                    customTitle = NULL, 
                    scaleData = TRUE)
        },
        
        
        height = 800, width = 800
        
        
        )
        
        observeEvent(input$doPlot, {
            updateTabsetPanel(session, "inTabset",
                                selected = "panel2"
            )
        })
        
        
    }
    
    shinyApp(ui, server)
    
    
    
    
}
