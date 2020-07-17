####################################################################
# Author    : M. Dubbelaar
# Date      : 17-jul-2020
# Purpose   : Shiny app to perform a small analysis.
####################################################################
#                     Install and Load libraries                   #
####################################################################
options(shiny.maxRequestSize=200*1024^2)
options(repos = BiocManager::repositories())
####################################################################
library(shiny)
library(shinyjs)
library(DESeq2)
library(shinycssloaders)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(lobstr)
library(edgeR)
library(DT)

## Use created functions
source("helpers/PCA.R")
source("helpers/geneExpression.R")
source("helpers/volcano.R")
source("helpers/generateRDS.R")
####################################################################
timeStamp <- format(Sys.time(), "%Y%m%d_%H%M%S_")
####################################################################
# Define UI
ui <- fluidPage(theme="style.css",
    ## Use shinyJs to support javaScript in the shiny app
    shinyjs::useShinyjs(),
    ## Application title
    titlePanel(" "),
    
    ## Sidebar that consists of a sidebar panel that where the content is replaced.
    ## This information is dependent on the analysis that is performed
    sidebarLayout(
        sidebarPanel(
            div(id="form", 
                ## The input panel that enables uploading of the (meta)data and conditions of interest
                div(id="inputPanel",       
                    helpText("Please upload a csv file with header and rownames.
                             The rows represent genes and columns different samples."),
                    fileInput(inputId = "datafile", label = "Upload a count matrix", multiple = FALSE, placeholder = "No file selected", accept = "csv"),
                    helpText("Please upload the metadata of the previous file as a csv file.
                             Provide two columns: \"Samples\" with samples names that correspond with the previous file, and \"Description\" with information about the samples.
                             Other information will not be used."),
                    fileInput(inputId = "metadatafile", label = "Upload a metadata", multiple = FALSE, placeholder = "No file selected", accept = "csv"),
                    br(),
                    helpText("This information is updated after uploading the (meta)data."),
                    selectInput(inputId = "con1", label = "Condition 1", choices = ""),
                    selectInput(inputId = "con2", label = "Condition 2", choices = ""),
                    actionButton(inputId = "analysis", label = "Analyse")),
                
                ## After calculation, the analysisoptions are revealed to perform the corresponding analysis
                shinyjs::hidden(radioButtons("analysisOptions", "Analysis", choices = c("Principal component"="pca", "Differentially expression"="DEA", "Quantitative expression"="QEA"), 
                                             selected = NULL,
                                             inline = FALSE, width = NULL, choiceNames = NULL,
                                             choiceValues = NULL)),
                ## Gene input that can be used to search the quantitaive expression of the gene in the data set
                shinyjs::hidden(div(id="geneInput", 
                                    br(),
                                    br(),
                                    selectizeInput(inputId = "genesymbol", label = "Gene symbol", choices = ""))),
                ## FDR and logFC slider for the DEG
                shinyjs::hidden(div(id="FDR slider", 
                                    br(),
                                    br(),
                                    sliderInput("FDR", "False discovery rate", value = 0.05, min = 0.0001, max = 0.1))),
                shinyjs::hidden(div(id="logfc slider", 
                                    br(),
                                    br(),sliderInput("logfc", "Log fold change", value = 3, min = 1, max = 10))),
                
                ## Reset button to go back to the initial state
                shinyjs::hidden(actionButton(inputId = "reset", label = "Reset"))
            )
        ),
           
        div(id="content", 
            # Show a plot that corresponds with the analysis type
            mainPanel(
                div(id="informativeText",
                    h1("Welcome to this shiny app!"),
                    br(),
                    br(),
                    p("This application provides some exploratory features for a bulk RNA-Seq analysis."),
                    p("Fill in your (meta)data and wait until the initial analysis is done."),
                    p("After this, you can choose to perform a PCA, DEA, or QEA on your data."),
                    p("Your data will be removed after the analysis."),
                    br(),
                    p("Note: code for the functions and this app are provided on Github."))
                ),
                
                shinyjs::hidden(
                    div(id="tabPanel",
                        tabsetPanel(id="tabs",
                                    tabPanel("PCA", plotOutput("pca") %>% withSpinner(color="#163459")),
                                    tabPanel("Gene expression", plotOutput("quant") %>% withSpinner(color="#163459")),
                                    tabPanel("Volcano and table",
                                             plotOutput("vol") %>% withSpinner(color="#163459"),
                                             DT::dataTableOutput("table"))
                        )
                    )
                )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    ## reactive part that loads the metadata to determine the available conditions for the DEA
    obtainCon <- reactive({
        inFile <- input$metadatafile
        if (is.null(inFile))
            return(NULL)
        read.csv(inFile$datapath)
    })
    
    ## Give the information, that was obtained by the reactive obtainCon, and fill this into the con1 and con2 ids
    observe({
        updateSelectInput(session, "con1", choices = unique(as.character(obtainCon()$Description)))
        updateSelectInput(session, "con2", choices = unique(as.character(obtainCon()$Description))[order(unique(as.character(obtainCon()$Description)), decreasing = T)])
    })
    
    ## If the analysis is started, continue to perform a variety of tasks
    observeEvent(input$analysis, {
        ## First, check if all of the information is filled in.
        if(!is.null(input$datafile) & !is.null(input$metadatafile) & input$con1 != input$con2){
            ## Disable the input forms
            shinyjs::disable("datafile")
            shinyjs::disable("metadatafile")
            shinyjs::disable("con1")
            shinyjs::disable("con2")
            shinyjs::disable("analysis")
            
            ## The first step is to create new files
            createFiles(input$datafile, input$metadatafile, input$con1, input$con2, timeStamp)
            
            ## Read the data that was created in the previous step
            DEAdat <- readRDS(paste0(timeStamp, "DEA.rds"))
            QEAdat <- readRDS(paste0(timeStamp, "QEA.rds"))
            PCAdat <- readRDS(paste0(timeStamp, "PCA-input.rds"))
            ## Hide the input panel
            shinyjs::hide("inputPanel")
            shinyjs::hide("informativeText")
            ## Show the tabPanel
            shinyjs::show(id="tabPanel")
            ## Show the analysis options (PCA, DEA and QEA)
            shinyjs::show("analysisOptions")
            shinyjs::show(id="reset")
            
            ## Check whether the user selects one of the analyses
            observe({
                ## If PCA is selected then only that tab is opened
                if (input$analysisOptions == "pca") {
                    showTab(inputId = "tabs", target = "PCA")
                } else {
                    hideTab(inputId = "tabs", target = "PCA")
                }
                ## If DEA is selected then only that tab is opened
                ## Furthermore, the slide for the FDR and logfc are shown
                if (input$analysisOptions == "DEA") {
                    showTab(inputId = "tabs", target = "Volcano and table")
                    shinyjs::show(id="FDR slider")
                    shinyjs::show(id="logfc slider")
                } else {
                    hideTab(inputId = "tabs", target = "Volcano and table")
                    shinyjs::hide(id="FDR slider")
                    shinyjs::hide(id="logfc slider")
                }
                ## If DEA is selected then only that tab is opened
                ## The selectizeInput is shown to provide a genesymbol
                if (input$analysisOptions == "QEA") {
                    showTab(inputId = "tabs", target = "Gene expression")
                    shinyjs::show(id="geneInput")
                } else {
                    hideTab(inputId = "tabs", target = "Gene expression")
                    shinyjs::hide(id="geneInput")
                }
                
            })
            
            ## Fill in the genesymbols for the gene search
            observe({
                updateSelectInput(session, "genesymbol", choices = rownames(QEAdat))
            })
            
            ## Create an PCA plot
            output$pca <- renderPlot({
                createCustomPCA(PCAdat, colData(QEAdat), "Description", colData(QEAdat)[,1])
            })
            
            ## Create a gene expression plot
            observe({
                output$quant <- renderCachedPlot({
                    if (input$genesymbol != "") {
                        suppressMessages(createCustomGeneExpressionPlot(QEAdat, input$genesymbol)) 
                    }
                }, cacheKeyExpr = {list(input$genesymbol)})
            })
            
            ## Create a volcano plot and the corresponding data table
            observe({
                output$table <- DT::renderDataTable(createDEGdf(DEAdat, input$logfc, input$FDR, input$con1, input$con2))
                output$vol <- renderCachedPlot({
                    createCustomVolcano(DEAdat, input$logfc, input$FDR, input$con1, input$con2)
                }, cacheKeyExpr = {list(input$FDR, input$logfc)})
            })    
        }
    })
    
    ## Remove the files and stop the app if the browser is closed
    session$onSessionEnded(function() {
        # Delete generated files
        unlink(paste0(timeStamp, "DEA.rds"))
        unlink(paste0(timeStamp, "QEA.rds"))
        unlink(paste0(timeStamp, "PCA-input.rds"))
    })
    
    observeEvent(input$reset, {
        # Delete generated files
        unlink(paste0(timeStamp, "DEA.rds"))
        unlink(paste0(timeStamp, "QEA.rds"))
        unlink(paste0(timeStamp, "PCA-input.rds"))
        reset("form")
        reset("content")
        
        ## Show the input panel
        shinyjs::show("inputPanel")
        shinyjs::show("informativeText")
        ## Hide the tabPanel
        shinyjs::hide(id="tabPanel")
        ## Hide the analysis options (PCA, DEA and QEA)
        shinyjs::hide("analysisOptions")
        shinyjs::hide(id="reset")
        
        shinyjs::enable("datafile")
        shinyjs::enable("metadatafile")
        shinyjs::enable("con1")
        shinyjs::enable("con2")
        shinyjs::enable("analysis")
        })  
}

# Run the application 
shinyApp(ui = ui, server = server)
