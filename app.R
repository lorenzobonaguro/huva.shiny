
# setwd("~/Documents/Uni Bonn/Lab Rotations/Schultze/package/huva 0.1.4")

library(shiny)
library(shinydashboard) # for the appearance
library(shinyBS)
library(huva)
library(huva.db)
library(useful)
library(shinybusy)    # for the loader
library(shinyjs)      # for the reset button
library(V8)           # also for the reset button 
library(plotly)       # for interactive plots
library(DT)           # for the DataTables

# Gene names
# create a vector containing all gene names
ImmVar_CD4T <- rownames(huva.db$ImmVar$data$CD4T)
ImmVar_CD14M <- rownames(huva.db$ImmVar$data$CD14M)
classifier_PBMC1 <- rownames(huva.db$classifier$data$PBMC.1)
classifier_PBMC2 <- rownames(huva.db$classifier$data$PBMC.2)
classifier_PBMC3 <- rownames(huva.db$classifier$data$PBMC.3)
FG500 <- rownames(huva.db$FG500$data$PBMC)
CEDAR_CD4T <- rownames(huva.db$CEDAR$data$CD4T)
CEDAR_CD8T <- rownames(huva.db$CEDAR$data$CD8T)
CEDAR_CD14M <- rownames(huva.db$CEDAR$data$CD14M)
CEDAR_CD19B <- rownames(huva.db$CEDAR$data$CD19B)
CEDAR_CD15G <- rownames(huva.db$CEDAR$data$CD15G)
CEDAR_PLA <- rownames(huva.db$CEDAR$data$PLA)

all_genes <- c(ImmVar_CD14M, ImmVar_CD4T, classifier_PBMC1, classifier_PBMC2, classifier_PBMC3, FG500, 
               CEDAR_CD4T, CEDAR_CD8T, CEDAR_CD14M, CEDAR_CD19B, CEDAR_CD15G, CEDAR_PLA)
# remove duplicates
genes_unsorted <- unique(all_genes)
# sort in alphabetical order
genes <- genes_unsorted[order(genes_unsorted)]


# Phenotype names
# only include FG500 study for phenotype
phenotype_cellcount <- sort(colnames(huva.db$FG500$metadata$cellcount))
phenotype_cytokines <- sort(colnames(huva.db$FG500$metadata$cytokines))


# Signature names 
signatures_unsorted <- names(hallmarks_V7.2)
signatures <- sort(signatures_unsorted)

# Define the js method that resets the page
jsResetCode <- "shinyjs.reset = function() {history.go(0)}" 

# HuVa logo
logo <- img(src = "huva_logo_V2.png", height = 50, width = 170)


ui <- dashboardPage(
  # skin = "green"
  dashboardHeader(title = logo, titleWidth = 180),
  
  dashboardSidebar(width = 180,
                   sidebarMenu(
                     menuItem("Start", tabName = "tab_start", icon = icon("chevron-right")),
                     menuItem("Expression", tabName = "tab_expr", icon = icon("chevron-right")),
                     menuItem("Annotation", tabName = "tab_anno", icon = icon("chevron-right")),
                     menuItem("DE genes", tabName = "tab_de", icon = icon("chevron-right")),
                     menuItem("GSEA", tabName = "tab_gsea", icon = icon("chevron-right")),
                     menuItem("Metadata", tabName = "tab_meta", icon = icon("chevron-right"))
                   )
  ),
  
  dashboardBody(
    # add a loader 
      add_busy_spinner(spin = "fading-circle", position = "bottom-right"), # ADD THIS TO THE FINAL CODE!!! 
    tabItems(
      # 1st tab
      tabItem(tabName = "tab_start",
              fluidRow(
                column(12, helpText("The", strong("HuVa app"), "uses the variances in gene expression to compare 
                                    2 groups – “low” and “high” – containing data from healthy individuals of 
                                    population-scale multiomics studies binned by 1 of the 3 following conditions:")),
                column(10, offset = 1, 
                       helpText("1)", strong("Gene of interest:"), "find individuals with “low” or “high” expression 
                                of this gene"),
                       helpText("2)", strong("Phenotype of interest:"), "find individuals with “low” or “high” values 
                                of that phenotype, e.g. # of monocytes in blood")),
                column(8, offset = 2, 
                       helpText("- find the 2 groups with ”low” and ”high” # of monocytes"),
                       helpText("- analyse the difference between these 2 groups at gene expression level")),
                column(10, offset = 1, 
                       helpText("3)", strong("Signature of interest:"), "a selected gene set will be used to 
                                stratify the individuals into the 2 groups."),
                       br(), br())
              ), 
              
              fluidRow(
                column(5, offset = 1,
                       box(width = NULL, status = "primary",
                           selectInput("select", "Select:",
                                       choices = c("", Gene = "goi", Phenotype = "poi", Signature = "soi")),
                           conditionalPanel(condition = "input.select == 'goi'",
                                            selectInput("select_goi", "Gene of interest:", 
                                                        c('Select or type' = "", genes))
                           ),
                           conditionalPanel(condition = "input.select == 'poi'",
                                            selectInput("select_meta", "Select Phenotype:",
                                                        c('Select or type' = "", Cellcount = "cellcount", Cytokines = "cytokines")),
                                            conditionalPanel(condition = "input.select_meta != ''",
                                                             selectInput("select_poi", "of:", 
                                                                         c('Select or type' = "", phenotype_cellcount)))
                           ),
                           conditionalPanel(condition = "input.select == 'soi'",
                                            selectInput("select_soi", "Signature of interest:", 
                                                        c('Select or type' = "", signatures))))
                ),
                column(4, offset = 1,
                       box(title = icon("cog"), width = NULL, 
                           selectInput("select_adjust.method", "Select an adjust method (optional):",
                                       choices = c('No adjustment' = "none", 'Benjamini & Hochberg' = "BH", 
                                                   'Benjamini & Yekutieli' ="BY", 'Holm' ="holm")),
                           sliderInput("quantile", "Quantile", 
                                       value = 10, min = 1, max = 50, post = " %"))
                )                                
              ),
              fluidRow(
                conditionalPanel(condition = "input.select_goi != '' || input.select_poi != '' || input.select_soi != ''",
                                 conditionalPanel(condition = "input.run == 0", 
                                                  column(5, offset = 5,
                                                         br(),
                                                         actionButton(inputId = "run", label = "RUN APP", 
                                                                      icon = icon("play-circle"), width = "110px"))),
                                 conditionalPanel(condition = "input.run > 0",
                                                  column(5, offset = 5,
                                                         br(),
                                                         useShinyjs(),                                           # Include shinyjs in the UI
                                                         extendShinyjs(text = jsResetCode),                      # Add the js code to the page
                                                         inlineCSS(list(.coral = "background: coral")),
                                                         actionButton("reset_button", "Reset App", icon("undo-alt"), 
                                                                      class = "coral")),
                                                  column(12,
                                                         br(),
                                                         helpText(em("Before rerunning the app with a new gene/phenotype/signature of 
                                                                     interest and/or different settings, hit the"), "RESET", em("button.")))))
              )
      ),
      #2nd tab
        # content of 2nd tab only appears after hitting the "RUN APP" button (on 1st tab)
      tabItem(tabName = "tab_expr",
              conditionalPanel(condition = "input.run == 0",
                               p(em("Press the"), "RUN APP", em("button on the Start tab.")),
                               helpText(em("Button appears after selecting a gene/phenotype/signature of interest."))),
              conditionalPanel(condition = "input.run > 0",
                               # enable to select a study and the cell type for goi & soi
                               conditionalPanel(condition = "input.select != 'poi'",
                                                fluidRow(
                                                  column(5, offset = 1, 
                                                         selectInput("select_study", "Select a study:",
                                                                     choices = c("CEDAR", "classifier", "FG500","ImmVar")),
                                                         selectInput("select_dataset", "Select a cell type:",
                                                                     choices = "")
                                                         # select_dataset will be updated (coded in server) depending on what is selected in select_study
                                                  ),
                                                  column(6,
                                                         br(), br(), br(),
                                                         useShinyjs(),# Set up shinyjs
                                                         inlineCSS(list(.powderblue = "background: powderblue")),
                                                         actionButton("go", "GO", 
                                                                      icon = icon("play-circle"), width = "100px"),
                                                         conditionalPanel(condition = "input.go == 0",
                                                                          helpText(em("Press the"), "GO", 
                                                                                   em("button to start the calculation."))),
                                                         conditionalPanel(condition = "input.go > 0",
                                                                          helpText(em("After selecting a new study and/or cell type,
                                                                                      press the"), "Update", em("button for recalculation."))))
                                                )
                               ),
                               tabsetPanel(type = "tabs",
                                           tabPanel("Plot", icon = icon("bar-chart-o"),
                                                    fluidRow(
                                                      column(width = 4,
                                                             conditionalPanel(condition = "input.select == 'goi'",
                                                                              box(width = NULL, solidHeader = TRUE,
                                                                                  uiOutput("more_genes_goi_UI")),
                                                                              conditionalPanel(condition = "input.go > 0",
                                                                                               helpText(em("Add more genes to visualize the expression levels 
                                                                                                           of the respective genes within the newly generated 
                                                                                                           2 groups – “low” and “high”.")),
                                                                                               helpText(em("After selecting the genes to be displayed in the plot, 
                                                                                                           press the"), "Update", em("button."))
                                                                              )
                                                             ),
                                                             conditionalPanel(condition = "input.select == 'poi'",
                                                                              box(width = NULL, solidHeader = TRUE,
                                                                                  uiOutput("more_genes_poi_UI")),
                                                                              conditionalPanel(condition = "input.more_genes_poi != ''",
                                                                                               box(width = NULL, solidHeader = TRUE,
                                                                                                   uiOutput("go_poi_UI")))),
                                                             conditionalPanel(condition = "input.select == 'soi'",
                                                                              box(width = NULL, solidHeader = TRUE,
                                                                                  uiOutput("more_genes_soi_UI")),
                                                                              conditionalPanel(condition = "input.go > 0",
                                                                                               helpText(em("Select the genes which shall be displayed in the plot.")),
                                                                                               helpText(em("After selecting the genes, press the"), "Update", 
                                                                                                        em("button."))))),
                                                      column(width = 8,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("expr_plot", height = 470)))
                                                    )
                                           ),
                                           tabPanel("Table", icon = icon("table"),
                                                    br(),
                                                    DT::dataTableOutput("expr_table"))
                                )
              )
      ),
      #3rd tab
      tabItem(tabName = "tab_anno",
              conditionalPanel(condition = "input.run == 0",
                               p(em("Press the"), "RUN APP", em("button on the Start tab.")),
                               helpText(em("Button appears after selecting a gene/phenotype/signature of interest."))),
              conditionalPanel(condition = "input.run > 0",
                               conditionalPanel(condition = "input.select != 'poi'",
                                                fluidRow(
                                                  column(5, offset = 1,
                                                         selectInput("t3_select_study", "Select a study:",
                                                                     choices = c("CEDAR", "classifier", "FG500", "ImmVar")),
                                                         selectInput("t3_select_dataset", "Select a cell type:",
                                                                     choices = "")),
                                                  column(6,
                                                         br(), br(), br(),
                                                         actionButton("t3_go", "GO", 
                                                                      icon = icon("play-circle"), width = "100px"),
                                                         conditionalPanel(condition = "input.t3_go == 0",
                                                                          helpText(em("Press the"), "GO", 
                                                                                   em("button to start the calculation."))),
                                                         conditionalPanel(condition = "input.t3_go > 0",
                                                                          helpText(em("After selecting a new study and/or cell type,
                                                                                      press the"), "Update", em("button for recalculation.")))))
                               ),
                               tabsetPanel(type = "tabs",
                                           tabPanel("Plot", icon = icon("bar-chart-o"), 
                                                    fluidRow(
                                                      column(width = 4, 
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 conditionalPanel(condition = "input.t3_go > 0",       # can only appear for GOI & SOI)
                                                                                  uiOutput("t3_plot_select_anno_UI")),
                                                                 conditionalPanel(condition = "input.select == 'poi'",
                                                                                  uiOutput("t3_plot_select_anno_poi_UI")))),
                                                      column(width = 8,
                                                             box(width = NULL, status = "primary",
                                                                 conditionalPanel(condition = "input.t3_go > 0",       # can only appear for GOI & SOI)
                                                                                  plotOutput("anno_plot", height = 470)
                                                                 ),
                                                                 conditionalPanel(condition = "input.select == 'poi'",
                                                                                  plotOutput("anno_plot_poi", height = 470)))))),
                                           tabPanel("Statistics", icon = icon("list-ul"), 
                                                    fluidRow(
                                                      column(width = 4,
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 conditionalPanel(condition = "input.t3_go > 0",
                                                                                  uiOutput("t3_stat_select_anno_UI")),
                                                                 conditionalPanel(condition = "input.select == 'poi'",
                                                                                  uiOutput("t3_stat_select_anno_poi_UI")))),
                                                      column(width = 8,
                                                             box(width = NULL, status = "primary",
                                                                 conditionalPanel(condition = "input.t3_go > 0",
                                                                                  tableOutput("anno_stat")),
                                                                 conditionalPanel(condition = "input.select == 'poi'",
                                                                                  tableOutput("anno_stat_poi")))))),
                                           tabPanel("Table", icon = icon("table"), 
                                                    br(),
                                                    conditionalPanel(condition = "input.t3_go > 0",
                                                                     DT::dataTableOutput("anno_table")),
                                                    conditionalPanel(condition = "input.select == 'poi'",
                                                                     DT::dataTableOutput("anno_table_poi")))
                               )
              )
      ),
      #4th tab
      tabItem(tabName = "tab_de",
              conditionalPanel(condition = "input.run == 0",
                               p(em("Press the"), "RUN APP", em("button on the Start tab.")),
                               helpText(em("Button appears after selecting a gene/phenotype/signature of interest."))),
              conditionalPanel(condition = "input.run > 0",
                               conditionalPanel(condition = "input.select != 'poi'",
                                                fluidRow(
                                                  column(5, offset = 1,
                                                         selectInput("t4_select_study", "Select a study:",
                                                                     choices = c("CEDAR", "classifier", "FG500", "ImmVar")),
                                                         selectInput("t4_select_dataset", "Select a cell type:",
                                                                     choices = "")),
                                                  column(6,
                                                         br(), br(), br(),
                                                         actionButton("t4_go", "GO", 
                                                                      icon = icon("play-circle"), width = "100px"),
                                                         conditionalPanel(condition = "input.t4_go == 0",
                                                                          helpText(em("Press the"), "GO", 
                                                                                   em("button to start the calculation."))),
                                                         conditionalPanel(condition = "input.t4_go > 0",
                                                                          helpText(em("After selecting a new study and/or cell type,
                                                                                      press the"), "Update", em("button for recalculation.")))))
                               ),
                               tabsetPanel(type = "tabs",
                                           tabPanel("PCA", icon = icon("braille"), 
                                                    fluidRow(
                                                      column(8, offset = 2,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("de_pca", height = 470))))),
                                           tabPanel("No. of DE genes", icon = icon("bar-chart-o"), 
                                                    fluidRow(
                                                      column(8, offset = 2,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("de_plot", height = 470))))),
                                           tabPanel("Heatmap", icon = icon("fire"), 
                                                    fluidRow(
                                                      column(12,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("de_hm", height = 470))))),
                                           tabPanel("Statistics", icon = icon("list-ul"),
                                                    br(), 
                                                    DT::dataTableOutput("de_stat")))
              )
      ),
      #5th tab
      tabItem(tabName = "tab_gsea",
              conditionalPanel(condition = "input.run == 0",
                               p(em("Press the"), "RUN APP", em("button on the Start tab.")),
                               helpText(em("Button appears after selecting a gene/phenotype/signature of interest."))),
              conditionalPanel(condition = "input.run > 0",
                               fluidRow(
                                 conditionalPanel(condition = "input.select != 'poi'",
                                                  column(5, offset = 1,
                                                         selectInput("t5_select_study", "Select a study:",
                                                                     choices = c("CEDAR", "classifier", "FG500", "ImmVar")),
                                                         selectInput("t5_select_dataset", "Select a cell type:",
                                                                     choices = ""))),
                                 column(6,
                                        conditionalPanel(condition = "input.select != 'poi'",
                                                         br(), br(), br()),
                                        actionButton("t5_go", "GO", 
                                                     icon = icon("play-circle"), width = "100px"),
                                        conditionalPanel(condition = "input.t5_go == 0",
                                                         helpText(em("Press the"), "GO", 
                                                                  em("button to start the calculation."))),
                                        conditionalPanel(condition = "input.t5_go > 0",
                                                         helpText(em("After selecting a new study and/or cell type,
                                                                     press the"), "Update", em("button for recalculation."))))
                               ),
                               tabsetPanel(type = "tabs",
                                           tabPanel("Plot of Ranked gene list", icon = icon("bar-chart-o"), 
                                                    fluidRow(
                                                      column(4,
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 sliderInput("top_genes_plot", "Nr. of the top genes:", 
                                                                             value = 5, min = 1, max = 20, ticks = FALSE)),          # --> how many max?!???))
                                                             conditionalPanel(condition = "input.t5_go > 0",
                                                                              helpText(em("After changing the number of top genes to be displayed,
                                                                                          press the"), "Update", em("button.")))),            
                                                      column(8,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("gsea_plot", height = 470))))),
                                           tabPanel("Ranked gene list", icon = icon("sort-amount-up"), 
                                                    fluidRow(
                                                      column(4,
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 sliderInput("top_genes_list", "Nr. of the top genes:", 
                                                                             value = 5, min = 1, max = 20, ticks = FALSE)),       # --> how many max?!???
                                                             conditionalPanel(condition = "input.t5_go > 0",
                                                                              helpText(em("After changing the number of top genes to be displayed,
                                                                                          press the"), "Update", em("button.")))),            
                                                      column(8,
                                                             box(width = NULL, status = "primary",
                                                                 tableOutput("gsea_ranked_list"))))),
                                           tabPanel("GSEA Volcano plot", icon = icon("chart-area"), 
                                                    fluidRow(
                                                      column(4,
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 checkboxInput("interactive_volcano", label = "Interactive Volcano Plot", 
                                                                               value = FALSE)),
                                                             conditionalPanel(condition = "input.t5_go > 0",
                                                                              helpText(em("For displaying an interactive Volcano Plot, 
                                                                              check the checkbox and press the"), "Update", em("button.")))),
                                                      column(8,
                                                             box(width = NULL, status = "primary",
                                                                 conditionalPanel(condition = "input.interactive_volcano == false",
                                                                                  plotOutput("gsea_volcano", height = 470)),
                                                                 conditionalPanel(condition = "input.interactive_volcano == true",
                                                                                  plotlyOutput("gsea_int_volcano", height = 470)))))),
                                           tabPanel("GSEA Table", icon = icon("table"), 
                                                    br(),
                                                    DT::dataTableOutput("gsea_table")))
              )
      ),
      #6th tab
      tabItem(tabName = "tab_meta",
              conditionalPanel(condition = "input.run == 0",
                               p(em("Press the"), "RUN APP", em("button on the Start tab.")),
                               helpText(em("Button appears after selecting a gene/phenotype/signature of interest."))),
              conditionalPanel(condition = "input.run > 0",
                               fluidRow(
                                 column(5, offset = 1,
                                        conditionalPanel(condition = "input.select != 'poi'",
                                                         selectInput("t6_select_study", "Select a study:",
                                                                     choices = c("CEDAR", "FG500"))), #"classifier" & "ImmVar" have no metadata
                                        selectInput("t6_select_metadata", "Select:",
                                                    choices = c(Cellcount = "FG500_PBMC_cellcount", Cytokines = "FG500_PBMC_cytokines"))),
                                 column(6,
                                        conditionalPanel(condition = "input.select != 'poi'",
                                                         selectInput("t6_select_dataset", "Select a cell type:",
                                                                     choices = "")),
                                        br(),
                                        actionButton("t6_go", "GO", 
                                                     icon = icon("play-circle"), width = "100px"),
                                        conditionalPanel(condition = "input.t6_go == 0",
                                                         helpText(em("Press the"), "GO", 
                                                                  em("button to start the calculation."))),
                                        conditionalPanel(condition = "input.t6_go > 0",
                                                         helpText(em("After selecting a new study and/or cell type,
                                                                    press the"), "Update", em("button for recalculation."))))
                               ),
                               tabsetPanel(type = "tabs",
                                           tabPanel("Plot", icon = icon("bar-chart-o"), 
                                                    fluidRow(
                                                      column(4,
                                                             box(width = NULL, solidHeader = TRUE,
                                                                 uiOutput("t6_select_metadata1_UI"))), 
                                                      column(8,
                                                             box(width = NULL, status = "primary",
                                                                 plotOutput("meta_plot"))))),
                                           tabPanel("Statistics", icon = icon("list-ul"),
                                                    br(),
                                                    DT::dataTableOutput("meta_stat")),
                                           tabPanel("Table", icon = icon("table"), 
                                                    br(),
                                                    DT::dataTableOutput("meta_table")))
              )) )))      #tabItems, dashboardBody, dashboardPage



server <- function(input, output, session) {
  
  # 1st tab 
    # update the selection of poi if cytokines is selected
  observe({
    input$select_meta
    if (input$select_meta == "cytokines") {
      updateSelectInput(session, "select_poi", 
                        choices = c('Select or type' = "", phenotype_cytokines))
    }
    else if (input$select_meta == "cellcount") {
      updateSelectInput(session, "select_poi", 
                        choices = c('Select or type' = "", phenotype_cellcount))
    }
  })
  
  # let HuVa run:  
  observeEvent(input$run, {            # calculation starts after hitting actionButton "run"
    showModal(modalDialog(title = "Calculation ...", 
                          "Please wait.", size = "l", footer = NULL))    # footer = NULL --> removes the "dismiss" button
        # removeModal (at end of code) is part of that
    if (input$select == "goi") {                                         #  for goi --> l. 150 of huva_workflow 
      binned_dataset <- run_huva_experiment(data = huva.db, 
                                            gene = input$select_goi, 
                                            quantiles = (input$quantile)/100, 
                                            gs_list = hallmarks_V7.2,
                                            summ = T, 
                                            datasets_list = NULL, 
                                            adjust.method = input$select_adjust.method)
      # 2nd tab
        # adjust the selectInputs "select_dataset" & "select_study" 
      huva_names <- names(binned_dataset)
      studies <- huva_names[huva_names != "summary"]
      updateSelectInput(session, "select_study",
                        choices = sort(studies))
      observe({
        input$select_study
        if (input$select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 3rd tab
      updateSelectInput(session, "t3_select_study",
                        choices = sort(studies))
      observe({
        input$t3_select_study
        if (input$t3_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t3_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t3_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t3_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t3_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 4th tab
      updateSelectInput(session, "t4_select_study",
                        choices = sort(studies))
      observe({
        input$t4_select_study
        if (input$t4_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t4_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t4_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t4_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t4_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 5th tab
      updateSelectInput(session, "t5_select_study",
                        choices = sort(studies))
      observe({
        input$t5_select_study
        if (input$t5_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t5_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t5_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t5_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t5_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 6th tab
        # here: classifier & ImmVar have no metadata --> exclude them here
      studies2 <- studies[studies != "classifier"]
      studies3 <- studies2[studies2 != "ImmVar"]
      updateSelectInput(session, "t6_select_study",
                        choices = sort(studies3))
      observe({
        input$t6_select_study
        if (input$t6_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t6_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t6_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t6_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        }
      })
      observe({
        input$t6_select_study
        if (input$t6_select_study == "CEDAR"){
          updateSelectInput(session, "t6_select_metadata",
                            choices = "Cellcount")
        } else if (input$t6_select_study == "FG500"){
          updateSelectInput(session, "t6_select_metadata",
                            choices = c(Cellcount = "FG500_PBMC_cellcount", Cytokines = "FG500_PBMC_cytokines"))
        }
      })
    } 
    
    else if (input$select == "poi") {           # --> l. 300 of huva_workflow
      binned_dataset <- run_huva_phenotype(data = huva.db,
                                           phenotype = input$select_poi,
                                           study = "FG500",
                                           metadata_table = input$select_meta,    
                                           quantiles = (input$quantile)/100, 
                                           adjust.method = input$select_adjust.method,
                                           gs_list = hallmarks_V7.2)
      # 2nd tab - Expression POI 
      expr_huva <- get_expr_huva(huva_exp = binned_dataset, 
                                 study = "FG500", dataset = "FG500_PBMC")
      output$expr_table = DT::renderDataTable(expr_huva, options = list(scrollX = T))
        # add selectInput to enable adding more genes to the boxplot 
      output$more_genes_poi_UI <- renderUI({
        selectInput("more_genes_poi", "Select genes:", 
                    choices = c('Multiple genes selectable' = "", sort(rownames(expr_huva))), multiple = TRUE)
      })
        # hit the actionbutton, to show the plot 
      output$go_poi_UI <- renderUI({
        actionButton("go_poi", label = "GO", 
                     icon = icon("play-circle"), width = "80px")
      })
      observeEvent(input$go_poi, {
        toggleClass("go_poi", "powderblue")
        updateActionButton(session, "go_poi", label = "Update", 
                           icon = character(0))
        plot_binned <- plot_binned_gene(goi = input$more_genes_poi, 
                                        huva_experiment = binned_dataset)
        output$expr_plot <- renderPlot(
          plot_binned$FG500_PBMC
        )
      })
      
      # 3rd tab - Annotation POI 
      anno_huva <- get_anno_huva(huva_exp = binned_dataset, study = "FG500", "FG500_PBMC")
      
      output$anno_table_poi = DT::renderDataTable(anno_huva, options = list(scrollX = T))
      
        # add selectInput for plot for the annotation parameters
      output$t3_plot_select_anno_poi_UI <- renderUI({
        selectInput("t3_plot_select_anno_poi", "Select an annotation parameter:",
                    choices = c('Select or type' = "", sort(colnames(anno_huva)[-1:-4])))
      })
        # generate the plot 
      anno.plot <- get_anno.plot_huva(huva_exp = binned_dataset, study = "FG500")
      output$anno_plot_poi <- renderPlot(
        anno.plot[["FG500_PBMC"]][input$t3_plot_select_anno_poi]
      )
        # add a SelectInput for statistics for the annotation parameters
      output$t3_stat_select_anno_poi_UI <- renderUI({
        selectInput("t3_stat_select_anno_poi", "Select an annotation parameter:",
                    choices = c('Select or type' = "", sort(colnames(anno_huva)[-1:-4])))
      })
      anno.stat <- get_anno.stat_huva(huva_exp = binned_dataset, study = "FG500")
      output$anno_stat_poi <- renderTable({
        anno.stat[["FG500_PBMC"]][input$t3_stat_select_anno_poi]},
        striped = TRUE, hover = TRUE
      )
      
      # 4th tab - DE Genes POI 
      DE_huva <- get_DE_huva(huva_exp = binned_dataset, 
                             study = "FG500", 
                             cluster_col = F, 
                             dataset = "FG500_PBMC")
        # PCA
      output$de_pca <- renderPlot(
        DE_huva["PCA_FG500_PBMC"]) 
        # plot
      output$de_plot <- renderPlot(
        DE_huva["plot_FG500_PBMC"])
        # HM
      output$de_hm <- renderPlot(
        plot_HM(DE_huva$HM_FG500_PBMC))
        # statistics
      output$de_stat = DT::renderDataTable(DE_huva$FG500_PBMC, options = list(scrollX = T))
      
      # 5th tab & 6th tab
       # have a go button --> will be implemented further down, where GOI & SOI are coded as well
    } 
    
    else if (input$select == "soi"){
      binned_dataset <- run_huva_signature(data = huva.db,   
                                           gene_set = hallmarks_V7.2[[input$select_soi]], 
                                           GSVA.method =  "gsva",
                                           quantiles = (input$quantile)/100, 
                                           gs_list = hallmarks_V7.2,
                                           summ = T, 
                                           datasets_list = NULL, 
                                           adjust.method = input$select_adjust.method)
      # 2nd tab
        # adjust the selectInputs "select_dataset" & "select_study" 
      huva_names <- names(binned_dataset)
      studies <- huva_names[huva_names!="summary"]
      updateSelectInput(session, "select_study",
                        choices = sort(studies))
      observe({
        input$select_study
        if (input$select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 3rd tab
      updateSelectInput(session, "t3_select_study",
                        choices = sort(studies))
      observe({
        input$t3_select_study
        if (input$t3_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t3_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t3_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t3_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t3_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t3_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 4th tab
      updateSelectInput(session, "t4_select_study",
                        choices = sort(studies))
      observe({
        input$t4_select_study
        if (input$t4_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t4_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t4_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t4_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t4_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t4_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 5th tab
      updateSelectInput(session, "t5_select_study",
                        choices = sort(studies))
      observe({
        input$t5_select_study
        if (input$t5_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t5_select_study == "classifier"){
          huva_dataset <- names(binned_dataset$classifier$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t5_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t5_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        } else if (input$t5_select_study == "ImmVar"){
          huva_dataset <- names(binned_dataset$ImmVar$data)
          updateSelectInput(session, "t5_select_dataset",
                            choices = sort(huva_dataset))
        }
      })
      # 6th tab
        # here: classifier & ImmVar have no metadata --> exclude them here
      studies2 <- studies[studies != "classifier"]
      studies3 <- studies2[studies2 != "ImmVar"]
      updateSelectInput(session, "t6_select_study",
                        choices = sort(studies3))
      observe({
        input$t6_select_study
        if (input$t6_select_study == "CEDAR"){
          huva_dataset <- names(binned_dataset$CEDAR$data)
          updateSelectInput(session, "t6_select_dataset",
                            choices = sort(huva_dataset))
        } else if (input$t6_select_study == "FG500"){      #don't exchange w/ huva_dataset, as there is only FG500_PBMC in FG500
          updateSelectInput(session, "t6_select_dataset",
                            choices = c(PBMC = "FG500_PBMC"))
        }
      })
      observe({
        input$t6_select_study
        if (input$t6_select_study == "CEDAR"){
          updateSelectInput(session, "t6_select_metadata",
                            choices = "Cellcount")
        } else if (input$t6_select_study == "FG500"){
          updateSelectInput(session, "t6_select_metadata",
                            choices = c(Cellcount = "FG500_PBMC_cellcount", Cytokines = "FG500_PBMC_cytokines"))
        }
      })
    }
    removeModal()
    
    # 2nd tab - Expression GOI & SOI 
    observeEvent(input$go, {     #actionbutton "go" only appears for GOI & SOI (written in ui)
      # actionbutton: switch from "go" to "update" after 1st click & change color
      toggleClass("go", "powderblue")
      updateActionButton(session, "go", label = "Update", 
                         icon = character(0))
      # expression table 
      expr_huva <- get_expr_huva(huva_exp = binned_dataset, 
                                 study = input$select_study, 
                                 dataset = input$select_dataset)
      output$expr_table = DT::renderDataTable(expr_huva, options = list(scrollX = T))
      
      if (input$select == "goi"){
        # add a selectInput to add more genes into the boxplot of expression plot for goi
        output$more_genes_goi_UI <- renderUI({
          selectInput("more_genes_goi", "Add more genes:", 
                      choices = c('Multiple genes selectable' = "", sort(rownames(expr_huva))), 
                      multiple = TRUE)
        })
        # generate the expression plot for goi
        plot_binned <- plot_binned_gene(goi = c(input$select_goi, input$more_genes_goi), 
                                        huva_experiment = binned_dataset)
        output$expr_plot <- renderPlot(
          plot_binned[input$select_dataset]
        )
      }
      else if (input$select == "soi"){
        # add SelectInput to enable adding more genes to the boxplot 
        output$more_genes_soi_UI <- renderUI({
          selectInput("more_genes_soi", "Select genes:", 
                      choices = c('Multiple genes selectable' = "", sort(rownames(expr_huva))), 
                      multiple = TRUE)
        })
        # hit the actionbutton, to show the plot 
        observeEvent(input$go, {
          plot_binned <- plot_binned_gene(goi = input$more_genes_soi, 
                                          huva_experiment = binned_dataset)
          output$expr_plot <- renderPlot(
            plot_binned[input$select_dataset]
          )
        })
      }
    })
    
    # 3rd tab - Annotation GOI & SOI 
    observeEvent(input$t3_go, {     #actionbutton t3_go only appear for GOI & SOI 
      toggleClass("t3_go", "powderblue")
      updateActionButton(session, "t3_go", label = "Update", 
                         icon = character(0))
      # annotation table
      anno_huva <- get_anno_huva(huva_exp = binned_dataset, 
                                 study = input$t3_select_study)
      output$anno_table = DT::renderDataTable(anno_huva[[input$t3_select_dataset]], options = list(scrollX = T))
      # annotation plot
        # add a SelectInput for the annotation parameters
      output$t3_plot_select_anno_UI <- renderUI({
        selectInput("t3_plot_select_anno", "Select an annotation parameter:",
                    choices = c('Select or type' = "", sort(colnames(anno_huva[[input$t3_select_dataset]])[-1:-4])))
      })
        # generate the plot 
      anno.plot <- get_anno.plot_huva(huva_exp = binned_dataset, 
                                      study = input$t3_select_study)
      output$anno_plot <- renderPlot(
        anno.plot[[input$t3_select_dataset]][input$t3_plot_select_anno] 
      )
      # statistics
        # add a SelectInput for the annotation parameters
      output$t3_stat_select_anno_UI <- renderUI({
        selectInput("t3_stat_select_anno", "Select an annotation parameter:",
                    choices = c('Select or type' = "", sort(colnames(anno_huva[[input$t3_select_dataset]])[-1:-4])))
      })
      anno.stat <- get_anno.stat_huva(huva_exp = binned_dataset, 
                                      study = input$t3_select_study)
      output$anno_stat <- renderTable({
        anno.stat[[input$t3_select_dataset]][input$t3_stat_select_anno]},
        striped = TRUE, hover = TRUE
      )
    })
    
    # 4th tab - DE genes GOI & SOI 
    observeEvent(input$t4_go, {     #actionbutton t4_go only appear for GOI & SOI 
      toggleClass("t4_go", "powderblue")
      updateActionButton(session, "t4_go", label = "Update", 
                         icon = character(0))
      DE_huva <- get_DE_huva(huva_exp = binned_dataset, 
                             study = input$t4_select_study, 
                             dataset = input$t4_select_dataset)
      # PCA
      de_PCA <- paste("PCA_", input$t4_select_dataset, sep = "")
      output$de_pca <- renderPlot(
        DE_huva[de_PCA]     
      )
      # plot
      de_PLOT <- paste("plot_", input$t4_select_dataset, sep = "")
      output$de_plot <- renderPlot(
        DE_huva[de_PLOT]           
      )
      # HM
      de_HM <- paste("HM_", input$t4_select_dataset, sep = "")
      output$de_hm <- renderPlot(
        plot_HM(DE_huva[[de_HM]])     
      )
      # statistics
      output$de_stat = DT::renderDataTable(DE_huva[[input$t4_select_dataset]], options = list(scrollX = T))
    })
    
    # 5th tab - GSEA GOI & SOI 
    observeEvent(input$t5_go, {   
      toggleClass("t5_go", "powderblue")
      updateActionButton(session, "t5_go", label = "Update", 
                         icon = character(0))
      if (input$select != "poi"){ 
        rank_huva_plot <- get_rank_huva(huva_exp = binned_dataset, 
                                        study = input$t5_select_study, 
                                        dataset = NULL, 
                                        n_top_genes = input$top_genes_plot)
        rank_huva_list <- get_rank_huva(huva_exp = binned_dataset, 
                                        study = input$t5_select_study, 
                                        dataset = NULL, 
                                        n_top_genes = input$top_genes_list)
        # Plot of Ranked gene list
        ranked_gene_plot <- paste("plot_", input$t5_select_dataset, sep = "")
        output$gsea_plot <- renderPlot(
          rank_huva_plot[ranked_gene_plot]
        )
        # Ranked gene list
        output$gsea_ranked_list <- renderTable({
          rank_huva_list[input$t5_select_dataset]},
          striped = TRUE, rownames = TRUE, hover = TRUE
        )
        gsea_huva <- get_gsea_huva(huva_exp = binned_dataset, 
                                   study = input$t5_select_study)
        # GSEA Volcano plot
        gsea_volc <- paste("plot_", input$t5_select_dataset, sep = "")
        gsea_volc_int <- paste("int_plot_", input$t5_select_dataset, sep = "")
        if(input$interactive_volcano){            
          output$gsea_int_volcano <- renderPlotly(
            gsea_huva[[gsea_volc_int]]
          )} else {
            output$gsea_volcano <- renderPlot(
              gsea_huva[gsea_volc]
            )}
        # GSEA table
        output$gsea_table = DT::renderDataTable(gsea_huva[[input$t5_select_dataset]], options = list(scrollX = T))
      }
      else if (input$select == "poi"){
        rank_huva_plot <- get_rank_huva(huva_exp = binned_dataset, 
                                        study = "FG500", 
                                        dataset = NULL, 
                                        n_top_genes = input$top_genes_plot)
        rank_huva_list <- get_rank_huva(huva_exp = binned_dataset, 
                                        study = "FG500", 
                                        dataset = NULL, 
                                        n_top_genes = input$top_genes_list)
        # Plot of Ranked gene list
        output$gsea_plot <- renderPlot(
          rank_huva_plot$plot_FG500_PBMC)
        # Ranked gene list
        output$gsea_ranked_list <- renderTable({
          rank_huva_list$FG500_PBMC},
          striped = TRUE, rownames = TRUE, hover = TRUE
        )
        gsea_huva <- get_gsea_huva(huva_exp = binned_dataset, 
                                   study = "FG500")
        # GSEA Volcano plot
        if(input$interactive_volcano){            
          output$gsea_int_volcano <- renderPlotly(
            gsea_huva$int_plot_FG500_PBMC
          )} else {
            output$gsea_volcano <- renderPlot(
              gsea_huva$plot_FG500_PBMC)
          }
        # GSEA table
        output$gsea_table = DT::renderDataTable(gsea_huva$FG500_PBMC, options = list(scrollX = T))
      }
    })
    
    #6th tab - Metadata GOI & POI & SOI 
    observeEvent(input$t6_go, {     
      toggleClass("t6_go", "powderblue")
      updateActionButton(session, "t6_go", label = "Update", 
                         icon = character(0))
      if (input$select != "poi"){      
        if (input$t6_select_study == "CEDAR"){
          # metadata here is always cellcount, only the dataset is changing
          meta1 <- paste(input$t6_select_dataset, "_cellcount_CEDAR", sep = "")
          # here: dataset always the same (PBMC), metadata either cellcount or cytokines
        } else if (input$t6_select_study == "FG500"){
          meta1 <- input$t6_select_metadata
        }
        # metadata plot
        meta.plot <- get_meta.plot_huva(huva_exp = binned_dataset, 
                                        study = input$t6_select_study)
        # metadata table
        meta_huva <- get_meta_huva(huva_exp = binned_dataset, 
                                   study = input$t6_select_study)
        # metadata statistics
        meta.stat <- get_meta.stat_huva(huva_exp = binned_dataset, 
                                        study = input$t6_select_study, 
                                        dataset = meta1)
        output$meta_stat = DT::renderDataTable(meta.stat, 
                                               options = list(scrollX = T))
      } else if (input$select == "poi"){
        meta1 <- input$t6_select_metadata
        # metadata plot
        meta.plot <- get_meta.plot_huva(huva_exp = binned_dataset, study = "FG500")
        # metadata table
        meta_huva <- get_meta_huva(huva_exp = binned_dataset, study = "FG500")
        # metadata statistics
        meta.stat <- get_meta.stat_huva(huva_exp = binned_dataset, study = "FG500")
        output$meta_stat = DT::renderDataTable(meta.stat[[meta1]], options = list(scrollX = T))
      }
      # plot 
        # add a selectInput for selecting metadata (cell, cytokines, bmi, ...)
      output$t6_select_metadata1_UI <- renderUI({
        selectInput("t6_select_metadata1", "Select the metadata:",
                    choices = c('Select or type' = "", names(meta.plot[[meta1]])[-1]))
      })
      output$meta_plot <- renderPlot(
        meta.plot[[meta1]][input$t6_select_metadata1])
      
      # table
      output$meta_table = DT::renderDataTable(meta_huva[[meta1]], options = list(scrollX = T))
    })
  })
  observeEvent(input$reset_button, {js$reset()})          # for the reset button
}


shinyApp(ui, server)
