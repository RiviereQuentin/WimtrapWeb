#' AdaptativeSelection UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' @importFrom tools file_ext
#' @importFrom data.table fread as.data.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom stats model.matrix cor predict
#' @importFrom caret findCorrelation confusionMatrix
#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train xgb.importance
#' @importFrom pROC roc auc
mod_AdaptativeSelection_ui <- function(id) {
  ns <- NS(id)
  tagList(navbarPage(
    theme = bs_theme(bootswatch = "cerulean"),
    title = "Wimtrap",
    tabPanel("Home",
                 img(src = "www/favicon.png", height = 175, width = 200, align = "right"),
                 h1("Welcome to Wimtrap!"),
                 h2("Introduction"),
                 tags$hr(),
             tags$div(
               "Wimtrap is a web-application for the prediction of ", tags$b("condition-specific"), 
               " transcription factor binding sites (TFBS) in plant species. The tool locates the potential TF binding-sites 
                   by pattern-matching and annotate them with genomic features that characterize their genomic context
                   (opening of the DNA, chromatin marks, DNA conservation, digital genomic footprints,...). Then, the potential TF
                   binding sites likely to be TFBS in the condition of interest are selected based on a \'decision rules` model. This model is obtained
                   by maching learning using ChIP-seq data as reference.", tags$br(), tags$br(), 
               "Prebuilt-models and genomic features integrated into the website allow to quickly
                   obtain prediction of TFBS and related potential gene targets in 10 different conditions for ", tags$i("Arabidopsis thaliana"), " ,
                   2 different conditions for ", tags$i("Solanum lycopersicum"), ", 2 different conditions for ", tags$i("Oryza sativa"),
               " , and 1 condition for ", tags$i("Zea mays"),". Go to the", tags$b("\`Query pre-built models\`"), " panel to take advantage of these ressources.",
               tags$br(), tags$br(), "The ", tags$b("\`Build custom models\`"), " and ",  tags$b("\`Use custom models\`"), " panels allow you to build your own models and use them for any other condition and organism"),
             h2("Tutorial"),
                tags$hr(),
                tags$div("Watch this tutorial video to take in hands Wimtrap. You will find all the necessary information related to the 
                         possibilities offered by the app, the data required and their format, the different options and the outputs."),
                HTML('<iframe width="560" height="315" src="https://www.youtube.com/watch?v=6371fN7dkak" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>'),
                h2("Further resources"),
                tags$hr(),
                tags$div("Please refer to the ", tags$a(href = "https://htmlpreview.github.io/?https://github.com/RiviereQuentin/Wimtrap/blob/main/vignettes/Wimtrap.html",
                "Wimtrap user guide"), "for further methodological details.
                Example data used in the tutorial and source data for ", tags$i("Arabidopsis thaliana"), " and ", tags$i("Solanum lycopersicum"),  " 
                can be found on the ", tags$a(href = "https://github.com/RiviereQuentin/carepat", "carepat github repository") ),
                h2("Authors"),
                tags$hr(),
                tags$div("This web-application was developed at the ", tags$a(href = "https://lpgmp.ulb.be/", "LPGMP"), " and ", tags$a(href = "https://mlg.ulb.ac.be/wordpress/", "MLG"), 
                         "groups of the Université Libre de Bruxelles by Quentin Rivière <qri@hotmail.be>, Madalina Ciortan, Massimiliano Corso, Grégoire Noël, Nathalie Verbruggen and Matthieu Defrance <Matthieu.Defrance@ulb.be>.")
                
                
                 ),
    tabPanel("Query pre-built models",
             
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(
                   inputId = ns("organism"),
                   label = "Organism",
                   choices = names(species_metadata)
                 ),
                 
                 selectizeInput(
                   inputId = ns("condition"),
                   label = "Condition",
                   choices = NULL
                 ),
                 
                 tags$hr(),
                 
                 selectizeInput(
                   inputId = ns("tf"),
                   label = "Transcription factor",
                   choices = NULL,
                   multiple = TRUE
                 ),
                 
                 
                  tags$div("OR"),
                 
                 tags$br(),
                   
                  WimtrapWeb:::fileWimtrap(
                     inputId = ns("pwm"),
                     label = tags$div("PWM file", "(.meme, .pfm, .jaspar,", tags$br(), ".transfac, .motif, or cis-bp format)"),
                     accept = c(".meme", ".pfm", ".jaspar", ".transfac",
                                ".motif", ".txt"),
                     multiple = TRUE
                   ),
                   
                 tags$hr(),
                 
                 sliderInput(
                   inputId = ns("score_threshold"),
                   label = "Prediction score threshold",
                   min = 0.5,
                   max = 1,
                   value = 0.86
                 ),
                 

                 actionButton(
                   inputId = ns("submit"),
                   label = "Submit",
                   icon = icon("leaf")
                 ),
                 
                 tags$h2("Prediction score"),
                 tags$div(
                   "The prediction score is the numeric value output by the model for each candidate TFBS. This value is comprised
                          between 0 and 1. A threshold needs to be set to convert the numeric predictions into binary predictions (transcription
                          factor binding sites or not). If the prediction score is superior to the threshold, then the candidate is predicted
                          as a TFBS. Otherwise, it is not. Higher is set the threshold, more specific are the predictions. Lower is set the
                          threshold, more sensitive are the predictions. The best balance between sensitivity and specificity for predicting
                          TFBS is obtained with a threshold of 0.5; for inferring target genes from the predicted TFBS, it is obtained with a
                          threshold of 0.86."),
             ),
               
               mainPanel(
                 downloadButton(ns("downloadData1"), "Download TFBS predictions"),
                 
                 downloadButton(ns("downloadData2"), "Download Target predictions"),
                 br(),
                 tags$div("N.B.: Genomic coordinates are expressed according to the following
                          assemblies: TAIR10 (", tags$i("Arabidopsis thaliana"), "), SL3.0 (",
                          tags$i("Solanum lycopersicum"), "), IRGSP-1.0 (", tags$i("Oryza sativa"),
                          "), and Zm-B73-REFERENCE-NAM-5.0 (", tags$i("Zea mays"), ")."),
                dataTableOutput(ns(
                 "CAREpredictions"
               )),
               
               dataTableOutput(ns(
                 "Targetpredictions"
               )))
             )),
    tabPanel("Build custom models",
             sidebarLayout(
               sidebarPanel(
                 tags$style(
                   "#fieldset {
          color: #373939;
          background-color: #eeeeee;
          border: solid;
                   }"
                 ),
                 
                 tags$div(id= "fieldset",
                          WimtrapWeb:::fileWimtrap(
                            inputId = ns("pwmb"),
                            label = tags$div("PWM file", "(.meme, .pfm, .jaspar,", tags$br(), ".transfac, .motif, or cis-bp format)"),
                            accept = c(".meme", ".pfm", ".jaspar", ".transfac",
                                       ".motif", ".txt"),
                            multiple = TRUE
                            )
                          ),
                 br(),
                 tags$div(id = "fieldset",
                          WimtrapWeb:::fileWimtrap(
                            inputId = ns("filesChIP"),
                            label = tags$div(
                              "TF ChIP-peaks data (.bed* or .narrowPeak)",
                              tags$br(),
                              tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "Press 'Ctrl' to select multiple files"),
                              tags$br(),
                              tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "*BED3 format")
                            ),
                            accept = ".bed",
                            placeholder = "No file(s) selected",
                            multiple = TRUE
                          ),
                          "Associate each TF encoded in the peak files with the name of its motif in the PWM file:",
                          br(),
                          tags$em("Tip: don't use twice the same motif name"),
                          br(),
                          br(),
                          tags$div(id = 'placeholder1')
                          ),
                 br(),
                 tags$div(id = "fieldset",
                   WimtrapWeb:::fileWimtrap(
                     inputId = ns("filesData"),
                     label = tags$div(
                       "Genomic data (.bed** or .gtf)",
                       tags$br(),
                       tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "Press 'Ctrl' to select multiple files"),
                       tags$br(),
                       tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "**Input genomic regions (such as conserved
                                 element, digital genomic footprint, open DNA region,...):"),
                       tags$br(),
                       tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "  - without a score => in BED3 format."),
                       tags$br(),
                       tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "  - with a score => in BED6 format. Fill the \'score\' field accordingly.
                                 The \'name\' and \'strand\' fields can be filled with \".\""),
                       tags$br(),
                       tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "  - with a category => in BED6 format. Fill the \'name\' field accordingly.
                                 The \'score\' and \'strand\' fields can be filled with \".\"")
                       
                     ),
                     accept = c(".bed", ".gtf"),
                     placeholder = "No file(s) selected",
                     multiple = TRUE
                   ),
                   
                   "Name the genomic features:",
                   br(),
                   br(),
                   tags$div(id = 'placeholder2')
                 ),
                 br(), 
                 
                 checkboxInput(ns("biomart"),
                               label = "Automatically download the genome sequence and the transcript models from Ensembl Plants and BioMart",
                               value = TRUE),
                 uiOutput(outputId = ns("biomartui")),

                 
                 
                 #checkboxGroupInput(inputId = ns("checkoptional"),
                 #                   label = "Check list",
                 #                   choices = c("TSSs are annotated with the name and orientation of the transcripts",
                 #                               "TTss are annotated with the name and orientation of the transcripts")),
                 #)
                 
                 tags$em(
                   "Tip: name consistently the chromosomes in all the files provided, following the Ensembl Plants notation"
                 ),
                 
                 br(),
                 
                 actionButton(inputId = ns("build"),
                              label = "Build!",
                              icon = icon("leaf")),
                 tags$h2("Demo"),
                 tags$div(
                   "Please browse to the Demo tab and find the instructions to run the demo example.")
                 
                 ),
               
               mainPanel(
                 downloadButton(ns("downloadModel"), "Download TFBS model"),

                 downloadButton(ns("downloadROC"), "Download ROC"),

                 downloadButton(
                   ns("downloadFeatureImportance"),
                   "Download plot of feature importance"
                 ),
                 verbatimTextOutput(ns("filepathtss")),
                 verbatimTextOutput(ns("filepathtts")),
                 h4("ROC curve (5x-cross validation):"),
                 tags$div("ROC curves were computed considering all the TFs that were input to build the model. ROC curve of the model, based
                 on the prediction score, is compared to that obtained when predicting TFBS based solely on the results of pattern-matching. The AUC of the model
                  is indicated on the y-axis."),
                 plotOutput(ns("roc")),
                 tags$hr(),
                 h4("Confusion matrix (prediction score threshold = 0.5):"),
                 tags$div("Confusion matrix and derived performance metrics obtained when applying the model to a randomly-selected test dataset (
                          comprising 20% of the TFBS validated by ChIP-seq for all the TFS that were input). The variable to predict takes 0 when
                          the candidate is not validated by ChIP-seq, 1 otherwise."),
                 verbatimTextOutput(ns("confusion")),
                 
                 tags$hr(),
                 h4("Feature Importance (Gain)"),
                 tags$div("The feature importance is assessed using the ", tags$a(href = "https://towardsdatascience.com/be-careful-when-interpreting-your-features-importance-in-xgboost-6e16132588e7", "Gain"),
                          ", a metrics specific to the algorithm of machine learning implemented by Wimtrap (XGBoost). The suffixes
                          \'_20bp\', \'_400bp\' and \'_1000bp\' refer to the length the window lenght on which the predictive features were extracted to characterize the regions surrounding the TFBS regions."),
                 plotOutput(ns("feature"))
               )
             )),
    
    tabPanel("Use custom models",
             sidebarLayout(
               sidebarPanel(
                 WimtrapWeb:::fileWimtrap(
                   inputId = ns("pwmc"),
                   label = "PWM file (.meme, .pfm, .jaspar, .transfac, .motif, or cis-bp format)",
                   accept = c(".meme", ".pfm", ".jaspar", ".transfac",
                              ".motif", ".txt"),
                   multiple = TRUE
                 ),
                 
                 textInput(
                   inputId = ns("motifsc"),
                   label = "Enter the name of the motifs of interest, as encoded in the PWM file",
                   placeholder = "Separate the name of the motifs with commas"
                 ),
                 
                 tags$div(
                   HTML("<label>"),
                   HTML("<input type=\"checkbox\" name=\"checkc\" onclick=\"onlyOne(this)\" value=1 id=\"AdaptativeSelection_ui_1-checkc1\" Checked>"),
                   HTML("%<%: Pipe directly from `build` panel. Reuse the data input previously so that predictions are made\n
                        in the same organism and condition considered to build the model</label>"),
                   HTML("<label>"),
                   HTML("<input type=\"checkbox\" name=\"checkc\" onclick=\"onlyOne(this)\" value=2 id=\"AdaptativeSelection_ui_1-checkc2\">"),
                   HTML("Automatically download the genome sequence and the transcript models from Ensembl Plants and BioMart</label>"),
                   HTML("<label>"),
                   HTML("<input type=\"checkbox\" name=\"checkc\" onclick=\"onlyOne(this)\" value=3 id=\"AdaptativeSelection_ui_1-checkc3\">"),
                   HTML("Manually input the data</label><br></br>")
                 ),
                 
                 conditionalPanel(
                   condition = "input.checkc3 == true",
                   tags$div(id = "fieldset",
                            WimtrapWeb:::fileWimtrap(
                              inputId = ns("TFBSmodel"),
                              label = "TFBS model (.Rdata,  .rda)",
                              accept = c(".RData", ".rda"),
                              multiple = FALSE
                            ),
                            
                            verbatimTextOutput(outputId = ns("PredictiveFeatures"))
                   ),
                   br(),
                   
                   tags$div(id = "fieldset",
                            tags$div(id = 'placeholder3')
                   ),
                   
                   br(),
                   
                   tags$div(id = "fieldset",
                            
                            WimtrapWeb:::fileWimtrap(
                              inputId = ns("genomec"),
                              label = "Genomic sequence (.fasta, .fa. fasta.gzip, fa.gzip)",
                              accept = c(".fa", ".fasta", ".fa.gzip", ".fasta.gzip")
                            ),
                            
                            tags$div(id = 'placeholder4'),
                            
                            
                            WimtrapWeb:::fileWimtrap(
                              inputId = ns("tssc"),
                              label = "Location of TSSs (.bed)",
                              accept = ".bed"
                            ),
                            
                            WimtrapWeb:::fileWimtrap(
                              inputId = ns("ttsc"),
                              label = "Location of TTSs (.bed)",
                              accept = ".bed"
                            ))
                   ,
                   ns = NS(id)
                 ),
                 
                 conditionalPanel(
                   condition = "input.checkc2 == true",
                   tags$div(id = "fieldset",
                            WimtrapWeb:::fileWimtrap(
                              inputId = ns("TFBSmodelc"),
                              label = "TFBS model (.RData,  .rda)",
                              accept = c(".RData", ".rda"),
                              multiple = FALSE
                            ),
                            
                            verbatimTextOutput(outputId = ns("PredictiveFeaturesc"))
                   ),
                   br(),
                   
                   tags$div(id = "fieldset",
                            tags$div(id = 'placeholder5')
                   ),
                   
                   br(),
                   
                   tags$div(
                     id = "fieldset",
                     textInput(inputId = ns("organismc"),
                               label = "Organism"),
                     numericInput(
                       inputId = ns("promoterc"),
                       label = "Promoter length (from TSS)",
                       value = 2000
                     ),
                     
                     numericInput(
                       inputId = ns("proximalec"),
                       label = "Proximal promoter length (from TSS)",
                       value = 500
                     ),
                     
                     numericInput(
                       inputId = ns("downstreamc"),
                       "Downstream region length (from TTS)",
                       value = 1000
                     )
                     ), 
                   ns = NS(id)
                 ),
                 
                 
                 
                 #checkboxGroupInput(inputId = ns("checkoptional"),
                 #                   label = "Check list",
                 #                   choices = c("TSSs are annotated with the name and orientation of the transcripts",
                 #                               "TTss are annotated with the name and orientation of the transcripts")),
                 #)

                 
                 sliderInput(
                   inputId = ns("score_thresholdc"),
                   label = "Prediction score threshold",
                   min = 0.5,
                   max = 1,
                   value = 0.86
                 ),
                 
                 tags$em(
                   "Tip: name consistently the chromosomes in all the files provided, following the Ensembl Plants notation"
                 ),
                 
                 br(),
                 
                 actionButton(inputId = ns("predict"),
                              label = "Predict!",
                              icon = icon("leaf")),
                 
                 tags$h2("Demo"),
                 tags$div(
                   "Please browse to the Demo tab and find the instructions to run the demo example."),
                 tags$br(),
                 tags$h2("Prediction score"),
                 tags$div(
                 "The prediction score is the numeric value output by the model for each candidate TFBS. This value is comprised
                          between 0 and 1. A threshold needs to be set to convert the numeric predictions into binary predictions (transcription
                          factor binding sites or not). If the prediction score is superior to the threshold, then the candidate is predicted
                          as a TFBS. Otherwise, it is not. Higher is set the threshold, more specific are the predictions. Lower is set the
                          threshold, more sensitive are the predictions. The best balance between sensitivity and specificity for predicting
                          TFBS is obtained with a threshold of 0.5; for inferring target genes from the predicted TFBS, it is obtained with a
                          threshold of 0.86.")
                  
               ),
               
               mainPanel(
                 
                 downloadButton(ns("downloadData3"), "Download TFBS predictions"),
                 
                 downloadButton(ns("downloadData4"), "Download Target predictions"),
                 
                 dataTableOutput(ns(
                 "CAREpredictions2"
               )),
               
               dataTableOutput(ns(
                 "Targetpredictions2"
               )))
               )),
    
    tabPanel(title = "Demo",
             mainPanel(
               tags$h2("Demo"),
               tags$br(),
               tags$div("Please download the ", 
                        tags$a(href = "https://github.com/RiviereQuentin/carepat", "carepat github repository"),
                        " (click Code -> Download ZIP) and browse to the \'demo\' folder. Then, we invite you to follow the instructions from the tutorial video.",
                        tags$br(),
                        tags$h4("Application cases considered:"),
                        tags$ul(
                          tags$li("->", tags$b("Training"), " from data obtained in the ", tags$b("condition 1"), " and predictions 
                                  in the ", tags$b("same condition"), " (in the organism 1)"),
                          tags$li("->", tags$b("Training"), " from data obtained in the ", tags$b("condition 1"), " and predictions 
                                  in the ", tags$b("condition2"), " (in the organism 1)"),
                        ),
                        tags$h4("TFs identified to build the model:"),
                        tags$b("TFexample1"), " and ", tags$b("TFexample2"), ", studied by ", 
                        tags$b("ChIP-seq"), " in the ", tags$b("condition1"), " => PWMs in PFM 
                          and location of ChIP-peaks in BED3 or narrowPeak",
                        tags$h4("TF considered for making predictions:"),
                        "-> ", tags$b("TFexample3"), " => PWM in PFM",
                        tags$h4("Predictive features:"),
                        tags$ul(
                          tags$li("-> ", tags$b("Pattern-macthing"), " results => Genome sequence in FASTA"),
                          tags$li("-> ", tags$b("Overlap"), " of the TFBS candidates with ", tags$b("promoter,
                                    proximal promoter, cds, intron, 5\'UTR, 3\'UTR or downstream regions"), " => BED6")
                        ),
                        "+ On windows of 20, 400 and 1000bp centered around the TFBS candidates:",
                        tags$ul(
                          tags$li("-> Average ", tags$b("density"), " of conserved elements (", tags$b("CE"),") -> BED3"),
                          tags$li("-> Average digital genomic footprint (", tags$b("DGF"), ") ", tags$b("score"), " -> BED6"),
                          tags$li("-> Average DNAseI hypersensitivity (", tags$b("DHS"), ") ", tags$b("score"), " -> BED6")
                        )
                        
                        
                        
               )
               ,
               width = 12
             ))
  ))
}

#' AdaptativeSelection Server Functions
#'
#' @noRd
mod_AdaptativeSelection_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    
    options(shiny.maxRequestSize = 200 * 1024 ^ 2)
    
    observeEvent(input$organism,
                 {
                   updateSelectizeInput(inputId = "condition",
                                        choices = species_metadata[[input$organism]][["conditions"]])
                 })
    
    observeEvent(input$organism,
                 {
                   updateSelectizeInput(inputId = "tf",
                                        choices = species_metadata[[input$organism]][["TFs"]])
                 })
    
    
    pwm_react <- reactiveValues(TFs = NULL)
    
    observeEvent(input$pwm, {
      if (is.null(input$pwm))
        return()
      pwmimp <- Wimtrap:::readPSFM(input$pwm$datapath)
      pwm_react$TFs <- names(pwmimp)
    })
    
    
    v <- reactiveValues(data = NULL)
    
    observeEvent(input$submit,
                 {
                   v$data <- rnorm(100)
                 })
    
    CARE_predictions <- reactive({
      if (is.null(v$data))
        return(data.table::data.table(NULL))
      show_modal_spinner(spin = "cube-grid",
                         color = "#428bca",
                         text = "Please wait 5-10 min...")
      
      results <- Wimtrap::carepat(
        organism = input$organism,
        condition = input$condition,
        pfm = input$pwm$datapath,
        TFnames = c(input$tf, pwm_react$TFs),
        score_threshold = input$score_threshold
      )
      remove_modal_spinner()
      return(results)
    })
    
    output$CAREpredictions <- renderDataTable({
      CARE_predictions()
    })
    
    Target_predictions <- reactive({
      if (is.null(v$data))
        return(data.table::data.table(NULL))
      show_modal_spinner(spin = "cube-grid",
                         color = "#428bca",
                         text = "Please wait...")
      transcripts <-
        split(x = CARE_predictions(), f = CARE_predictions()$TF)
      genes <- lapply(transcripts, function(pred) {
        targets <-
          unlist(strsplit(x = as.character(pred$transcript), split = "[.]"))[seq(1, 2 *
                                                                                   nrow(pred), 2)]
        targets <- targets[!duplicated(targets)]
        targets <- data.table::data.table(Gene = targets,
                                          TF = pred$TF[1])
        return(targets)
      })
      genes <- do.call(rbind, genes)
      remove_modal_spinner()
      return(genes)
    })
    
    output$TargetPredictions <- renderDataTable({
      Target_predictions()
    })
    
    output$downloadData1 <- downloadHandler(
      filename = function() {
        paste("TFBS_", Sys.Date(), "_", Sys.time(), ".tsv", sep = "")
      },
      content = function(file) {
        write.table(CARE_predictions(),
                    file,
                    sep = "\t",
                    row.names = FALSE)
      }
    )
    
    output$downloadData2 <- downloadHandler(
      filename = function() {
        paste("Target_", Sys.Date(), "_", Sys.time(), ".tsv", sep = "")
      },
      content = function(file) {
        write.table(Target_predictions(),
                    file,
                    sep = "\t",
                    row.names = FALSE)
      }
    )
    
    output$biomartui <- renderUI({
      if (input$biomart){
        tags$div(
          id = "fieldset",
          textInput(inputId = ns("organismb"),
                    label = "Organism"),
          numericInput(
            inputId = ns("promoterb"),
            label = "Promoter length (from TSS)",
            value = 2000
          ),
          
          numericInput(
            inputId = ns("proximalb"),
            label = "Proximal promoter length (from TSS)",
            value = 500
          ),
          
          numericInput(
            inputId = ns("downstreamb"),
            "Downstream region length (from TTS)",
            value = 1000
          )
        )
        
        } else {
        tags$div(
          id = "fieldset",
          WimtrapWeb:::fileWimtrap(
            inputId = ns("genome"),
            label = "Genomic sequence (.fasta, .fa, .fasta.gzip, .fa.gzip)",
            accept = c(".fa", ".fasta", ".fa.gzip", ".fasta.gzip")
          ),
          
          numericInput(
            inputId = ns("promoter"),
            label = "Promoter length (from TSS)",
            value = 2000
          ),
          
          numericInput(
            inputId = ns("proximal"),
            label = "Proximal promoter length (from TSS)",
            value = 500
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("x5utr"),
            label = "Location of 5'UTRs (.bed***) (facultative)",
            accept = ".bed"
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("cds"),
            label = "Location of coding sequences (.bed***) (facultative)",
            accept = ".bed"
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("intron"),
            "Location of introns (.bed***) (facultative)",
            accept = ".bed"
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("x3Utr"),
            label = "Location of 3'UTRs (.bed***) (facultative)",
            accept = ".bed"
          ),
          
          numericInput(
            inputId = ns("downstream"),
            "Downstream region length (from TTS)",
            value = 1000
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("tss"),
            label = "Location of TSSs (.bed***)",
            accept = ".bed"
          ),
          
          WimtrapWeb:::fileWimtrap(
            inputId = ns("tts"),
            label = "Location of TTSs (.bed***)",
            accept = ".bed"
          ),
          tags$br(),
          tags$span(style = "font-weight: 500; font-size: 14px; color: grey", "***BED6 format. Fill the \'name\' field accordingly to the transcript name. For TSS and TTS, 
                    specify the orientation of the transcripts in the \'strand\' field. Otherwise, the \'strand'\ field can be filled with '\".\"")
          )
        
      }
    })
  
    
    
    output$filepathtss <- renderText({
      if (is.null(input$tss)){
        return()
        } else {
          if (length(input$tss$filepath) == 0) {return()}
          else {
      tsstest <-
        rtracklayer::import(as.character(input$tss$filepath))
      print(tsstest)
      if ((as.character(GenomicRanges::strand(tsstest)[1]) != "*") &&
          ("name" %in% names(GenomicRanges::mcols(tsstest))))
        return()
      "Strand and name fields in the TSSs files need to be fullfilled according to the orientation and name of the
      transcript models"
          }
          }
    })
    
    
    output$filepathtts <- renderText({
      if (is.null(input$tts)){
        return()
      } else {
        if (length(input$tts$filepath) == 0) {return()}
        else {
      tsstest <-
        rtracklayer::import(as.character(input$tts$filepath))
      if ((as.character(GenomicRanges::strand(ttstest)[1]) != "*") &&
          ("name" %in% names(GenomicRanges::mcols(ttstest))))
        return()
      "Strand and name fields in the TTSs files need to be fullfilled according to the orientation and name of the
      transcript models"
        }
      }
    })
    
    w <- reactiveValues(data = NULL)
    
    observeEvent(input$build,
                 {
                   w$data <- rnorm(100)
                 })
    
    Wimtrap_model <- eventReactive(input$build, {
      shinybusy::show_modal_spinner(spin = "cube-grid",
                                    color = "#428bca",
                                    text = "Please wait 20-60 min...")
      bed_data <- reactiveValuesToList(input)
      bed_data <-
        bed_data[grep(pattern = "bed", x = names(bed_data))]
      bed_data <- unlist(bed_data)
      TFnames.carepat <- reactiveValuesToList(input)
      TFnames.carepat <-
        TFnames.carepat[grep(pattern = "txt", x = names(TFnames.carepat))]
      TFnames.carepat <- unlist(TFnames.carepat)
      genomic_data.carepat <- as.character(input$filesData$datapath)
      names(genomic_data.carepat) <- bed_data
      if (input$biomart) {
        imported_genomic_data.carepat <- tryCatch(
          Wimtrap::importGenomicData(organism = as.character(input$organismb),
                                     genomic_data = genomic_data.carepat,
                                     promoter_length = as.numeric(input$promoterb),
                                     downstream_length = as.numeric(input$downstreamb),
                                     proximal_length = as.numeric(input$proximalb)),
          error = function(e){return(NA)},
          finally = message("Interrogating Ensembl and Biomart..."))
        if(is.na(imported_genomic_data.carepat)){
          imported_genomic_data.carepat <- 
            Wimtrap::importGenomicData(organism = as.character(input$organismb),
                                       genomic_data = genomic_data.carepat,
                                       promoter_length = as.numeric(input$promoterb),
                                       downstream_length = as.numeric(input$downstreamb),
                                       proximal_length = as.numeric(input$proximalb))
        } else {}
        TFBSdata.carepat <-
          Wimtrap::getTFBSdata(
            pfm = as.character(input$pwmb$datapath),
            TFnames = TFnames.carepat,
            organism = as.character(input$organismb),
            imported_genomic_data = imported_genomic_data.carepat
          )
      
      } else {
      transcript_data.carepat <-
        c(
          X5UTR = as.character(input$x5utr$datapath),
          X3UTR = as.character(input$x3utr$datapath),
          CDS = as.character(input$cds$datapath),
          Intron = as.character(input$intron$datapath)
        )
      genomic_data.carepat <-
        c(transcript_data.carepat, genomic_data.carepat)
      imported_genomic_data.carepat <-
        Wimtrap::importGenomicData(
          biomart = FALSE,
          genomic_data = genomic_data.carepat,
          tss = as.character(input$tss$datapath),
          tts = as.character(input$tts$datapath),
          promoter_length = as.numeric(input$promoter),
          proximal_length = as.numeric(input$proximal),
          downstream_length = as.numeric(input$downstream)
        )
      TFBSdata.carepat <-
        Wimtrap::getTFBSdata(
          pfm = as.character(input$pwmb$datapath),
          TFnames = TFnames.carepat,
          genome_sequence = as.character(input$genome$datapath),
          imported_genomic_data = imported_genomic_data.carepat
        )
      
      }
      ChIPpeaks.carepat <- as.character(input$filesChIP$datapath)
      names(ChIPpeaks.carepat) <- TFnames.carepat
      model.carepat <-
        WimtrapWeb:::buildTFBSmodel(TFBSdata = TFBSdata.carepat,
                                      ChIPpeaks = ChIPpeaks.carepat)
      shinybusy::remove_modal_spinner()
      return(model.carepat)
    })
    
    ww <- reactiveValues(data = NULL)
    observeEvent(input$build,
                 {
                   ww$data <- rnorm(100)
                 })

    
    output$roc <- renderPlot({
      if (is.null(ww$data))
        return()
      plotROC(Wimtrap_model())
    })
    
    output$downloadROC <- downloadHandler({
      filename = function() {
        paste("ROC_", Sys.Date(), "_", Sys.time(), ".png", sep = "")
      }
    },
    content = function(file) {
      png(file)
      print(plotROC(Wimtrap_model()))
      dev.off()
    })
    
    output$confusion <- renderPrint({
      if (is.null(ww$data))
        return()
      printConfusionMatrix(Wimtrap_model())
    })
    
    output$feature <- renderPlot({
      if (is.null(ww$data))
        return()
      plotFeatureImportance(Wimtrap_model())
    })
    output$downloadFeatureImportance <- downloadHandler({
      filename = function() {
        paste("FeatureImportance_",
              Sys.Date(),
              "_",
              Sys.time(),
              ".png",
              sep = "")
      }
    },
    content = function(file) {
      png(file)
      print(plotFeatureImportance(Wimtrap_model()))
      dev.off()
    })
    
    observe({
      if (!(is.null(w$data))){
        Wimtrap_model2 <- reactiveValues()
        isolate(Wimtrap_model2 <<- Wimtrap_model()$model)
        w$data <- NULL
      }
    })
    
    output$downloadModel <- downloadHandler(
      filename = function() {
        paste("TFBSmodel_", Sys.Date(), "_", Sys.time(), ".RData", sep = "")
      },
      content = function(file) {
        save(Wimtrap_model2, file = file)
      }
    )
    
    z <- reactiveValues(data = NULL)
    
    observeEvent(input$TFBSmodel,
                 {
                   z$data <- rnorm(100)
                 })
    
    PredFeatures <- reactive({
      if (is.null(z$data))
        return ("The predictive features are: ")
      predictiveModel <- get(load(file = input$TFBSmodel$datapath))
      predictiveFeatures <- predictiveModel$feature_names
      predictiveFeatures <-
        gsub(pattern = "_20bp",
             replacement = "",
             x = predictiveFeatures)
      predictiveFeatures <-
        gsub(pattern = "_400bp",
             replacement = "",
             x = predictiveFeatures)
      predictiveFeatures <-
        gsub(pattern = "_1000bp",
             replacement = "",
             x = predictiveFeatures)
      if(length(grep(pattern = "matchLogPval", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "matchLogPval", x = predictiveFeatures)]
      }
      if(length(grep(pattern = "Matches", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "Matches", x = predictiveFeatures)]
      }
      if (length(grep(pattern = "DistToClosestTSS", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
         predictiveFeatures[-grep(pattern = "DistToClosestTSS", x = predictiveFeatures)]
      }
      if (length(grep(pattern = "DistToClosestTTS", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "DistToClosestTTS", x = predictiveFeatures)]
      }
      predictiveFeatures <-
        predictiveFeatures[!duplicated(predictiveFeatures)]
      predictiveFeatures <-
        predictiveFeatures[order(predictiveFeatures)]
      return(paste0(
        "The predictive features are: " ,
        paste(predictiveFeatures, collapse = ", ")
      ))
    })
    
    output$PredictiveFeatures <- renderPrint({
      PredFeatures()
    })
    
    zc <- reactiveValues(data = NULL)
    
    observeEvent(input$TFBSmodelc,
                 {
                   zc$data <- rnorm(100)
                 })
    
    PredFeaturesc <- reactive({
      if (is.null(zc$data))
        return ("The predictive features are: ")
      predictiveModel <- get(load(file = input$TFBSmodelc$datapath))
      predictiveFeatures <- predictiveModel$feature_names
      predictiveFeatures <-
        gsub(pattern = "_20bp",
             replacement = "",
             x = predictiveFeatures)
      predictiveFeatures <-
        gsub(pattern = "_400bp",
             replacement = "",
             x = predictiveFeatures)
      predictiveFeatures <-
        gsub(pattern = "_1000bp",
             replacement = "",
             x = predictiveFeatures)
      if(length(grep(pattern = "matchLogPval", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "matchLogPval", x = predictiveFeatures)]
      }
      if(length(grep(pattern = "Matches", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "Matches", x = predictiveFeatures)]
      }
      if (length(grep(pattern = "DistToClosestTSS", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "DistToClosestTSS", x = predictiveFeatures)]
      }
      if (length(grep(pattern = "DistToClosestTTS", x = predictiveFeatures))==0){
        
      } else {
        predictiveFeatures <-
          predictiveFeatures[-grep(pattern = "DistToClosestTTS", x = predictiveFeatures)]
      }
      predictiveFeatures <-
        predictiveFeatures[!duplicated(predictiveFeatures)]
      predictiveFeatures <-
        predictiveFeatures[order(predictiveFeatures)]
      return(paste0(
        "The predictive features are: " ,
        paste(predictiveFeatures, collapse = ", ")
      ))
    })
    
    output$PredictiveFeaturesc <- renderPrint({
      PredFeaturesc()
    })
    
    zz <- reactiveValues(data = NULL)
    
    observeEvent(input$predict,
                 {
                 zz$data <- rnorm(100)
                 })
    
    CARE_predictions2 <- reactiveValues()
    
    observeEvent(input$predict,{
      show_modal_spinner(spin = "cube-grid",
                         color = "#428bca",
                         text = "Please wait 5-10 min...")
      
      if (input$checkc1) {
        message("Alexandra")
        if (is.null(input$organismb)){
          predictiveModel <- Wimtrap_model2
          bed_data <- reactiveValuesToList(input)
          bed_data <-
            bed_data[grep(pattern = "bed", x = names(bed_data))]
          bed_data <- unlist(bed_data)
          genomic_data.carepat <- as.character(input$filesData$datapath)
          names(genomic_data.carepat) <- bed_data
          imported_genomic_data2 <- Wimtrap::importGenomicData(
                                    biomart = FALSE,
                                    genomic_data = genomic_data.carepat,
                                    tss = as.character(input$tss$datapath),
                                    tts = as.character(input$tts$datapath),
                                    promoter_length = as.numeric(input$promoter),
                                    proximal_length = as.numeric(input$proximal),
                                    downstream_length = as.numeric(input$downstream)
                                  )
          genome2 <- input$genome$datapath
        } else {
          predictiveModel <- Wimtrap_model2
          bed_data <- reactiveValuesToList(input)
          bed_data <-
            bed_data[grep(pattern = "bed", x = names(bed_data))]
          bed_data <- unlist(bed_data)
          genomic_data.carepat <- as.character(input$filesData$datapath)
          names(genomic_data.carepat) <- bed_data
          imported_genomic_data2 <- tryCatch(
            Wimtrap::importGenomicData(organism = as.character(input$organismb),
                                       genomic_data = genomic_data.carepat,
                                       promoter_length = as.numeric(input$promoterb),
                                       downstream_length = as.numeric(input$downstreamb),
                                       proximal_length = as.numeric(input$proximalb)),
            error = function(e){return(NA)},
            finally = message("Interrogating Ensembl and Biomart..."))
          if(is.na(imported_genomic_data2)){
            imported_genomic_data2 <- 
              Wimtrap::importGenomicData(organism = as.character(input$organismb),
                                         genomic_data = genomic_data.carepat,
                                         promoter_length = as.numeric(input$promoterb),
                                         downstream_length = as.numeric(input$downstreamb),
                                         proximal_length = as.numeric(input$proximalb))
          }
        }
      } else { 
        if (input$checkc2){
          genomicFeatures <-
            unlist(strsplit(x = as.character(PredFeaturesc()), split = "The predictive features are: "))
          genomicFeatures <-
            unlist(strsplit(x = genomicFeatures, split = ", "))
          genomicData <- reactiveValuesToList(input)
          genomicData <-
            genomicData[paste0("bdc", seq(0, (length(
              genomicFeatures
            ) - 1)))]
          names(genomicData) <- genomicFeatures
          if (length((names(genomicData) %in% c("Promoter", "Downstream", "ProximalPromoter",
                                                "CDS", "X5UTR", "X3UTR", "Intron"))) > 0){
            genomicData <-
              genomicData[!(names(genomicData) %in% c("Promoter", "Downstream", "ProximalPromoter",
                                                      "CDS", "X5UTR", "X3UTR", "Intron"))]
          } else {
            
          }
          genomicFeatures <- names(genomicData)
          genomicData <- do.call(rbind, genomicData)
          genomicData <- genomicData$datapath
          names(genomicData) <- genomicFeatures
          imported_genomic_data2 <- tryCatch(
            Wimtrap::importGenomicData(organism = as.character(input$organismc),
                                       genomic_data = genomicData,
                                       promoter_length = as.numeric(input$promoterc),
                                       proximal_length = as.numeric(input$proximalc),
                                       downstream_length = as.numeric(input$downstreamc)),
            error = function(e){return(NA)},
            finally = message("Interrogating Ensembl and Biomart..."))
          if(is.na(imported_genomic_data2)){
            imported_genomic_data2 <- 
              Wimtrap::importGenomicData(organism = as.character(input$organismc),
                                         genomic_data = genomicData,
                                         promoter_length = as.numeric(input$promoterc),
                                         proximal_length = as.numeric(input$proximalc),
                                         downstream_length = as.numeric(input$downstreamc)
                                         )
          } else {}

      } else {
          genomicFeatures <-
            unlist(strsplit(x = as.character(PredFeatures()), split = "The predictive features are: "))
          genomicFeatures <-
            unlist(strsplit(x = genomicFeatures, split = ", "))
          genomicData <- reactiveValuesToList(input)
          genomicData <-
            genomicData[paste0("bdc", seq(0, (length(
              genomicFeatures
            ) - 1)))]
          names(genomicData) <- genomicFeatures
          promoterLength <- genomicData[["Promoter"]]
          downstreamLength <- genomicData[["Downstream"]]
          proximalLength <- genomicData[["ProximalPromoter"]]
          genomicData <-
            genomicData[!(names(genomicData) %in% c("Promoter", "Downstream", "ProximalPromoter"))]
          genomicFeatures <- names(genomicData)
          genomicData <- do.call(rbind, genomicData)
          genomicData <- genomicData$datapath
          names(genomicData) <- genomicFeatures
          imported_genomic_data2 <-
            Wimtrap::importGenomicData(
              genomic_data = genomicData,
              biomart = FALSE,
              promoter_length = promoterLength,
              downstream_length = downstreamLength,
              proximal_length = proximalLength,
              tss = input$tssc$datapath,
              tts = input$ttsc$datapath
            )
          predictiveModel <- get(load(file = input$TFBSmodel$datapath))
          genome2 <- input$genomec$datapath
        }
      }
      
      TFnames2 <-
        gsub(
          pattern = "\\s",
          replacement = "",
          x = as.character(input$motifsc)
        )
      TFnames2 <- unlist(strsplit(x = TFnames2, split = ","))
      
      if (input$checkc2){
        TFBSdata2 <- Wimtrap::getTFBSdata(
          pfm = input$pwmc$datapath,
          TFnames = TFnames2,
          organism = input$organismc,
          imported_genomic_data = imported_genomic_data2
        )
        
      } else {
        if (input$checkc1 & input$biomart){
          TFBSdata2 <- Wimtrap::getTFBSdata(
            pfm = input$pwmc$datapath,
            TFnames = TFnames2,
            organism = as.character(input$organismb),
            imported_genomic_data = imported_genomic_data2
          )
          
        } else { 
          TFBSdata2 <- Wimtrap::getTFBSdata(
            pfm = input$pwmc$datapath,
            TFnames = TFnames2,
            genome_sequence = genome2,
            imported_genomic_data = imported_genomic_data2
          )
        }
      }
      
      results2 <- Wimtrap::predictTFBS(
        TFBSmodel = predictiveModel,
        TFBSdata = TFBSdata2,
        studiedTFs = TFnames2,
        score_threshold = input$score_thresholdc
      )
      CARE_predictions2$data <- results2
      remove_modal_spinner()
    })
    
    output$CAREpredictions2 <- renderDataTable({
      CARE_predictions2$data
    })
    
    Target_predictions2 <- reactive({
      if (is.null(zz$data))
        return(data.table::data.table(NULL))
      show_modal_spinner(spin = "cube-grid",
                         color = "#428bca",
                         text = "Please wait...")
      transcripts2 <-
        split(x = CARE_predictions2$data, f = CARE_predictions2$data$TF)
      genes2 <- lapply(transcripts2, function(pred) {
        targets <-
          unlist(strsplit(x = as.character(pred$transcript), split = "[.]"))[seq(1, 2 *
                                                                                   nrow(pred), 2)]
        targets <- targets[!duplicated(targets)]
        targets <- data.table::data.table(Gene = targets,
                                          TF = pred$TF[1])
        return(targets)
      })
      genes2 <- do.call(rbind, genes2)
      remove_modal_spinner()
      return(genes2)
    })
    
    output$TargetPredictions2 <- renderDataTable({
      Target_predictions2()
    })
    
    output$downloadData3 <- downloadHandler(
      filename = function() {
        paste("TFBS_", Sys.Date(), "_", Sys.time(), ".tsv", sep = "")
      },
      content = function(file) {
        write.table(CARE_predictions2$data,
                    file,
                    sep = "\t",
                    row.names = FALSE)
      }
    )
    
    output$downloadData4 <- downloadHandler(
      filename = function() {
        paste("Target_", Sys.Date(), "_", Sys.time(), ".tsv", sep = "")
      },
      content = function(file) {
        write.table(Target_predictions2(),
                    file,
                    sep = "\t",
                    row.names = FALSE)
      }
    )
    
    })
  }

## To be copied in the UI
# mod_AdaptativeSelection_ui("AdaptativeSelection_ui_1")

## To be copied in the server
# mod_AdaptativeSelection_server("AdaptativeSelection_ui_1")
