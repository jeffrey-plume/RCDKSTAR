options(shiny.maxRequestSize=500*1024^2) 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx12288m"))
library(shiny)
library(rpubchem)
library(Rcpi)
library(caret)
library(doParallel)
library(ChemmineR)
library(rcdk)
library(tictoc)
library(gbm)
library(ChemmineR)
library(caret)
library(tidyverse)
library(org.Hs.eg.db)

descriptors <- map(readxl::excel_sheets("Descriptors.xls")[4:7], ~ readxl::read_excel("Descriptors.xls",    sheet = .x)) %>%
  setNames(str_replace(readxl::excel_sheets("Descriptors.xls")[4:7], 'Fingerprinter', ""))

descriptors[['PubChem']] <- ChemmineR::pubchemFPencoding %>%
  dplyr::mutate(SMARTS = make.unique(as.character(Bit_Substructure)))

ui <- fluidPage(
  
  titlePanel("RCDKSTAR"),
  shinybusy::add_busy_spinner(spin = "fading-circle"),
  shinyjs::useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        inputId = 'aid',
        label = 'Target', 
        choices = ""),
      
      selectizeInput(
        inputId = 'fingerprints_in',
        label = 'Fingerprinting Method',
        choices = names(descriptors),
        selected = 'PubChem'),
      
      selectInput(inputId = 'properties_in',
                  label = 'Chemical Properties',
                  choices = ls("package:Rcpi", pattern = 'extractDrug')[!ls("package:Rcpi", pattern = 'extractDrug') %in% c(str_subset(ls("package:Rcpi"), 'Complete'), str_sub(str_subset(ls("package:Rcpi"), 'Complete'), 1, -9))],
                  selected = "",
                  multiple=TRUE,
                  size = 10,
                  selectize=FALSE),
      
      actionButton(
        inputId = 'load_properties', 
        label = 'Load Properties', 
        style = 'background-color:blue; font-weight:bold; color:white;'),
    
      HTML('<br>'),
      HTML('<br>'),
      
       checkboxGroupInput(
         inputId = 'preprocess',
         label = 'Pre-Process Dataset',
         inline = T,
         choices = c("BoxCox", "YeoJohnson", "expoTrans", "center", "scale",
                     "range", "knnImpute", "bagImpute", "medianImpute", 
                     "pca", "ica", "spatialSign", "corr", "zv", "nzv", "conditionalX" ),
         selected = c('zv', 'medianImpute')),
      
      selectInput(
        inputId = 'algorithm',
        label = 'Learning Method',
        choices = c('glmnet',
                    'kknn',
                    'naive_bayes',
                    'nnet',
                    'kernelpls',
                    'rf',
                    'gbm',
                    'svmPoly',
                    'svmRadial',
                    'svmRadialWeights'),
        selected = 'rf'),
      
      sliderInput(
        inputId = 'split',
        label = 'Percentage of Dataset for Training',
        min = 25,
        value = 75, 
        max = 100,
        step = 5),
      
      shinyjs::hidden(
        div(id = "advanced",
            
            HTML('<br>'),
            checkboxGroupInput(
              inputId = 'sampling',
              label = 'Sampling',
              choices = c("UpSample",
                          "DownSample")),
            
            sliderInput(
              inputId = 'repeats',
              label = 'Number of k-nodes',
              min = 1,
              value = 3, 
              max = 10,
              step = 1),
            
            sliderInput(
              inputId = 'number',
              label = 'Resample Iterations',
              min = 1,
              value = 5, 
              max = 10,
              step = 1),
            
            sliderInput(
              inputId = 'tune_length',
              label = 'Set Tune Length',
              min = 1,
              value = 5, 
              max = 10,
              step = 1)
        )
      ),
      
      actionButton('train', 'Train Model', 
                   style = 'background-color:blue; font-weight:bold; color:white'),
      HTML('<br>'),
      
      a(id = "toggleAdvanced", "Advanced Options"),
      
      HTML('<br>'),
      HTML('<br>'),
      
      fileInput(
        inputId = 'unknowns_in', 
        label = 'Upload Unknowns in SDF Format',
        placeholder = "repurposing.sdf"),
      
      actionButton('unknowns_pred', 'Test Unknowns', 
                   style = 'background-color:blue; font-weight:bold; color:white'),
      
      actionButton(
        inputId = 'unknown_help', 
        label = '',
        icon = icon("info-circle"),
        style = 'background-color:#f5f5f5; border:none; font-weight:bold; color:blue; font-size: 18px; margin: 0px; padding: 2px;'),
      HTML('<br>')
      
    ),
    
    mainPanel(
      
      uiOutput('genename'),
      
      div(id='prop1',
          textOutput('summary_small'),
          a(id = "toggleSmall", "Show More")
      ),
      
      shinyjs::hidden(
        div(id = "prop2",
            textOutput('summary_big'),
            a(id = "toggleBig", "Show Less")
        )
      ),
      
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'overview'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      uiOutput('error'),
      
      HTML('<br>'),
      
      div(DT::dataTableOutput(
        outputId = 'prop'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'corr'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'prop_corr'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'summary'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      plotOutput(outputId = 'roc'),
      HTML('<br>'),
      
      plotOutput(outputId = 'imp_out'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'imp_table'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'prediction'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>'),
      
      div(DT::dataTableOutput(outputId = 'unknowns_out'), style = 'overflow-x: auto; white-space: nowrap; text-overflow: ellipsis;'),
      HTML('<br>')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  rv <- reactiveValues()
  
  updateSelectizeInput(session, 'aid', 'Target', server = T, 
                       selected = 'BRAF',
                       choices = sort(AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db), columns = "SYMBOL")$SYMBOL))
  
  output$genename <- renderUI(h2(input$aid))
  
  # Retreives Gene descriptions from the NCBI
  summaryRV <- reactive({
    if(is.null(input$aid) | length(input$aid)==0 | input$aid == "") return()

    rentrez::entrez_summary('gene', AnnotationDbi::select(org.Hs.eg.db, keys = input$aid, columns = "ENTREZID", keytype="SYMBOL")$ENTREZID)$summary
  })
  
  # The code below expands retracts the gene description summary
  output$summary_small <- renderText({
    
    if(is.null(summaryRV())|length(summaryRV())==0) return()
    
    str_extract(summaryRV(), '.*?[a-z0-9][.?!](?= )')
  })
  
  output$summary_big <- renderText({
    if(is.null(summaryRV())|length(summaryRV())==0) return()
    
    summaryRV()
  })
  
  shinyjs::onclick("toggleSmall",{
    shinyjs::hide("prop1")
    shinyjs::show('prop2')
  })
  
  shinyjs::onclick("toggleBig",{
    shinyjs::show("prop1")
    shinyjs::hide('prop2')
  })
  
  # Expand hidden options in the sidebar
  shinyjs::onclick("toggleAdvanced",
                   shinyjs::toggle(id = "advanced", anim = TRUE))    
  
  # Aggreagate bioassay results from protein target
  aid <- reactive({
    
    if(is.null(input$aid)|length(input$aid)==0 | input$aid == "") return()
    
    a <- map(AnnotationDbi::select(org.Hs.eg.db, keys = input$aid, columns = "UNIPROT", keytype="SYMBOL")$UNIPROT,
             ~read_csv(RCurl::getURL(paste0("https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22protacxn%22:%22notnull%22},{%22cid%22:%22notnull%22},{%22repacxn%22:%22", 
                                            .x,
                                            "%22}]},%22order%22:[%22activity,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22{PROTACXN_", 
                                            .x,
                                            "}_bioactivity_protein%22}"))))
    
    bind_rows(a[sapply(a, nrow) > 0])
  })
  
  output$overview <- DT::renderDataTable(server = TRUE, {
    Sys.sleep(1)
    
    if(is.null(aid())|length(aid())==0) return()
    
    DT::datatable(
      mutate_all(aid(), as.factor),
      extensions = c('Buttons', 'AutoFill', 'KeyTable'), 
      filter = 'top', 
      editable = TRUE,
      options = list(
        searching = TRUE,
        autoWidth = TRUE,
        dom = 'Bfrtlp',
        searchCols = list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, list(search = '["Confirmatory"]')), 
        buttons = c('copy', 'csv', 'excel'),
        keys = TRUE
      ),
      rownames = FALSE)
  })
  
  # Warnings/hints regarding filtered AID data
  
  error <- reactive({
    if(is.null(aid())|length(aid())==0) return()
    
    error <- list()
    
    if(any(summarise(group_by(aid()[input$overview_rows_all, ], cid), n=n_distinct(activity))$n > 1)){
      
      error[['act']] <- '<div style="color: red; font-weight:bold">Conflicting Activitiy Reports.  Consider filtering data.</div>'
    }
    
    aid()$cid[input$overview_rows_all][duplicated(aid()$cid[input$overview_rows_all])]
    if(any(duplicated(aid()$cid[input$overview_rows_all]))){
      error[['dup']] <- '<div style="color: GoldenRod; font-weight:bold">Duplicates Found. First value will be used.</div>'
    }
    
    if(n_distinct(aid()$activity[input$overview_rows_all])==1){
      error[['out']] <- '<div style="color: red; font-weight:bold">Multiple Activity Outcomes Required.</div>'
    }
    
    return(error)
  })
  
  output$error <- renderUI({
    if(is.null(error())|length(error())==0) return()
    
    HTML(str_c(error(), sep = "", collapse = ''))
  })
  
  url <- reactive(paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/', 
                         map(split(aid()$cid[input$overview_rows_all], ceiling(seq(aid()$cid[input$overview_rows_all])/400)), ~str_c(.x, collapse = ',')), '/SDF'))
  
  observeEvent(
    input$load_properties, {
      
      gc(reset = T, full = T)
      
      
      tictoc::tic()
      mol <-  unlist(map(isolate(url()), purrr::insistently(~load.molecules(normalizePath(.x), typing = TRUE))))
      names(mol) <- aid()$cid[input$overview_rows_all]

      properties <- eval(parse(text = paste0('extractDrug', input$fingerprints_in, 'Complete')))(mol) %>%
        as.data.frame() %>%
        setNames(descriptors[[input$fingerprints_in]]$SMARTS)
      
      
      if(!is.null(input$properties_in) & length(input$properties_in)>0){
        
        desc <- do.call(cbind, map(input$properties_in, ~eval(parse(text = str_c('Rcpi::', .x, '(mol)'))))) %>%
          data.frame() 
        
        desc <- bind_rows(map(desc, ~replace_na(.x, median(.x, na.rm = T))))
        
        properties <- cbind(properties, desc)
        
      }
      
      properties <- properties %>%
        mutate(cid = as.numeric(names(mol))) 
      
      rv$properties <- inner_join(properties, isolate(aid()[input$overview_rows_all, c('activity', 'cid')]))  %>%
        group_by(cid) %>%
        slice_head(n=1) %>%
        ungroup() %>%
        tibble::column_to_rownames('cid')
      
     
      tictoc::toc()
      
    })
  
  # Data is indexed into training and testing sets
  index <- reactive(createDataPartition(
    y = rv$properties$activity,
    p = input$split/100,
    list = FALSE))
  
  # The pre-processioning model is created using only data in the training set
  prepro <- reactive({
      caret::preProcess(rv$properties[index(), ], input$preprocess)
     
  })
    
  # Pre-Processed training data is displayed
  
    output$prop <- DT::renderDataTable(server = TRUE, {
      if(is.null(rv$properties) | length(rv$properties)==0) return()
      
      DT::datatable(predict(prepro(), rv$properties)[isolate(index()), ],
                    extensions = 'Buttons', 
                    filter = 'top',
                    options = list(
                      searching = TRUE,
                      autoWidth = TRUE,
                      dom = 'Bfrtlp',
                      buttons = c('copy', 'csv', 'excel')), 
                    rownames = TRUE)
    })
  
    
    observeEvent(
      input$train, {
        
        properties <- predict(prepro(), rv$properties)
        
        tictoc::tic()
        
        cl <- parallel::makeCluster(parallel::detectCores()-1)
        doParallel::registerDoParallel(cl)
        
        training <- predict(prepro(), rv$properties[index(), ])
        
        if(is.null(input$sampling)|length(input$sampling)==0){
          
          training <- training
          
        } else if('UpSample' %in% input$sampling){
          training <- caret::downSample(x = training %>%
                                      dplyr::select(-activity),
                                    y = as.factor(training$activity), list = FALSE,
                                    yname = 'activity')
        }
        
        if('DownSample' %in% input$sampling){
          training <- caret::downSample(x = training %>%
                                      dplyr::select(-activity),
                                    y = as.factor(training$activity), list = FALSE,
                                    yname = 'activity')
        }
        
        ctrl <- trainControl(
          method = 'repeatedcv', 
          repeats = input$repeats,
          number = input$number,
          allowParallel = TRUE,
          classProbs = TRUE, 
          savePredictions = 'all',
          verboseIter = TRUE)
        
        fit <- caret::train(as.factor(activity) ~ ., data = training, 
                            method = input$algorithm,
                            trControl = ctrl, 
                            metric = 'Accuracy',
                            tuneLength = input$tune_length,
                            na.action=na.exclude)
        
        parallel::stopCluster(cl)
        
        tictoc::toc()
        
        importance <- caret::varImp(fit)

        output$summary <- DT::renderDataTable(server = TRUE, {
          if(is.null(fit$results)) return()
          
          DT::datatable(fit$results,
                        extensions = 'Buttons', filter = 'top',
                        options = list(
                          searching = TRUE,
                         # autoWidth = TRUE,
                          dom = 'Bfrtlp',
                          buttons = c('copy', 'csv', 'excel')
                        ), rownames = T)
        })
        
        
        output$roc <- renderPlot({
          if (is.null(fit)) return()
          plot(fit, type='b')
        })
        
        output$imp_out <- renderPlot({
          if(is.null(fit)) return()
          plot(importance, top = 20)
        })

        # proprocessing is applied to testing set
        testing <- predict(prepro(), rv$properties[-index(), ])
        
        pred <- predict(fit, testing, na.action = na.omit, type = 'prob') %>%
          setNames(paste('Probs', colnames(.))) %>%
          tibble::rownames_to_column('cid') %>%
          mutate(cid = as.numeric(cid)) %>%
          inner_join(isolate(aid())) 
        

        output$prediction <- DT::renderDataTable(server = TRUE, {
          if(is.null(pred)) return()
          DT::datatable(pred,
                        extensions = 'Buttons', 
                        options = list(
                          searching = TRUE,
                          autoWidth = TRUE,
                          dom = 'Bfrtlp',
                          buttons = c('copy', 'csv', 'excel'),
                          pageLength = 10), 
                        rownames = T)
        })
        
        rv$fit <- fit
        
      })
    
    
    observeEvent(
      input$unknowns_pred, {
        
        # if no SDF is uploaded, drugs the Drug Repurposing Hub from Broad Institute will be used
        if(is.null(input$unknowns_in)){
          unknown_sdf <- "repurposing.sdf"
        }
        else{
          unknown_sdf <- input$unknowns_in$datapath
        }
        
        unknown_mol <- load.molecules(unknown_sdf)
        names(unknown_mol) <- sdfid(read.SDFset(unknown_sdf))
        
        
        unknown_properties <- eval(parse(text = paste0('extractDrug', input$fingerprints_in, 'Complete')))(unknown_mol) %>%
          as.data.frame() %>%
          setNames(descriptors[[input$fingerprints_in]]$SMARTS) 
        
        if(!is.null(input$properties_in) & length(input$properties_in)>0){
          
          unknown_desc <- do.call(cbind, map(input$properties_in, ~eval(parse(text = str_c('Rcpi::', .x, '(unknown_mol)'))))) %>%
            data.frame() 
          
          unknown_desc <- bind_rows(map(desc, ~replace_na(.x, median(.x, na.rm = T))))
          
          unknown_properties <- cbind(unknown_properties, unknown_desc)
          
        }
        
        unknown_properties <- predict(prepro(), unknown_properties)
        
        fit <- isolate(rv$fit)
        
        unknown_pred <- predict(fit, unknown_properties, na.action = na.pass, type = 'prob') %>%
          setNames(paste('Probs', colnames(.)))
          
        # Drug Identifiers are pulled from the datablock of the SDF file
        unknown_pred <- bind_cols(bind_rows(datablock(read.SDFset(unknown_sdf))), unknown_pred)
        
        
        
        output$unknowns_out <- DT::renderDataTable(server = TRUE, {
          if (is.null(unknown_pred)) return()
          DT::datatable(
            unknown_pred, 
            extensions = 'Buttons', 
            options = list(
              searching = TRUE,
              autoWidth = TRUE,
              dom = 'Bfrtlp',
              buttons = c('copy', 'csv', 'excel'),
              pageLength = 10),
            rownames = T)
        })
        
      })
    
  
}

# Run the application 
shinyApp(ui = ui, server = server)
