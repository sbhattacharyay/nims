#### Classification Shiny App ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay and Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(devtools)
if (!require(lolR))
  install_github('neurodata/lol',
                 build_vignettes = TRUE,
                 force = TRUE)
library(lolR)
library(R.matlab)
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(readxl)
library(plotly)
library(naniar)
library(MASS)
library(glmnet)
library(caret)
library(kernlab)
library(rlist)
library(nnet)
library(e1071)
library(randomForest)
library(foreach)
library(gridExtra)
library(ggplotify)
library(precrec)
library(cvAUC)
library(grid)


setwd("..")

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source("./functions/cross_val_splits.R")
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')
source('./functions/classification_function_shiny_dis.R')
source('./functions/classification_function_shiny_12mo.R')
source('./functions/train_caret_models.R')
source('./functions/predict_caret_models.R')
source('./functions/get_auc_info.R')
source('./functions/get_auc_plots.R')
source('./functions/get_precrec_info.R')
source('./functions/get_precrec_plots.R')

ui <- fluidPage(
  fluidRow(
    column(12,
           h1("ICU Accelerometry Classification with 5-fold Cross-Validation")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      cellArgs = list(style='white-space: normal;'),
      
      checkboxGroupInput("classifier_choice", label="Choice of Classifier",
                         choices = c(
                           "SVM Radial Weights" = "svmRadialWeights",
                           "k-Nearest Neighbors" = "knn",
                           "Logistic Regression (Elastic Net)" = "glmnet",
                           "Linear Discriminant Analysis" = "lda",
                           "Model Averaged Neural Network" = "avNNet",
                           "Parallel Random Forests" = "parRF"
                         ), 
                         selected=c("svmRadialWeights","knn","lda","glmnet")),
      radioButtons("time_choice", label="Measure Time Interval by:",
                   choices=c("Time of Day (TOD)"="tod","Time from Recording (TFR)"="tfr"), 
                   selected="tfr"),
      uiOutput("slider_type"),
      numericInput("r",label="# of Dimensions to Reduce to",
                   min=1,max=10,step=1,value=3),
      checkboxGroupInput("mf_choice",label="Choice of Motion Features",
                         choices = c(
                           "Bandpower (0.3 - 3.5 Hz)" = "band_powerFeats",
                           "Frequency-Domain Entropy" = "freq_entropyFeats",
                           "High-Frequency, Low-Frequency Time-Domain Median Pairs" = "freq_pairsFeats",
                           "Median Frequency" = "med_freqFeats",
                           "Signal Magnitude Area" = "smaFeats",
                           "db5 Wavelet (L2-6) Detail Coefficients" = "waveletsFeats"
                         ),
                         selected="band_powerFeats",
                         inline = FALSE),
      checkboxGroupInput("clinicalVars", label="Choice of Clinical Variables",
                         choices = c(
                           "Age" = "Age",
                           "Sex" = "Sex",
                           "APACHE" = "APACHE",
                           "GCS at Enrollment" = "GCS_en",
                           "Diagnosis" = "diag"
                         ),
                         selected = c("APACHE", "Sex", "diag"), 
                         inline=TRUE),
      checkboxGroupInput("sensor_loc", label="Choice of Sensor Locations",
                         choices = c(
                           "Left Ankle" = "left_ank",
                           "Left Elbow" = "left_el",
                           "Left Wrist" = "left_wr",
                           "Right Ankle" = "right_ank",
                           "Right Elbow" = "right_el",
                           "Right Wrist" = "right_wr"
                         ),
                         selected=c("left_wr","right_wr"),
                         inline=TRUE),
      br(),
      actionButton("button","Click Here to Train Models")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("ROC Discharge", 
                 plotOutput("plot_roc_dis", height = "1000px")
                 ), 
        tabPanel("ROC 12mo", 
                 plotOutput("plot_roc_12mo", height = "1000px")
        ),         
        tabPanel("Precision-Recall Discharge", 
                 plotOutput("plot_precrec_dis", height = "1000px")
                 ), 
        tabPanel("Precision-Recall 12mo", 
                 plotOutput("plot_precrec_12mo", height = "1000px")
        ),         
        tabPanel("Calibration", "contents")
      )
    )
  )
)

server <- function(input, output) {
  
  output$slider_type <- renderUI({
    if (input$time_choice == 'tod') {
      sliderInput("time_slide", label="TOD Time Interval Between 6 PM and 12 PM (+1)", timeFormat = "%I:%M %p",
                  step = 1800,
                  min = as.POSIXct("2020-05-05 18:00:00"),
                  max = as.POSIXct("2020-05-06 12:00:00"), 
                  value = c(as.POSIXct("2020-05-05 18:00:00"),as.POSIXct("2020-05-06 12:00:00")))
    } else {
      sliderInput("time_slide", label="TFR in hours",  
                  min = 0, max = 8, step=0.5, value = c(0,8))
    }
  })
  
  preds_dis <- eventReactive(input$button,{
    classification_function_shiny_dis(input$time_choice,input$time_slide,input$classifier_choice,input$r,
                                      input$mf_choice,input$clinicalVars,input$sensor_loc)
  })
  output$plot_roc_dis <- renderPlot({
    get_auc <- get_auc_info(preds_dis())
    get_plots <- get_auc_plots(get_auc)    
    do.call(grid.arrange, c(unlist(get_plots, recursive = F), nrow=length(get_auc[[1]])))
  })
  output$plot_precrec_dis <- renderPlot({
    get_precrec <- get_precrec_info(preds_dis())
    get_plots <- get_precrec_plots(get_precrec)    
    do.call(grid.arrange, c(unlist(get_plots, recursive = F), nrow=length(get_precrec[[1]])))
  })
  
  preds_12mo <- eventReactive(input$button,{
     classification_function_shiny_12mo(input$time_choice,input$time_slide,input$classifier_choice,input$r,
                                       input$mf_choice,input$clinicalVars,input$sensor_loc)
  })
  output$plot_roc_12mo <- renderPlot({
    get_auc <- get_auc_info(preds_12mo())
    get_plots <- get_auc_plots(get_auc)    
    do.call(grid.arrange, c(unlist(get_plots, recursive = F), nrow=length(get_auc[[1]])))
  })
  output$plot_precrec_12mo <- renderPlot({
    get_precrec <- get_precrec_info(preds_12mo())
    get_plots <- get_precrec_plots(get_precrec)    
    do.call(grid.arrange, c(unlist(get_plots, recursive = F), nrow=length(get_precrec[[1]])))
  })
}  


shinyApp(ui, server)