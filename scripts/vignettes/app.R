library(shiny)

ui <- fluidPage(
  fluidRow(
    column(12,
           h1("A JOOBY HOOBY COLLEGE PRESENTATION"),
    )
  ),
    sidebarLayout(
      sidebarPanel(
          cellArgs = list(style='white-space: normal;'),
          
          submitButton("Click Here to Train Models"),
          br(),
          checkboxGroupInput("classifier_choice", label="Choice of Classifier",
                             choices=c("SVM Radial Weights"="svm", "k-Nearest Neighbors"="knn", "Logistic (Ridge) Regression"="lrr", 
                                       "Linear Discriminant Analysis"="lda", "Model Averaged Neural Network"="mann","Parallel Random Forests"="prf"),
                             selected=c("svm","knn","lda")),
          radioButtons("time_choice", label="Measure Time Interval by:",
                       choices=c("Time of Day (TOD)"="tod","Time from Recording (TFR)"="tfr"), 
                       selected="tfr"),
          verticalLayout( 
           sliderInput("tod_slice", label="TOD Time Interval Between 6 PM and 12 PM (+1)", timeFormat = "%T%p",
                       step = 30,
                       ticks=FALSE,
                       min = as.POSIXct("2020-05-05 18:00:00"),
                       max = as.POSIXct("2020-05-06 12:00:00"), 
                       value = c(as.POSIXct("2020-05-06 0:00:00"),as.POSIXct("2020-05-06 6:00:00"))),
           sliderInput("tfr_slide", label="TFR in hours",  
                        min = 0, max = 8, step=1, value = c(2,6))
            ),
         
         numericInput("dim_red",label="# of Dimensions to Reduce to",
                      min=1,max=3,step=0.5,value=3),

         checkboxGroupInput("mf_choice",label="Choice of Motion Features",
                            choices=c("BP"="bp","FE"="fe","FP"="fp","MF"="mf","SM"="sm","WV"="wv"),
                            selected="bp",
                            inline = TRUE),
         checkboxGroupInput("clinvar_choice", label="Choice of Clinical Variables",
                            choices=c("Age"="age","Sex"="sex","APACHE"="apache","GCS_Enrollment"="gcs_en","Diagnosis"="diag"),
                            selected=c("apache","sex","diag"),
                            inline=TRUE),
         checkboxGroupInput("sensor_loc", label="Choice of Sensor Locations",
                            choices=c("Left Elbow"="left_el","Left Wrist"="left_wr","Left Ankle"="left_ank",
                                      "Right Elbow"="right_el","Right Wrist"="right_wr","Right Ankle"="right_ank"),
                            selected=c("left_el","left_wr","left_ank"),
                            inline=TRUE)
      ),
      mainPanel(
             tabsetPanel(
               tabPanel("ROC", "contents"), 
               tabPanel("Precision-Recall", "contents"), 
               tabPanel("Calibration", "contents")
        )
        )
      )
    )

server <- function(input, output) {
  
  
  
}


shinyApp(ui, server)