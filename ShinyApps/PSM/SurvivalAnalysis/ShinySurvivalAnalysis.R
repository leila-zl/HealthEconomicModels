library(shiny)
library(survival)
library(survminer)
library(ggpubr)
library(ggthemes)
library(flexsurv)
library(asaur)
library(ggplot2)
library(ggsurvfit)

# Load the publicly available dataset for pancreatic cancer
pancreatic2 <- asaur::pancreatic2

# Define UI
ui <- fluidPage(
  titlePanel("Survival Analysis for Pancreatic Cancer"),
  
  mainPanel(
    # set which plots to display, each of them in a separate tab
    tabsetPanel(
      tabPanel("Progression-Free Survival (PFS)", plotOutput("pfs_plot")),
      tabPanel("Overall Survival (OS)", plotOutput("os_plot")),
      tabPanel("Stage Comparison (PFS)", plotOutput("pfs_stage_plot")),
      tabPanel("Stage Comparison (OS)", plotOutput("os_stage_plot"))
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive function to perform survival analysis
  analysis_results <- reactive({
    
    # Create survival objects for PFS and OS
    #  define time to event, and event status (0 censored and 1 event happened)
    pfs_surv <- Surv(time = pancreatic2$pfs, event = pancreatic2$status) # PFS survival object
    os_surv <- Surv(time = pancreatic2$os, event = pancreatic2$status)   # OS survival object
    
    # Fit Kaplan-Meier models for PFS and OS
    # use survival object, and still need to refer back to the original data
    # ~ 1 means there is no subgroup analysis
    pfs_fit <- survfit(pfs_surv ~ 1, data = pancreatic2) # PFS Kaplan-Meier fit
    os_fit <- survfit(os_surv ~ 1, data = pancreatic2)   # OS Kaplan-Meier fit
    
    # Generate Kaplan-Meier plots
    pfs_plot <- ggsurvfit(pfs_fit) +
      add_confidence_interval()+                     # PFS plot with confidence intervals
      labs(x = "Time (days)", y = "PFS Probability")
    
    os_plot <- ggsurvfit(os_fit) +
      add_confidence_interval() +                       # OS plot with confidence intervals
      labs(x = "Time (days)", y = "OS Probability")
    
    # Fit Kaplan-Meier models by cancer stage
    # after the ~ symbol the variable that is used for the subgroup analysis is specified
    pfs_fit_stage <- survfit(pfs_surv ~ stage, data = pancreatic2) # PFS by stage
    os_fit_stage <- survfit(os_surv ~ stage, data = pancreatic2)   # OS by stage
    
    # Generate stage comparison plots
    pfs_stage_plot <- ggsurvfit(pfs_fit_stage) +
      add_confidence_interval() +         # PFS by stage plot
      labs(x = "Time (days)", y = "PFS Probability by Stage")
    
    os_stage_plot <- ggsurvfit(os_fit_stage) +
      add_confidence_interval() +           # OS by stage plot
      labs(x = "Time (days)", y = "OS Probability by Stage")
    
    # Return all results as a list
    list(
      pfs_plot = pfs_plot,
      os_plot = os_plot,
      pfs_stage_plot = pfs_stage_plot,
      os_stage_plot = os_stage_plot
    )
  })
  
  # Render PFS and OS plots
  output$pfs_plot <- renderPlot({ analysis_results()$pfs_plot })
  output$os_plot <- renderPlot({ analysis_results()$os_plot })
  
  # Render stage comparison plots
  output$pfs_stage_plot <- renderPlot({ analysis_results()$pfs_stage_plot })
  output$os_stage_plot <- renderPlot({ analysis_results()$os_stage_plot })
}

# Run the application
shinyApp(ui = ui, server = server)
