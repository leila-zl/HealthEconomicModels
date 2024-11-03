library(shiny)
library(survival)
library(survminer)
library(ggpubr)
library(ggthemes)
library(flexsurv)
library(asaur)
library(ggplot2)
library(ggsurvfit)


# Load the dataset
pancreatic2 <- asaur::pancreatic2

# Define UI
ui <- fluidPage(
  titlePanel("Survival Analysis with Different Fitted Models"),
  sidebarLayout(
    sidebarPanel(
      selectInput("model_type", "Choose Model Type:",
                  choices = list("Weibull" = "weibull", 
                                 "Exponential" = "exp", 
                                 "Gamma" = "gamma", 
                                 "Generalized Gamma" = "gengamma", 
                                 "Generalized F" = "genf", 
                                 "Log-Normal" = "lnorm", 
                                 "Gompertz" = "gompertz", 
                                 "Log-Logistic" = "llogis")),
      actionButton("fit_model", "Fit Model")
    ),
    mainPanel(
      plotOutput("combined_plot"),
      verbatimTextOutput("model_summary")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  fit_results <- eventReactive(input$fit_model, {
    # Create the survival object
    os_surv <- Surv(time = pancreatic2$os, event = pancreatic2$status)
    
    # Fit Kaplan-Meier model
    km_fit <- survfit(os_surv ~ 1, data = pancreatic2)
    
    # Fit selected model based on input
    selected_model <- input$model_type
    fit <- flexsurvreg(os_surv ~ 1, data = pancreatic2, dist = selected_model)
    
    # Predict survival probabilities for a range of time points
    time_seq <- seq(0, max(pancreatic2$os), by = 1)
    fitted_summary <- summary(fit, type = "survival", t = time_seq)
    
    # Convert to a dataframe
    fitted_data <- do.call(rbind, lapply(fitted_summary, function(x) {
      data.frame(time = x$time, surv = x$est)
    }))
    
    # Ensure valid data only
    fitted_data <- fitted_data[complete.cases(fitted_data), ]
    
    # Create combined plot
    km_plot <- ggsurvfit(km_fit) +
      labs(x = "Time (days)", y = "Survival Probability", title = "Kaplan-Meier and Fitted Model") +
      add_confidence_interval() +
      theme_minimal() +
      geom_line(data = fitted_data, aes(x = time, y = surv), color = "red", size = 1)
    
    list(combined_plot = km_plot, summary = summary(fit))
  })
  
  output$combined_plot <- renderPlot({
    fit_results()$combined_plot
  })
  
  output$model_summary <- renderPrint({
    fit_results()$summary
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
