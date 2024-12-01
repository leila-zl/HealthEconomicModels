# Unpacking Diagnostic Metrics: What Matters Most
#
# This Shiny application enables users to explore the relationships between
# diagnostic test metrics - sensitivity, specificity, and prevalence — and their
# impact on Positive Predictive Value (PPV) and Negative Predictive Value (NPV).
# Users can adjust these parameters to visualize how changes affect diagnostic
# outcomes, thereby gaining a deeper understanding of the importance of each metric.

# Author: Leila Zakka
# Date: 01/12/2024


library(shiny)

# Define UI
ui <- fluidPage(
  titlePanel("Unpacking Diagnostic Metrics: What Matters Most"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("prevalence", "Select Disease Prevalence (%)", min = 0.1, max = 50, value = 5, step = 0.1),
      sliderInput("sensitivity", "Select Sensitivity (%)", min = 50, max = 100, value = 90, step = 1),
      sliderInput("specificity", "Select Specificity (%)", min = 50, max = 100, value = 95, step = 1),
      numericInput("population", "Enter Population Size", value = 1000, min = 1)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Text Summary", verbatimTextOutput("summary")),
        tabPanel("Bar Chart 1", plotOutput("barPlot1")),
        tabPanel("Bar Chart 2", plotOutput("barPlot")),
        tabPanel("PPV and Prevalence", plotOutput("lineGraph"))
      )
    )
  )
)

# Define Server Logic
server <- function(input, output) {
  # Reactive expressions to calculate metrics
  reactive_metrics <- reactive({
    prevalence <- input$prevalence / 100
    sensitivity <- input$sensitivity / 100
    specificity <- input$specificity / 100
    population <- input$population
    
    true_positives <- prevalence * sensitivity * population
    false_positives <- (1 - prevalence) * (1 - specificity) * population
    true_negatives <- (1 - prevalence) * specificity * population
    false_negatives <- prevalence * (1 - sensitivity) * population
    
    ppv <- true_positives / (true_positives + false_positives)
    npv <- true_negatives / (true_negatives + false_negatives)
    
    list(
      true_positives = true_positives,
      false_positives = false_positives,
      true_negatives = true_negatives,
      false_negatives = false_negatives,
      ppv = ppv,
      npv = npv
    )
  })
  
  # Text Summary
  output$summary <- renderText({
    metrics <- reactive_metrics()
    paste0(
      "With a population of ", input$population, ":\n",
      "- True Positives: ", round(metrics$true_positives), "\n",
      "- False Positives: ", round(metrics$false_positives), "\n",
      "- True Negatives: ", round(metrics$true_negatives), "\n",
      "- False Negatives: ", round(metrics$false_negatives), "\n",
      "Positive Predictive Value (PPV): ", round(metrics$ppv * 100, 2), "%\n",
      "Negative Predictive Value (NPV): ", round(metrics$npv * 100, 2), "%"
    )
  })
  
  # Bar Chart
  output$barPlot1 <- renderPlot({
    metrics <- reactive_metrics()
    bar_data <- c(metrics$true_positives, metrics$false_positives, metrics$false_negatives)
    bar_labels <- c("True Positives", "False Positives", "False Negatives")
    barplot(bar_data, names.arg = bar_labels, main = "Diagnostic Outcomes", ylab = "Number of Cases")
  })
  
  output$barPlot <- renderPlot({
    metrics <- reactive_metrics()
    bar_data <- c(metrics$true_positives, metrics$false_positives, metrics$true_negatives, metrics$false_negatives)
    bar_labels <- c("True Positives", "False Positives", "True Negatives", "False Negatives")
    barplot(bar_data, names.arg = bar_labels, main = "Diagnostic Outcomes", ylab = "Number of Cases")
  })
  
  # Line Graph (PPV vs Prevalence)
  output$lineGraph <- renderPlot({
    prevalence_seq <- seq(0.1, 50, by = 0.1) / 100
    ppv_seq <- prevalence_seq * (input$sensitivity / 100) /
      (prevalence_seq * (input$sensitivity / 100) + (1 - prevalence_seq) * (1 - input$specificity / 100))
    plot(prevalence_seq * 100, ppv_seq * 100, type = "l", col = "blue", lwd = 2,
         xlab = "Prevalence (%)", ylab = "Positive Predictive Value (%)",
         main = "PPV vs Prevalence")
  })
}

# Run the App
shinyApp(ui = ui, server = server)
