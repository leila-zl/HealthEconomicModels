# Cancer Diagnostic Metrics Explorer
#
# This Shiny application allows users to explore diagnostic metrics such as
# sensitivity, specificity, and prevalence for various cancer types across different demographics. 
# Users can select specific cancer types, genders,and age groups
# to visualize and understand the impact of these metrics on diagnostic outcomes.
#
# The prevalence data utilized in this application was extracted from the
# National Disease Registration Service (NDRS) 2021 cancer statistics, provided by NHS Digital.

# Author: Leila Zakka
# Date: 01/12/2024


library(shiny)
library(dplyr)

cancer_data <- read.csv("cancer_prev.csv")

cancer_data <- cancer_data %>%
  mutate(
    cancer_type = as.factor(cancer_type),
    gender = as.factor(gender),
    age_group = as.factor(age_group)
  )

# Define UI
ui <- fluidPage(
  titlePanel("Cancer Prevalence and Diagnostic Metrics Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("cancer_type", "Select Cancer Type", choices = levels(cancer_data$cancer_type), selected = "Prostate"),
      selectInput("gender", "Select Gender", choices = levels(cancer_data$gender), selected = "Male"),
      selectInput("age_group", "Select Age Group", choices = levels(cancer_data$age_group), selected = "65-69"),
      sliderInput("sensitivity", "Test Sensitivity (%)", min = 10, max = 100, value = 90, step = 1),
      sliderInput("specificity", "Test Specificity (%)", min = 10, max = 100, value = 95, step = 1),
      numericInput("population", "Enter Population Size", value = 1000, min = 1)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary Metrics", verbatimTextOutput("summary")),
        tabPanel("Bar Chart 1", plotOutput("barPlot1")),
        tabPanel("Bar Chart 2", plotOutput("barPlot")),
        tabPanel("PPV and Prevalence", plotOutput("lineGraph"))
      )
    )
  )
)


# Define Server Logic
server <- function(input, output) {
  # Reactive expression to filter data based on user input
  filtered_data <- reactive({
    cancer_data %>%
      filter(
        cancer_type == input$cancer_type,
        gender == input$gender,
        age_group == input$age_group
      )
  })
  
  # Reactive expression to calculate metrics
  reactive_metrics <- reactive({
    prevalence <- filtered_data()$prevalence / 100 
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
      npv = npv,
      prevalence = prevalence
    )
  })
  
  # Text Summary
  output$summary <- renderText({
    metrics <- reactive_metrics()
    paste0(
      "For ", input$cancer_type, " cancer in ", input$gender, "s between age ", input$age_group, " the test would produce:\n",
      "- True Positives: ", round(metrics$true_positives), "\n",
      "- False Positives: ", round(metrics$false_positives), "\n",
      "- True Negatives: ", round(metrics$true_negatives), "\n",
      "- False Negatives: ", round(metrics$false_negatives), "\n",
      "Positive Predictive Value (PPV): ", round(metrics$ppv * 100, 2), "%\n",
      "Negative Predictive Value (NPV): ", round(metrics$npv * 100, 2), "%"
    )
  })
  
  # Bar Charts
  output$barPlot1 <- renderPlot({
    metrics <- reactive_metrics()
    bar_data <- c(metrics$true_positives, metrics$false_positives, metrics$false_negatives)
    bar_labels <- c("True Positives", "False Positives", "False Negatives")
    bp <- barplot(bar_data, names.arg = bar_labels, main = "Diagnostic Outcomes", ylab = "Number of Cases", ylim = c(0, max(bar_data) * 1.2))
    text(x = bp, y = bar_data, label = round(bar_data), pos = 1, cex = 1.2, col = "black")
    prevalence <- metrics$prevalence
    prevalence_text <- paste("at", round(prevalence * 100, 2), "% disease prevalence")
    mtext(prevalence_text, side = 3, line = 0.3, cex = 1.1, col = "blue")
  })
  
  output$barPlot <- renderPlot({
    metrics <- reactive_metrics()
    bar_data <- c(metrics$true_positives, metrics$false_positives, metrics$true_negatives, metrics$false_negatives)
    bar_labels <- c("True Positives", "False Positives", "True Negatives", "False Negatives")
    barplot(bar_data, names.arg = bar_labels, main = "Diagnostic Outcomes", ylab = "Number of Cases", ylim = c(0, max(bar_data) * 1.2))
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

# Run the Shiny App
shinyApp(ui = ui, server = server)
