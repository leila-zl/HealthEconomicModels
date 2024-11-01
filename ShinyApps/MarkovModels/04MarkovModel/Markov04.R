library(shiny)
library(shinythemes)
library(ggplot2)

# Create a 5-state Markov model with tunnel states for a hypothetical drug
# drug is taken by healthy as well not just the diseased
# costs are simplified: cost of healthy state = drug only, cost of diseased = drug + hospitalization

# Define constants for the model
n_states <- 5 # Number of health states
state_names <- c("Healthy", "Diseased", "Progressed", "Progressed2", "Dead") # Names for health states
HR_progression <- 0.95 # hazard ratio for probability of developing the disease on new drug
HR_death <- 0.9 # hazard ratio for probability of dying on new drug (same for healthy and diseased)

# Define UI (User Interface)
ui <- fluidPage(
  theme = shinytheme("simplex"),

  titlePanel("Leila's Markov Model App_04"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("cohort_size", "Enter Cohort size:", value = 1000, min = 10, max = 100000),
      sliderInput("psa_iterations", "Number of PSA iterations:", min = 1, max = 1000, value = 500),
      sliderInput("cycle_number", "Select Number of Cycles:", min = 1, max = 100, value = 40),
      sliderInput("cost_drug", "Select Base Case Drug Cost:", min = 1, max = 50000, value = 2000)
    ),
    
    mainPanel(
      plotOutput("plot1"),
      tags$h4("Base Case ICER", style = "font-size: 16px; font-weight: bold;"),
      verbatimTextOutput("ICER")
      
    )
  )
)

server <- function(input, output) {
  
  model_results <- reactive({
    
    # Base case parameters for New Drug and SoC
    params_basecase_new <- data.frame(
      cost_drug = input$cost_drug, # cost of new drug
      cost_hosp = 15000, # cost of hospitalisation
      
      u_healthy = 0.9, # Utility value for healthy state
      u_diseased = 0.6, # Utility value for diseased state
      u_progressed = 0.5,  # Utility value for progressed state       
      u_progressed2 = 0.45 # Utility value for progressed2 state   
    )
    
    params_basecase_soc <- data.frame(
      cost_drug = 1000,            # SoC drug cost
      cost_hosp = 15000,        # Same hospitalization cost
      
      u_healthy = 0.9,         
      u_diseased = 0.6,         
      u_progressed = 0.5,         
      u_progressed2 = 0.45         
    )
    
    # PSA parameters (probabilistic sampling for sensitivity analysis)
    params_psa_samples_new <- data.frame(
      cost_drug = rgamma(input$psa_iterations, shape = 25, scale = input$cost_drug / 25),
      cost_hosp = rgamma(input$psa_iterations, shape = 25, scale = 15000 / 25),
      
      u_healthy = rbeta(input$psa_iterations, 90, 10),
      u_diseased = rbeta(input$psa_iterations, 60, 40),
      u_progressed = rbeta(input$psa_iterations, 50, 60),
      u_progressed2 = rbeta(input$psa_iterations, 45, 65)
    )
    
    params_psa_samples_soc <- data.frame(
      cost_drug = rgamma(input$psa_iterations, shape = 25, scale = 1000 / 25), 
      cost_hosp = rgamma(input$psa_iterations, shape = 25, scale = 15000 / 25),
      
      u_healthy = rbeta(input$psa_iterations, 90, 10),
      u_diseased = rbeta(input$psa_iterations, 60, 40),
      u_progressed = rbeta(input$psa_iterations, 50, 60),
      u_progressed2 = rbeta(input$psa_iterations, 45, 65)
    )
    
    m_probability_base <- matrix(c(
      0.96, 0.03, 0.00, 0.00, 0.01, # from healthy
      0.00, 0.85, 0.12, 0.00, 0.03, # from diseased
      0.00, 0.00, 0.70, 0.20, 0.10, # from progressed
      0.00, 0.00, 0.00, 0.75, 0.25, # from progressed2
      0.00, 0.00, 0.00, 0.00, 1.00  # from dead
    ), nrow = 5, ncol = 5, byrow = TRUE, 
    dimnames = list(from = state_names, to = state_names))
    
    # Function to apply hazard ratios for the New Drug matrix
    apply_hazard_ratios <- function(base_matrix, HR_progression, HR_death) {
      m_probability_drug <- base_matrix
      
      # Apply hazard ratios to progression transitions
      m_probability_drug["Healthy", "Diseased"] <- base_matrix["Healthy", "Diseased"] * HR_progression
      m_probability_drug["Diseased", "Progressed"] <- base_matrix["Diseased", "Progressed"] * HR_progression
      m_probability_drug["Progressed", "Progressed2"] <- base_matrix["Progressed", "Progressed2"] * HR_progression
      
      # Apply hazard ratios to death transitions
      m_probability_drug["Healthy", "Dead"] <- base_matrix["Healthy", "Dead"] * HR_death
      m_probability_drug["Diseased", "Dead"] <- base_matrix["Diseased", "Dead"] * HR_death
      m_probability_drug["Progressed", "Dead"] <- base_matrix["Progressed", "Dead"] * HR_death
      m_probability_drug["Progressed2", "Dead"] <- base_matrix["Progressed2", "Dead"] * HR_death
      
      # Adjust self-transitions to sum rows to 1
      m_probability_drug["Healthy", "Healthy"] <- 1 - sum(m_probability_drug["Healthy", -1])
      m_probability_drug["Diseased", "Diseased"] <- 1 - sum(m_probability_drug["Diseased", -2])
      m_probability_drug["Progressed", "Progressed"] <- 1 - sum(m_probability_drug["Progressed", -3])
      m_probability_drug["Progressed2", "Progressed2"] <- 1 - sum(m_probability_drug["Progressed2", -4])
      
      return(m_probability_drug)
    }
    
    # Create a Markov model function to run simulations depending on input parameters, and HR
    model <- function(.params, apply_hr = FALSE) {
      with(.params, {
        # Select the appropriate transition matrix
        m_probability <- if (apply_hr) {
          apply_hazard_ratios(m_probability_base, HR_progression, HR_death)
        } else {
          m_probability_base
        }
        
        # Define the state membership matrix for each cycle
        state_membership <- array(NA_real_, dim = c(input$cycle_number, n_states), 
                                  dimnames = list(cycle = 1:input$cycle_number, state = state_names))
        
        # Initialize with full cohort in the healthy state at the start
        state_membership[1,] <- c(input$cohort_size, 0, 0, 0, 0)
        
        # Calculate the state memberships for the rest of the cycles
        for (i in 2:input$cycle_number) {
          state_membership[i, ] <- state_membership[i - 1, ] %*% m_probability
        }
        
        # Define payoff matrix for costs and QALYs in each state
        m_payoffs <- array(0, dim = c(n_states, 2), 
                           dimnames = list(states = state_names, payoffs = c("cost", "qaly")))
        m_payoffs["Healthy", "cost"] <- cost_drug
        m_payoffs["Healthy", "qaly"] <- u_healthy
        m_payoffs["Diseased", "cost"] <- cost_drug + cost_hosp
        m_payoffs["Diseased", "qaly"] <- u_diseased
        m_payoffs["Progressed", "cost"] <- cost_drug + cost_hosp * 1.2
        m_payoffs["Progressed", "qaly"] <- u_progressed 
        m_payoffs["Progressed2", "cost"] <- cost_drug + cost_hosp * 1.5
        m_payoffs["Progressed2", "qaly"] <- u_progressed2
        
        # Calculate costs and QALYs per cycle
        payoff_trace <- state_membership %*% m_payoffs
        colSums(payoff_trace) / input$cohort_size
      })
    }
    
    # Run the model for the base case and PSA for both New Drug and SoC
    basecase_new <- model(params_basecase_new, apply_hr = TRUE)
    basecase_soc <- model(params_basecase_soc, apply_hr = FALSE)
    
    # Run the model for each set of PSA parameters and compile results
    psa_new <- t(sapply(split(params_psa_samples_new, 1:input$psa_iterations), model, apply_hr = TRUE))
    psa_soc <- t(sapply(split(params_psa_samples_soc, 1:input$psa_iterations), model, apply_hr = FALSE))
    
    
    # Calculate incremental costs and QALYs for base case and PSA
    incremental_basecase <- basecase_new - basecase_soc
    ICER <- incremental_basecase[[1]] / incremental_basecase[[2]]
    incremental_psa <- as.data.frame(psa_new - psa_soc)
    
    # Plot the PSA results in a cost-effectiveness plane
    plot1 <- ggplot(incremental_psa, aes(x = qaly, y = cost)) +
      geom_point(color = "blue", size = 1) +
      # Highlight the base case result in red
      annotate("point", x = incremental_basecase["qaly"], y = incremental_basecase["cost"], color = "red", size = 4) +
      labs(title = "Cost-Effectiveness Plane",
           x = "Incremental QALYs", 
           y = "Incremental Cost") +
      theme(
        plot.title = element_text(size = 16, face = "bold"), 
        axis.title = element_text(size = 14) 
      ) +
      # Add cost-effectiveness threshold line
      geom_abline(slope = 30000, intercept = 0, color = "brown", linetype = "dashed") +
      # Label for cost-effectiveness threshold line
      annotate("text", x = max(incremental_psa$qaly) * 0.4, y = 30000 * max(incremental_psa$qaly) * 1,
               label = "Cost-Effectiveness Threshold", color = "brown", hjust = -0.1, vjust = -0.5, size = 4)
    
    # Return the plot and base case results
    list(plot1 = plot1, ICER = ICER)
  })
      
  
  # Render ICER
  output$ICER <- renderPrint({
    model_results()[["ICER"]]
  })
  
  # Render the PSA plot in the main panel
  output$plot1 <- renderPlot({
    model_results()[["plot1"]]
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
