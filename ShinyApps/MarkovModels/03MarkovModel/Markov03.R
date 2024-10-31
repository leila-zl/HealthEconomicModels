library(shiny)
library(shinythemes)
library(ggplot2)

# Create a simple 3-state Markov model for a hypothetical preventative drug
# drug is taken by healthy as well not just the diseased
# costs are simplified: cost of healthy state = drug only, cost of diseased = drug + hospitalization
# cost of the drug is adjustable in the app to see the impact of lowering or increasing its price

# Define constants for the model
n_states <- 3 # Number of health states
state_names <- c("Healthy", "Diseased", "Dead") # Names for health states
HR_disease <- 0.8 # hazard ratio for probability of developing the disease on new drug
HR_dead <- 0.8 # hazard ratio for probability of dying on new drug (same for healthy and diseased)

# Define UI (User Interface)
ui <- fluidPage(
  theme = shinytheme("simplex"),

  titlePanel("Leila's Markov Model App"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("cohort_size", "Enter Cohort size:", value = 1000, min = 10, max = 100000),
      sliderInput("psa_iterations", "Number of PSA iterations:", min = 1, max = 1000, value = 500),
      sliderInput("cycle_number", "Select Number of Cycles:", min = 1, max = 100, value = 40),
      sliderInput("cost_drug", "Select Base Case Drug Cost:", min = 1, max = 50000, value = 20000)
    ),
    
    mainPanel(
      # display incremental costs and qalys plot
      plotOutput("plot1"),
      # display ICER
      tags$h4("Base Case ICER", style = "font-size: 16px; font-weight: bold;"),
      verbatimTextOutput("ICER")
      
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive function to calculate results for both New Drug and Standard of Care
  model_results <- reactive({
    
    # Base case parameters for New Drug and SoC
    params_basecase_new <- data.frame(
      p_healthy_diseased = 0.03 * HR_disease, # Probability of transitioning from healthy to diseased
      p_healthy_dead = 0.01 * HR_dead, # Probability of transitioning from healthy to dead
      p_diseased_dead = 0.05, # Probability of transitioning from diseased to dead
      cost_drug = input$cost_drug, # cost of new drug
      cost_hosp = 15000, # cost of hospitalisation
      u_healthy = 0.9, # Utility value for healthy state
      u_diseased = 0.6 # Utility value for diseased state
    )
    
    params_basecase_soc <- data.frame(
      p_healthy_diseased = 0.03,  
      p_healthy_dead = 0.01,
      p_diseased_dead = 0.05,
      cost_drug = 1000,            # SoC drug cost
      cost_hosp = 15000,        # Same hospitalization cost
      u_healthy = 0.9,         
      u_diseased = 0.6         
    )
    
    # PSA parameters (probabilistic sampling for sensitivity analysis)
    params_psa_samples_new <- data.frame(
      
      # Transition probabilities with sampling from beta distribution
      # alpha (shape1) = mean * 100, beta (shape2) = 100 - alpha
      p_healthy_diseased = rbeta(input$psa_iterations, 3, 97) * HR_disease,
      p_healthy_dead = rbeta(input$psa_iterations, 1, 99) * HR_dead,
      p_diseased_dead = rbeta(input$psa_iterations, 5, 95),
      
      # costs and utilities generated using random sampling from gamma and beta distributions
      # set shape parameter for gamma distribution to desired value based on prior knowledge
      # set scale parameter for gamma distribution = mean / shape
      cost_drug = rgamma(input$psa_iterations, shape = 25, scale = input$cost_drug / 25),
      cost_hosp = rgamma(input$psa_iterations, shape = 25, scale = 15000 / 25),
      u_healthy = rbeta(input$psa_iterations, 90, 10),
      u_diseased = rbeta(input$psa_iterations, 60, 40)
    )
    
    params_psa_samples_soc <- data.frame(
      p_healthy_diseased = rbeta(input$psa_iterations, 3, 97),
      p_healthy_dead = rbeta(input$psa_iterations, 1, 99),
      p_diseased_dead = rbeta(input$psa_iterations, 5, 95),
      cost_drug = rgamma(input$psa_iterations, shape = 25, scale = 1000 / 25), 
      cost_hosp = rgamma(input$psa_iterations, shape = 25, scale = 15000 / 25),
      u_healthy = rbeta(input$psa_iterations, 90, 10),
      u_diseased = rbeta(input$psa_iterations, 60, 40)
    )
    
    # Markov model function
    model <- function(.params) {
      with(.params, {
        # Transition probability matrix (3x3), defining transition probabilities for each health state
        m_probability <- matrix(0, nrow = n_states, ncol = n_states, 
                                dimnames = list(from = state_names, to = state_names))
        
        # Populate matrix with transition probabilities
        m_probability["Healthy", "Healthy"] <- 1 - p_healthy_diseased - p_healthy_dead
        m_probability["Healthy", "Diseased"] <- p_healthy_diseased
        m_probability["Healthy", "Dead"] <- p_healthy_dead
        m_probability["Diseased", "Diseased"] <- 1 - p_diseased_dead
        m_probability["Diseased", "Dead"] <- p_diseased_dead
        m_probability["Dead", "Dead"] <- 1
        
        # Define the state membership matrix for each cycle
        state_membership <- array(NA_real_, dim = c(input$cycle_number, n_states), 
                                  dimnames = list(cycle = 1:input$cycle_number, state = state_names))
        # Initialize with full cohort in the healthy state at the start
        state_membership[1,] <- c(input$cohort_size, 0, 0)
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
        
        # Calculate costs and QALYs per cycle by multiplying state membership and payoff matrix
        payoff_trace <- state_membership %*% m_payoffs
        colSums(payoff_trace) / input$cohort_size
      })
    }
    
    # Run the model for the base case and PSA for both New Drug and SoC
    basecase_new <- model(params_basecase_new)
    basecase_soc <- model(params_basecase_soc)
    
    # Run the model for each set of PSA parameters and compile results
    psa_new <- t(sapply(split(params_psa_samples_new, 1:input$psa_iterations), model))
    psa_soc <- t(sapply(split(params_psa_samples_soc, 1:input$psa_iterations), model))
    
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
      geom_abline(slope = 20000, intercept = 0, color = "brown", linetype = "dashed")
    
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
