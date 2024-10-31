library(shiny)
library(ggplot2)

# Define constants for the model
n_states <- 3 # Number of health states
state_names <- c("Healthy", "Diseased", "Dead") # Names for health states

# Define UI (User Interface)
ui <- fluidPage(
  titlePanel("Markov Model - Adjustable Costs and Utilities"),
  
  sidebarLayout(
    sidebarPanel(
      # Input for initial cohort size
      numericInput("cohort_size", "Enter Cohort size:", value = 1000, min = 10, max = 100000),
      
      # Slider inputs for number of PSA iterations and cycle length
      sliderInput("psa_iterations", "Number of PSA iterations:", min = 1, max = 1000, value = 500),
      sliderInput("cycle_length", "Choose Cycle Length:", min = 1, max = 100, value = 40),
      
      # Adjustable sliders for health state costs and utility values
      sliderInput("cost_h", "Base case cost of healthy:", min = 1, max = 10000, value = 1000),
      sliderInput("cost_d", "Base case cost of diseased:", min = 1, max = 50000, value = 20000),
      sliderInput("utility_h", "Base case utility - Healthy:", min = 0.85, max = 0.99, value = 0.9),
      sliderInput("utility_d", "Base case utility - Diseased:", min = 0.5, max = 0.75, value = 0.6)
    ),
    
    mainPanel(
      # Output for displaying the payoffs and the PSA plot
      verbatimTextOutput("payoffs"),
      plotOutput("plot1")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive function to calculate Markov model results based on user inputs
  model_results <- reactive({
    
    # Base case parameters
    params_basecase <- data.frame(
      p_healthy_diseased = 0.03,  # Probability of transitioning from healthy to diseased
      p_healthy_dead     = 0.01,  # Probability of transitioning from healthy to dead
      p_diseased_dead    = 0.05,  # Probability of transitioning from diseased to dead
      
      cost_healthy  = input$cost_h,   # Cost for healthy state
      cost_diseased = input$cost_d,   # Cost for diseased state
      u_healthy  = input$utility_h, # Utility value for healthy state
      u_diseased = input$utility_d  # Utility value for diseased state
    )
    
    # PSA parameters (probabilistic sampling for sensitivity analysis)
    params_psa_samples <- data.frame(
      # Transition probabilities with sampling from beta distribution
      p_healthy_diseased = rbeta(input$psa_iterations, 3, 97),
      p_healthy_dead     = rbeta(input$psa_iterations, 1, 99),
      p_diseased_dead    = rbeta(input$psa_iterations, 5, 95),
      
      # Variable costs and utilities generated using gamma and beta distributions
      cost_healthy  = rgamma(input$psa_iterations, shape = 25, scale = input$cost_h / 25),
      cost_diseased = rgamma(input$psa_iterations, shape = 25, scale = input$cost_d / 25),
      u_healthy  = rbeta(input$psa_iterations, (input$utility_h) * 100, 10),
      u_diseased = rbeta(input$psa_iterations, (input$utility_d) * 100, 40)
    )
    
    # Markov model function to easily calculate base case and PSA payoffs from parameters
    model <- function(.params) {
      with(.params, {
        # Transition probability matrix (3x3), defining transition probabilities for each health state
        m_probability <- matrix(0, 
                                nrow = n_states, ncol = n_states, 
                                dimnames = list(from = state_names, to = state_names))
        # Populate matrix with transition probabilities
        m_probability["Healthy",  "Healthy"]  <- 1 - p_healthy_diseased - p_healthy_dead
        m_probability["Healthy",  "Diseased"] <- p_healthy_diseased
        m_probability["Healthy",  "Dead"]     <- p_healthy_dead
        m_probability["Diseased", "Diseased"] <- 1 - p_diseased_dead
        m_probability["Diseased", "Dead"]     <- p_diseased_dead
        m_probability["Dead",     "Dead"]     <- 1
        
        # Define the state membership matrix for each cycle
        state_membership <- array(NA_real_, 
                                  dim = c(input$cycle_length, n_states), 
                                  dimnames = list(cycle = 1:input$cycle_length, state = state_names))
        # Initialize with full cohort in the healthy state at the start
        state_membership[1,] <- c(input$cohort_size, 0, 0)
        for (i in 2:input$cycle_length) {
          state_membership[i, ] <- state_membership[i - 1, ] %*% m_probability
        }
        
        # Define payoff matrix for costs and QALYs in each state
        m_payoffs <- array(0,
                           dim = c(n_states, 2),
                           dimnames = list(states = state_names, payoffs = c("cost", "qaly")))
        # Populate with base case costs and QALYs
        m_payoffs["Healthy","cost"] <- cost_healthy
        m_payoffs["Healthy","qaly"] <- u_healthy
        m_payoffs["Diseased","cost"] <- cost_diseased
        m_payoffs["Diseased","qaly"] <- u_diseased
        
        # Calculate costs and QALYs per cycle by multiplying state membership and payoff matrix
        payoff_trace <- state_membership %*% m_payoffs
        summary_results <- colSums(payoff_trace) / input$cohort_size # Average costs and QALYs
        
        return(summary_results)
      })
    }
    
    # Run the model for the base case
    basecase_results <- model(params_basecase)
    
    # Run the model for each set of PSA parameters and compile results
    psa_results <- t(sapply(
      X = split(params_psa_samples, 1:input$psa_iterations), 
      FUN = model, 
      simplify = TRUE))
    
    # Plot the PSA results in a cost-effectiveness plane
    plot1 <- ggplot(as.data.frame(psa_results),
                    aes(x = qaly, y = cost)) +
      geom_point(color = "blue", size = 1) + 
      # Highlight the base case result in red
      annotate("point", x = basecase_results["qaly"], y = basecase_results["cost"], color = "red", size = 4) +
      labs(title = "Cost-Effectiveness Plane",
           x = "QALYs", 
           y = "Cost") +
      theme_minimal() +
      # Add a line representing a willingness-to-pay threshold (£20,000/QALY)
      geom_abline(slope = 20000, intercept = 0, color = "brown", linetype = "dashed")
    
    # Return the plot and base case results
    list(plot1 = plot1, basecase_results = basecase_results)
  })
  
  # Render the payoffs as text output in the main panel
  output$payoffs <- renderPrint({
    model_results()[["basecase_results"]]
  })
  
  # Render the PSA plot in the main panel
  output$plot1 <- renderPlot({
    model_results()[["plot1"]]
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
