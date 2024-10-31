# Load the Shiny library to create the app
library(shiny)

# Define model constants outside the reactive block for clarity and easy adjustment
num_states <- 3 # Number of health states
state_names <- c("healthy", "disease", "death") # Names for each health state

# Define UI (User Interface) for the app
ui <- fluidPage(
  titlePanel("Basic Markov Model"), # Title for the app
  
  sidebarLayout(
    sidebarPanel(
      # Input: Numeric field to set the size of the initial cohort/population
      # value: starting value
      numericInput("cohort_size", "Enter Cohort size:", value = 1000, min = 10, max = 100000),
      
      # Input: Slider to set the number of cycles (steps) for the model simulation
      sliderInput("cycle_length", "Choose Cycle Length:", min = 1, max = 100, value = 40)
    ),
    
    mainPanel(
      # Output: Display the calculated average costs and QALYs in a table format
      tableOutput("payoffs")
    )
  )
)

# Define the server logic where the main calculations happen
server <- function(input, output) {
  
  # Reactive function to calculate the model results when inputs change
  model_results <- reactive({
    
    # Set up the transition probability matrix
    # Rows represent the "from" state; columns represent the "to" state
    # Each row should add up to 1, as it represents a full probability distribution
    transition_matrix <- matrix(c(0.96, 0.03, 0.01,  # Probabilities of moving from "healthy"
                                  0.00, 0.95, 0.05,  # Probabilities of moving from "disease"
                                  0.00, 0.00, 1.00), # Probabilities of moving from "death"
                                nrow = num_states, ncol = num_states, byrow = TRUE, 
                                dimnames = list(from = state_names, to = state_names))
    
    # Initialize an array to store state membership for each cycle
    # Each row represents a cycle, each column represents a state
    # cycle lenght is being fed from the input cycle_length, updated at each change
    state_membership <- array(NA_real_, dim = c(input$cycle_length, num_states), 
                              dimnames = list(cycle = 1:input$cycle_length, state = state_names))
    
    # Set the initial state (all individuals start in the "healthy" state)
    # cohort size is defined by the input cohort_size
    state_membership[1, ] <- c(input$cohort_size, 0, 0) # Full cohort in "healthy" at cycle 1
    
    # Calculate state memberships for each cycle by applying the transition matrix
    # for each loop we take one row (the previous row) of the state_membership matrix
    # and we multiply that row (look at it as a 1x3 matrix) by the transition matrix (3x3 matrix)
    # the result matrix of this matrix multiplication will give you a 1x3 matrix
    # where the columns represent the 'to' health state, which in the next loop will be used as a 'from' health state
    for (i in 2:input$cycle_length) {
      state_membership[i, ] <- state_membership[i - 1, ] %*% transition_matrix
    }
    
    # Define the payoff matrix, specifying costs and QALYs for each state
    # rows are health states, columns are payoffs (costs and QALYs)
    payoff_matrix <- matrix(c(50, 1000, 0,   # Costs for states
                              0.9, 0.6, 0),  # QALY values for states
                            byrow = FALSE, nrow = num_states, ncol = 2, 
                            dimnames = list(states = state_names, payoff = c("cost", "QALY")))
    
    # Calculate the total payoffs by multiplying state memberships with payoffs
    # the result matrix will have cycles as rows and total payoffs (from all health states) as columns
    payoff_trace <- state_membership %*% payoff_matrix
    payoffs <- colSums(payoff_trace) / input$cohort_size # Average cost and QALY
    
    # Return the payoffs as a data frame for clear tabular output
    return(data.frame(Payoff = names(payoffs), Value = payoffs))
  })
  
  # Render the payoffs table as output in the main panel
  output$payoffs <- renderTable({
    model_results()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
