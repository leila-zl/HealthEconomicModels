library(shiny)

# Define UI (User Interface)
ui <- fluidPage(
  titlePanel("Basic Markov Model"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("cohort_size", "Enter Cohort size:", value = 1000, min = 10, max = 100000),
      
      # Slider input
      sliderInput("cycle_length", "Choose Cycle Length:", min = 1, max = 100, value = 40)
      
    ),
    
    mainPanel(
      verbatimTextOutput("payoffs")
      
    )
  )
)

# Define server logic
server <- function(input, output) {
  

  
  model_results <- reactive({
    
    n_s <- 3 # number of health states
    v_state_names <- c("healthy", "disease", "death")
    #matrix for transition probabilities
    m_p <- matrix(c(0.96, 0.03, 0.01, 0.00, 0.95, 0.05, 0.00, 0.00, 1), 
                  nrow = 3, ncol = 3, byrow = TRUE, 
                  dimnames = list(from = v_state_names, to = v_state_names))
    m_p
    rowSums(m_p)
    # create an empty array with NAs, 
    # then fill it up by multiplying cycle rows one by with each probability matrix column
    # each of the 3 columns multiplied by the whole row creating a new element in the cycle row
    state_membership <- array(NA_real_, dim = c(input$cycle_length,n_s), 
                              dimnames = list(cycle=1:input$cycle_length, state= v_state_names))
    state_membership
    state_membership [1,] <- c(input$cohort_size,0,0)
    for (i in 2:input$cycle_length) {state_membership[i, ] <- state_membership[i-1, ] %*% m_p}
    
    # calculate payoffs
    m_payoff <- matrix(c(50, 1000, 0, 0.9, 0.6, 0), byrow = FALSE, nrow = 3, ncol = 2, 
                       dimnames = list(states = v_state_names, payoff = c("cost", "QALY")))
    # calculate costs and QALYs for each cycle
    payoff_trace <- state_membership %*% m_payoff
    payoffs <- colSums(payoff_trace)/input$cohort_size #average cost and qualys
    
    return(payoffs)
    
  })

    
  # Render the payoffs table
  output$payoffs <- renderPrint({
    model_results()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
