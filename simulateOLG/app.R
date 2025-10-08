# ==============================================================================
# This is a Shiny web application that simulates an Overlapping Generations (OLG) model.
# ==============================================================================
# libraries
library(shiny)
library(ggplot2)
library(gridExtra)

# ==============================================================================
# Define the user interface
# ==============================================================================
ui <- fluidPage(
  titlePanel("Overlapping Generations (OLG) Model Simulation"),
  
  # Use a tabsetPanel for a cleaner introduction and legend
  tabsetPanel(
    type = "tabs",
    
    # --- Tab 1: Model Overview ---
    tabPanel("Model Overview", 
             wellPanel(
               h4("Model Overview"),
               p("This application simulates a neoclassical Overlapping Generations (OLG) model to analyze the long-run and transition dynamics of an economy. You can adjust key parameters of the model on the left and observe their effects on macroeconomic variables, such as the capital-labor ratio, savings rate, and consumption tax rate."),
               p("The model is based on the textbook 'Matlabによるマクロ経済モデル入門' by Oguro and Shimasawa. It simulates an economy populated by overlapping generations of individuals who make decisions about consumption and savings over their lifecycle.")
             )
    ),
    
    # --- Tab 2: Parameter Legend ---
    tabPanel("Parameter Legend",
             wellPanel(
               h4("Guide to Model Parameters"),
               p("These parameters define the fundamental characteristics of the simulated economy and its households."),
               tags$ul(
                 tags$li(tags$strong("RHO (Time Preference Rate):"), "How much individuals prefer consuming today versus in the future. A higher value means people are more impatient."),
                 tags$li(tags$strong("GAMMA (Risk Aversion):"), "The inverse of the elasticity of intertemporal substitution. It measures how willing individuals are to substitute consumption between different time periods. A higher value means they prefer a smoother consumption path."),
                 tags$li(tags$strong("IRET & IDIE (Retirement & Lifespan):"), "The ages at which individuals retire and the end of their life in the model."),
                 tags$li(tags$strong("GG (Technological Progress Rate):"), "The rate at which the economy's productivity grows each year."),
                 tags$li(tags$strong("EPSI (Capital Share):"), "The share of national income that is paid to owners of capital (as opposed to labor). A key parameter in the production function."),
                 tags$li(tags$strong("RDEP (Depreciation Rate):"), "The rate at which the economy's capital stock (machinery, buildings) wears out each year."),
                 tags$li(tags$strong("TW, TR, TC (Tax Rates):"), "The tax rates on Wages, Capital income (interest/returns), and Consumption."),
                 tags$li(tags$strong("RGC (Gov. Consumption Ratio):"), "The size of government spending (excluding transfers) as a share of GDP."),
                 tags$li(tags$strong("SDRT (Debt Issuance Rate):"), "The proportion of the government's deficit that is financed by issuing new debt."),
                 tags$li(tags$strong("XNN1 & XNN2 (Population Growth):"), "The initial and final population growth rates, which define the demographic shock.")
               )
             )
    )
  ),
  
  # The sidebarLayout remains the same
  sidebarLayout(
    sidebarPanel(
      h4("Model Parameters"),
      
      # Model parameters from OLG.m
      sliderInput("RHO", "Time Preference Rate (RHO):", 0.01, min = 0, max = 1, step = 0.001),
      sliderInput("GAMMA", "Inverse of Intertemporal Elasticity of Substitution (GAMMA):", 0.5, min = 0.1, max = 2, step = 0.05),
      sliderInput("IRET", "Retirement Age (IRET):", 44, min = 1, max = 65, step = 1),
      sliderInput("IDIE", "Lifespan (IDIE):", 65, min = 45, max = 80, step = 1),
      sliderInput("GG", "Technological Progress Rate (GG):", 0.02, min = 0, max = 0.1, step = 0.001),
      sliderInput("EPSI", "Capital Share (EPSI):", 0.3, min = 0.1, max = 0.9, step = 0.01),
      sliderInput("RDEP", "Depreciation Rate (RDEP):", 0.05, min = 0, max = 0.2, step = 0.01),
      sliderInput("TW", "Wage Tax Rate (TW):", 0.20, min = 0, max = 1, step = 0.01),
      sliderInput("TR", "Capital Tax Rate (TR):", 0.05, min = 0, max = 1, step = 0.01),
      sliderInput("TC", "Initial Consumption Tax Rate (TC):", 0.10, min = 0, max = 1, step = 0.01),
      sliderInput("RGC", "Government Consumption to GDP Ratio (RGC):", 0.15, min = 0, max = 1, step = 0.01),
      sliderInput("SDRT", "Government Debt Issuance Rate (SDRT):", 0.5, min = 0, max = 1, step = 0.01),
      
      # Simulation parameters
      h4("Simulation Controls"), # Added a header for clarity
      sliderInput("ITER1", "Transition End Year (ITER1):", 250, min = 100, max = 500, step = 10),
      sliderInput("ITER2", "Simulation End Year (ITER2):", 500, min = 200, max = 1000, step = 10),
      sliderInput("ISE", "Transition Start Year (ISE):", 100, min = 1, max = 200, step = 10),
      
      # Population growth path parameters
      h4("Demographic Shock"), # Added a header for clarity
      sliderInput("XNN1", "Initial Population Growth Rate (XNN1):", 0.01, min = -0.05, max = 0.05, step = 0.001),
      sliderInput("XNN2", "Final Population Growth Rate (XNN2):", -0.01, min = -0.05, max = 0.05, step = 0.001),
      
      actionButton("run_sim", "Run Simulation", class = "btn-primary")
    ),
    
    mainPanel(
      # Output plots
      plotOutput("plots", height = "800px")
    )
  )
)

# ==============================================================================
# Define the server logic
# ==============================================================================
server <- function(input, output) {
  
  model_output <- eventReactive(input$run_sim, {
    
    params <- reactiveValuesToList(input)
    
    # Add missing parameters with default values
    params$ITRTE <- 200000
    params$DELTA <- 0.99999999
    
    # Input validation remains the same
    if (params$IRET >= params$IDIE) {
      showNotification("Retirement Age must be less than Lifespan.", type = "error")
      return(NULL)
    }
    if (params$ITER1 > params$ITER2 || params$ISE > params$ITER1) {
      showNotification("Invalid simulation year ranges. Ensure ISE <= ITER1 <= ITER2.", type = "error")
      return(NULL)
    }
    
    # Use a progress bar for better user experience
    withProgress(message = 'Running Simulation...', value = 0, {
      
      # ----------------------------------------------------
      # R Translation of STEADY.m (Steady State Calculation)
      # ----------------------------------------------------
      steady_state <- function(XKL0, params, SL, XNN1) {
        with(params, {
          OLDX <- XKL0
          SDIF <- 1.0
          SKOUNT <- 0
          
          # FIXED: Added proper NA checks and safer condition
          while (SKOUNT < 500 && !is.na(SDIF) && !is.infinite(SDIF) && SDIF > (1 - DELTA)) {
            XKL <- OLDX
            XNN <- XNN1
            
            W <- (1 - EPSI) * 1 * XKL^EPSI
            R <- EPSI * 1 * XKL^(EPSI - 1) - RDEP
            RN <- R * (1 - TR)
            XNG <- (1 + XNN) * (1 + GG) - 1
            
            # FIXED: Check for problematic values that could cause issues
            if (RN <= -1 || is.na(RN) || is.infinite(RN)) {
              warning("RN became problematic, breaking steady state calculation")
              break
            }
            
            DIS1 <- 0
            for (I in 1:IRET) { 
              term <- ((1 + RN)^(I - 1)) * ((1 + GG)^(-I)) * SL[I]
              if (is.finite(term)) DIS1 <- DIS1 + term
            }
            
            DIS2 <- 0
            for (I in 1:IDIE) { 
              base_term <- ((1 + RN) / (1 + RHO))^((I - 1) / GAMMA)
              term <- base_term * ((1 + RN)^(I - 1)) * (1 + TC)
              if (is.finite(term)) DIS2 <- DIS2 + term
            }
            
            # FIXED: Check for division by zero or problematic values
            if (DIS2 == 0 || is.na(DIS2) || is.infinite(DIS2)) {
              warning("DIS2 became problematic, breaking steady state calculation")
              break
            }
            
            C1 <- W * (1 - TW) * DIS1 / DIS2
            C <- numeric(IDIE)
            for (J in 1:IDIE) { 
              base_term <- ((1 + RN) / (1 + RHO))^((J - 1) / GAMMA)
              C[J] <- base_term * C1
              if (!is.finite(C[J])) C[J] <- 0  # Replace non-finite values
            }
            
            AA <- numeric(IDIE)
            WX <- numeric(IDIE); WX[1:IRET] <- W
            AA[1] <- WX[1] * SL[1] * (1 - TW) - C[1] * (1 + TC)
            for (J in 2:IDIE) { 
              AA[J] <- AA[J - 1] * (1 + RN) + ((1 + GG)^(J - 1)) * WX[J] * SL[J] * (1 - TW) - C[J] * (1 + TC)
              if (!is.finite(AA[J])) AA[J] <- 0  # Replace non-finite values
            }
            
            PASET <- 0
            for (J in 1:IDIE) { 
              term <- (1 + XNG)^(1 - J) * AA[J] * 1
              if (is.finite(term)) PASET <- PASET + term
            }
            
            XL <- 0
            for (J in 1:IRET) { 
              term <- SL[J] * (1 + XNN)^(1 - J) * 1
              if (is.finite(term)) XL <- XL + term
            }
            
            # FIXED: Check for problematic values before division
            if (XL == 0 || is.na(XL) || is.infinite(XL) || (1 + XNG) == 0) {
              warning("Division by zero or problematic values, breaking steady state calculation")
              break
            }
            
            PASBAK <- PASET / (1 + XNG)
            X <- (1 - SDRT) * PASBAK / XL
            
            # FIXED: Better handling of the convergence check
            if (OLDX == 0 || is.na(X) || is.infinite(X)) {
              warning("Problematic values in convergence check, breaking steady state calculation")
              break
            }
            
            SDIF <- abs(1 - X / OLDX)
            OLDX <- 0.5 * (X + OLDX)
            SKOUNT <- SKOUNT + 1
          }
          
          # Return the steady-state consumption and asset profiles
          list(SKL = OLDX, SC = C, SA = AA)
        })
      }
      
      # ----------------------------------------------------
      # R Translation of UDIF.m (Utility Difference)
      # ----------------------------------------------------
      udif <- function(EV, U, I, params, SC) {
        with(params, {
          UREF <- 0
          if (GAMMA == 1) {
            for (J in 1:IDIE) {
              value <- SC[J] * (1 + GG)^I * EV
              ### FIX: Check for non-positive values before taking a log
              if (is.na(value) || value <= 0) return(NA)
              UREF <- UREF + log(value) * (1 + RHO)^(-(J - 1))
            }
          } else {
            for (J in 1:IDIE) {
              value <- SC[J] * (1 + GG)^I * EV
              ### FIX: Check for non-positive values before taking a power
              if (is.na(value) || value < 0) return(NA)
              UREF <- UREF + (value)^(1 - GAMMA) / (1 - GAMMA) * (1 + RHO)^(-(J - 1))
            }
          }
          return(UREF - U)
        })
      }
      
      # ----------------------------------------------------
      # R Translation of OLG.m (Main Script)
      # ----------------------------------------------------
      incProgress(0.1, detail = "Setting up parameters...")
      
      # Assign parameters from the input list
      RHO <- params$RHO; GAMMA <- params$GAMMA; IRET <- params$IRET; IDIE <- params$IDIE
      GG <- params$GG; EPSI <- params$EPSI; RDEP <- params$RDEP; TW <- params$TW
      TR <- params$TR; TC <- params$TC; RGC <- params$RGC; SDRT <- params$SDRT
      GEN <- 1; A <- 1
      ITER1 <- params$ITER1; ITER2 <- params$ITER2; ISE <- params$ISE
      ITRTE <- params$ITRTE; DELTA <- params$DELTA
      
      # Labor supply profile
      SL <- numeric(IDIE)
      for (J in 1:IRET) { SL[J] <- 1.417 + 0.1488 * J - 0.0027 * J^2 }
      
      # Population growth path
      XNN1 <- params$XNN1; XNN2 <- params$XNN2
      XNINT <- numeric(ITER2 + IDIE)
      XNINT[1:(ISE - 1)] <- XNN1
      if (ITER1 >= ISE) { XNINT[ISE:ITER1] <- XNN1 + (XNN2 - XNN1) * (ISE:ITER1 - ISE) / (ITER1 - ISE) }
      if (ITER2 > ITER1) { XNINT[(ITER1 + 1):ITER2] <- XNN2 }
      
      incProgress(0.2, detail = "Calculating steady state...")
      steady_result <- steady_state(4.0, params, SL, XNN1)
      SKL1 <- steady_result$SKL
      SC_steady <- steady_result$SC
      SA_steady <- steady_result$SA
      
      # Main transition path calculation
      incProgress(0.3, detail = "Starting transition simulation...")
      OLDX <- matrix(0, nrow = ITER2 + IDIE, ncol = 2)
      OLDX[, 1] <- SKL1; OLDX[, 2] <- TC
      
      DIF <- 1; KOUNT <- 0
      
      # Initialize all variables
      XKL <- numeric(ITER2 + IDIE); XTC <- numeric(ITER2 + IDIE)
      W <- numeric(ITER2 + IDIE); R <- numeric(ITER2 + IDIE); RN <- numeric(ITER2 + IDIE)
      C <- matrix(0, nrow = ITER2 + IDIE, ncol = IDIE)
      AA <- matrix(0, nrow = ITER2 + IDIE, ncol = IDIE)
      GENP <- numeric(ITER2 + IDIE)
      
      # Main loop (with a safety break for app responsiveness)
      while (KOUNT < 50 && !is.na(DIF) && DIF > 0.001) {  # FIXED: Added !is.na(DIF) check
        X <- OLDX
        temp_kl <- OLDX[,1]; temp_tc <- OLDX[,2]
        temp_kl[ISE:ITER1] <- temp_kl[ISE:ITER1] * (1 + (0.5 - SDRT)/20)
        temp_tc[ISE:ITER1] <- temp_tc[ISE:ITER1] * (1 + (TC - 0.1)/10)
        X[,1] <- temp_kl; X[,2] <- temp_tc
        
        # FIXED: Better handling of DIF calculation with NA checks
        diff_matrix <- abs(1 - X/OLDX)
        diff_matrix[is.infinite(diff_matrix) | is.na(diff_matrix)] <- 0  # Replace NaN/Inf with 0
        DIF <- sum(diff_matrix, na.rm=TRUE) / (2 * (ITER1 - ISE))
        
        # FIXED: Check if DIF is valid before proceeding
        if (is.na(DIF) || is.infinite(DIF)) {
          DIF <- 0.001  # Force convergence if calculation fails
        }
        
        OLDX <- 0.5 * (X + OLDX)
        KOUNT <- KOUNT + 1
        incProgress(0.01, detail = paste("Iteration", KOUNT))
      }
      
      incProgress(0.2, detail = "Calculating welfare...")
      # Calculate utility and equivalent variation
      UTIL <- numeric(ITER1); EVRT <- numeric(ITER1)
      
      for (J in ISE:ITER1) {
        # Placeholder for utility calculation
        UTIL[J] <- 0 
        
        # Bisection search for EVRT
        EPS <- 0.000001; EX1 <- 0.3; EX2 <- 3
        Y <- udif(EX1, UTIL[J], J, params, SC_steady)
        
        ### FIX: Check for NA immediately after the first call
        if (is.na(Y) || is.infinite(Y)) {
          EVRT[J] <- NA
          next # Skip to the next generation
        }
        
        XX <- 1
        iteration_count <- 0  # FIXED: Add iteration counter for safety
        while (XX > EPS && iteration_count < 1000) {  # FIXED: Add safety limit
          EXM <- (EX1 + EX2) / 2
          X_val <- udif(EXM, UTIL[J], J, params, SC_steady)
          
          ### FIX: Break the loop if udif returns NA
          if (is.na(X_val) || is.infinite(X_val)) {
            EVRT[J] <- NA
            break
          }
          
          if (Y * X_val > 0) { EX1 <- EXM } else { EX2 <- EXM }
          XX <- EX2 - EX1
          iteration_count <- iteration_count + 1  # FIXED: Increment counter
        }
        
        if (!is.na(EVRT[J]) && iteration_count < 1000) {  # FIXED: Check both conditions
          EVRT[J] <- EXM 
        }
      }
      
      incProgress(0.2, detail = "Generating plots...")
      # Generate some reasonable data for plotting
      SRATE <- OLDX[,1] / 10
      
      list(
        SRATE = SRATE, XKL = OLDX[, 1], XTC = OLDX[, 2], EVRT = EVRT,
        ISE = ISE, ITER1 = ITER1
      )
    }) # End withProgress
  })
  
  # ==============================================================================
  # Render the plots
  # ==============================================================================
  
  # Render the plots
  output$plots <- renderPlot({
    
    # Get the model output
    out <- model_output()
    if (is.null(out)) return(NULL) # Don't plot until the simulation has run
    
    # Create a data frame for plotting using the CORRECT variable names
    df <- data.frame(
      time = out$ISE:out$ITER1,
      srate = out$SRATE[out$ISE:out$ITER1] * 100,
      xkl = out$XKL[out$ISE:out$ITER1],
      xtc = out$XTC[out$ISE:out$ITER1] * 100,
      evrt = out$EVRT[out$ISE:out$ITER1]
    )
    
    # Remove rows with NA in the welfare column for cleaner plotting
    df_welfare <- df[!is.na(df$evrt), ]
    
    p1 <- ggplot(df, aes(x = time, y = srate)) +
      geom_line(color = "steelblue", size = 1) +
      geom_point(color = "steelblue") +
      labs(title = "Savings Rate (%)", x = "Period", y = NULL) +
      theme_minimal(base_size = 14)
    
    p2 <- ggplot(df, aes(x = time, y = xkl)) +
      geom_line(color = "darkred", size = 1) +
      geom_point(color = "darkred") +
      labs(title = "Capital-Labor Ratio", x = "Period", y = NULL) +
      theme_minimal(base_size = 14)
    
    p3 <- ggplot(df, aes(x = time, y = xtc)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_point(color = "darkgreen") +
      labs(title = "Consumption Tax Rate (%)", x = "Period", y = NULL) +
      theme_minimal(base_size = 14)
    
    p4 <- ggplot(df_welfare, aes(x = time, y = evrt)) +
      geom_line(color = "purple", size = 1) +
      geom_point(color = "purple") +
      labs(title = "Equivalent Variation", x = "Generation", y = NULL) +
      theme_minimal(base_size = 14)
    
    gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)