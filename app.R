#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(plotly)
library(bslib)


# Example functions
examples <- list(
  "Absolute value" = "abs(x)",
  "Smooth function" = "sin(x)",
  "Vertical tangent at 0" = "x^(1/3)",
  "Piecewise corner" = "ifelse(x < 0, 3, x)"
)

# Function evaluation
evaluate_function <- function(expr, x) {
  tryCatch({
    result <- eval(
      parse(text = expr),
      envir = list(
        x = x,
        sin = sin,
        cos = cos,
        tan = tan,
        exp = exp,
        log = log,
        sqrt = sqrt,
        abs = abs,
        pi = pi,
        ifelse = ifelse
      )
    )
    
    as.numeric(result)
  }, error = function(e) {
    rep(NA_real_, length(x))
  })
}

# UI
ui <- page_sidebar(
  title = "Interactive Derivative App",
  theme = bs_theme(version = 5),
  
  sidebar = sidebar(
    width = 320,
    
    h4("Function"),
    
    selectInput(
      "example",
      "Choose an example:",
      choices = names(examples)
    ),
    
    textInput(
      "function_expr",
      "Enter f(x):",
      value = "abs(x)"
    ),
    
    hr(),
    
    h4("Graph settings"),
    
    sliderInput(
      "domain",
      "Domain:",
      min = -10,
      max = 10,
      value = c(-5, 5),
      step = 0.5
    ),
    
    sliderInput(
      "a",
      "Point a:",
      min = -5,
      max = 5,
      value = 0,
      step = 0.05
    ),
    
    sliderInput(
      "h",
      "Increment h:",
      min = 0.001,
      max = 2,
      value = 0.5,
      step = 0.001
    ),
    
    hr(),
    
    h5("Instructions"),
    
    p(
      "Move a towards the suspect point; ",
      "shrink h; compare left vs right estimates."
    ),
    
    hr(),
    
    hr(),
    
    actionButton(
      "reset",
      "Reset example",
      class = "btn-primary"
    )
  ),
  
  layout_columns(
    col_widths = c(8, 4),
    
    # Main graph
    card(
      card_header("Function and derivative estimates"),
      plotlyOutput("graph", height = "500px")
    ),
    
    # Calculations
    card(
      card_header("Numerical estimates"),
      
      h5("Function value"),
      verbatimTextOutput("function_value"),
      
      hr(),
      
      h5("Derivative estimates"),
      verbatimTextOutput("estimates"),
      
      hr(),
      
      h5("Conclusion"),
      uiOutput("conclusion")
    ),
    
    # Additional information
    card(
      card_header("How to interpret the graph"),
      p(
        "The secant line connects the function values ",
        "at a and a+h."
      ),
      p(
        "The tangent-line estimate uses the numerical ",
        "derivative at a."
      ),
      p(
        "Compare the left and right derivative estimates ",
        "as h becomes smaller."
      )
    ),
    
    card(
      card_header("Theory"),
      h5("Definition of Differentiability"),
      p(
        "A function is differentiable at a point a, if the limit of
        (f(a+h)-f(a))/h, as h approaches 0, exists and is finite. Similarly,
        the left and right-hand difference quotients approach the same finite value."
      ),
      h5("Numerical Limitations"),
      p(
        "Numerical estimates provide evidence, not a proof. ",
        "The app compares left and right difference quotients ",
        "at progressively smaller step sizes and may return ",
        "'Inconclusive' when the results are unstable."),
        
      h5("Example of a Differentiable Function"),
      p("f(x) = x"),
      h5("Example of a Continuous BUT Non-Differentiable Function"),
      p("f(x) = x^1/3")
    )
  )
)


# Classify differentiability using multiple step sizes
classify_differentiability <- function(expr, a, h) {
  
  # Use progressively smaller positive step sizes
  hs <- h / (2^(0:7))
  
  f0 <- evaluate_function(expr, a)
  
  left_est <- sapply(hs, function(step) {
    (f0 - evaluate_function(expr, a - step)) / step
  })
  
  right_est <- sapply(hs, function(step) {
    (evaluate_function(expr, a + step) - f0) / step
  })
  
  # Use the three smallest step sizes for the final comparison
  left_tail <- tail(left_est, 3)
  right_tail <- tail(right_est, 3)
  
  # Check whether values are finite
  finite_values <- all(is.finite(c(left_tail, right_tail)))
  
  # Relative variation among the last three estimates
  relative_variation <- function(values) {
    max(abs(values - median(values))) /
      max(1, abs(median(values)))
  }
  
  # Detect estimates that grow as h decreases
  growing <- function(values) {
    z <- abs(values)
    
    all(is.finite(z)) &&
      all(diff(z) > 0) &&
      z[3] > 1.25 * max(z[1], .Machine$double.eps)
  }
  
  if (!is.finite(f0)) {
    status <- "Inconclusive"
    reason <- "The function value at a could not be evaluated."
    derivative <- NA_real_
    
  } else if (
    growing(left_tail) || growing(right_tail)
  ) {
    status <- "Likely not differentiable"
    reason <- "Both derivative estimates grow as h decreases."
    derivative <- NA_real_
    
  } else if (finite_values) {
    
    left_center <- median(left_tail)
    right_center <- median(right_tail)
    
    difference <- abs(left_center - right_center) /
      max(1, abs(left_center), abs(right_center))
    
    left_stable <- relative_variation(left_tail) < 0.02
    right_stable <- relative_variation(right_tail) < 0.02
    
    if (left_stable && right_stable && difference < 0.02) {
      status <- "Likely differentiable"
      reason <- "Both estimates stabilize near the same finite value."
      derivative <- mean(c(left_center, right_center))
      
    } else if (
      left_stable && right_stable && difference >= 0.02
    ) {
      status <- "Likely not differentiable"
      reason <- "The left and right estimates stabilize at different values."
      derivative <- NA_real_
      
    } else {
      status <- "Inconclusive"
      reason <- "The estimates have not stabilized sufficiently."
      derivative <- NA_real_
    }
    
  } else {
    status <- "Inconclusive"
    reason <- "Some estimates are non-finite."
    derivative <- NA_real_
  }
  
  list(
    status = status,
    reason = reason,
    derivative = derivative,
    hs = hs,
    left_est = left_est,
    right_est = right_est
  )
}

# Server
server <- function(input, output, session) {
  
  # Update function when an example is selected
  observeEvent(input$example, {
    updateTextInput(
      session,
      "function_expr",
      value = examples[[input$example]]
    )
  })
  
  # Reset sliders
  observeEvent(input$reset, {
    updateSliderInput(
      session,
      "domain",
      value = c(-5, 5)
    )
    
    updateSliderInput(
      session,
      "a",
      value = 0
    )
    
    updateSliderInput(
      session,
      "h",
      value = 0.5
    )
  })
  
  # Function
  f <- reactive({
    req(input$function_expr)
    input$function_expr
  })
  
  
  calculations <- reactive({
    
    expr <- f()
    a <- input$a
    h <- input$h
    
    fa <- evaluate_function(expr, a)
    f_right <- evaluate_function(expr, a + h)
    f_left <- evaluate_function(expr, a - h)
    
    # Forward and backward differences
    right <- (f_right - fa) / h
    left <- (fa - f_left) / h
    
    # Central difference
    tangent <- (f_right - f_left) / (2 * h)
    
    # Multi-step differentiability test
    classification <- classify_differentiability(
      expr = expr,
      a = a,
      h = h
    )
    
    list(
      fa = fa,
      f_right = f_right,
      f_left = f_left,
      forward = right,
      left = left,
      tangent = tangent,
      classification = classification
    )
  })
  

  # Plot
  output$graph <- renderPlotly({
    
    expr <- f()
    
    a <- input$a
    h <- input$h
    
    domain <- input$domain
    
    calc <- calculations()
    
    x <- seq(
      domain[1],
      domain[2],
      length.out = 500
    )
    
    y <- evaluate_function(expr, x)
    
    data <- data.frame(x = x, y = y)
    
    # Function plot
    p <- ggplot(data, aes(x = x, y = y)) +
      geom_line(
        color = "black",
        linewidth = 1
      ) +
      
      # Secant line
      geom_abline(
        intercept = calc$fa - calc$forward * a,
        slope = calc$forward,
        color = "blue",
        linewidth = 0.9,
        linetype = "dashed"
      ) +
      
      # Tangent-line estimate
      geom_abline(
        intercept = calc$fa - calc$tangent * a,
        slope = calc$tangent,
        color = "red",
        linewidth = 0.9
      ) +
      
      # Point a
      geom_point(
        aes(x = a, y = calc$fa),
        color = "red",
        size = 3
      ) +
      
      # Point a+h
      geom_point(
        aes(
          x = a + h,
          y = calc$f_right
        ),
        color = "blue",
        size = 3
      ) +
      
      coord_cartesian(
        xlim = domain
      ) +
      
      labs(
        x = "x",
        y = "f(x)",
        title = paste("f(x) =", expr)
      ) +
      
      theme_minimal()
    
    ggplotly(p)
  })

  # Function value
  output$function_value <- renderText({
    
    calc <- calculations()
    
    paste0(
      "f(a) = ",
      signif(calc$fa, 6)
    )
  })

  # Derivative estimates
  output$estimates <- renderText({
    
    calc <- calculations()
    
    paste0(
      "Forward difference:\n",
      signif(calc$forward, 6),
      "\n\nLeft derivative estimate:\n",
      signif(calc$left, 6),
      "\n\nRight derivative estimate:\n",
      signif(calc$forward, 6),
      "\n\nTangent slope estimate:\n",
      signif(calc$tangent, 6)
    )
  })
  
  # Differentiability conclusion
  output$conclusion <- renderUI({
    
    calc <- calculations()
    result <- calc$classification
    
    tagList(
      h5(result$status),
      p(result$reason),
      
      if (is.finite(result$derivative)) {
        p(
          paste(
            "Estimated derivative:",
            signif(result$derivative, 6)
          )
        )
      }
    )
  })
  
}

# Run app
shinyApp(ui = ui, server = server)

