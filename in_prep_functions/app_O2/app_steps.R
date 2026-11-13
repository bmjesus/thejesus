library(shiny)

# Define UI
ui <- fluidPage(
  titlePanel("Dynamic Input App with Ranges"),
  sidebarLayout(
    sidebarPanel(
      fileInput("fileIn", 'O2 file', buttonLabel = "Upload file...", accept = ".txt"),
      fileInput("file", "Upload a text file (comma-separated vector):"),
      uiOutput("dynamicInputs"),
      actionButton("compile", "Retrieve Data from Ranges")
    ),
    mainPanel(
      h3("Retrieved Data from 'oxygen':"),
      verbatimTextOutput("outputData"),
      h3("Scatter Plots:"),
      uiOutput("scatterPlots")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  `%>%` <- dplyr::`%>%`
  app.data <- shiny::reactiveValues(data = NULL)
  
  # Load the O2 file and process it
  shiny::observeEvent(input$fileIn, {
    my.file <- input$fileIn
    matriz_O2 <- read.table(my.file$datapath, skip = 14, header = FALSE, sep = "\t")
    names(matriz_O2) <- c("date", "time", "time_s", "comment", "ch1_O2", "ch2_O2", "ch3_O2", "ch4_O2",
                          "ch1_temp", "ch2_temp", "ch3_temp", "ch4_temp", "pressure", "humidity",
                          "temp_probe", "int_temp", "analogue_in", "ch1_raw", "ch2_raw", "ch3_raw", "ch4_raw",
                          "ch1_signal", "ch2_signal", "ch3_signal", "ch4_signal",
                          "ch1_light", "ch2_light", "ch3_light", "ch4_light")
    matriz_O2 <- matriz_O2[, -30]
    matriz_O2$ch1_O2 <- as.numeric(gsub(",", ".", matriz_O2$ch1_O2))
    matriz_O2$time_s <- as.numeric(gsub(",", ".", matriz_O2$time_s))
    matriz_O2 <- matriz_O2[, c("ch1_O2", "time", "time_s")]
    app.data$data <- matriz_O2
  })
  
  oxygen <- reactive(app.data$data$ch1_O2)
  time <- reactive(app.data$data$time_s)
  
  # Reactive value to store the uploaded vector
  vectorData <- reactiveVal(NULL)
  
  # Observe the uploaded file
  observeEvent(input$file, {
    req(input$file)
    data <- readLines(input$file$datapath)
    vector <- unlist(strsplit(data, ","))
    vectorData(vector)
  })
  
  # Generate dynamic input fields
  output$dynamicInputs <- renderUI({
    req(vectorData())
    vector <- vectorData()
    lapply(seq_along(vector), function(i) {
      textInput(paste0("input_", i), label = paste("Input range for item", i), value = "")
    })
  })
  
  # Extract data from 'oxygen' and 'time' based on input ranges (by position)
  retrievedData <- eventReactive(input$compile, {
    req(vectorData())
    vector <- vectorData()
    result <- lapply(seq_along(vector), function(i) {
      rangeInput <- input[[paste0("input_", i)]]
      if (!is.null(rangeInput) && grepl("-", rangeInput)) {
        range <- as.numeric(unlist(strsplit(rangeInput, "-")))
        if (length(range) == 2 && range[1] > 0 && range[2] <= length(oxygen())) {
          list(
            oxygen = oxygen()[range[1]:range[2]],
            time = time()[range[1]:range[2]]
          )
        } else {
          NULL
        }
      } else {
        NULL
      }
    })
    result
  })
  
  # Output the retrieved data
  output$outputData <- renderText({
    req(retrievedData())
    sapply(seq_along(retrievedData()), function(i) {
      paste0(
        "Range ", i, ": ",
        "Oxygen - ", paste(retrievedData()[[i]]$oxygen, collapse = ", "),
        "; Time - ", paste(retrievedData()[[i]]$time, collapse = ", ")
      )
    }) %>% paste(collapse = "\n")
  })
  
  # Generate scatter plots dynamically
  output$scatterPlots <- renderUI({
    req(retrievedData())
    lapply(seq_along(retrievedData()), function(i) {
      plotOutput(outputId = paste0("scatterPlot_", i))
    }) %>% tagList()
  })
  
  observe({
    req(retrievedData())
    lapply(seq_along(retrievedData()), function(i) {
      local({
        index <- i
        output[[paste0("scatterPlot_", index)]] <- renderPlot({
          data <- retrievedData()[[index]]
          req(data)
          plot(
            x = data$time,
            y = data$oxygen,
            xlab = "Time (s)",
            ylab = "Oxygen",
            main = paste("Scatter Plot for Range", index),
            pch = 19
          )
        })
      })
    })
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)