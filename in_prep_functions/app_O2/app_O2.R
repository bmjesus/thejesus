ui <- fluidPage(

  # Application title
  titlePanel("O2 analysis - pyroscience"),

shiny::fluidRow(
  shiny::column(3,
          shiny::br(),
          fileInput("fileIn", 'O2 file', buttonLabel = "Upload file...",accept = ".txt"),
          shiny::fileInput('fileLight', 'Light file',
                           multiple = F, accept = c(".txt",".csv")),
          shiny::actionButton('confirm', 'Confirm', width='100%'),
          DT::dataTableOutput('lightSteps')
    ),
  shiny::column(6,
      plotly::plotlyOutput('O2_all', width = '100%', height='400px')
          ),
  shiny::column(3,
                DT::dataTableOutput('all_data')
  )
)
)


server <- function(input, output, session) {

`%>%` <- dplyr::`%>%`

#creating a reactive object to store data

app.data <- shiny::reactiveValues(data = NULL)



#Importing and parsing the O2 file

shiny::observeEvent(input$fileIn, {
  #get file details - also makes this reactive on input$fileIn
  my.file <- input$fileIn

  #load the data
    #import and skip the 13 lines of comments, data separated by Tab
  matriz_O2 <- read.table(my.file$datapath,skip = 14, header=F, sep="\t")

    #replace names for something more obvious
    names(matriz_O2) <- c("date","time","time_s","comment","ch1_O2","ch2_O2","ch3_O2","ch4_O2",
                           "ch1_temp","ch2_temp","ch3_temp","ch4_temp","pressure","humidity",
                           "temp_probe","int_temp","analogue_in","ch1_raw","ch2_raw","ch3_raw","ch4_raw",
                           "ch1_signal","ch2_signal","ch3_signal","ch4_signal",
                           "ch1_light","ch2_light","ch3_light","ch4_light")

    #remove last column
    matriz_O2<-matriz_O2[,-30]


    #if the decimal separator is a ',' (e.g. French computer) then the data is not numeric
    # added a step to replace ',' by '.' and convert it to numeric format
    matriz_O2$ch1_O2<-gsub(",",".",  matriz_O2$ch1_O2)
    matriz_O2$time<-gsub(",",".",  matriz_O2$time)
    matriz_O2$time_s<-gsub(",",".",  matriz_O2$time_s)


    matriz_O2$ch1_O2 <- as.numeric(matriz_O2$ch1_O2)
    matriz_O2$time <- matriz_O2$time
    matriz_O2$time_s <- as.numeric(matriz_O2$time_s)

    #make a smaller table for just one channel
    matriz_O2 <- matriz_O2[,c("ch1_O2","time","time_s")]


    print(matriz_O2)

    #push to app.data
  app.data$data <- matriz_O2
  app.data$O2 <- matriz_O2


})



##Importing light levels

shiny::observeEvent(input$fileLight, {

#get light file details - also makes this reactive on input$fileLight
light_file <- input$fileLight

#for the moment it only works with light levels separated by commas
#maybe there are more flexible functions to read this file?

light <- tryCatch({
  as.numeric(utils::read.csv(light_file$datapath, header = FALSE))
},error=function(e){NA}
)


#storing light
app.data$light <- light

#function to print the data table using light levels from the input file
output$lightSteps <- file_light()


})


#function to interatively add textInput with an ID and initial value to dataframe
shinyInput = function(FUN, len, id, int.value, ...) {
  #validate(need(character(len)>0,message=paste("")))
  inputs = character(len)
  for (i in seq_len(len)) {
    inputs[i] = as.character(FUN(paste0(id, i), label = NULL,
                                 value = int.value[i], width = '100px',...))
  }
  inputs
}

#a reactive function to fill the light data table when the user imports a light file

file_light <- reactive({

  req(is.null(app.data$light) == FALSE)

  DT::renderDataTable({

    #establish number of light levels to have (counts the items in the input file)
    no.light.steps <- length(app.data$light)

    #create empty slots that will store the time
    my.int.par <- rep('NA', no.light.steps)

    #make datatable for rendering
    DT::datatable(
      data.frame('Step' = 1:no.light.steps,
                 'PAR' = shinyInput(shiny::textInput,no.light.steps,"ls_", my.int.par)),
      selection="multiple",
      rownames = FALSE,
      escape = FALSE,
      #filter = list(position = 'top', clear = FALSE),
      #extensions = list("ColReorder" = NULL, "Buttons" = NULL),
      options = list(
        dom = 't',
        ColReorder = TRUE,
        #extensions = list("ColReorder" = NULL, "Buttons" = NULL),
        autoWidth=TRUE,
        lengthMenu = list(c(no.light.steps), c(as.character(no.light.steps))),
        rowCallback = DT::JS("function(r,d) {$(r).attr('height', '50px')}"),
        preDrawCallback = DT::JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = DT::JS('function() { Shiny.bindAll(this.api().table().node()); } '),
        columnDefs = list(list(targets=c(0), visible = T, width = '100px'))
      ))

  })


})

#an observe event that keeps storing the time intervals the user inputs
#once the user clicks on button confirm
shiny::observeEvent(input$confirm, {

  no.light.steps <- length(app.data$light)
  my.par <- vector()
  for(i in 1:no.light.steps){
    my.in <- paste('ls_',i,sep='')
    my.par[i] <- input[[my.in]] %>% as.numeric

  }

  #if the user has input light levels they get directed to app.data$par.levels
  app.data$par.levels <-  app.data$light
  if(is.null(app.data$light)==TRUE){
    #app.data$par.levels <- my.int.par
    app.data$par.levels <- my.par
  }

  test_output <<- app.data$par.levels



})

#Plotting section

plotO2_all <- shiny::reactive({
  plotly::plot_ly(x = app.data$data$time_s, y = app.data$data$ch1_O2,
                  type='scatter', mode='line') %>%
    plotly::layout(xaxis = list(title = 'Time (sec)'),
                   yaxis = list(title = 'O2 (umol/L)')
    )
})

output$O2_all <- plotly::renderPlotly({
  plotO2_all()
})





output$all_data <- DT::renderDataTable({app.data$O2})


}


shinyApp(ui,server)

