#' @title Function to lunch a shiny App that fits RLC models to fluorescence data
#' @description Function to lunch a shiny App that fits RLC models to fluorescence data
#' @return It does not return any R object
#' @keywords external
#' @export


rlc_fit<- function(){

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("RLC model fit"),

   sidebarLayout(

     sidebarPanel(

   fileInput("file", NULL, buttonLabel = "Upload file...",accept = ".csv"),

   selectInput(
     "rlc",
     "Select RLC",
     choices =  'No choices here yet'),

   selectInput(
     "f",
     "Select F channel",
     choices =  'No choices here yet'
     ),

   selectInput(
     "fm",
     "Select Fm channel",
     choices =  'No choices here yet'
   ),

   textInput("data_range", "Data range"),

   actionButton("calculate","Fit the models"),

   radioButtons("npq_type", "Chose NPQ type:",
                c("NPQ" = "npq",
                  "YNPQ" = "ynpq")),


   tableOutput("table_etr"),

   tableOutput("table_slope_parameters"),

   tableOutput("npq_table")

   #end of slidebar
     ),

   mainPanel(

   plotOutput("plot_etr")
   ,
   plotOutput("plot_npq")

   )
)
)


###############################################################################
#Server function
###############################################################################



#
server <- function(input, output, session) {

#importing a data file
  data <- reactive({
  req(input$file)

#checking the extension to see if it is a csv file
  ext <- tools::file_ext(input$file$name)
  switch(ext,
         csv = vroom::vroom(input$file$datapath),
         validate("Invalid file type")
         )

  })

#observes the event of loading the file so that the column names are loaded and the selection menus are filled
observeEvent(input$file,{

  updateSelectInput(session, "rlc",
                      label = "Select RLC",
                      choices = unique(data()$Label) )

  updateSelectInput(session, "f",
                    label = "Select F channel",
                    choices = colnames(data()),
                    selected = 'F')

  updateSelectInput(session, "fm",
                    label = "Select Fm channel",
                    choices = colnames(data()),
                    selected = 'Fm')

#calculating the number of light steps to provide the first
#range values


initial_range <- reactive(
    {
      req(input$file)
      #output <- length(data()$PAR[data()$label == input$rlc])
      output <- length(data()$PAR[data()$Label==unique(data()$Label[1])])
      #print(input$rlc)
      return(output)

    })


observeEvent(input$rlc,{

new_range <- length(data()$PAR[data()$Label == input$rlc])

#updates the box with the initial range (range of the first dataset)
updateTextInput(session, "data_range",
                  label = "Data range",
                  value = paste("1:",new_range, sep = ""))
                  #value = paste("1:",initial_range(), sep = ""))
  })
})

#almost everything is linked to clicking on the button "input$calculate"
#like that the calculations are only done once the right columns
#and range are selected

range <- reactive({
  eval(parse(text=input$data_range))
}
)

#range <- eventReactive(input$calculate,{
#  eval(parse(text=input$data_range))
#}
#)



#Calculating the different parameters

f <- eventReactive(input$calculate,{
unlist(data()[data()$Label == as.character(input$rlc), input$f])[range()]}
)



light <- eventReactive(input$calculate,{

  #unlist(unique(data()$PAR))[range()]
  unlist(data()$PAR[ data()$Label == input$rlc])[range()]

}
)


fm <- eventReactive(input$calculate,{
  unlist(data()[data()$Label == as.character(input$rlc), input$fm])[range()]
}
)

eff <- eventReactive(input$calculate,{
  (fm()-f())/fm()
}
)

etr <- eventReactive(input$calculate,{

  light()*eff()
}
  )

npq <- eventReactive(input$calculate,{

  thejesus::npq(fm = max(fm()), fmp = fm())
}
)


ynpq <- eventReactive(input$calculate,{

  thejesus::ynpq(f = f(), fm = max(fm()), fmp = fm())
}
)

#section to fit an ETR model
fit_etr <- eventReactive(input$calculate,{

  output <- thejesus::fit_etr(light = light(), etr = etr(),
                    model = 'P', plots = FALSE)


  return(output)
}
  )

#section to estimate alpha and Ek using the slope of the first 3 points

parameters_lm <- eventReactive(input$calculate,{

a_lm <- lm(etr()[1:3]~light()[1:3])

b_lm <- summary(a_lm)$coefficient[2]

return(b_lm)

}
)

#section to fit a NPQ model to the variable NPQ
fit_npq <- eventReactive(input$calculate,{

  thejesus::fit_npq_2011 (light = light(), npq = npq(),plots = FALSE)
}
)

#section to fit a NPQ model to the variable YNPQ
fit_ynpq <- eventReactive(input$calculate,{

  thejesus::fit_npq_2011 (light = light(), npq = ynpq(),plots = FALSE)
}
)


#ETR plot
output$plot_etr <- renderPlot({

  # Make sure requirements are met
  req(input$calculate)

  #eff <- unlist(data()[data()$salinity == as.character(input$replicate

  plot(light(),
       #unlist(data()[data()$salinity == as.character(input$rlc),5]),
       etr(),
       ylab = "rETR", xlab = "Light", pch =21, bg = 1, las = 1)

  #print(fit_etr())

  if (is.null(fit_etr()$pred.par) == FALSE){
    points(fit_etr()$pred.par,fit_etr()$pred,
           type = 'l', col = 2, lwd = 2)
    }


},res = 96)

#NPQ plot


output$plot_npq <- renderPlot({

  req(input$calculate)

if (input$npq_type=="npq"){

plot(light(),
       npq(),
       ylab = "NPQ", xlab = "Light", pch =21, bg = 1, las = 1)

  points(fit_npq()$predicted$par,fit_npq()$predicted$npq,
         type = 'l', col = 2, lwd = 2)
}else{

plot(light(),
       ynpq(),
       ylab = "YNPQ", xlab = "Light", pch =21, bg = 1, las = 1)

  points(fit_ynpq()$predicted$par,fit_ynpq()$predicted$npq,
         type = 'l', col = 2, lwd = 2)

}




  },res = 96)



#Constructing a table for the NPQ outputs


output$npq_table <- renderTable({

  req(input$calculate)

  if (input$npq_type=="npq"){
  fit_npq()$parameters
  }else{
    fit_ynpq()$parameters
    }
  })

#Constructing a table for the ETR Platt model outputs
output$table_etr <- renderTable({

  req(input$calculate)

table_etr <- as.data.frame(

  cbind(
  fit_etr()$alpha,
  fit_etr()$etrmax,
  fit_etr()$Ek ,
  fit_etr()$beta
))

  names(table_etr) <- c("Alpha","rETRmax","Ek", "Beta")
  return(table_etr)

})

output$table_slope_parameters <- renderTable({

  Ek_slope <- fit_etr()$etrmax/parameters_lm()

  alpha_lm_temp <- as.data.frame(cbind(parameters_lm(), Ek_slope))

  names(alpha_lm_temp) <- c("Alpha slope", "Ek slope")

  return(alpha_lm_temp)

})

}


shinyApp(ui,server)


#TO DO list


#Add other EP model

#Add the possibility of fitting using YII?

}




