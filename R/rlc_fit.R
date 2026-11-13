###############################################################################
#Internal helpers for rlc_fit() (not exported)
###############################################################################

#guess the column separator from the first lines of the file
#tab and ";" are tested before "," so that files with a decimal comma
#(e.g. French computers: "0,52;1,3") are not split on the decimals
rlc_guess_delim <- function(path){

  lines <- readLines(path, n = 50, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  for (d in c("\t", ";", ",")){
    if (length(lines) > 0 && all(grepl(d, lines, fixed = TRUE))){
      return(d)
    }
  }

  ","
}


#read a delimited text file: guesses the separator, reads everything as
#text and then converts to numeric every column that only contains numbers,
#accepting both "." and "," as decimal separator
rlc_read_file <- function(path){

  delim <- rlc_guess_delim(path)

  d <- vroom::vroom(path, delim = delim,
                    col_types = vroom::cols(.default = vroom::col_character()),
                    trim_ws = TRUE, show_col_types = FALSE, progress = FALSE)
  d <- as.data.frame(d, stringsAsFactors = FALSE)

  decimal_comma <- FALSE

  for (i in seq_along(d)){
    x <- trimws(d[[i]])
    empty <- is.na(x) | x == ""
    num <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
    if (any(!empty) && all(!is.na(num[!empty]))){
      if (any(grepl(",", x[!empty], fixed = TRUE))) decimal_comma <- TRUE
      d[[i]] <- num
    }
  }

  list(data = d, delim = delim, decimal_comma = decimal_comma)
}


#simplify a column name to compare it with the list of known names:
#lower case, prime (') becomes "p" (Fm' -> fmp), other symbols removed
rlc_normalize_name <- function(x){
  x <- tolower(x)
  x <- gsub("['\u2019\u2032]", "p", x)
  gsub("[^a-z0-9]", "", x)
}


#guess which column corresponds to a variable from a list of known names
#first tries the whole column name ("Fm'"), then its first word
#("PAR (umol m-2 s-1)" -> "PAR"); returns "" if nothing matches
rlc_guess_column <- function(candidates, aliases){

  if (length(candidates) == 0) return("")

  full <- rlc_normalize_name(candidates)
  first <- rlc_normalize_name(sub("^\\s*([^ (\\[_.-]+).*$", "\\1", candidates))

  for (key in list(full, first)){
    for (a in aliases){
      hit <- candidates[key == a]
      if (length(hit) > 0) return(hit[1])
    }
  }

  ""
}


#known names for each column (already normalized)
rlc_aliases <- list(
  label = c("label", "rlc", "curve", "sample", "sampleid", "id", "name", "group"),
  par = c("par", "light", "ppfd", "pfd", "irradiance", "e"),
  f = c("f", "fp", "ft", "fs"),
  fm = c("fm", "fmp", "fmprime")
)


###############################################################################
#Shiny app
###############################################################################

#' @title Function to lunch a shiny App that fits RLC models to fluorescence data
#' @description Function to lunch a shiny App that fits RLC models to
#' fluorescence data. The data file (.csv, .txt or .tsv; separated by ",", ";"
#' or tabs; "." or "," as decimal separator) needs one column with the light
#' (PAR), one with F and one with Fm, and optionally one column identifying
#' each light curve. The app tries to recognise these columns from their names
#' and the user can change the selection in the app.
#' @return It does not return any R object
#' @keywords external
#' @export


rlc_fit<- function(){

no_choice <- c("Choose a column" = "")

ui <- shiny::fluidPage(

   # Application title
   shiny::titlePanel("RLC model fit"),

   shiny::sidebarLayout(

     shiny::sidebarPanel(

   shiny::fileInput("file", NULL, buttonLabel = "Upload file...",
                    accept = c(".csv", ".txt", ".tsv")),

   #column selection, filled in (with a best guess) when a file is loaded
   shiny::selectInput("label_col", "RLC column", choices = no_choice),

   shiny::selectInput("par_col", "PAR column", choices = no_choice),

   shiny::selectInput("f", "F column", choices = no_choice),

   shiny::selectInput("fm", "Fm column", choices = no_choice),

   shiny::selectInput(
     "rlc",
     "Select RLC",
     choices =  'No choices here yet'),

   #one tickbox per light step of the selected RLC
   #(filled in when a RLC is selected, all ticked by default)
   shiny::checkboxGroupInput("light_steps", "Light levels to use (PAR)",
                             choices = character(0)),

   shiny::radioButtons("etr_model", "Chose fitting model:",
                c("Platt et al" = "P",
                  "Eilers Petters" = "EP")),

   shiny::actionButton("calculate","Fit the models"),

   shiny::radioButtons("npq_type", "Chose NPQ type:",
                c("NPQ" = "npq",
                  "YNPQ" = "ynpq")),


   shiny::tableOutput("table_etr"),

   shiny::tableOutput("table_slope_parameters"),

   shiny::tableOutput("npq_table")

   #end of slidebar
     ),

   shiny::mainPanel(

     shiny::tabsetPanel(

       shiny::tabPanel("Plots",
         shiny::plotOutput("plot_etr"),
         shiny::plotOutput("plot_npq")
       ),

       shiny::tabPanel("Data preview",
         shiny::verbatimTextOutput("preview_info"),
         shiny::tableOutput("preview")
       )
     )

   )
)
)


###############################################################################
#Server function
###############################################################################

server <- function(input, output, session) {

#importing a data file
#if there is no file, or the file is not valid, the plots show a message
#instead (shiny::validate) and nothing else is calculated
  imported <- shiny::reactive({

  shiny::validate(shiny::need(!is.null(input$file),
    "No valid dataset: please upload a .csv, .txt or .tsv file with the light (PAR), F and Fm columns"))

#checking the extension
  ext <- tolower(tools::file_ext(input$file$name))
  shiny::validate(shiny::need(ext %in% c("csv", "txt", "tsv"),
    "No valid dataset: the file must be a .csv, .txt or .tsv file"))

  imp <- tryCatch(rlc_read_file(input$file$datapath), error = function(e) NULL)

  shiny::validate(
    shiny::need(!is.null(imp) && nrow(imp$data) > 0,
                "No valid dataset: the file could not be read or is empty"),
    shiny::need(sum(vapply(imp$data, is.numeric, logical(1))) >= 3,
                "No valid dataset: the file needs at least 3 numeric columns (PAR, F and Fm)")
  )

  imp

  })

  data <- shiny::reactive(imported()$data)


#when a file is loaded, fill the column menus and pre-select the most likely
#column for each variable (the user can change it)
shiny::observeEvent(input$file,{

  d <- tryCatch(data(), error = function(e) NULL)
  shiny::req(d)

  all_cols <- names(d)
  num_cols <- all_cols[vapply(d, is.numeric, logical(1))]

  par_guess <- rlc_guess_column(num_cols, rlc_aliases$par)
  f_guess <- rlc_guess_column(num_cols, rlc_aliases$f)
  fm_guess <- rlc_guess_column(setdiff(num_cols, f_guess), rlc_aliases$fm)

  #RLC column: if no name matches and there is no text column, the file
  #is probably a single light curve; otherwise the user has to choose
  label_guess <- rlc_guess_column(all_cols, rlc_aliases$label)
  if (label_guess == "" && length(num_cols) == length(all_cols)){
    label_guess <- ".none"
  }

  shiny::updateSelectInput(session, "label_col",
    choices = c(no_choice,
                "(none: whole file is one curve)" = ".none",
                all_cols),
    selected = label_guess)

  shiny::updateSelectInput(session, "par_col",
    choices = c(no_choice, num_cols), selected = par_guess)

  shiny::updateSelectInput(session, "f",
    choices = c(no_choice, num_cols), selected = f_guess)

  shiny::updateSelectInput(session, "fm",
    choices = c(no_choice, num_cols), selected = fm_guess)
})


#light curve identifier of each row ("all" if the whole file is one curve)
labels <- shiny::reactive({

  d <- data()

  shiny::validate(shiny::need(
    input$label_col == ".none" || input$label_col %in% names(d),
    "Please choose the RLC column (or 'none' if the file has only one curve)"))

  if (input$label_col == ".none"){
    rep(".all", nrow(d))
  }else{
    as.character(d[[input$label_col]])
  }
})


#fill the RLC menu when the RLC column is chosen
shiny::observeEvent(list(input$file, input$label_col),{

  l <- tryCatch(labels(), error = function(e) NULL)
  shiny::req(l)

  if (input$label_col == ".none"){
    choices <- c("Whole file" = ".all")
  }else{
    choices <- unique(l[!is.na(l) & l != ""])
  }

  shiny::updateSelectInput(session, "rlc", label = "Select RLC",
                           choices = choices)
})


#when a RLC (or the PAR column) is selected, build one tickbox per light step
#the value of each box is the position of the light step in the RLC
#(so that repeated PAR values are still treated as separate steps)
#and the label is the PAR value
shiny::observeEvent(list(input$rlc, input$par_col, input$label_col),{

  l <- tryCatch(labels(), error = function(e) NULL)
  shiny::req(l, input$rlc %in% l, input$par_col %in% names(data()))

  par_values <- data()[[input$par_col]][l == input$rlc]

  shiny::updateCheckboxGroupInput(session, "light_steps",
                    choiceNames = as.character(par_values),
                    choiceValues = as.character(seq_along(par_values)),
                    selected = as.character(seq_along(par_values)))
  })


###############################################################################
#data preview
###############################################################################

output$preview_info <- shiny::renderText({

  imp <- imported()

  sep_name <- c("\t" = "tab", ";" = "semicolon (;)", "," = "comma (,)")[imp$delim]

  paste0(nrow(imp$data), " rows, ", ncol(imp$data), " columns\n",
         "Column separator: ", sep_name, "\n",
         "Decimal separator: ",
         if (imp$decimal_comma) "comma (converted to point)" else "point", "\n",
         "Numeric columns: ",
         paste(names(imp$data)[vapply(imp$data, is.numeric, logical(1))],
               collapse = ", "))
})

output$preview <- shiny::renderTable({
  utils::head(data(), 10)
})


###############################################################################
#data of the selected light curve (updated live, no fitting)
###############################################################################

curve <- shiny::reactive({

  d <- data()
  l <- labels()

  shiny::validate(
    shiny::need(input$par_col %in% names(d), "Please choose the PAR column"),
    shiny::need(input$f %in% names(d), "Please choose the F column"),
    shiny::need(input$fm %in% names(d), "Please choose the Fm column"),
    shiny::need(input$rlc %in% l, "Please select a RLC")
  )

  rows <- l == input$rlc

  light <- as.numeric(d[[input$par_col]][rows])
  f <- as.numeric(d[[input$f]][rows])
  fm <- as.numeric(d[[input$fm]][rows])

  shiny::validate(
    shiny::need(!all(is.na(light)) && !all(is.na(f)) && !all(is.na(fm)),
                "No valid dataset: PAR, F and Fm must be numeric")
  )

  #ticked light steps (TRUE/FALSE for each step of the curve)
  used <- seq_along(light) %in% as.integer(input$light_steps)

  #maximum Fm of the ticked steps (all steps if none is ticked)
  fm_max <- max(fm[if (any(used)) used else TRUE], na.rm = TRUE)

  eff <- (fm - f) / fm

  list(light = light,
       f = f,
       fm = fm,
       used = used,
       etr = light * eff,
       npq = npq(fm = fm_max, fmp = fm),
       ynpq = ynpq(f = f, fm = fm_max, fmp = fm))

})


###############################################################################
#fitting the models (only when the button is clicked)
###############################################################################

#fits are stored here; NULL means "no model fitted"
fits <- shiny::reactiveVal(NULL)

#any change to the data selection makes the previous fit invalid
shiny::observeEvent(list(input$file, input$label_col, input$par_col,
                         input$rlc, input$f, input$fm, input$light_steps),{
  fits(NULL)
}, ignoreInit = TRUE)


#helpers to keep only the fits that converged
#fit_etr returns a text message when it does not converge
etr_ok <- function(x) is.list(x) && !is.null(x$pred.par)
#fit_npq_2011 returns NULL when it does not converge
npq_ok <- function(x) is.list(x) && !is.null(x$predicted)

safe_fit <- function(expr) tryCatch(expr, error = function(e) NULL)


shiny::observeEvent(input$calculate,{

  cv <- tryCatch(curve(), error = function(e) NULL)
  shiny::req(cv)

  if (sum(cv$used) < 3){
    fits(NULL)
    shiny::showNotification("Please select at least 3 light levels",
                            type = "warning")
    return()
  }

  light <- cv$light[cv$used]
  etr <- cv$etr[cv$used]

  fit_P <- safe_fit(fit_etr(light = light, etr = etr,
                                      model = 'P', plots = FALSE))
  fit_EP <- safe_fit(fit_etr(light = light, etr = etr,
                                       model = 'EP', plots = FALSE))
  fit_npq <- safe_fit(fit_npq_2011(light = light,
                                             npq = cv$npq[cv$used],
                                             plots = FALSE))
  fit_ynpq <- safe_fit(fit_npq_2011(light = light,
                                              npq = cv$ynpq[cv$used],
                                              plots = FALSE))

  #alpha estimated using the slope of the first 3 selected points
  alpha_slope <- safe_fit(
    summary(stats::lm(etr[1:3] ~ light[1:3]))$coefficients[2])

  fits(list(
    P = if (etr_ok(fit_P)) fit_P else NULL,
    EP = if (etr_ok(fit_EP)) fit_EP else NULL,
    npq = if (npq_ok(fit_npq)) fit_npq else NULL,
    ynpq = if (npq_ok(fit_ynpq)) fit_ynpq else NULL,
    alpha_slope = alpha_slope
  ))

})

#fit currently selected with the radio buttons (NULL if not fitted)
etr_fit <- shiny::reactive({
  if (is.null(fits())) NULL else fits()[[input$etr_model]]
})

npq_fit <- shiny::reactive({
  if (is.null(fits())) NULL else fits()[[input$npq_type]]
})


###############################################################################
#plotting section
###############################################################################

#plots the light curve: ticked steps as black points, unticked steps as
#hollow grey points, and the model line or a "Model not fitted" message
plot_curve <- function(light, y, used, ylab, fit_x = NULL, fit_y = NULL){

  graphics::plot(light, y, type = "n", ylab = ylab, xlab = "Light", las = 1)

  graphics::points(light[used], y[used], pch = 21, bg = 1)

  if (any(!used)){
    graphics::points(light[!used], y[!used], pch = 21, bg = "white",
                     col = "grey50")
    graphics::legend("bottomright", legend = c("used", "not used"),
                     pch = 21, pt.bg = c(1, "white"), col = c(1, "grey50"),
                     bty = "n")
  }

  if (is.null(fit_x)){
    graphics::mtext("Model not fitted", side = 3, line = -1.5,
                    col = 2, font = 2)
  }else{
    graphics::lines(fit_x, fit_y, col = 2, lwd = 2)
  }
}


#ETR plot
output$plot_etr <- shiny::renderPlot({

  cv <- curve()
  fit <- etr_fit()

  plot_curve(cv$light, cv$etr, cv$used, ylab = "rETR",
             fit_x = fit$pred.par, fit_y = fit$pred)

},res = 96)


#NPQ plot
output$plot_npq <- shiny::renderPlot({

  cv <- curve()
  fit <- npq_fit()

  if (input$npq_type == "npq"){
    y <- cv$npq
    ylab <- "NPQ"
  }else{
    y <- cv$ynpq
    ylab <- "YNPQ"
  }

  plot_curve(cv$light, y, cv$used, ylab = ylab,
             fit_x = fit$predicted$par, fit_y = fit$predicted$npq)

},res = 96)


###############################################################################
#tables (only shown when the models are fitted)
###############################################################################

#Constructing a table for the NPQ outputs
output$npq_table <- shiny::renderTable({

  shiny::req(npq_fit())

  npq_fit()$parameters

  })


#Constructing a table for the ETR model outputs
output$table_etr <- shiny::renderTable({

  fit <- etr_fit()
  shiny::req(fit)

    # choosing the model
    if (input$etr_model == 'P'){

      table_etr <- data.frame(fit$alpha, fit$etrmax, fit$Ek, fit$beta, fit$ps)
      names(table_etr) <- c("Alpha","rETRmax","Ek", "Beta", "Ps")

    }else{

      table_etr <- data.frame(fit$alpha, fit$etrmax, fit$Ek, fit$eopt)
      names(table_etr) <- c("Alpha","rETRmax","Ek","Eopt")

    }

  table_etr

})


#Output table slope parameters
#Ek estimated as rETRmax of the selected model / alpha from the first 3 points
output$table_slope_parameters <- shiny::renderTable({

  fit <- etr_fit()
  shiny::req(fit, fits()$alpha_slope)

  alpha_lm_temp <- data.frame(fits()$alpha_slope,
                              fit$etrmax / fits()$alpha_slope)

  names(alpha_lm_temp) <- c("Alpha slope", "Ek slope")

  alpha_lm_temp

})

}


shiny::shinyApp(ui,server)


#TO DO list


#Add other EP model

#Add the possibility of fitting using YII?

}
