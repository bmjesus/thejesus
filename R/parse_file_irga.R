#' @title Function to import IRGA files
#' @description Function to import and parse files created by the EGM-4 IRGA machine
#' @param file_path Path to the file
#' @return The function returns a dataframe with only the columns used in the calculations.
#' @keywords external
#' @export

#Function to import the CO2 IRGA file
parse_file_irga <- function(file_path) {
  # read all lines
  lines <- readLines(file_path, warn = FALSE)

  # remove empty lines
  lines <- lines[nzchar(trimws(lines))]

  # remove the BEGIN RECORD line if present
  lines <- lines[!grepl("^BEGIN RECORD$", trimws(lines))]

  # the first remaining line is now the true header
  header <- lines[1]

  # remaining lines are data/event lines
  data_lines <- lines[-1]

  # keep only true measurement rows
  data_lines <- data_lines[grepl(",\\s*M5,", data_lines)]

  # rebuild text
  txt <- paste(c(header, data_lines), collapse = "\n")

  # read into data frame
  df <- read.csv(
    text = txt,
    header = TRUE,
    stringsAsFactors = FALSE,
    strip.white = TRUE,
    check.names = FALSE
  )

  #keeping only the useful columns
  df <- df[,c(1,7,8,11)]

  #converting the first column in time format
  df[,1] <- as.POSIXct(df[,1], format = "%d/%m/%Y %H:%M:%S")

  #adding a column with absolute time in seconds
  df$cum_time <- as.numeric(df[,1] - df[1,1])

  # cleaning up a bit the column names

  names(df) <- c('date_time','co2','pressure','temperature','cum_time')

  return(df)
}
