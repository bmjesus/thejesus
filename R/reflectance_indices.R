#' @title Function to calculate normalized difference vegetation index (NDVI)
#' @description Function to calculate NDVI
#' @param red path to the file to import (text vector)
#' @param nir parameter to activate plotting (default = TRUE)

#' @return a list with the following objects: filename- name of the imported file, metadata- all the metadata in the file header, datetime- time and date of the aquisition, all_data- a dataframe with the wavelengths and respective reflectance values.
#' @keywords external
#' @export
ndvi<-function(red,nir){

ndvi<-(nir-red)/(nir+red)

return(ndvi)

}
