#' @title Function to read reflectance spectra from the Ocean Optics machine (text format)
#' @description Function to read reflectance spectra from the Ocean Optics machine (text format)
#' @param filename path to the file to import (text vector)
#' @param plot parameter to activate plotting (default = TRUE)
#' @param indices parameter to activating calculating indices (default = FALSE)
#' @return a list with the following objects: filename- name of the imported file, metadata- all the metadata in the file header, datetime- time and date of the aquisition, all_data- a dataframe with the wavelengths and respective reflectance values.
#' @keywords external
#' @export


read_ooptics<-function(filename="./",plot=TRUE, indices=FALSE){

#read file
#the files can have different row numbers depending on the machine used
#so the function will look for markers to know how many lines to skip and when to stop

data_raw<-scan(filename,sep='\t',character(0),quiet=TRUE)

#a copy for the metadata
data_raw2<-data_raw

#sometimes the computer operating system uses commas insteat of points as a
#decimal separator, the next line replaces commas by points to fix that
data_raw<-gsub(",",".",data_raw)

#data_sep <- which(data_raw=="^>>>>>")
data_start <- grep(pattern="^>>>>>Begin",x=data_raw)+1

#some files have a marker at the end. others not...
data_end <- grep(pattern="^>>>>>End",x=data_raw)-1

#extracting just the wavelengths and reflectance
if (length(data_end)!=0){
  data_raw<-data_raw[data_start:data_end]
}else{ data_raw<-data_raw[data_start:length(data_raw)]}


#convert to numeric
data_raw<-as.numeric(data_raw)

#store it in a dataframe
all_data<-as.data.frame(matrix(data_raw,length(data_raw)/2,2,byrow = TRUE))
names(all_data)<-c("wl","reflect")


#group reflectance by wavelength
all_data$wl<-round(all_data$wl,0)
all_data<-aggregate(all_data$reflect,list(all_data$wl),mean)
all_data$wl<-as.numeric(as.character(all_data$Group.1))
all_data<-as.data.frame(cbind(all_data$wl,all_data$x))
names(all_data)<-c("wl","reflect")


if (plot==TRUE){
#the values at extremes are usually very noisy so it's useful to have
#  the min and max between 400 and 900
min_y<-min(all_data$reflect[all_data$wl>=400&all_data$wl<=900])
max_y<-max(all_data$reflect[all_data$wl>=400&all_data$wl<=900])
  plot(all_data$wl,all_data$reflect,xlab='wavelength',ylab = 'reflectance',
       ylim=c(min_y,max_y),pch=21,cex=0.7,bg=1)
}



#store filename
filename<-basename(filename)


#store all the metadata inside a dataframe

data_raw2<-data_raw2[1:(data_start-1)]

if (data_raw2[2]=="++++++++++++++++++++++++++++++++++++"){
  data_raw2<-data_raw2[3:(length(data_raw2)-1)]
}else{
  data_raw2<-data_raw2[2:(length(data_raw2)-1)]
}

data_raw2<-as.data.frame(data_raw2)

metadata<-tidyr::separate(data = data_raw2, col = data_raw2, into = c("parameter", "value"), sep = ": ")


#print(metadata)

#store date and time in separate object to make easier to produce time-series
#this seems to work in my MacBook but not on my Linux, confirm and debug
datetime<-lubridate::parse_date_time(metadata[1,2],orders = "abdTY")


###############################################################################
#optional steps for reflectance indices
###############################################################################

#use an external NDVI function
#if (indices=TRUE){
#  thejesus::ndvi(red = ,nir = )
#}


#use an external NDVI diatom function


###############################################################################
###############################################################################



#return list with results

output<-list()
output<-list(filename,metadata,datetime,all_data)
names(output)<-c("filename","metadata","datetime","all_data")
return(output)


}
