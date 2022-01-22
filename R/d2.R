#' @title Function to calculate second derivatives from reflectance spectra
#' @description Function to calculate second derivatives from reflectance spectra
#' @param wl wavelengths
#' @param w reflectance values
#' @param knots smoothing parameter for the loess function (default = 0.08)
#' @param reflectance optional parameter to preview the effect of smoothing on the reflectance data (default=FALSE)
#' @param smooth_reflectance (default = 0.08)
#' @param plots parameter to turn off plotting, default TRUE
#' @return a list with two dataframes: one with the original data and another with the derivative values
#' @keywords external
#' @export

d2<-function (wl,w,knots=0.08,reflectance=FALSE, smooth_reflectance, plots=TRUE)
{

#Add an optional smoothing parameter for the reflectance data
#Add a way of looking just at the reflectance data to test the reflectance smoothing parameters

if (reflectance==TRUE){
  if (plots==TRUE){
plot(wl,w,ylab="reflectance",xlab="Wavelengths",
     pch=21,bg=1,type = 'l')
}
alisa<-loess(w~wl,span = smooth_reflectance)
pred<-predict(alisa,wl)

if (plots==TRUE){
points(wl,pred,type='l',col=2)
}

#assigning the smoothed reflectance values the new w
w<-pred

#calculating the derivatives

#calculating 1st derivative
di<-length(wl)
passo1<-(w[2:di]-w[1:(di-1)])/(wl[2:di]-wl[1:(di-1)])
pmax<-max(pretty(w))
pmin<-min(pretty(w))

#calculating 2nd derivative
x<-(wl[2:di]+wl[1:(di-1)])/2
y<-passo1
di<-length(x)
passo2<-(y[2:di]-y[1:(di-1)])/(x[2:di]-x[1:(di-1)])

#new x limits
x<-(x[2:di]+x[1:(di-1)])/2


#alisa<-smooth.spline(x,passo2,nknot=knots)
alisa<-loess(passo2~x,span=knots)
pred<-predict(alisa,x)

norm<-pred/w[-c(1,length(w))]


}


if (reflectance==FALSE){

#calculating 1st derivative
di<-length(wl)
passo1<-(w[2:di]-w[1:(di-1)])/(wl[2:di]-wl[1:(di-1)])
pmax<-max(pretty(w))
pmin<-min(pretty(w))

#calculating 2nd derivative
x<-(wl[2:di]+wl[1:(di-1)])/2
y<-passo1
di<-length(x)
passo2<-(y[2:di]-y[1:(di-1)])/(x[2:di]-x[1:(di-1)])

#new x limits
x<-(x[2:di]+x[1:(di-1)])/2


#alisa<-smooth.spline(x,passo2,nknot=knots)
alisa<-loess(passo2~x,span=knots)
pred<-predict(alisa,x)

norm<-pred/w[-c(1,length(w))]


}



###############################################################################
if (plots==TRUE){

  #plotting section
  plot(wl,w,type="l",ylim=c(pmin,pmax),
       xlab="Wavelengths",ylab="Relectance")
  limd<-c(min(pretty(norm)),max(pretty(norm)))
  smax<-max(pretty(norm))
  smin<-min(pretty(norm))
  lines(wl[-c(1,length(w))],(pmax-pmin)*((norm-smin)/(smax-smin))+pmin,col="red")

  #adding lines and peaks

  for (v in 1:(length(norm)-2))
  {
    if(norm[v+1]>norm[v]&norm[v+2]<norm[v+1])
    #if((norm[v+1]>norm[v]&norm[v+2]<norm[v+1]) || (norm[v+1]<norm[v]&norm[v+2]>norm[v+1]))

    {
      #abline(v=x[v+1])
      #segments(x[v],0,x[v],(pmax-pmin)*((norm[v+1]-smin)/(smax-smin))+pmin)
      text(x[v+1],(pmax-pmin)*((norm[v+1]-smin)/(smax-smin))+pmin,labels=round(x[v+1],0),pos=3,cex=0.8)

    }
  }

}
###############################################################################



###############################################################################
  #output section
  der2y<-norm
  der2x<-x

output<-list(der2x,der2y)

  names(output)<-c("der2x","der2y")

  return(output)

###############################################################################


}
