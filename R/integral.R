#' @title Function to calculate the integral of derivative peaks
#' @description Function to calculate the integral of derivative peaks
#' @param x x variable, numerical vector
#' @param y y variable, numerical vector
#' @return  a number with the integral value
#' @keywords external
#' @export
integral<-function(x,y)

{

#to know how many steps to integrate
n<-length(x)

#to create the object
integral<-0

plot(x,y)

for (w in 1:(n-1))
{

integral<-integral+ ((x[w+1]-x[w])*y[w]+0.5*(y[w+1]-y[w])*(x[w+1]-x[w]))

polygon( c(x[w],x[w+1],x[w+1],x[w]),c(y[w],y[w+1],0,0),col="grey" )

}

return(integral)


}
