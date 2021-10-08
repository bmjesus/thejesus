#' @title Fit growth curve model to cell counts
#' @description Fits growth curve model (Gompertz, 1825) to cell per day and returns 
#' the parameters: A, umax and l. Where A is the maximum concentration of the biomass (e.g. cells . ml-1),
#'  μmax is maximum growth rate (e.g. cells . day-1) and l is the latency time (e.g. days)  
#' @param days numerical vector with the time information
#' @param cells numerical vector with the cell growth per time information
#' @param A number with an estimation of "A" starting value for the model fit
#' @param umax number with an estimation of "umax" starting value for the model fit
#' @param l number with an estimation of "l" starting value for the model fit
#' @return A list with 3 dataframes: parameters, model and original. The dataframe "parameters" contains the 
#' parameters estimated by the model; the dataframe "model" contains the growth data estimated by the model; 
#' and the dataframe "original" contains the values used to fit the model. 
#' @keywords external
#' @export
gompertz<-function(days,cells,A,umax,l){

  cells<-cells
  days<-days

temp<-minpack.lm::nlsLM(cells ~ A * exp(-exp(umax * (exp(1)/A) * days * (l-days)+1)), start = list("A" = A,"umax" = umax,"l" = l)
                  ,lower = c(0,0,0))  

out_A<-coefficients(temp)[1]
out_umax<-coefficients(temp)[2]
out_l<-coefficients(temp)[3]

parameters<-coefficients(temp)

new_x<-seq(0,max(days),by=0.1)
predicted<- out_A * exp(-exp(out_umax * (exp(1)/out_A) * new_x * (out_l-new_x)+1))

model<-data.frame(cbind(new_x,predicted))
names(model)<-c("time","predicted")

original<-data.frame(cbind(days,cells))
names(original)<-c("time","cells")

output<-list(parameters,model,original)
names(output)<-c("parameters","model","original")

################################################
#plotting section
################################################

plot(output$original$time,output$original$cells,ylab="cells",xlab="Time (days)",pch=21,bg=1)
points(output$model$time,output$model$predicted,col=2, type="l",lwd=2)
legend("topleft",legend=c(paste("A: ",round(out_A,0)),paste("umax: ", round(out_umax,0)), paste("Latency: ",round(out_l,2))), bty="n")


################################################


return(output)
                 
  
}


