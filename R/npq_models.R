#' @title Function to fit NPQ model
#' @description Function to fit NPQ model using the Serodio & Lavaud 2011
#' @param light Light values (numeric vector)
#' @param npq npq values
#' @param starting_values The starting values of the parameters to fit (list).
#' @param num_obs Number of observations to fit the model (optional), if nothing is defined it will use the number of light steps
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @param limit_npq_max when deling with YNPQ data the NPQmax parameter should be limited to 1 by turning this parameter to TRUE (default is FALSE)
#' @return The function returns a list with: 1- par and NPQ values of the fitted model; 2 - fitted parameters.
#' @keywords external
#' @export
fit_npq_2011<-function(light,
                    npq,
                    starting_values=NULL,
                    num_obs=NULL,
                    plots=TRUE,
                    limit_npq_max=FALSE){

  app.data<-list()

  if (is.null(num_obs)==FALSE){
    light<-light[1:num_obs]
    npq<-npq[1:num_obs]
  }

  if (is.null(starting_values)==FALSE){
    starting_values=starting_values
  }

  if (is.null(starting_values)==TRUE){
      starting_values=list(NPQmax = 2,E50 = 200, hill = 1)
    }

  if (limit_npq_max==TRUE){
    my_upper <- c(1,Inf,Inf)
  }else{
    my_upper <- c(Inf,Inf,Inf)
  }


    my.res<-tryCatch({
      minpack.lm::  nlsLM(npq ~ NPQmax*(light^hill/((E50^hill + light^hill))) ,
                          start=starting_values, algorithm="port", trace=F,
                          lower = c(0,0,0), upper = my_upper, control=stats::nls.control(maxiter=1024))

    },error=function(e){NaN}
    )


    if(!is.na(my.res[1])){
      coefs <- stats::coef(my.res)
      my.npqmax <- as.numeric(coefs[1])
      my.E50 <- as.numeric(coefs[2])
      my.hill <- as.numeric(coefs[3])


      #run predict to get fit line
      new.dat <- data.frame(light=seq(0,max(light), by=1))
      #print(new.dat$light)
      pred <- stats::predict(my.res, newdata=new.dat)
      #print(pred)

      #make available
      app.data$npqmax <- my.npqmax
      app.data$E50 <- my.E50
      app.data$hill <- my.hill
      app.data$pred.par <- new.dat$light
      app.data$pred <- pred
    }


################################################################################

    if(is.na(my.res[1])==TRUE){
      print("The model did not converge, please choose different starting values,, or add more points to your dataset")
    }


################################################################################
##Plotting section
################################################################################

if (plots==TRUE){
      plot(light,npq,pch=21,bg=1,las=1,ylab = "NPQ", xlab = "Light")
      points(app.data$pred.par,app.data$pred,type='l',col=2,lwd=2)
    }

    #add the parameters to the plot if the model converged

    if(is.na(my.res[1])==FALSE){
      legend("bottomright",legend = c(paste("NPQmax = ",round(app.data$npqmax,2)),
                                      paste("E50 = ",round(app.data$E50,0)),
                                      paste("Sigmoidicity = ",round(app.data$hill,2)))
             ,bty="n",cex=0.9)
    }else{
      legend("bottomright",legend ="The model did not converge, please choose different starting values, or add more points to your dataset"
             ,bty="n",cex=0.9)
    }


    ################################################################################
    ##Return results
    ################################################################################
    if(is.na(my.res[1])==FALSE){

    predicted<-data.frame(new.dat$light,pred)
    names(predicted)<-c("par","npq")

    parameters<-data.frame(cbind(my.npqmax,my.E50,my.hill))
    names(parameters)<-c("NPQmax","E50","hill")

    output<-list(predicted,parameters)

    names(output)<-c("predicted","parameters")

    return(output)
    }



}


################################################################################
################################################################################

#' @title Function to fit NPQ model
#' @description Function to fit NPQ model using the Serodio & Lavaud 2021.
#' NPQ = NPQo * e^(-kd*E) + NPQm * (E_n/(E_n_50 + E_n)) + C
#NPQo measured at dark E=0
#Kd  coefficient of NPQ dissipation when exposed to light
#' @param light Light values (numeric vector)
#' @param npq npq values
#' @param starting_values The starting values of the parameters to fit (list with the starting values for NPQo=,NPQmax,E50,hill,Kd,C).
#' @param num_obs Number of observations to fit the model (optional), if nothing is defined it will use the number of light steps
#' @param limit_npq_max when deling with YNPQ data the NPQmax parameter should be limited to 1 by turning this parameter to TRUE (default is FALSE)
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @return The function returns a list with: 1- par and NPQ values of the fitted model; 2 - fitted parameters
#' @keywords external
#' @export

fit_npq_2021<-function(light,npq,
                       starting_values = NULL,
                       num_obs = NULL,
                       limit_npq_max = FALSE,
                       plots = TRUE
                       ){

  npq<-npq
  light<-light

  if (is.null(num_obs)==FALSE){
    light<-light[1:num_obs]
    npq<-npq[1:num_obs]
  }


  df<-data.frame(na.omit(cbind(light,npq)),stringsAsFactors=FALSE)

  if (is.null(starting_values)==TRUE){
    starting_values = list(NPQo=npq[1],NPQmax=max(na.omit(npq)),E50=max(light)/2,hill=1,Kd=0.01,C=0)
  }

  if (limit_npq_max == TRUE){
    my_upper <- c(Inf,1,Inf,Inf,Inf,Inf)
  }else{
    my_upper <- c(Inf,Inf,Inf,Inf,Inf,Inf)
  }

  npq_sim<-tryCatch({minpack.lm::nlsLM(npq~(NPQo*exp(-Kd*light) + NPQmax*(light^hill/(E50^hill + light^hill)) - C)
                             ,start=starting_values
                             , data=df
                             ,lower = c(0,0,0,0,0,-Inf),upper = my_upper)},error=function(e){NaN}
  )

  ################################################################################

  if(is.na(npq_sim[1])==TRUE){
    print("The model did not converge, please choose different starting values, or add more points to your dataset")
  }


  if(!is.na(npq_sim[1])){

  NPQo<-coef(npq_sim)[1]
  NPQmax<-coef(npq_sim)[2]
  E50<-coef(npq_sim)[3]
  hill<-coef(npq_sim)[4]
  Kd<-coef(npq_sim)[5]
  C<-coef(npq_sim)[6]


  x2<-seq(0,max(light),by=1)
  y2<- NPQo*exp(-Kd*x2) + NPQmax*(x2^hill/(E50^hill + x2^hill)) - C

  }

  ################################################################################

  if(is.na(npq_sim[1])==TRUE){
    print("No plot possible")
  }

  if(is.na(npq_sim[1]) == FALSE){

    if (plots == TRUE){
      plot(light,npq,xlab="Light",ylab="NPQ",las=1,pch=21,bg=1,cex=1,cex.lab=1,cex.axis=1)

      lines(x2,y2,col="blue",lwd=2)

      legend("bottomright",legend = c(paste("NPQo = ",round(NPQo,3)),
                                      paste("NPQmax = ",round(NPQmax,2)),
                                      paste("E50 = ",round(E50,0)),
                                      paste("Sigmoidicity = ",round(hill,2)),
                                      paste("Kd = ",round(Kd,4)),
                                      paste("C = ",round(C,4)))
             ,bty="n",cex=0.9)

    }


#Output section

  predicted<-data.frame(x2,y2)
  names(predicted)<-c("par","npq")

  parameters<-data.frame(cbind(NPQo,NPQmax,E50,hill,Kd,C))

  output<-list(predicted,parameters)

  names(output)<-c("predicted","parameters")

  }

  if(is.na(npq_sim[1])==FALSE){
    return(output)
  }




}








