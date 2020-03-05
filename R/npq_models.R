#' @title Function to fit NPQ model
#' @description Function to fit NPQ model using the MODEL
#' @param light Light values (numeric vector)
#' @param npq npq values
#' @param starting_values the starting values of the parameters to fit (list).
#' @param num_obs Number of observations to fit the model (optional), if nothing is defined it will use the number of light steps
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @return The function returns a list with: 1- fitted parameters and respective; 2 - modeled light levels (numeric vector); 3- predicted NPQ values (numeric vector).
#' @keywords external
#' @export
fit_npq<-function(light,
                    npq,
                    starting_values=NULL,
                    num_obs=NULL,
                    plots=TRUE){

  app.data<-list()

  if (is.null(num_obs)==FALSE){
    light<-light[1:num_obs]
    npq<-npq[1:num_obs]
  }

  if (is.null(starting_values)==FALSE){
    starting_values=starting_values
  }



    if (is.null(starting_values)==TRUE){
      starting_values=list(npqmax = 2,E50 = 200, hill = 1)
    }
    my.res<-tryCatch({
      minpack.lm::  nlsLM(npq ~ npqmax*(light^hill/((E50^hill + light^hill))) ,
                          start=starting_values, algorithm="port", trace=F,
                          lower = c(0,0,0), control=stats::nls.control(maxiter=1024))

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
      print("The model did not converge, please choose different starting values")
    }


################################################################################
##Plotting section
################################################################################

if (plots==TRUE){
      plot(light,npq,pch=21,bg=1,las=1)
      points(app.data$pred.par,app.data$pred,type='l',col=2,lwd=2)
    }

    #add the parameters to the plot if the model converged

    if(is.na(my.res[1])==FALSE){
      legend("bottomright",legend = c(paste("NPQmax = ",round(app.data$npqmax,2)),
                                      paste("E50 = ",round(app.data$E50,0)),
                                      paste("Sigmoidicity = ",round(app.data$hill,2)))
             ,bty="n")
    }else{
      legend("bottomright",legend ="The model did not converge, please choose different starting values"
             ,bty="n",cex=0.7)
    }


    ################################################################################
    ##Return results
    ################################################################################

    if(is.na(my.res[1])==FALSE){
      return(app.data)
    }



}

