#' @title Function to fit PI models using ETR data and extract fitting parameters
#' @description Function to fit PI models using ETR data and extract fitting parameters. The available models are: DESCRIPTION of the models and respective options
#' @param etr ETR values (numeric vector)
#' @param light Light values (numeric vector)
#' @param model parameter to select the model ("JP", "P", "EP")
#' @param starting_values the starting values of the parameters to fit (list).
#' @param num_obs Number of observations to fit the model (optional), if nothing is defined it will use the number of light steps
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @param  alpha_lm Optional to estimate alpha from the first 2 points of the light curve. Default=FALSE
#' @return The function returns a list with: 1- fitted parameters and respective; 2 - modeled light levels (numeric vector); 3- predicted ETR values (numeric vector).
#' @keywords external
#' @export

fit_etr<-function(light,etr,
                  model,
                  starting_values=NULL,
                  num_obs=NULL,
                  plots=TRUE,
                  alpha_lm=FALSE){

app.data<-list()

if (is.null(num_obs)==FALSE){
  light<-light[1:num_obs]
  etr<-etr[1:num_obs]
}

if (is.null(starting_values)==FALSE){
  starting_values=starting_values
}


#Jassy & Platt model
  #Jassby & Platt model selected
  if(model == 'JP'){
    if (is.null(starting_values)==TRUE){
      starting_values=list(alpha = 0.1, etrmax = 40)
    }
    my.res<-tryCatch({
      minpack.lm::nlsLM(etr ~ etrmax * tanh((alpha * light) / etrmax),
                        start=starting_values, algorithm="port", trace=F, lower = c(0,0),
                        control=stats::nls.control(maxiter=1024))
    },error=function(e){NaN}
    )


    if(!is.na(my.res[1])){
      coefs <- stats::coef(my.res)
      my.alpha <- as.numeric(coefs[1])
      my.etrmax <- as.numeric(coefs[2])

      #run predict to get fit line
      new.dat <- data.frame(light=seq(0,max(light), by=1))
      #print(new.dat$light)
      pred <- stats::predict(my.res, newdata=new.dat)
      #print(pred)

      #make available
      app.data$alpha <- my.alpha
      app.data$etrmax <- my.etrmax
      app.data$Ek <- my.etrmax/my.alpha
      app.data$pred.par <- new.dat$light
      app.data$pred <- pred
    }
  }#end of if


################################################################################

#Platt et al model
if(model == 'P'){
  if (is.null(starting_values)==TRUE){
    starting_values=list(alpha = 0.2, Ps = 2000, beta=0)
  }

  my.res <- tryCatch({
    minpack.lm::nlsLM(etr ~ Ps*(1-exp(-alpha*light/Ps))*exp(-beta*light/Ps),
                      start = starting_values, algorithm = "port",
                      trace = F, control = stats::nls.control(maxiter=1024),
                      lower = c(0,0,0))
  }, error = function(e) {NaN} )


  #extract parameters
  if(!is.na(my.res[1])){
    coefs <- stats::coef(my.res)

    my.alpha <- as.numeric(coefs[1])
    my.ps <- as.numeric(coefs[2])
    my.beta <- as.numeric(coefs[3])

    #run predict to get fit line
    new.dat <- data.frame(light=seq(0,max(light), by=1))
    #print(new.dat$light)
    pred <- stats::predict(my.res, newdata=new.dat)
    #print(pred)

    #make available
    app.data$alpha <- my.alpha
    print(my.alpha)
    app.data$ps <- my.ps
    app.data$beta <- my.beta
    app.data$etrmax <- my.ps*(my.alpha/(my.alpha+my.beta))*(my.beta/(my.alpha+my.beta))^(my.beta/my.alpha)
    app.data$Ek <- app.data$etrmax/my.alpha

    app.data$pred.par <- new.dat$light
    app.data$pred <- pred
  }

}
#end of platt model if

################################################################################
#Eilers Petters model

#based on alpha, ETRmax,
if(model == 'EP'){

  if (is.null(starting_values)==TRUE){
    starting_values=list(alpha = 0.4, etrmax = 40, Eopt=150)
  }

  my.res <- tryCatch({
    minpack.lm::nlsLM(etr ~ light/(light^2*(1/(alpha*Eopt^2))+(light/etrmax)-((2*light)/(alpha*Eopt))+(1/alpha)),
                      start=starting_values, algorithm="port", trace=F,
                      control=stats::nls.control(maxiter=1024),lower=c(0,0,0))
  },error=function(e){NaN}
  )



  if(!is.na(my.res[1])){
    coefs <- stats::coef(my.res)

    my.alpha <- as.numeric(coefs[1])
    my.etrmax <- as.numeric(coefs[2])
    my.eopt <- as.numeric(coefs[3])

    #run predict to get fit line
    new.dat <- data.frame(light=seq(0,max(light), by=1))
    #print(new.dat$light)
    pred <- stats::predict(my.res, newdata=new.dat)
    #print(pred)

    #make available
    app.data$alpha <- my.alpha
    app.data$etrmax <- my.etrmax
    app.data$eopt <- my.eopt
    app.data$Ek <- my.etrmax/my.alpha
    app.data$pred.par <- new.dat$light
    app.data$pred <- pred
  }
}

#based on a,b,c


################################################################################

#Webber?





################################################################################

if(is.na(my.res[1])==TRUE){
  print("The model did not converge, please choose different starting values")
}

#app.data<<-app.data



################################################################################
##Plotting section
################################################################################
if (plots==TRUE){
  plot(light,etr,pch=21,bg=1,las=1)
  points(app.data$pred.par,app.data$pred,type='l',col=2,lwd=2)
}

#add the parameters to the plot if the model converged

if(is.na(my.res[1])==FALSE){
  legend("bottomright",legend = c(paste("alpha = ",round(app.data$alpha,2)),
                                  paste("ETRmax = ",round(app.data$etrmax,0)),
                                  paste("Ek = ",round(app.data$Ek,0)))
         ,bty="n")
  print("The model converged")
}else{
  legend("bottomright",legend ="The model did not converge, please choose different starting values"
         ,bty="n",cex=0.7)
  }


################################################################################
##Return results
################################################################################

#adding original ETR and PAR values as a dataframe to the output list

original_values<-as.data.frame(cbind(light,etr))
app.data$original_values<-original_values

if(is.na(my.res[1])==FALSE){
return(app.data)
}else{return("The model did not converge, please choose different starting values")}

}
