#' @title Function to fit PI models using PSII efficiency data and extract fitting parameters
#' @description Function to fit PI models using PSII efficiency data and extract fitting parameters. They are based on the publication Modeling the irradiance dependency of the quantum efficiency of photosynthesis by Greg M. Silsbe and Jacco C. Kromkamp (Limnol. Oceanogr.: Methods 10, 2012, 645–652). The available models are: DESCRIPTION of the models and respective options
#' @param light Light values (numeric vector)
#' @param eff PSII quantum efficiency values (numeric vector)
#' @param model parameter to select the model ("JP", "P", "EP")
#' @param starting_values the starting values of the parameters to fit (list).
#' @param num_obs Number of observations to fit the model (optional), if nothing is defined it will use the number of light steps
#' @param  plots Turns on and off plotting of data and fitted curves (TRUE, FALSE). Default is TRUE.
#' @param  alpha_lm Optional to estimate alpha from the first 2 points of the light curve. Default=FALSE
#' @return The function returns a list with: 1- fitted parameters and respective; 2 - modeled light levels (numeric vector); 3- predicted ETR values (numeric vector).
#' @keywords external
#' @export
fit_eff<-function(light,
                  eff,
                  model,
                  starting_values=NULL,
                  num_obs=NULL,
                  plots=TRUE,
                  alpha_lm=FALSE)
{


  #Jassy & Platt model
  #Jassby & Platt model selected
  if(model == 'JP'){
    if (is.null(starting_values)==TRUE){
      starting_values=list(alpha = 0.3, Ek = 40)
    }
    my.res<-tryCatch({
      minpack.lm::nlsLM(eff ~ alpha * Ek * tanh(light * Ek^{-1}) * light^{-1},
                        start=starting_values, trace=F,  algorithm = "port", lower = c(0,0),
                        control=stats::nls.control(maxiter=1024))
    },error=function(e){NaN}
    )

   # alpha * Ek * tanh(light * Ek^{-1}) * light^{-1}




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
      minpack.lm::nlsLM(eff~((Ps/light)*(1-exp((-light*alpha)/Ps))*exp((-light*beta)/Ps)),
                        start = starting_values, algorithm = "port",
                        trace = F, control = stats::nls.control(maxiter=1024),
                        lower = c(0,0,0))
    }, error = function(e) {NaN} )


    #(Ps/light)*(1-exp((-light*alpha)/Ps))*exp((-light*beta)/Ps)

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

  #based on a, b, c
  if(model == 'EP'){

    if (is.null(starting_values)==TRUE){
      starting_values=list(a = 0.0001, b = 0.01, c = 1)
    }

    my.res <- tryCatch({
      minpack.lm::nlsLM(eff~1/(a*light^2+b*light+c),start=starting_values,
                        algorithm="port", trace=T,
                        control=stats::nls.control(maxiter=1024),lower=c(0,0,0))
    },error=function(e){NaN}
    )



    if(!is.na(my.res[1])){
      coefs <- stats::coef(my.res)

      my.a <- as.numeric(coefs[1])
      my.b <- as.numeric(coefs[2])
      my.c <- as.numeric(coefs[3])

      #transform the parameters a,b and c in alpha, etrmax, eopt
      my.alpha<-(1/my.c)
      my.etrmax<-1/(my.b+2*(my.a*my.c)^0.5)
      my.eopt<-(my.c/my.a)^0.5



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
    plot(light,eff,pch=21,bg=1,las=1)
    points(app.data$pred.par,app.data$pred,type='l',col=2,lwd=2)
  }

  #add the parameters to the plot if the model converged

  if(is.na(my.res[1])==FALSE){
    legend("topright",legend = c(paste("alpha = ",round(app.data$alpha,2)),
                                    paste("ETRmax = ",round(app.data$etrmax,0)),
                                    paste("Ek = ",round(app.data$Ek,0)))
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
