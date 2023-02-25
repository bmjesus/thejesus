#' @title Function to fit PI models using PSII efficiency data and extract fitting parameters
#' @description Function to fit PI models using PSII efficiency data and extract fitting parameters. They are based on the publication Modeling the irradiance dependency of the quantum efficiency of photosynthesis by Greg M. Silsbe and Jacco C. Kromkamp (Limnol. Oceanogr.: Methods 10, 2012, 645–652). The available models are: DESCRIPTION of the models and respective options
#' @param light Light values (numeric vector)
#' @param eff PSII quantum efficiency values (numeric vector)
#' @param model parameter to select the model ("JP", "P", "EP")
#' @return The function returns a list with: 1- fitted parameters and respective; 2 - modeled light levels (numeric vector); 3- predicted ETR values (numeric vector).
#' @keywords external
#' @export
fit_eff<-function(light,
                  eff,
                  model
                  )
{

#adding a small value if the curve starts at zero so that the fitting algorithm
#doesn't fail

  if(light[1] == 0){
    light <- light + 0.01
  }else{
    light <- light
  }

#constructing a df with the data
  df <- as.data.frame(cbind(light,eff))
  names(df) <- c('light','eff')

#creating a new x-axis for the prediction line
  x2 <- seq(min(light),max(light), by = 10)


#####################################################################
#Jassy & Platt model
  #Jassby & Platt model selected
  if(model == 'JP'){

    fit <- minpack.lm::nlsLM(eff ~ alpha * Ek * tanh(light * Ek^{-1}) * light^{-1},start = list(alpha = 0.5, Ek = 40),
                             trace=F, data = df)

    alpha <- summary(fit)$coefficients[1]
    Ek <- summary(fit)$coefficients[2]
    ETRmax <- Ek/alpha

#plotting section

    nf<-graphics::layout(matrix(1:2,nrow=2,ncol=1,byrow=TRUE))

    graphics::par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

    plot(light,eff,pch=21,xlab="Light",ylab="PSII",las=1,bg=1)
    lines(x2,alpha * Ek * tanh(x2 * Ek^{-1}) * x2^{-1},col="red",lwd=2)

    legend("topright", legend = c('model = JP' ,
                                  paste("alpha = ", round(alpha,2)),
                                  paste("ETRmax = ", round(ETRmax,0)),
                                  paste("Ek = ",round(Ek,0))))



    parameters <- as.data.frame(cbind(    round(alpha,2),
                                          round(ETRmax,0),
                                          round(Ek,0)))

    #plotting the ETR estimation based on the fitted model
    rETR <- light*eff
    plot(light,rETR,xlab="",ylab="rETR",las=1, pch=21, bg=1,
         xaxt ='n')

    #pred_etr <- ETRmax * tanh((alpha * x2) / ETRmax)

    pred_etr <- alpha * Ek * tanh(x2/Ek)

    lines(x2,pred_etr,col="red",lwd=2)




    # constructing the output
    names(parameters) <- c('alpha', 'ETRmax', 'Ek')
    model <- fit


  }#end of if JP


#####################################################################

  #Platt et al model
  else if (model == 'P'){


    # fit <- minpack.lm::nlsLM(eff~((Ps/light)*(1-exp((-light*alpha)/Ps))*exp((-light*beta)/Ps)),
    #                          start = list(alpha = 0.2, Ps = 1500, beta=0),trace = F, lower = c(0,0,0), data =df)

    fit <- minpack.lm::nlsLM(eff~ Ps * (1 - exp(-alpha * light * Ps^{-1}) * exp(-beta * light * Ps^{-1})) * light^{-1},
                             start = list(alpha = 0.2, Ps = 15000, beta=0),trace = F, lower = c(0,0,0), data =df)


    alpha <- summary(fit)$coefficients[1]
    Ps <- summary(fit)$coefficients[2]
    beta <- summary(fit)$coefficients[3]
    ETRmax <- Ps*(alpha/(alpha + beta))*(beta/(alpha + beta))^(beta/alpha)
    ETRmax <- Ek/alpha


#plotting section
nf<-graphics::layout(matrix(1:2,nrow=2,ncol=1,byrow=TRUE))

    graphics::par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

    plot(light,eff,pch=21,xlab="",ylab="PSII",las=1,bg=1, xaxt = 'n')
    #lines(x2,((Ps/x2)*(1-exp((-x2*alpha)/Ps))*exp((-x2*beta)))/x2,col="red",lwd=2)
    lines(x2,Ps * (1 - exp(-alpha * x2 * Ps^{-1}) * exp(-beta * x2 * Ps^{-1})) * x2^{-1},col="red",lwd=2)

    legend("topright", legend = c('model = P' ,
                                  paste("alpha = ", round(alpha,2)),
                                  paste("ETRmax = ", round(ETRmax,0)),
                                  paste("Ek = ",round(Ek,0))),
                                  paste("beta = ",round(beta,0)))



    parameters <- as.data.frame(cbind(    round(alpha,2),
                                          round(ETRmax,0),
                                          round(beta,0),
                                          round(Ek,0)))
    #plotting the ETR estimation based on the fitted model
    rETR <- light*eff
    plot(light,rETR,xlab="Ligth",ylab="rETR",las=1, pch=21, bg=1)

    pred_etr <- Ps * (1 - exp(-alpha * x2 * Ps^{-1}) * exp(-beta * x2 * Ps^{-1}))

    lines(x2,pred_etr,col="red",lwd=2)



    # constructing the output
    names(parameters) <- c('alpha', 'ETRmax' ,'beta', 'Ek')

    model <- fit


  }
  #end of if platt model

#####################################################################

#Eilers Petters model

  else if(model == 'EP'){

    fit <- minpack.lm::nlsLM(eff~1/(((alpha * Eopt^2)^{-1})*light^2+ (ETRmax^{-1} - 2 * (alpha * Eopt)^{-1})*light+(alpha^{-1})),data = df,
                             start=list(alpha = 0.4,Eopt = 90,ETRmax = 30),trace = F)

    alpha <- summary(fit)$coefficients[1]
    Eopt <- summary(fit)$coefficients[2]
    ETRmax <- summary(fit)$coefficients[3]
    Ek <- ETRmax/alpha

#Plotting section

 nf<-graphics::layout(matrix(1:2,nrow=2,ncol=1,byrow=TRUE))

    graphics::par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

    plot(light,eff,pch=21,xlab = "",ylab="PSII",las=1,bg=1,xaxt = 'n')

    lines(x2,1/(((alpha * Eopt^2)^{-1})*x2^2+ (ETRmax^{-1} - 2 * (alpha * Eopt)^{-1})*x2+(alpha^{-1})),col="red",lwd=2)

    legend("topright", legend = c('model = EP' ,
                                  paste("alpha = ", round(alpha,2)),
                                  paste("ETRmax = ", round(ETRmax,0)),
                                  paste("Eopt = ",round(Eopt,0)),
                                  paste("Ek = ",round(Ek,0))))


    parameters <- as.data.frame(cbind(    round(alpha,2),
                                          round(ETRmax,0),
                                          round(Eopt,0),
                                          round(Ek,0))
                                )

#plotting the ETR estimation based on the fitted model
    rETR <- light*eff
    plot(light,rETR,xlab="Light",ylab="rETR",las=1, pch=21, bg=1)

    pred_etr <- x2/(x2^2*(1/(alpha*Eopt^2))+(x2/ETRmax)-((2*x2)/(alpha*Eopt))+(1/alpha))

    lines(x2,pred_etr,col="red",lwd=2)



#output section
    names(parameters) <- c("alpha","ETRmax","Eopt","Ek")
    model <- fit
  }



#####################################################################

  #Webber?
  else if(model == 'W'){

    fit <- minpack.lm::nlsLM(eff ~ alpha * Ek * (1 - exp(-light * Ek)) * light^{-1},data = df,start=list(alpha = 0.4,Ek = 500),trace = F)


    alpha <- summary(fit)$coefficients[1]
    Ek <- summary(fit)$coefficients[2]
    ETRmax <- Ek*alpha

#plotting section
    nf<-graphics::layout(matrix(1:2,nrow=2,ncol=1,byrow=TRUE))

    graphics::par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)
    #
    plot(light,eff,pch=21,xlab="Light",ylab="PSII",las=1,bg=1)


    lines(x2, alpha * Ek * (1 - exp(-x2 * Ek)) * x2^{-1},col="red",lwd=2)

    legend("topright", legend = c('model = W' ,
                                  paste("alpha = ", round(alpha,2)),
                                  paste("ETRmax = ", round(ETRmax,0)),
                                  paste("Ek = ",round(Ek,0))))


    #plotting the ETR estimation based on the fitted model
     rETR <- light*eff
     plot(light,rETR,xlab="Light",ylab="rETR",las=1, pch=21, bg=1)

pred_etr <- alpha * Ek * (1 - exp(-x2 * Ek))


lines(x2,pred_etr,col="red",lwd=2)


#output
    parameters <- as.data.frame(cbind(    round(alpha,2),
                                          round(ETRmax,0),
                                          round(Ek,0))
    )
    model <- fit

}
#####################################################################

# if(is.na(my.res[1])==TRUE){
#     print("The model did not converge, please choose different starting values")
#   }

######################################################################
##Return results
######################################################################



  output <- list(parameters, model)

  names(output) <- c('parameters', 'model')

  return(output)





}
