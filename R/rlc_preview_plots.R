#' @title function to fit several PI models
#' @description  Preview of RLC curves and associated parameters
#' @param x rlc object from the process_psi() function
#' @param starting_values list of starting values to fit the rETR model, needs to be changed to auto values
#' @param model not implemented yet
#' @return the fitted values for several PI models
#' @export
#' @keywords external
rlc_analysis<-function(x,starting_values=list(alpha = 0.5, Eopt = 500, etrmax = 100),
                       model = NULL
                       ){

###############################################################################
#variables and parameters section
###############################################################################

#calculating the variables

#light
par<-x$data$measuring_steps$par

#number of light steps
n_light<-length(x$data$measuring_steps$par)

#measurement type (light or after 2 seconds of dark, i.e., light, dark)
type<-c(rep("light",n_light),rep("dark",n_light))


#PSII quantum efficiency values
eff<-x$sti_parameters$psII_eff_sti

#relative ETR
retr <- par * eff


#NPQ
#(Fm-Fm')/Fm'

fm_light<-x$sti_parameters$fm_sti[1]
fm_dark<-x$sti_parameters$fm_sti[n_light + 1]
fmp_light<-x$sti_parameters$fm_sti[1:n_light]
fmp_dark<-x$sti_parameters$fm_sti[(1 + n_light):(n_light * 2)]

npq_light<-(fm_light-fmp_light)/fmp_light
npq_dark<-(fm_dark-fmp_dark)/fmp_dark

npq<-c( npq_light, npq_dark )


#same as before but using the maximum Fm observed in the series instead of
#the same value
#NPQm

fmm_light<-max(na.omit(fmp_light))
fmm_dark<-max(na.omit(fmp_dark))

npq_light_m<-(fmm_light-fmp_light)/fmp_light
npq_dark_m<-(fmm_dark-fmp_dark)/fmp_dark

npq_m<-c( npq_light_m, npq_dark_m )


#YNPQ
#(F/Fm')-(F/Fm)

f_light<-x$sti_parameters$fo_sti[1:n_light]
f_dark<-x$sti_parameters$fo_sti[(1 + n_light):(n_light * 2)]

ynpq_light<-(f_light/fmp_light)-(f_light/fm_light)
ynpq_dark<-(f_dark/fmp_dark)-(f_dark/fm_dark)

ynpq<-c( ynpq_light, ynpq_dark )

#same as before but using the maximum Fm observed in the series instead of
#the same value
#YNPQm

ynpqm_light<-(f_light/fmp_light)-(f_light/fmm_light)
ynpqm_dark<-(f_dark/fmp_dark)-(f_dark/fmm_dark)

ynpqm<-c( ynpqm_light, ynpqm_dark )

#parameters for the ETR equations

#Fv_Fm_light<-(fmm_light-f_light)/fmm_light
#Fv_Fm_dark<-(fmm_dark-f_dark)/fmm_dark

Fv_Fm_light<-eff[1:n_light]

Fv_Fm_dark<-eff[(1 + n_light):(n_light * 2)]


sigma_p_light<-x$sti_parameters$sigma_sti[1:n_light]
sigma_light<-sigma_p_light[1]

sigma_p_dark<-x$sti_parameters$sigma_sti[(1 + n_light):(n_light * 2)]
sigma_dark<-sigma_p_dark[1]

sigma<-c(sigma_p_light,sigma_p_dark)

sigma_se_light<-x$sti_parameters$sigma_se_sti[1:n_light]
sigma_se_dark<-x$sti_parameters$sigma_se_sti[(1 + n_light):(n_light * 2)]


tau_light<-x$str_parameters$tau1_str[1:n_light]
tau_dark<-x$str_parameters$tau1_str[(1 + n_light):(n_light * 2)]

tau<-c(tau_light,tau_dark)

tau_se_light<-x$str_parameters$tau1_se_str[1:n_light]
tau_se_dark<-x$str_parameters$tau1_se_str[(1 + n_light):(n_light * 2)]


rho_light<-x$sti_parameters$rho_sti[1:n_light]
rho_dark<-x$sti_parameters$rho_sti[(1 + n_light):(n_light * 2)]

rho_se_light<-x$sti_parameters$rho_se_sti[1:n_light]
rho_se_dark<-x$sti_parameters$rho_se_sti[(1 + n_light):(n_light * 2)]


###############################################################################
#fitting the models
###############################################################################

#relative ETR, EP model

retr_ep_light<-tryCatch({thejesus::fit_etr(par,retr[1:n_light],model="EP",
                                 starting_values = starting_values)},error=function(e){
  print("Error: could not fit model rETR light");
  return(NULL)
})



#print(retr_ep_light)


retr_ep_dark<-tryCatch({thejesus::fit_etr(par,retr[(1 + n_light):(n_light * 2)],model="EP",
                                  starting_values = starting_values)},error=function(e){
  print("Error: could not fit model rETR dark");
  return(NULL)
})




#section to select if the NPQ model should include "dark" NPQ or not
#the 2021 model does not fit well when NPQ starts at zero


if (ynpqm_light[1]==0){
  ynpq_m_2011_light<-tryCatch({thejesus::fit_npq_2011(par,ynpqm_light,
                                              starting_values = list("NPQmax"=2,"E50"=700,"hill"= 1))},error=function(e){
                                                print("Error: could not fit model YNPQm 2011 light");
                                                return(NULL)
                                              })
  ynpq_m_2011_dark<-tryCatch({thejesus::fit_npq_2011(par,ynpqm_dark,
                                           starting_values = list( "NPQmax"=2,"E50"=700,"hill"= 1))},error=function(e){
                                               print("Error: could not fit model YNPQm 2011 dark");
                                               return(NULL)
                                             })
ynpq_m_2021_light<-NULL

}

if (ynpqm_light[1]!=0){

ynpq_m_2011_light<-NULL

  #YNPQ_m

  ynpq_m_2021_light<-tryCatch({thejesus::fit_npq_2021(par,ynpqm_light
                                              )},error=function(e){
                                                print("Error: could not fit model YNPQm light");
                                                return(NULL)
                                              })
  ynpq_m_2021_dark<-tryCatch({thejesus::fit_npq_2021(par,ynpqm_dark)},error=function(e){
                                               print("Error: could not fit model YNPQm dark");
                                               return(NULL)
                                             })

}




###############################################################################
#plotting section
###############################################################################

par(mfrow=c(3,2),oma=c(4,4,3,4),mar=c(0,0,0,0),las=1)

####################
#PSII efficiency


plot(par,Fv_Fm_light,
     ylab = "", xaxt = 'n',xlab = "",type = 'b',xaxt = "n",
     pch = 21, bg = 1, col = 1,ylim=c(min(na.omit(c(Fv_Fm_light,Fv_Fm_dark))),max(na.omit(c(Fv_Fm_light,Fv_Fm_dark)))))
points(par,Fv_Fm_dark,
     ylab = "", xaxt = 'n',xlab = "",type = 'b',
     pch = 21, bg = 2, col = 2)
legend("bottomleft",legend=expression(phi[II]),bty = "n", cex=1.5)
legend("topright",legend = c("Light","Dark 2s"),
       pch = 21, pt.bg = c(1,2),col=c(1,2),bty="n",lty=1)

####################
#Sigma PSII



thejesus::my_plot(par,sigma_p_light,sd = sigma_se_light,
     ylab = "", xlab = "",type = 'b',xlim=c(min(par),max(par)),
     pch = 21, bg = 1, ylim=c(na.omit(min(c(sigma_p_light-na.omit(max(sigma_se_light)),sigma_p_dark-na.omit(max(sigma_se_dark))))),
                                      na.omit(max(c(sigma_p_light+na.omit(max(sigma_se_light)),sigma_p_dark+na.omit(max(sigma_se_dark)))))),
     xaxt = "n", yaxt = "n")
thejesus::my_plot(par,sigma_p_dark,sd = sigma_se_dark,
       ylab = "", xlab = "",type = 'b',
       pch = 21, bg = 2,overlay=TRUE)
axis(4)

legend("bottomleft",legend=expression(sigma[II]),bty = "n", cex=1.5)

#relative ETR
#light step
plot(par,retr[1:n_light],pch=21,bg=1,xaxt="n",xlab="",ylab = "rETR (a.u.)")
points(retr_ep_light$pred.par,retr_ep_light$pred,type="l",col=2,xaxt="n")
legend("bottomright",legend=c(paste("Alpha = ",round(retr_ep_light$alpha,2))
                              ,paste("ETRmax = ",round(retr_ep_light$etrmax,0))
                              ,paste("Ek = ",round(retr_ep_light$Ek,0))
                              ,paste("Eopt = ",round(retr_ep_light$eopt,0))
                              )
       ,bty="n")
legend("topleft",legend=expression("rETR"),bty = "n", cex=1)

#YNPQ model
if (ynpqm_light[1]==0){

  #print(ynpq_m_2021_light)
  if (is.null(ynpq_m_2011_light)==FALSE){
    #NPQ, model Serodio & Lavaud 2021
    #light step
    plot(par,ynpqm_light,pch=21,bg=1, ylab="YNPQ",xaxt = "n", yaxt = "n")
    points(ynpq_m_2011_light$predicted$par,ynpq_m_2011_light$predicted$npq
           ,type="l",col=2)
    axis(4)
    legend("topleft",legend=c(paste("NPQmax = ",round(ynpq_m_2011_light$parameters$NPQmax,2))
                              ,paste("E50 = ",round(ynpq_m_2011_light$parameters$E50,0))

                              ,paste("h = ",round(ynpq_m_2011_light$parameters$hill,2))

    )
    ,bty="n")

  }


}else{

  #print(ynpq_m_2021_light)
  if (is.null(ynpq_m_2021_light)==FALSE){
    #NPQ, model Serodio & Lavaud 2021
    #light step
    plot(par,ynpqm_light,pch=21,bg=1, ylab="YNPQ",xaxt = "n", yaxt = "n",ylim = c(0,max(na.omit(ynpqm_light))))
    points(ynpq_m_2021_light$predicted$par,ynpq_m_2021_light$predicted$npq
           ,type="l",col=2)
    axis(4)
    legend("topleft",legend=c(paste("NPQmax = ",round(ynpq_m_2021_light$parameters$NPQmax,2))
                              ,paste("E50 = ",round(ynpq_m_2021_light$parameters$E50,0))
                              ,paste("NPQo = ",round(ynpq_m_2021_light$parameters$NPQo,2))
                              ,paste("h = ",round(ynpq_m_2021_light$parameters$hill,2))
                              ,paste("Kd = ",round(ynpq_m_2021_light$parameters$Kd,3))
    )
    ,bty="n")

  }
}

legend("bottomright",legend=expression("YNPQ"),bty = "n", cex=1)

#########
#tau
thejesus::my_plot(x = par,y = tau_light, sd = tau_se_light,type = "b",
        ylim=c(min(c(tau_light-max(tau_se_light),tau_dark-max(tau_se_dark))),
                                                                  max(c(tau_light+max(tau_se_light),tau_dark+max(tau_se_dark))))
)
thejesus::my_plot(x = par,y = tau_dark,overlay = TRUE, sd = tau_se_dark,type="b",bg=2)

legend("bottomleft",legend=expression(tau[1]),bty = "n", cex=1.5)


#########
#rho
thejesus::my_plot(x = par,y = rho_light, sd = rho_se_light,type = "b",yaxt = "n",
        ylim=c(min(na.omit(c(rho_light-max(na.omit(rho_se_light)),rho_dark-max(na.omit(rho_se_dark))))),
               max(na.omit(c(rho_light+max(na.omit(rho_se_light)),rho_dark+max(na.omit(rho_se_dark))))))
        ,xlim=c(min(par),max(par)))

thejesus::my_plot(x = par,y = rho_dark,overlay = TRUE, sd = rho_se_dark,type="b",bg=2)
axis(4)

legend("bottomleft",legend=expression(rho),bty = "n", cex=1.5)








###############################################################################
#output section
###############################################################################


#bulding the output dataframe



if (is.null(ynpq_m_2011_light)==FALSE){

      npq_model_values<-as.data.frame(cbind(ynpq_m_2011_light$predicted$par,ynpq_m_2011_light$predicted$npq,ynpq_m_2011_dark$predicted$npq),stringsAsFactors=FALSE)

    names(npq_model_values)<-c("par","predicted_light","predicted_dark")

    model_param<-cbind.data.frame(ynpq_m_2011_light$parameters$NPQmax,
                                  ynpq_m_2011_light$parameters$E50,
                                  ynpq_m_2011_light$parameters$hill

                                 )

    names(model_param)<-c("NPQmax", "E50","h")
    sti_values<-data.frame(type, par, retr, npq, npq_m, ynpq,ynpqm,eff,sigma, tau,stringsAsFactors=FALSE)
    output<-list(sti_values,model_param,npq_model_values)
    names(output)<-c("sti_values","npq_model_param","npq_model_values")
}



if (is.null(ynpq_m_2021_light)==FALSE){

    npq_model_values<-as.data.frame(cbind(ynpq_m_2021_light$predicted$par,ynpq_m_2021_light$predicted$npq,ynpq_m_2021_dark$predicted$npq),stringsAsFactors=FALSE)

    names(npq_model_values)<-c("par","predicted_light","predicted_dark")

    model_param<-cbind.data.frame(ynpq_m_2021_light$parameters$NPQmax,
                                  ynpq_m_2021_light$parameters$E50,
                                  ynpq_m_2021_light$parameters$NPQo,
                                  ynpq_m_2021_light$parameters$hill,
                                  ynpq_m_2021_light$parameters$Kd,
                                  ynpq_m_2021_light$parameters$C)

    names(model_param)<-c("NPQmax", "E50","NPQo","h","Kd","C")

   # sti_values<-as.data.frame(cbind(type, par, retr, npq, npq_m, ynpq,ynpqm,eff,sigma, tau),stringsAsFactors=FALSE)
    sti_values<-data.frame(type, par, retr, npq, npq_m, ynpq,ynpqm,eff,sigma, tau,stringsAsFactors=FALSE)
    sti_values$par<-as.numeric(sti_values$par)

    output<-list(sti_values,model_param,npq_model_values)
    names(output)<-c("sti_values","npq_model_param","npq_model_values")

  }


if (is.null(ynpq_m_2021_light)==TRUE&is.null(ynpq_m_2011_light)==TRUE){

  sti_values<-data.frame(type, par, retr, npq, npq_m, ynpq,ynpqm,eff,sigma, tau,stringsAsFactors=FALSE)

  output<-sti_values

}


return(output)


}
