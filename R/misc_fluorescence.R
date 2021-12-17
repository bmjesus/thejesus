#series of functions to calculate secondary fluorescence derived parameters

######################################
#ETR related
######################################

#' @title ETR calculation from equation 1 Schuback et al. 2021
#' @description JPII = E x sigma PII x (Fv/Fm)−1 x ( Fq'/Fm')
#' @param light light values (numerical vector)
#' @param Fm maximum fluorescence yield in dark adapted state (single value)
#' @param Fmp maximum fluorescence yield in light adapted state (numerical vector)
#' @param Fo minimum fluorescence yield in dark adapted state (single value)
#' @param Fp F', minimum fluorescence yield in light adapted state (numerical vector)
#' @param sigmaPSII Sigma PSII (dark adapted) (single value)
#' @return numeric vector with the ETR values
#' @keywords external
#' @export
ETR_eq1<-function(light, sigmaPSII, Fo, Fm, Fp, Fmp){

  Fv <- (Fm - Fo)/Fm

  Fqp <- (Fmp - Fp)/Fmp

  JPII <- light * sigmaPSII * (Fv/Fm)^{-1} * (Fqp/Fmp)

  return(JPII)

}


#' @title ETR calculation from equation 2 Schuback et al. 2021
#' @description JSPII = E x Sigma PSII' x (Fq'/Fv'), Fq'/Fv' , calculated as (Fm' -F')/(Fm' - Fo')
#' @param light light values (numerical vector)
#' @param Fmp maximum fluorescence yield in light adapted state (numerical vector)
#' @param Fop minimum fluorescence yield in light adapted state with open reaction centers (numerical vector)
#' @param Fp F', minimum fluorescence yield in light adapted state (numerical vector)
#' @param sigmaPSIIp Sigma' PSII (light adapted) (numerical vector)
#' @return numeric vector with the ETR values
#' @keywords external
#' @export
ETR_eq2<-function(light, sigmaPSIIp, Fop, Fmp, Fp){

  Fvp <- (Fmp - Fp) / (Fmp - Fop)
  Fqp <- (Fmp - Fp)/Fmp
  JPII <- light * sigmaPSIIp * (Fqp/Fvp)
  return(JPII)

}



#' @title ETR calculation from equation 3 Schuback et al. 2021
#' @description JPSII = E x Sigma PSII' x (1 / (1 + (Sigma PSII' x E x tau1 )))
#' @param light light values (numerical vector)
#' @param tau1 tau1 reopening rate (numerical vector)
#' @param sigmaPSIIp Sigma' PSII (light adapted) (numerical vector)
#' @return numeric vector with the ETR values
#' @keywords external
#' @export
ETR_eq3<-function(light, sigmaPSIIp, tau1){

  JPII <- light * sigmaPSIIp * ((1/(1 + (sigmaPSIIp * light * tau1 ))))

  return(JPII)

}


#to be implemented, I don't understand the parameter Emax
#Equation 4 Schuback et al. 2021
#JPSII = 1/tau1 * [(E * (Fq'/Fm') / (Emax * (Fq'/Fm'_Emax)))                           ]




#Equation 5 Schuback et al. 2021
#JVII = E * Fq'/Fm' * [Ka * (Fm * Fo)/Fv ]

#' @title ETR calculation from equation 5 Schuback et al. 2021
#' @description JVII = E x Fq'/Fm' x (Ka x (Fm x Fo) / Fv)
#' @param light light values (numerical vector)
#' @param Fm maximum fluorescence yield in dark adapted state (single value)
#' @param Fmp maximum fluorescence yield in light adapted state (numerical vector)
#' @param Fo minimum fluorescence yield in dark adapted state (single value)
#' @param Fp F', minimum fluorescence yield in light adapted state (numerical vector)
#' @param Ka machine specific calibration factor (single value)
#' @return numeric vector with the ETR values
#' @keywords external
#' @export
ETR_eq5<-function(light, Fmp, Fo, Fm, Fp, Ka){

  Fqp <- (Fmp - Fp)/Fmp
  Fv <- (Fm - Fo) / Fm

  JPII <- light * Fqp/Fmp * ( Ka * (Fm * Fo)/Fv )

  return(JPII)

}















#' @title calculate NPQ
#' @description calculate NPQ
#' @param fm maximum fluorescence yield in dark adapted state
#' @param fmp maximum fluorescence yield in dark adapted state
#' @return numeric vector with the NPQ values
#' @keywords external
#' @export
npq<-function(fm,fmp){

npq <- (fm-fmp)/fmp
return(npq)

}

#' @title calculate YNPQ
#' @description calculate YNPQ
#' @param f fluorescence yield
#' @param fm maximum fluorescence yield in dark adapted state
#' @param fmp maximum fluorescence yield in dark adapted state
#' @return numeric vector with the YNPQ values
#' @keywords external
#' @export
ynpq <- function(f,fm,fmp){

ynpq <- (f/fmp) - (f/fm)

return(ynpq)

}


#' @title calculate YNO
#' @description calculate YNO
#' @param f fluorescence yield
#' @param fm maximum fluorescence yield in dark adapted state
#' @return numeric vector with the YNO values
#' @keywords external
#' @export
yno <- function(f,fm){

  yno <- (f/fm)

  return(yno)

}


#' @title calculate PSII quantum efficiency
#' @description calculate PSII quantum efficiency
#' @param f fluorescence yield
#' @param fm maximum fluorescence yield
#' @return numeric vector with the PSII quantum efficiency values
#' @keywords external
#' @export
eff <- function(f,fm){

  eff <- (fm-f)/fm

  return(eff)

}


