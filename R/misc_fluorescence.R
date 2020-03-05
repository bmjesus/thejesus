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


