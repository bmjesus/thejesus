#' @title Function to apply the MPBOM
#' @description Applies the MPBOM and calculates the associated indices. It's based on
#' the paper: "Launeau, Patrick; Méléder, Vona; Verpoorter, Charles; Barillé, Laurent; Kazemipour-Ricci, Farzaneh; Giraud, Manuel; Jesus, Bruno; Le Menn, Erwan. 2018. "Microphytobenthos Biomass and Diversity Mapping at Different Spatial Scales with a Hyperspectral Optical Model" Remote Sens. 10, no. 5: 716. https://doi.org/10.3390/rs10050716 "
#' @param wl vector with the wavelength values
#' @param reflectance vector with the reflectance values
#' @param xmin value with the minimum wl
#' @param xmax value with the maximum wl
#' @return a list with the following objects:
#' @keywords external
#' @export

mpbom<-function(wl,reflectance,xmin,xmax){

  reflectance<-as.data.frame(cbind(wl,reflectance))
  names(reflectance)<-c("wl","reflectance")

  par(mfrow=c(3,1),oma=c(0,0,0,0))

  par(mar=c(0,4,0,1))
  plot(reflectance$wl,reflectance$reflectance,type='l',
       xlab="",ylab='reflectance',las=1,xaxt='n')
  lm_mpb<-lm(reflectance$reflectance[reflectance$wl>xmin&reflectance$wl<xmax]~
               reflectance$wl[reflectance$wl>xmin&reflectance$wl<xmax])
  Ra<-reflectance$reflectance
  #Rb<-lm_mpb$coefficients[2]
  Rb<<-lm_mpb$coefficients[1]+lm_mpb$coefficients[2]*reflectance$wl
  points(reflectance$wl,Rb,type='l',col=1,lty=2)
  points(reflectance$wl[reflectance$wl>xmin&reflectance$wl<xmax],
         predict(lm_mpb),type='l',col=2,lty=1)
  if (any((Ra/Rb<0)==TRUE)){
    trans<-rep(NA,length(reflectance$wl))
    alpha<-rep(NA,length(reflectance$wl))
    par(mar=c(0,4,0,1))
    plot(reflectance$wl, reflectance$reflectance,type='n')
    par(mar=c(4,4,0,1))
    plot(reflectance$wl,reflectance$reflectance,type='n')
    print("Impossible to fit MPBOM")
  }else{
    trans<-sqrt(Ra/Rb)
    par(mar=c(0,4,0,1))
    plot(reflectance$wl,trans,col=1,type='l',
         ylab='Transmitance',xlab='',las=1,xaxt='n')
    par(mar=c(4,4,0,1))
    alpha <- -log(trans^(1/3))
    plot(reflectance$wl,alpha,col=1,type='l',
         ylab='Alpha',xlab='wl',las=1)
    abline(h=0)

  }



  #building a dataframe with the output data

  output<- as.data.frame(cbind(reflectance$wl,reflectance$reflectance,trans,alpha))

  names(output)<-c("wl","refl","trans","abs")

#extracting absorption value at 675 nm

  abs675<-mean(output$abs[output$wl>=674&output$wl<=676])

output<-list(output,abs675)
return(output)


  }
