#' @title Function to fit PI models to O2 data from the FireSting (Pyroscience)
#' @description Function to fit PI models to O2 data from the FireSting (Pyroscience)
#' @param a File with the O2 data
#' @param light_time Dataframe with one column with light data, one column with the starting times and one column with the finishing times
#' @return
#' @keywords external
#' @export


light_curvesO2<-function(a="./",light_time,num_rows=13){

  #######################################################################
  #read the O2 file
  #######################################################################

  #import and skip the 13 lines of comments, data separated by Tab
  matriz_O2<-read.table(a,skip=num_rows,header=T,sep="\t")

  #replace names for something more obvious
  names (matriz_O2)<-c("date","time","time_s","comment","ch1_O2","ch2_O2","ch3_O2","ch4_O2",
                       "ch1_temp","ch2_temp","ch3_temp","ch4_temp","pressure","humidity",
                       "temp_probe","int_temp","analogue_in","ch1_raw","ch2_raw","ch3_raw","ch4_raw",
                       "ch1_signal","ch2_signal","ch3_signal","ch4_signal",
                       "ch1_light","ch2_light","ch3_light","ch4_light")
  #remove last column
  matriz_O2<-matriz_O2[,-30]

  #if the decimal separator is a ',' (e.g. French computer) then the data is not numeric
  # added a step to replace ',' by '.' and convert it to numeric format

  matriz_O2$ch1_O2<-gsub(",",".",  matriz_O2$ch1_O2)
  matriz_O2$time_s<-gsub(",",".",  matriz_O2$time_s)

  matriz_O2$ch1_O2<-as.numeric(matriz_O2$ch1_O2)
  matriz_O2$time_s<-as.numeric(matriz_O2$time_s)

  #str(matriz_O2)

  #######################################################################
  #plotting the data and the selected time ranges
  #######################################################################
  #closing ALL open devices
  grDevices::graphics.off()
  grDevices::dev.set()
  grDevices::dev.new()

  plot(matriz_O2$time_s,matriz_O2$ch1_O2,typ="l",ylab="O2 (umol/L)", xlab="Time (sec)",las=1 )

  #######################################################################
  #calculating the slopes and storing the rates
  #######################################################################

  my_list_lm<-list()
  my_list_x<-list()
  my_list_y<-list()
  my_list_predict<-list()
  my_list_rates<-list()

  for (i in 1:length(light_time[,1])){

    y<-with(matriz_O2,ch1_O2[time_s>=light_time[i,2]&time_s<=light_time[i,3]])
    x<-with(matriz_O2,time_s[time_s>=light_time[i,2]&time_s<=light_time[i,3]])

    lm_O2<-lm(y~x)


    points(x,predict(lm_O2),col="red",typ="l",lwd=2)

    my_list_lm[[i]]<-lm_O2
    my_list_x[[i]]<-x
    my_list_y[[i]]<-y
    my_list_predict[[i]]<-as.numeric(predict(lm_O2))
    my_list_rates[[i]]<-as.numeric(round(lm_O2$coefficients[2],5)*60)
  }

  #my_list_rates<<-my_list_rates


  #######################################################################
  #plot all in a matrix by counting the number of light steps
  #######################################################################

  #number of plots
  num_plots<-length(light_time[,1])





  grDevices::dev.new()
  nf<-graphics::layout(matrix(c(1:(ceiling(num_plots/3)*3)),nrow=ceiling(num_plots/3),ncol=3,byrow=TRUE))

  # nf<-graphics::layout(matrix(c(1:(ceiling(num_plots/4)*4)),nrow=ceiling(num_plots/4),ncol=4,byrow=TRUE))

  #nf<-layout(matrix(c(1:(ceiling(no_light_steps/4)*4)),nrow=ceiling(no_light_steps/4),ncol=4,byrow=TRUE))

  graphics::par(oma=c(4,4,3,4),mar=c(0,0,0,0), xpd=NA,tcl=-0.3,bg="white",cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

  for (i in 1:num_plots){
    graphics::plot(my_list_x[[i]],my_list_y[[i]],xlab="", ylab="", bty="l",xaxt='n',yaxt='n',pch=21,col="blue")
    graphics::points(my_list_x[[i]],my_list_predict[[i]],type="l",col="red",lwd=2)
    legend("bottomright",legend = my_list_rates[[i]])
  }

  #######################################################################
  #plot the rates vs light intensity and return it as an object
  #######################################################################

  my_list_rates<-as.numeric(my_list_rates)

  grDevices::dev.new()
  plot(light_time[,1],my_list_rates)


  #creating a dataframe object to return

  output<-as.data.frame(cbind(light_time[,1],my_list_rates))
  names(output)<-c("time","O2_rate")
  return(output)


  #fitting the Eilers Peters model
  #print(max(as.numeric(my_list_rates)))


  #this equation needs to have a component to inclued the respiration rate also. i.e. it needs to be adapted to these type of data
  #
  # fit_O2<-function(light,etr,
  #                  starting_values=NULL,
  #                  num_obs=NULL,
  #                  plots=TRUE,
  #                  alpha_lm=FALSE){
  #
  #   app.data<-list()
  #
  #   my.res <- tryCatch({
  #     minpack.lm::nlsLM(etr ~ Ps*(1-exp(-alpha*light/Ps))*exp(-beta*light/Ps)+Rs,
  #                       start = starting_values, algorithm = "port",
  #                       trace = FALSE, control = stats::nls.control(maxiter=1024))
  #   }, error = function(e) {NaN} )
  #
  #
  #   #,lower = c(0,0,0,Inf)
  #
  #   #extract parameters
  #   if(!is.na(my.res[1])){
  #     coefs <- stats::coef(my.res)
  #
  #     my.alpha <- as.numeric(coefs[1])
  #     my.ps <- as.numeric(coefs[2])
  #     my.beta <- as.numeric(coefs[3])
  #     my.Rs <- as.numeric(coefs[4])
  #
  #     #run predict to get fit line
  #     new.dat <- data.frame(light=seq(0,max(light), by=1))
  #     #print(new.dat$light)
  #     pred <- stats::predict(my.res, newdata=new.dat)
  #     #print(pred)
  #
  #     #make available
  #     app.data$alpha <- my.alpha
  #     #print(my.alpha)
  #     app.data$ps <- my.ps
  #     #print(my.ps)
  #     app.data$Rs <- my.Rs
  #     app.data$beta <- my.beta
  #     app.data$etrmax <- my.ps*(my.alpha/(my.alpha+my.beta))*(my.beta/(my.alpha+my.beta))^(my.beta/my.alpha)
  #     app.data$Ek <- app.data$etrmax/my.alpha
  #
  #     app.data$pred.par <- new.dat$light
  #     app.data$pred <- pred
  #   }
  #
  #
  #   plot(light,etr,pch=21,bg=1,las=1,xlab="light",ylab="Prod O2 (unit?)")
  #   points(app.data$pred.par,app.data$pred,type='l',col=2,lwd=2)
  #
  #   #add the parameters to the plot if the model converged
  #
  #   if(is.na(my.res[1])==FALSE){
  #     legend("bottomright",legend = c(paste("alpha = ",round(app.data$alpha,4)),
  #                                     paste("Pmax = ",round(app.data$etrmax,1)),
  #                                     paste("Ek = ",round(app.data$Ek,0)),
  #                                     paste("Rs = ",round(app.data$Rs,4)))
  #            ,bty="n")
  #   }else{
  #     legend("bottomright",legend ="The model did not converge, please choose different starting values"
  #            ,bty="n",cex=0.7)
  #   }
  #
  #
  #   return(app.data)
  #
  # }
  #
  # fit_O2(light = light,etr = my_list_rates,
  #        starting_values = list(
  #          alpha=(my_list_rates[2]/my_list_rates[1])/100,
  #          Ps=max(my_list_rates)*100,
  #          beta=0,Rs=my_list_rates[1]))
  #
  #


}


