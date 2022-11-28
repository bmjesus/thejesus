#' @title plot with standard deviation bars
#' @description a customized plot function to plot x,y and SD bars. With the parameter "overlay = TRUE" it uses the function points to overall a data series.
#' @param x x
#' @param y y
#' @param sd numerical vector with the dispersion measurement to plot
#' @param sd_x numerical vector with the dispersion measurement in the x-axis to plot
#' @param overlay logical, if TRUE then the function "points" is used instead of "plot"
#' @param ylab ylab
#' @param xlab xlab
#' @param ylim  ylim
#' @param xlim xlim
#' @param pch type of symbol
#' @param bg color of the symbol background
#' @param col color of the symbol
#' @param lty type of line
#' @param type type of plot
#' @param xaxt xaxt
#' @param yaxt yaxt
#' @param cex cex
#' @param cex.axis cex.axis
#' @param cex.lab cex.lab
#' @param log log
#' @return plots a x,y plot with sd bars
#' @keywords external
#' @export
my_plot<-function(x,y,sd = NULL, sd_x = NULL, overlay = NULL,
                  ylab = NULL,xlab = NULL, ylim = NULL,xlim = NULL,
                  pch = NULL,bg = NULL, type = NULL, xaxt = NULL,
                  yaxt = NULL,lty = NULL, col = NULL,
                  cex = NULL, cex.axis = NULL, cex.lab = NULL, log = NULL){

  if (is.null(sd)==TRUE){
    sd <- 0
  }else{
    sd <- sd
  }

  if (is.null(sd_x)==TRUE){
    sd_x <- 0
  }else{
    sd_x <- sd_x
  }

  if (is.null(ylim)==TRUE){
    y_min<-min(c(y - sd))
    y_max<-max(c(y + sd))
    ylim = c(y_min,y_max)
  }else{
    ylim<-ylim
  }

  if (is.null(xlim)==TRUE){
    x_min<-min(c(x - sd))
    x_max<-max(c(x + sd))
    xlim = c(x_min,x_max)
  }else{
    xlim<-xlim
  }

  if (is.null(pch)==TRUE){
    pch<-21
  }else{
    pch<-pch
  }

  if (is.null(bg)==TRUE){
    bg<-1
  }else{
    bg<-bg
  }

  if (is.null(col)==TRUE){
    col <- 1
  }else{
    col <- col
  }



  if (is.null(lty)==TRUE){
    lty <- 1
  }else{
    lty <- lty
  }

  if (is.null(type)==TRUE){
    type<-"p"
  }else{
    type <- type
  }

  if (is.null(log)==FALSE){
    log <- log
  }


if (is.null(overlay)==TRUE){
      plot(x,y,ylab=ylab,xlab=xlab,las=1,ylim = ylim, xlim = xlim,
           pch = pch, bg = bg,type = type, col = col, lty = lty,
           xaxt = xaxt, yaxt = yaxt, cex = cex, cex.axis = cex.axis,
           cex.lab = cex.lab)
      segments(x,y-sd,x,y+sd,col = bg)
      segments(x-sd_x,y,x+sd_x,y,col = bg)
    }else{
      points(x,y,ylab=ylab,xlab=xlab,las=1,ylim = ylim, xlim = xlim,
             pch = pch, bg = bg, type = type, col = col, lty = lty,
             xaxt = xaxt, yaxt = yaxt, cex = cex , cex.axis = cex.axis,
             cex.lab = cex.lab)
      segments(x,y-sd,x,y+sd,col = bg)
      segments(x-sd_x,y,x+sd_x,y,col = bg)
    }





}

