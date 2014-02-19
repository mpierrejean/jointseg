# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, min, max, ab=NULL,mar=c(4, 5, 2.5, 2), ...){
  
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
     ## check for additional function arguments
     if( length(list(...)) ){
       Lst <- list(...)
       if( !is.null(Lst$zlim) ){
         min <- Lst$zlim[1]
         max <- Lst$zlim[2]
       }
       if( !is.null(Lst$yLabels) ){
         yLabels <- c(Lst$yLabels)
       }
       if( !is.null(Lst$xLabels) ){
         xLabels <- c(Lst$xLabels)
       }
       if( !is.null(Lst$title) ){
         title <- Lst$title
       }
     }
     ## check for null values
     if( is.null(xLabels) ){
       xLabels <- c(1:ncol(x))
     }
     if( is.null(yLabels) ){
       yLabels <- c(1:nrow(x))
     }
     
     layout(matrix(data=c(1, 2), nrow=1, ncol=2), widths=c(5, 1.2), heights=c(1, 1))

     ## white and blue range 
     ColorRamp <- colorRampPalette(c("white", "darkblue"))(199)
     ColorLevels <- seq(min, max, length=length(ColorRamp))
     
     ## Reverse Y axis
     reverse <- nrow(x) : 1
     yLabels <- yLabels[reverse]
     x <- x[reverse, ]

     ## Data Map
     par(mar=mar)
     image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="", 
           ylab="", axes=FALSE, zlim=c(min, max))
     if( !is.null(title) ){
       title(main=title)
     }
     if(!is.null(ab)){abline(h=ab, col="black")}
     axis(BELOW <-1, at=1:length(xLabels), labels=xLabels, cex.axis=1.3, col.ticks=NULL, line=1, tick=FALSE)
     axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las=
          HORIZONTAL<-1, 
          cex.axis=1.5)

     ## Color Scale
     par(mar=c(3, 2.5, 2.5, 2))
     image(1,  ColorLevels, 
           matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1), 
           col=ColorRamp, 
           xlab="", ylab="",
           xaxt="n")
     
     layout(1)
   }
## ----- END plot function ----- #
