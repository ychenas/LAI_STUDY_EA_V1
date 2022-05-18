
#----------------------------------------------------------------------
# This script is designed to plot the result from lai_change_analysis.R
# First Date: 2020-01-21
# Author: Yi-Ying Chen yiyingchen@gate.sinica.edu.tw
#----------------------------------------------------------------------

fun_plot_eff_size <- function( in_file = "WP_20100110to20181231_table.combw_5_p_200.rda",
                               ld_pdf = FALSE, min_size= 10000 ) {
# in_file = "WP_20010110to20150120_table.combw_3_p_200_125D.rda"
# ld_pdf = FALSE
  load(file= in_file)
   
# subset table by removing small affected area cases 
  #min_size=50000
  table.comb <- subset(table.comb,table.comb$eve.for.pix.ref>min_size) 

  library("plyr")
  library("fields")
  library("rgdal")
  library("maptools")
  library("raster")
  #load the coastlines data by readOGR function from sp package
  coastlines <- readOGR("/data1/home/ychen/R/ne_110m_coastline/ne_110m_coastline.shp")  
  #call the lai_dynamic function 

  if (ld_pdf) {
    pdf(file = paste("effect_size_w",iwind,"_p",irain,".pdf",sep=""),
        width = 7, # The width of the plot in inches
        height = 7 ) # The height of the plot in inches
  }  
  #set the plot layout
  layout(matrix(data=seq(1,4,1),nrow=2, ncol=2, byrow=TRUE))
  #plot the geolocation of TCs 
  plot(y=table.comb$eve.lat.med.aff,x=table.comb$eve.lon.med.aff,
     pch=16, col="red", xlim=c(100,150), ylim=c(10,50),
     xlab="longitude", ylab="latitude",main="Geolocation")
  plot(coastlines, col="black", add=TRUE)  

  #plot zonal mean of the effect size
  plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$eve.lat.med.aff)),
     horizontal=TRUE, col="red", ylim=c(-1,1.0), outline=FALSE,
     xlab="latitude", ylab="TC effect size",main="Zonal average")

  #plot seasonal change of effect size
  plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$date.mon)),
     col="red", outline=FALSE,ylim=c(-1.,0.5),main="Seasonal change",
     xlab="months",ylab="TC effect size" )

  #plot histogram of effect size 
  hist(table.comb$eff.size.for,nclass=35,  
     col="red", xlim=c(-1.,1.),main="Histogram",
     xlab="TC effect size", ylab="event frequency")


  # Get the density estimate
  #dens=density(table.comb$eff.size.for)
  # Plot y-values scaled by number of observations against x values
  #lines(x=dens$x,
  #     y=length(table.comb$eff.size.for)*dens$y*dens$bw,
  #     type="l", col="red",
  #     xlab="TC effect size",ylab="Count estimate")

   
  # copy current plot to pdf 
  #dev.copy2pdf(file=paste("effect_size_w",iwind,"_p",irain,".pdf",sep=""), width = 7, height = 5)
  #close the pdf device
  if (ld_pdf) { dev.off() } 
}#end of function
