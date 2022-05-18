
#----------------------------------------------------------------------
# This script is designed to plot the result from lai_change_analysis.R
# First Date: 2020-01-21
# Author: Yi-Ying Chen yiyingchen@gate.sinica.edu.tw
#----------------------------------------------------------------------
 
wind.s <- c(3)
rain.s <- c(100)
#post days after TC event
pdays <-  c(60)
#offset days for LAI observation after TC event
offdays <- c(50)
#rain.s <- c(40, 80, 100, 200) 
# test cases 

#it1=90   #20010630
#it2=90
#it1=275 #20060820
#it2=275
#it1=205 #20040910 Japan Songda
#it2=205
#it1=203
#it2=203
#it1=382 #20090810 Morakot
#it2=382
##it1=528 #20130831
#it2=528
#
  library("plyr")
  library("fields")
  library("rgdal")
  library("maptools")
  library("raster")
  #load the coastlines data by readOGR function from sp package
  #coastlines <- readOGR("/data1/home/ychen/R/ne_10m_coastline/ne_10m_coastline.shp")
  #call the lai_dynamic function 
  source("lai_change_analysis.R")

for (iwind in wind.s) {

  for (irain in rain.s) {

  #call the lai_change_analyis.R 
  table.comb <- fun.dyn.lai(Ref.wind=iwind, Ref.rain=irain, 
                            target=c("WP"),track.dir="/TRACK_DATA_TS_125D/",wrk.dir="/work/ychen/LAI_STUDY/",
                            it1=36*3+1, it2=36*19,
                            mon.str=1, mon.end=12, pdays=pdays, offdays=offdays ) 

  #load the data
  #load("WP_20100110to20181231_table.combw_4_p_40.rda")
  pdf(file = paste("effect_size_w",iwind,"_p",irain,".pdf",sep=""),
      width = 7, # The width of the plot in inches
      height = 7) # The height of the plot in inches
  #set the plot layout
  layout(matrix(data=seq(1,4,1),nrow=2, ncol=2, byrow=TRUE))
  #plot the geolocation of TCs 
  plot(y=table.comb$eve.lat.med.aff,x=table.comb$eve.lon.med.aff,
     pch=16, col="red", xlim=c(100,150), ylim=c(10,50),
     xlab="longitude", ylab="latitude",main="Geolocation")
  #plot(coastlines, col="gray", add=TRUE)  

  #plot zonal mean of the effect size
  plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$eve.lat.med.aff)),
     horizontal=TRUE, col="red", ylim=c(-1.5,1.0), outline=FALSE,
     xlab="latitude", ylab="TC effect size",main="Zonal average")

  #plot seasonal change of effect size
  plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$date.mon)),
     col="red", outline=FALSE,ylim=c(-1.5,0.5),main="Seasonal change",
     xlab="months",ylab="TC effect size" )

  #plot histogram of effect size 
  hist(table.comb$eff.size.for, nclass=30, 
     col="red", xlim=c(-1.5,1.5),main="Histogram",
     xlab="TC effect size", ylab="event frequency")
   
  # copy current plot to pdf 
  #dev.copy2pdf(file=paste("effect_size_w",iwind,"_p",irain,".pdf",sep=""), width = 7, height = 5)
  dev.off()
  }# end of irain
}# end of iwind 


