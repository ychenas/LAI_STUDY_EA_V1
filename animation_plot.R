go.pdf <- FALSE
#copy default par()
par.def <- par()      

library("rgdal")
library("maptools")
library("raster")

# Generate a lancover map  overlap with TC tracks and effect size 
wrk.dir   <- c("/lfs/work/ychen/LAI_STUDY/")
track.dir <- c("/TRACK_DATA/")
lc.dir    <- c("LANDCOVER_DATA/")

# load the land cover map 
load(paste(wrk.dir,lc.dir,"2015.esa.landcover.rda",sep=""))

# variable name: esa.lc
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30
# set water(210) to NA
#generate the color palette 
lc.leg   <- read.csv("ESACCI-LC-Legend.csv",sep=";")
lc.col    <- rgb(lc.leg$R/255, lc.leg$G/255, lc.leg$B/255)
lc.breaks <- c(0, lc.leg$NB_LAB+0.1)

#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/data1/home/ychen/R/ne_110m_coastline/ne_110m_coastline.shp")

#convert arrary to raster obj for ploting
xymn <- c(105.0089, 150, 5, 49.99107)
lc.raster <- raster( x=t(esa.lc),
                         xmn=xymn[1], xmx=xymn[2], ymn=xymn[3], ymx=xymn[4],
                         crs=CRS("+proj=longlat +datum=WGS84"))


# convert the land cover array to raster object 
# load the tc track data 

track.data <- data.frame()

for (iyr in 2001:2017) {
    
   tmp.data <- read.csv(file=paste(wrk.dir,track.dir,"/hourly_data/wpb_",iyr,"_hourly.txt",sep=""))
   track.data <- rbind(tmp.data,track.data)

}  



# load the effect size data
library("fields")
# load data table 
load("./Rda//WP_20010110to20181231_table.combw_3_p_350.rda")

#load("WP_20010110to20140710_table.combw_3_p_350_TY_125D.rda")

#  Generate a figure print plot with x-axis: months; y-axis: latitude


#library(maptools)
library("raster")
library("rgdal")
# start to plot the figure print
#library("fields")

#dev.new()
#reset plot parameters
#par(par.def)
#set graph parameter
par(omi = c(1.5, 1.0, 0., 2 ))
par(mai = c(1, 1, 0, 2 ))
#asign lengend classes and texts
x.lab   <- c(0,10,20,30,40,50,60,70,80,90,
             100,110,120,130,140,150,160,170,180,190,200,210)
txt.lab <- c("No-Obs","Agr-Rf","Agr-Ir","AF-M1","AF-M2","For-BE","For-BD","For-NE","For-ND","For-Mix",
             "Mos-Tr","Mos-Hb","Shrub","Grass","Mosses","Spas-V","Fld-T1","Fld-T2","Fld-Hb","Urban","Bare","Water")

#create variables for plotting from table.comb
   x_pos <- as.numeric(table.comb$eve.lon.med.aff)   
   y_pos <- as.numeric(table.comb$eve.lat.med.aff)
eff_cex <-  sqrt(table.comb$eve.for.pix.aff/3.1416)/100 

#my.col    <- colorRampPalette(c("darkblue","blue","cyan","lightgray","yellow","pink","red","brown","brown"))(25)
eff.col    <- colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(16)
eff.breaks <- c( seq(-0.8, 0.8, length.out = 15))
eff.breaks <- round(eff.breaks,2)
#assign the color base on effect size
table.comb$col <-  colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(15)[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]

print(table.comb$col)


  #add event-base TC tracks
  for (ieve in 1:length(table.comb$date.time) ) {
#    for (ieve in 1:5 ) {
   
   # create a png device for generating the plots within an animation
   png(file=paste("./animation/",formatC(ieve, width = 4, format = "d", flag = "0"),".png",sep=""),units="px",width=1000,height=970)    
   #subset track.data to event table 
      #plot  raster land cover layer
      plot(lc.raster, col=lc.col, breaks=lc.breaks, legend=F, ylim=c(5,45),
           ylab="Latitude", xlab="Longitude", cex.lab=1.5, cex.axis=1.5,
           cex.main=2.0,main=paste("TC date: ",table.comb$date.time[ieve],sep=""))

      #add past events 
      for (ii in 1:ieve) { 
           eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ii]) &
                                          (track.data$CY == table.comb$eve.id[ii]) ))
           eve.table  <- subset (eve.table, eve.table$VMAX >=50) 
           # plot the lines 
           lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="longdash",col="white",lwd=0.5 )
           #add event-based TC location 
           points(x=x_pos[ii] , y=y_pos[ii], xlab="", ylab="",
                  col="gray", bg="NA", cex=eff_cex[ii], pch=21,lwd=0.5)
      }

      #add single event information  

      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      eve.table  <- subset (eve.table, eve.table$VMAX >=50) 
      # plot the tc track line
      lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="solid",col=table.comb$col[ieve],lwd=3.0 )
      #add TC date in the very first geo-location 
      text(x=eve.table$LonEW[1],y=eve.table$LatNS[1], labels=table.comb$date.time[ieve],pos=4,cex=1.0,font=2) 
      #add event-based TC location 
      points(x=x_pos[ieve] , y=y_pos[ieve], xlab="", ylab="",
             col=table.comb$col[ieve], bg="NA", cex=eff_cex[ieve]*2.0, pch=21,lwd=3.0)

      plot(coastlines, add=T) 
      par(bty = "n")
      #add legend
      plot(lc.raster, col=lc.col, breaks=lc.breaks, legend.only=T, horizontal=F, ylim=c(5,45),
           add=T, smallplot=c(0.93,0.95,0.2,0.9), legend.width=0.25, legend.shrick=0.5,
           axis.arg=list(at= x.lab- 2,
                        labels= txt.lab ,
                        las=2, cex.axis=0.8, mgp=c(5,0.5,0)),  
         legend.arg= list(text="Land cover type",las=3, side=2, font=2, line=0.2, cex=0.8) )
    #close png device 
   #close png device 
   dev.off() 
  print(paste("finished event:",ieve, "plots",sep=""))
   }
   #add coastlines
   #plot(coastlines, add=T) 
   #par(bty = "n")
   #add legend
   #plot(lc.raster, col=lc.col, breaks=lc.breaks, legend.only=T, horizontal=F,
   #     add=T, smallplot=c(0.93,0.96,0.2,0.9), legend.width=0.25, legend.shrick=0.5,
   #     axis.arg=list(at= x.lab- 2,
   #                   labels= txt.lab ,
   #                    las=2, cex.axis=0.8, mgp=c(5,0.5,0)),  
   #     legend.arg= list(text="Land cover type",las=3, side=2, font=2, line=0.2, cex=0.8) )
#reset plot parameters
par(par.def)
if(go.pdf) {dev.off()}
#par(oma = c(0.5, 0.5, 0.5, 0.5))
