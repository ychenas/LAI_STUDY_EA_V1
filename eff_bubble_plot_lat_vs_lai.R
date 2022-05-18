go.pdf <- F
#copy default par()
par.def <- par()      

library("rgdal")
library("maptools")
library("raster")

# Generate a lancover map  overlap with TC tracks and effect size 
wrk.dir   <- c("/work/ychen/LAI_STUDY/")
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
#lc.col    <- rgb(lc.leg$R/255, lc.leg$G/255, lc.leg$B/255)
#lc.breaks <- c(0, lc.leg$NB_LAB+0.1)


#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/data1/home/ychen/R/ne_10m_coastline/ne_10m_coastline.shp")

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
#load("./Rda//WP_20010110to20181231_table.combw_3_p_350.rda")

load("./Rda_Allruns/WP_20010110to20171231_table.combw_4_p_150_S_125Dpos_050_off_000_cut_0.rda")

#  Generate a figure print plot with x-axis: months; y-axis: latitude


# define array 

nlai = 12
nlat = 16
fig.arr <- array(NA, dim=c(nlai,nlat))

lai.txt <- c("2.5","2.75","3.0","3.25","3.5","3.75","4.0","4.25","4.5","4.75","5.0","5.25")
lat.txt <- c(seq(9,37,2))

#create day column
table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)

ini.lat = 8.
del.lat = 2.0 

for (ix in 1:12) {
    for (iy in 1:16) {

      lat1 <- ini.lat + del.lat * as.numeric(iy-1) 
      lat2 <- ini.lat + del.lat * as.numeric(iy)
      
      mon1 <- as.integer(ix)
      # subset the table to the condiction
      tmp.table <-  subset(table.comb, (table.comb$eve.lat.med.aff>lat1 ) & (table.comb$eve.lat.med.aff <=lat2) 
                                      & (table.comb$date.mon ==mon1 )   ) 
 
      #calculate the effect size in the table 
      tmp.eff.size <- mean(tmp.table$eff.size.for, na.rm=T)
 
      fig.arr[ix,iy] <- tmp.eff.size
                    
      print(paste("lat:",lat1,"to",lat2,"; ", "month:",mon1))
      print(tmp.table$eve.lat.med.aff)
      print(tmp.table$date.mon)
   }
}

#library(maptools)
library("raster")
library("rgdal")
# start to plot the figure print
#library("fields")

#my.col    <- colorRampPalette(c("darkblue","blue","cyan","lightgray","yellow","pink","red","brown","brown"))(25)
eff.col    <- colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(16)
eff.breaks <- c( seq(-1.25, 1.25, length.out = 15))
eff.breaks <- round(eff.breaks,2)

#assign the color base on effect size 
table.comb$col <-  colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(15)[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]

#convert array  to raster obj for ploting 
r_p <-raster( t(fig.arr[,nlat:1]) ,
               xmn=2, xmx=6,
               ymn=8, ymx=38, 
               crs=CRS("+proj=longlat +datum=WGS84") )

#plot(r_p, xlab="month", ylab="Latitude",cex.lab=1.5,
#    xlim=c(1,24),ylim=c(8,38),axes=FALSE,
#    legend.args=list(text="Eff. Size", side=3, font=2, line=0.2, cex=1.0),legend.width=2,
#    col=my.col,breaks=my.breaks )

x_pos <- as.numeric(table.comb$bef.for.lai.aff)
y_pos <- as.numeric(table.comb$eve.lat.med.aff)
eff_cex <- 1.5* sqrt(table.comb$eve.for.pix.aff/3.1416)/100 
#points(x=x_pos , y=y_pos, col=table.comb$col, bg="NA",  cex=my_cex, pch=21)
#mtext(text=c(paste(mon.txt[1:12])), side=1, line=-1.5, at=seq(1,24,2)+0.5, las=2, cex=0.8)
#mtext(text=c(paste(lat.txt[1:15])), side=2, line=-2, at=seq(9,37,2), las=2, cex=0.8)
#box()

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=9, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(x=c(2,6), y=c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='')
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

if (go.pdf) { pdf(file="acc.plot_b.pdf",width=8, height=8) }
#layoutsetting for combine plot 
par(oma = c(2, 0.5, 0.1, 0.1))
layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1, 1), c(1, 3), TRUE)
op = par(no.readonly = TRUE)

out_col <- colorRampPalette(c("lightgray","gray","black"))(length(x_pos))
 
#plot 2 scater plot of the effect size
par(mar = c(4, 4, 0, 0))
   plot(x=x_pos , y=y_pos, xlab="LAI_bef", ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
        col="black", bg=table.comb$col, xlim=c(2,6), ylim=c(8,41.5), cex=eff_cex, pch=21,lwd=1.5)
   grid()
#plot 4  averaged LAI pre-condition  vs  the effect size
par(mar = c(0, 4, 3, 0))
   #plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$date.mon)),
   #     col="lightgray", outline=FALSE,ylim=c(-0.8,0.8),main="",
   #     xlab="",ylab="Effect size",  xaxt = "n" )
   # 
   
   plot(y=table.comb$eff.size.for, x=table.comb$bef.for.lai.aff,

         col="white",ylim=c(-1.5,1.5),xlim=c(2,6),
         xlab="", ylab="Effect size",main="", xaxt = "n" )
 
   #get file-list 
   in.files <- list.files(path="./Rda_Allruns/", pattern="WP_20010110to20171231_table.combw_*")
   for (id in 1:length(in.files)  ) { 
   #for (id in 1:1) { 
   load(paste("./Rda_Allruns/",in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   
   # subset datataset (this is super critical!! we develop this QC flags)
   table.comb <- subset(table.comb, ((table.comb$eve.for.pix.aff)>5000 & as.integer(as.character(table.comb$date.time)) <= 20180101  )  )
   table.comb <- subset(table.comb, (abs(table.comb$qc1.score)<=0.2 &  abs(table.comb$qc2.score)>=0.2 ))

   sm.ln <- smooth.spline(y=table.comb$eff.size.for, x=table.comb$bef.for.lai.aff, df=10)
   lines(x=sm.ln$x, y=sm.ln$y, col = "darkgray",lwd=1,lty="dotted")
   text(x=sm.ln$x[1], y=sm.ln$y[1],id)
 
   #scatter.smooth(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30 ),
   #               col="lightgray", span=0.5, degree=1, 
   #               lpars=list(lwd=3,col="gray"), ylim=c(-0.5,0.5), xlab="",ylab="Effect size",  xaxt = "n")
   }
   #add sselected dataset 
   load("./Rda_Allruns/WP_20010110to20171231_table.combw_4_p_150_S_125Dpos_050_off_000_cut_0.rda")
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   #subset datataset
   table.comb <- subset(table.comb, ((table.comb$eve.for.pix.aff)>5000 & as.integer(as.character(table.comb$date.time)) <= 20180101  )  )
   table.comb <- subset(table.comb, (abs(table.comb$qc1.score)<=0.2 &  abs(table.comb$qc2.score)>=0.2 ))

   sm.ln <- smooth.spline(y=table.comb$eff.size.for,x=table.comb$bef.for.lai.aff, df=10)
   #lines(x=sm.ln$x, y=sm.ln$y, col = "black",lwd=3)
  #grid(col="gray")
   abline(h=0,col="black",lty="dashed")
#plot 1 zonal average of the effect size
par(mar = c(4, 0, 0, 0))
   #plot(table.comb$eff.size.for ~ as.factor(as.integer(table.comb$eve.lat.med.aff)),
   #     horizontal=TRUE, col="orange", ylim=c(-0.5,0.5), outline=FALSE,
   #     xlab="", ylab="Effect size",main="", yaxt = "n")
   
   plot(x=table.comb$eff.size.for, y=table.comb$eve.lat.med.aff,ylim=c(10,43), 
        col="white",xlim=c(-1.5,1.5) ,ylab="", xlab="Effect size",main="", yaxt = "n" )
   #scatter.smooth(x=table.comb$eff.size.for, y=table.comb$eve.lat.med.aff,
   #               col="gray", span=0.5, degree=1, xlab="", ylab="Effect size", 
   #               lpars=list(lwd=2,col="black"), xlim=c(-0.3,0.3))
   #grid(col="gray")
   for (id in  1:length(in.files)  ) {
   load(paste("./Rda_Allruns/",in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   #subset datataset
   table.comb <- subset(table.comb, ((table.comb$eve.for.pix.aff)>5000 & as.integer(as.character(table.comb$date.time)) <= 20180101  )  )
   table.comb <- subset(table.comb, (abs(table.comb$qc1.score)<=0.2 &  abs(table.comb$qc2.score)>=0.2 ))

  
   sm.ln <- smooth.spline(y=table.comb$eff.size.for,x=table.comb$eve.lat.med.aff, df=10)
   lines(x=sm.ln$y, y=sm.ln$x, col = "darkgray",lwd=1, lty="dotted")
   text(x=sm.ln$y[1], y=sm.ln$x[1],id)
   }


   #add the selected dataset 
   load("./Rda_Allruns/WP_20010110to20171231_table.combw_4_p_150_S_125Dpos_050_off_000_cut_0.rda")
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   #subset datataset
   table.comb <- subset(table.comb, ((table.comb$eve.for.pix.aff)>5000 & as.integer(as.character(table.comb$date.time)) <= 20180101  )  )
   table.comb <- subset(table.comb, (abs(table.comb$qc1.score)<=0.2 &  abs(table.comb$qc2.score)>=0.2 ))

   sm.ln <- smooth.spline(y=table.comb$eff.size.for,x=table.comb$eve.lat.med.aff, df=10)
   #lines(x=sm.ln$y, y=sm.ln$x, col = "black",lwd=3)
  
   abline(v=0,col="black",lty="dashed")
   
 
# plot 3 legend on the top-right 
par(mar = c(0.5, 5, 3.5, 2), xaxt = "n")
    color.bar(lut=eff.col,min=-1.5,title="Effect size scale")
#close devic()
if(go.pdf){dev.off()}
ld.go <- FALSE
if (ld.go) {
#open new device
#dev.new()
if(go.pdf) {pdf(file="eff_plot_sub_a.pdf",width=8, height=8) } else{dev.new() }
#reset plot parameters
par(par.def)
#set graph parameter
par(oma = c(1, 2.5, 4.0, 0.5))
par(mar = c(4, 4, 0, 0))
#asign lengend classes and texts
#x.lab   <- c(0,10,20,30,40,50,60,70,80,90,
#             100,110,120,130,140,150,160,170,180,190,200,210)
#txt.lab <- c("No-Obs","Agr-Rf","Agr-Ir","AF-M1","AF-M2","For-BE","For-BD","For-NE","For-ND","For-Mix",
#             "Mos-Tr","Mos-Hb","Shrub","Grass","Mosses","Spas-V","Fld-T1","Fld-T2","Fld-Hb","Urban","Bare","Water")
 x.lab   <- c(0, 40,90,120,180,190,200,210,220)
txt.lab  <- c("No-Obs","Agr-all","For-all","Grass","Bare","Urban","Water")
 
lc.col  <- c("gray","yellow","forestgreen","darkolivegreen","orange","brown","orange","white")

lc.breaks <- c(0, x.lab+0.5)

  #replace x_pos as longitude
   x_pos <- as.numeric(table.comb$eve.lon.med.aff)   

   #plot  raster land cover layer
   plot(lc.raster, col=lc.col, breaks=lc.breaks, legend=F, ylim=c(7,50),xlim=c(105,157), 
        ylab="Latitude", xlab="Longitude", cex.lab=1.5, cex.axis=1.5)
   #add event-base TC tracks
   for (ieve in 1:length(table.comb$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      eve.table  <- subset (eve.table, eve.table$VMAX >=70) 
      # plot the tracks
      lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="dashed",col=table.comb$col[ieve],lwd=0.5 )
      # print(eve.table)
      #add event-based TC location 
      points(x=x_pos[ieve] , y=y_pos[ieve], xlab="", ylab="",
         bg=table.comb$col[ieve], col="black", cex=0.8, pch=21,lwd=0.5)
 
      # plot the tracks
      #lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="dashed",col="gray",lwd=0.5 )
      # print(eve.table)
      #add event-based TC location 
      #points(x=x_pos[ieve] , y=y_pos[ieve], xlab="", ylab="",
      #   bg=table.comb$col[ieve], col="black", cex=0.8, pch=21,lwd=0.5)
  
 } 
 
   #add coastlines
   plot(coastlines,lwd=1.2, add=T) 
   par(bty = "n")
   #add legend
   #plot(lc.raster, col=lc.col, breaks=lc.breaks, legend.only=T, horizontal=F,
   #     add=T, smallplot=c(0.93,0.96,0.2,0.9), legend.width=0.25, legend.shrick=0.5,
   #     axis.arg=list(at= x.lab,
   #                   labels= txt.lab ,
   #                    las=2, cex.axis=0.8, mgp=c(5,0.5,0)),  
   #     legend.arg= list(text="Land cover type",las=3, side=2, font=2, line=0.2, cex=0.8) )
#reset plot parameters
par(par.def)
if(go.pdf) {dev.off()}
}
#par(oma = c(0.5, 0.5, 0.5, 0.5))
