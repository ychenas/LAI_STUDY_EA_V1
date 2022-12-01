go.pdf <- T



#copy default par()
 go.pdf = T
 if(go.pdf) {pdf(file=paste("Fig2.pdf",sep=""),width=16, height=16) } else{ print("go x-windows") }
 #reset plot parameters
 #par(mfrow=c(2,2))
 par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5)      


# load library 
library("rgdal")
library("maptools")
library("mapplots")
library("raster")
library("plotrix")
library("fields")
#library(maptools)
all_data <- data.frame()
#set directories
wrk.dir   <- c("/lfs/home/ychen/LAI_STUDY_EAsia/")
track.dir <- c("/TRACK_DATA/")
lc.dir    <- c("/LANDCOVER_DATA/")
#table.dir <- c("/lfs/home/ychen/scripts/R/Rscripts/LAI_STUDY_EA/Rda_60_tc.occ/")
table.dir <- c("/lfs/home/ychen/scripts/R/Rscripts/LAI_STUDY_EA/spei02/BG_R1/")


# load the land cover map 
load("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda")

#convert arrary to raster obj for ploting
xymn <- c(90, 150, 0, 60)

#set NA = 999 for the ocean 
esa.lc[is.na(esa.lc)] <- 250
#apply the tc occurence map as mask
#esa.lc[arr.for.plot < 0.05] <- 0 

lc.raster <- raster( x=t(esa.lc),
                         xmn=xymn[1], xmx=xymn[2], ymn=xymn[3], ymx=xymn[4],
                         crs=CRS("+proj=longlat +datum=WGS84"))

# variable name: esa.lc  
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30 
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- NA 
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 40 ) ] <- 1
esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3 
# create the LC mask
lc.mask <- esa.lc
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- NA
lc.for.mask[lc.mask==2] <- 1


#only for tc accected area 
load("../TC_track_mask/tc.occ.1999to2018.rda")
#only accounted the tc affected area 
lc.for.mask[arr.for.plot<=0.1] <- NA 

tc.occ <- arr.for.plot
tc.occ[tc.occ==0] <- NA
tc.occ[lc.for.mask !=1] <- NA   

# variable name: esa.lc
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30
# set water(210) to NA
#generate the color palette 
lc.leg   <- read.csv("ESACCI-LC-Legend.csv",sep=";")
#lc.col    <- rgb(lc.leg$R/255, lc.leg$G/255, lc.leg$B/255)
#lc.breaks <- c(0, lc.leg$NB_LAB+0.1)

source("read_ncdf4.R")
nx=6722
ny=6722 

eff.acc.pos <- array( 0, dim=c( nx,ny) )  
eff.acc.neg <- array( 0, dim=c( nx,ny) )  
eff.acc.neu <- array( 0, dim=c( nx,ny) )  
#nav.lat <- array(NA,dim=c(nx,ny)) 
#nav.lon <- array(NA,dim=c(nx,ny)) 
    
lai.1 <- fun_read_nc( arg1="/lfs/home/ychen/LAI_STUDY_EAsia/LAI_DATA/c_gls_LAI_202003200000_PROBAV_V2_EAST_ASIA.nc", var_st=2 )
#for (ix in 1:nx) {
#      for (iy in 1:ny) {
#         nav.lat[ix,iy] <- lai.1$lat[iy]
#         nav.lon[ix,iy] <- lai.1$lon[ix] 
#      }
#}

#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_10m_coastline/ne_10m_coastline.shp")
#coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")
# Read this shape file with the rgdal library. 
#library(rgdal)
#land_mask <- readOGR( 
#  dsn= paste0("/lfs/home/ychen/GIS/world_border/") , 
#  layer="TM_WORLD_BORDERS_SIMPL-0.3",
#  verbose=FALSE
#)
#load river
#river <- readOGR("/lfs/home/ychen/GIS/ne_10m_rivers_lake/ne_110m_reviers_lake.shp")

#set work rda file
wrk.rda <-c(
	"EA_19990220to20181231_table.combw_8_p_0_2D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_10_p_0_3D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_12_p_0_4D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_0_p_60_2D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_0_p_80_3D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_0_p_100_4D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_12_p_100_4D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_10_p_80_3D_pos_060_off_000_cut_0_tc.occ.rda",
	"EA_19990220to20181231_table.combw_8_p_60_2D_pos_060_off_000_cut_0_tc.occ.rda")

# convert the land cover array to raster object 
# load the tc track data 

track.data <- data.frame()

for (iyr in 1999:2018) {
    
   tmp.data <- read.csv(file=paste(wrk.dir,track.dir,"/hourly_data/wpb_",iyr,"_hourly.txt",sep=""))
   track.data <- rbind(tmp.data,track.data)

}  

#load TC occurence
load("../TC_track_mask/tc.occ.1999to2018.rda")


#creat an array to calculate the forest area over each 10 by 10 degree 
xygd <- c(10)
nx10d <- (xymn[2]-xymn[1])/xygd
ny10d <- (xymn[4]-xymn[3])/xygd 
for.acc.10d <- array( 0, dim=c( nx10d,ny10d) )  
for.navlon.10d <- array( 0, dim=c( nx10d,ny10d) )  
for.navlat.10d <- array( 0, dim=c( nx10d,ny10d) )  

#array for accounting for positive or negtive effect size over 10 by 10 degree
eff.acc.10d.pos <- array( 0, dim=c( nx10d,ny10d) )  
eff.acc.10d.neg <- array( 0, dim=c( nx10d,ny10d) )  
eff.acc.10d.neu <- array( 0, dim=c( nx10d,ny10d) )  

tc.occ.10d.avg <- array( 0, dim=c( nx10d,ny10d) )

# find forest pixel for each grid for 10 by 10 degree

for (ix in nx10d:1) {
    for (iy in ny10d:1) {
       # find pixel for the seach area 
        cn.lat <-  xymn[4]- (xygd/2) - ((ix-1)*xygd)
        cn.lon <-  xymn[2]- (xygd/2) - ((iy-1)*xygd)
        dl.lat <- xygd/2
        dl.lon <- xygd/2
        id.x.10d <- which ( (lai.1$lat > cn.lat-dl.lat) & (lai.1$lat <= cn.lat+dl.lat) )
        id.y.10d <- which ( (lai.1$lon > cn.lon-dl.lon) & (lai.1$lon <= cn.lon+dl.lon) )
        #count total pixels 
        for.navlon.10d[ix,iy] <- cn.lon 
        for.navlat.10d[ix,iy] <- cn.lat 
        for.acc.10d[ix,iy] <-  sum(lc.for.mask[id.x.10d,id.y.10d],na.rm=T)
   } # 

} #  
for.acc.10d <- t(for.acc.10d[nx10d:1,ny10d:1])


#QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
qc1.set <- 0.5
qc2.set <- 0.1
area.set <- 0

fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.1, area.set=0) {
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      #QC1 difference of mean LAI
      #QC2 measurement uncertainty of ES  
      wrk.table <- subset(wrk.table, ( (abs(wrk.table$qc1.score)<=qc1.set) &  (abs(wrk.table$eff.size.for) >= abs(wrk.table$qc2.score))  ))
      return(wrk.table)
      }

get_table.qc1 <- function(wrk.table, qc1.set=0.5, qc2.set=0.1,area.set=0) {
      #get big difference
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      wrk.table <- subset(wrk.table, ( abs(wrk.table$qc1.score)>qc1.set)  )
      return(wrk.table)  
}

get_table.qc2 <- function(wrk.table, qc1.set=0.5, qc2.set=0.1, area.set=0) {
      #get neutral event
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      wrk.table <- subset(wrk.table, (abs(wrk.table$eff.size.for)<qc2.set) & ( abs(wrk.table$qc1.score)<qc1.set) )  
      return(wrk.table)  
}



#===================================
runs <- c(1,2,3,4,5,6,7,8,9)
run.txt <-c("a","b","c","d","e","f","g")
 
# loop for different run cases
for (irun in 2:2 ) {
# get run-case information 
run.case <- sub("_pos.*","",sub(".*table.comb","",wrk.rda[irun])) 
# get the data from Rda files 
load(paste(table.dir,"/",wrk.rda[irun], sep=""))

#remove no available pixels in affected area 
table.comb <- subset(table.comb, table.comb$eve.for.mws.aff > 0) 
 
# get how many events 
nstp <- length(table.comb$date.time)
#calculate how many events
neve <- 0
for (i in 2:nstp) {

   if (table.comb$eve.id[i] != table.comb$eve.id[i-1]) {
      neve =neve +1 
   }
}     


#get table neutral table  
table.neutral <- get_table.qc2( wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)

neve.neu <- length(table.neutral$eve.id) 

#gwt table for big LAI difference 
table.bigdiff <- get_table.qc1( wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)


#table quality check
table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 

neve.neq <- length(table.comb$eff.size.for[table.comb$eff.size.for<=0])
neve.pos <- length(table.comb$eff.size.for[table.comb$eff.size.for>0])


# Generate a figure print plot with x-axis: months; y-axis: latitude
# define array 

nmon = 12
nlat = 20
fig.arr <- array(NA, dim=c(nmon,nlat))

mon.txt <- c("2.5","2.75","3.0","3.25","3.5","3.75","4.0","4.25","4.5","4.75","5.0","5.25")
lat.txt <- c(seq(10,50,2))

#create day column
table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)

ini.lat = 10.
del.lat = 2. 

for (ix in 1:12) {
    for (iy in 1:20) {

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

#my.col    <- colorRampPalette(c("darkblue","blue","cyan","lightgray","yellow","pink","red","brown","brown"))(25)
#eff.col    <- colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(16)
#eff.col    <- colorRampPalette(c("forestgreen","green","lightgreen","white","gray","orange","brown"))(16)
#eff.col    <- rev(colorRampPalette(c("forestgreen","green","lightgreen","white","pink","red","brown"))(20) ) 
eff.col    <- colorRampPalette(c("orange4","orange2","lightblue","green","forestgreen"))(20)



eff.breaks <- c( seq(-1., 1., length.out = 19))
eff.breaks <- round(eff.breaks,2)

#assign the color base on effect size 
table.comb$col <-  eff.col[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]


#assign the maximum radius to the event
  for (ieve in 1:length(table.comb$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      table.comb$RAD[ieve] <- max(eve.table$RAD)
  }  


#convert array  to raster obj for ploting 
r_p <-raster( t(fig.arr[,nlat:1]) ,
               xmn=2, xmx=13,
               ymn=10, ymx=50, 
               crs=CRS("+proj=longlat +datum=WGS84") )

#plot(r_p, xlab="month", ylab="Latitude",cex.lab=1.5,
#    xlim=c(1,24),ylim=c(8,38),axes=FALSE,
#    legend.args=list(text="Eff. Size", side=3, font=2, line=0.2, cex=1.0),legend.width=2,
#    col=my.col,breaks=my.breaks )

#x_pos <- as.numeric(table.comb$bef.for.lai.aff)
x_pos <- (table.comb$date.mon+ as.integer(table.comb$date.day)/30 )
y_pos <- as.numeric(table.comb$eve.lat.med.aff)
#eff_cex <- 2.0* sqrt(table.comb$eve.for.pix.aff/3.1416)/100 
#eff_cex <- 2.5 
eff_cex <-  table.comb$RAD

# calculate the accumulative effect size for the grid of 10 by 10 degree 
for (ieve in 1:length(table.comb$RAD) ) {
     #allocate array
     tmp.arr <- array(0,dim=c(nx,ny))
     #for positive event
     if (table.comb$eff.size.for[ieve] > 0.18) {
        #find pixel for the events
        x.id <- which(  (lai.1$lon <= (table.comb$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.comb$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.comb$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.comb$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.pos <- tmp.arr + eff.acc.pos
     }#end if
      tmp.arr <- array(0,dim=c(nx,ny))
     #for negtive event 
     if (table.comb$eff.size.for[ieve] < -0.18) {
        x.id <- which(  (lai.1$lon <= (table.comb$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.comb$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.comb$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.comb$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.neg <- tmp.arr + eff.acc.neg
 
     }#end if

     #for neutral event abs(ES) <=0.2 
     if ( abs(table.comb$eff.size.for[ieve]) <= 0.18) {
        x.id <- which(  (lai.1$lon <= (table.comb$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.comb$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.comb$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.comb$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.neu <- tmp.arr + eff.acc.neu
 
     }#end if


}#end for

#work on neutral table by adding the event didn't pass the flag-2
#----> table.neutral <------ 
for (ieve in 1:length(table.neutral$eff.size.for) ) {

        tmp.arr <- array(0,dim=c(nx,ny))
         #for neutral event 
        x.id <- which(  (lai.1$lon <= (table.neutral$eve.lon.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lon >= (table.neutral$eve.lon.med.aff[ieve]-eff_cex[ieve]))  ) 
 
        y.id <- which(  (lai.1$lat <= (table.neutral$eve.lat.med.aff[ieve]+eff_cex[ieve])) &
                        (lai.1$lat >= (table.neutral$eve.lat.med.aff[ieve]-eff_cex[ieve]))  ) 
        tmp.arr[c(x.id),c(y.id)] <- 1    
        eff.acc.neu <- tmp.arr + eff.acc.neu

}


# aggregate to 10 degree by 10 degree
for (ix in 1:nx10d) {
   for (iy in 1:ny10d) {
    #set spacinf for 10 degree
    x.del =xygd/2
    y.del =xygd/2
    #find xid
    x.id <- which (  (lai.1$lon >= for.navlon.10d[ix,iy]- x.del) &
                     (lai.1$lon <  for.navlon.10d[ix,iy]+ x.del) )  
    #find yid    
    y.id <- which (  (lai.1$lat >= for.navlat.10d[ix,iy]- y.del) &
                     (lai.1$lat <  for.navlat.10d[ix,iy]+ y.del) )  
    #aggregate by sum 
    eff.acc.10d.neg[ix,iy] <- mean(eff.acc.neg[c(x.id),c(y.id)],na.rm=T) 
    eff.acc.10d.pos[ix,iy] <- mean(eff.acc.pos[c(x.id),c(y.id)],na.rm=T) 
    eff.acc.10d.neu[ix,iy] <- mean(eff.acc.neu[c(x.id),c(y.id)],na.rm=T) 
    #count average tc occurence at 10 by 10 
    tc.occ.10d.avg[ix,iy] <- mean( tc.occ[c(x.id), c(y.id)],na.rm=T)
   }#end for
}#end for 


#calculate the tc return frequency and Negtive ES LAI signal(damage)  


#library(fields)
#image.plot(t(tc.occ.10d.avg[6:1,6:1]),main="TC OCC AVG")
#dev.new()
#image.plot(t(eff.acc.10d.neg[6:1,6:1]), main="ES_LAI AVG")
#dev.new()


#eff.acc.10d.neg <- t(eff.acc.10d.neg[nx10d:1,ny10d:1])
#eff.acc.10d.pos <- t(eff.acc.10d.pos[nx10d:1,ny10d:1])
#eff.acc.10d.neu <- t(eff.acc.10d.neu[nx10d:1,ny10d:1])
#
#points(x=x_pos , y=y_pos, col=table.comb$col, bg="NA",  cex=my_cex, pch=21)
#mtext(text=c(paste(mon.txt[1:12])), side=1, line=-1.5, at=seq(1,24,2)+0.5, las=2, cex=0.8)
#mtext(text=c(paste(lat.txt[1:15])), side=2, line=-2, at=seq(9,37,2), las=2, cex=0.8)
#box()

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(x=c(2.,12.9), y=c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='')
    axis(2, round(ticks,digits=1), las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

ld.go <- T
if (ld.go) {
 #open new device
#par(mar = c(0, 0, 0, 0), mgp=c(3,1,0),xpd=FALSE)
 #set graph parameter
 #layout(matrix(c(1),1 , 1, byrow = TRUE), widths=c(1), heights=c(1), TRUE)
 par(oma = c(1,1,1,1), mar = c(5,6.5,1,1),mgp=c(4,2,0), xaxs="i", yaxs="i")
 # assign the table for plot 
 load(paste(table.dir,"/",wrk.rda[irun],sep=""))
 table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
 #call table quality check
 table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set) 
 #assign the color base on effect size 
 table.comb$col <-  eff.col[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]
 load("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda")
 #
  #asign lengend classes and texts
  x.lab   <- c(0, 40,90,120,180,190,200,210,230,250)
  txt.lab  <- c("No-Obs","Agr-all","For-all","Grass","Bare","Urban","Water","Sea")
  #lc.col  <- c("gray","yellow","forestgreen","darkolivegreen","orange","brown","orange","white")
  lc.col  <- c("NA","NA","NA","NA","NA","NA","NA","NA")
  lc.breaks <- c(0, x.lab+0.5)
  #replace x_pos as longitude
   x_pos <- as.numeric(table.comb$eve.lon.med.aff)   
  # set subplot margin
  par( usr=c(100,150,0,60)) #no buffer in x-y axis 
  plot(coastlines,col="white",lwd=0,ylim=c(0,60),xlim=c(100,150), ylab="", xlab="")
  #plot  raster land cover layer
  #plot legend
  #add coastaline
  par( usr=c(100,150,0,60)) #no buffer in x-y axis 

  plot(lc.raster, col=lc.col, breaks=lc.breaks, legend=F, ylim=c(0.,60.),xlim=c(100,150),add=T)
  par( usr=c(100,150,0,60), xpd=FALSE)  
  plot(coastlines, lwd=1.0,xlim=c(100,150), col="black", add=T) 

  grid(lwd=1,col="black")
 
# add Ticks Labels
  degs.N <- seq(0,60, length=7)
  axis(side = 2, at = degs.N, srt=0, las=1,
       labels = paste0(degs.N,"°","N") , tck = -0.02)
  mtext("Latitude" , side=2, line= 5.0, cex=2.5)
  deg.E <- c("","100","110","120","130","140","150")
  axis(side = 1, at = seq(90,150,length.out=7),
       labels = paste0(deg.E,"°","E"), tck = -0.02)
  mtext("Longitude" , side=1, line=4.0, cex=2.5)
  box() 
   #add land sea mask
   #add event-base TC tracks
   t.blue=rgb(red=0.1, green=.2, blue=.8, alpha=.5,  maxColorValue = 1)
   t.red=rgb(red=0.8, green=.2, blue=.1, alpha=.5,  maxColorValue = 1)
   t.green=rgb(red=0.1, green=.8, blue=.1, alpha=.5,  maxColorValue = 1)

   #write out track data table
   write.table(track.data, file="all_tc_track_table.txt",sep=",")
   write.table(table.comb, file="table_comb_3a_run.txt",sep=",")
   write.table(table.neutral, file="table_neutral_3a_run.txt",sep=",")



   #plot neutral-track lines
   for (ieve in 1:length(table.neutral$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.neutral$date.yr[ieve]) &
                                          (track.data$CY == table.neutral$eve.id[ieve]) ))
      eve.table  <- subset (eve.table, ((eve.table$VMAX) >=40 & (eve.table$LonEW <= 148 )& (eve.table$LonEW >= 102 )) ) 
      eve.stp <- length(eve.table$VMAX)
      eve.rad <- max(eve.table$RAD)
      print(paste("TC max RAD:",eve.rad,", ","TC_stps:",eve.stp,sep=""))
      # plot the tracks
      lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="solid",
             col="lightblue",
             lwd=ifelse( abs(table.neutral$eff.size.for[ieve])>0.4  ,0.5, 0.5))
   }
   #plot color-track lines
   for (ieve in 1:length(table.comb$date.time) ) {
      #subset track.data to event table 
      eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == table.comb$date.yr[ieve]) &
                                          (track.data$CY == table.comb$eve.id[ieve]) ))
      eve.table  <- subset (eve.table, ((eve.table$VMAX) >=40 & (eve.table$LonEW <= 148 )& (eve.table$LonEW >= 102 ))  ) 
      eve.stp <- length(eve.table$VMAX)
      eve.rad <- max(eve.table$RAD)
      print(paste("TC max RAD:",eve.rad,", ","TC_stps:",eve.stp,sep=""))
      # plot the tracks
      lines( x=eve.table$LonEW, y=eve.table$LatNS, lty="solid",
             col=table.comb$col[ieve],
             lwd=ifelse( abs(table.comb$eff.size.for[ieve])>0.4  ,1, 1))
   }
   #add coastline
   par( usr=c(100,150,0,60), xpd=FALSE)  
   plot(coastlines, lwd=0.5,xlim=c(100,150), col="black", add=T) 
 
   #add forest area that affected by TCs
    eff.acc.10d.neg[eff.acc.10d.neg==0]<- 0
    eff.acc.10d.pos[eff.acc.10d.pos==0]<- 0
    eff.acc.10d.neu[eff.acc.10d.neu==0]<- 0
    tot.a  <- 0.
    lab.text <- c("NA","NA","NA")
    for (ix in 1:nx10d) {
      for (iy in 1:ny10d) {
          # add pie chart 
            if ( (eff.acc.10d.neu[ix,iy] >= 0) & (for.acc.10d[ix,iy]/10000 > 0.1)) {
     
              print(paste("lon:",for.navlon.10d[ix,iy],"lat:",for.navlat.10d[ix,iy],sep=" "))
              print(paste("pos:",eff.acc.10d.pos[ix,iy]/10000,"neu:",eff.acc.10d.neu[ix,iy]/10000,
                            "neg:",eff.acc.10d.neg[ix,iy]/10000, sep=""))
               tot.a <- eff.acc.10d.pos[ix,iy] + eff.acc.10d.neu[ix,iy] +eff.acc.10d.neg[ix,iy]
               print( paste("Total Area:", tot.a,sep=""))
               if (tot.a > 0.005 ) {     
               lab.text[1] <- paste(as.character(format( eff.acc.10d.pos[ix,iy]/tot.a*100, digits=0,nsmall=0)),"%",sep="")  
               lab.text[2] <- format( eff.acc.10d.neu[ix,iy]/tot.a*100, digits=0,nsmall=0)
               if (lab.text[2] < 1) {lab.text[2] <- c("")} 
               else { lab.text[2] <- paste(lab.text[2],"%",sep="")}
  
               lab.text[3] <- paste(as.character(format( eff.acc.10d.neg[ix,iy]/tot.a*100, digits=0,nsmall=0)),"%",sep="")

                add.pie(z=c(eff.acc.10d.pos[ix,iy], eff.acc.10d.neu[ix,iy], eff.acc.10d.neg[ix,iy]),init.angle=90,
                        x=for.navlon.10d[ix,iy], y=for.navlat.10d[ix,iy],label.dist=1.2,cex=0.6,bg="red", 
                        col=c("forestgreen","lightblue","orange2"), 
                        radius=2.5, label=NA)
                area.txt <- round(for.acc.10d[ix,iy]/10000, digits=1) 
                #add information of percentage 
                #boxed.labels(bg= "white",border=NA,x=for.navlon.10d[ix,iy]+2.2, y=for.navlat.10d[ix,iy]-4.0,
                #             labels=paste(area.txt,"Mha",sep=""),cex=0.8,font=2)
                #add information of total tc events
                #boxed.labels(bg= "white",border=NA,x=for.navlon.10d[ix,iy]-2.5, y=for.navlat.10d[ix,iy]+4.0,
                #        labels=paste("TC:",format(tc.occ.10d.avg[ix,iy],digits=1,nsmall=1),sep=""),cex=0.8,font=2)
                } #if tot.a > 0
           } #if forest and signal are avariable   
       }# iy loop  
    }#ix loop

#  plot(coastlines,lwd=1.2, add=T) 
   add.pie(z=c(33,34,33), init.angle=90,
           x=109, y=56,cex=2.,label.dist=1.4, 
           col=c("forestgreen","lightblue","orange2"),
           radius=2.5, labels=c("Positive","Neutral","Negative"))
 
 # add color bar
   colorbar.plot(x=102.,y=49,adj.y=.5,adj.x=0.0, strip=seq(-1,1,0.1), strip.width=.02, strip.length=0.02*15,
                horizontal=T,col=eff.col)
   text(x=110,y=47.5,srt=0,"Effect Size",cex=2.5)
   text(x=110,y=51.0,srt=0,"-1.0    -0.5    0     0.5     1.0 ",cex=2.0) 
 
  if(go.pdf) {dev.off()}
} #end of ld_go


#library(ggpubr)
#x <- as.vector(tc.occ.10d.avg)
#y <- as.vector(eff.acc.10d.neu/(eff.acc.10d.neu+eff.acc.10d.neg+eff.acc.10d.pos))
#y[1] = x[1]
#y[y==0] <- NaN
#y[y>=0.8] <- NaN

#xy_data <- data.frame(occ=x, es.neu=y)
#all_data <- rbind(xy_data,all_data)
# make a plot 
#pdf(file=paste(run.case,"TCC_vs_ES_LAI(neu).pdf",sep=""),width=8, height=8)

#aa <- ggscatter(xy_data, x = "occ", y = "es.neu",ylim=c(0.,1.0), xlim=c(0,3.5), 
#          add = "reg.line", conf.int = TRUE,cex=3.0, 
#          cor.coef =F , cor.method = "pearson",
#          xlab = "TC Occurence", ylab = "Shared ES propotion (Neutral effects)")

#aa
#aa + stat_cor(method = "pearson", label.x = 1, label.y = 2)

#dev.off()


print(paste("Total events: ", neve,sep=""))
print(paste("Neutral events: ",neve.neu,sep=""))
print(paste("Positive events:",neve.pos,sep=""))
print(paste("Negative events:",neve.neq,sep=""))
 
} #end of irun case loop 
#===============================================================
