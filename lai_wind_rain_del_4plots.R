# -----------------------------------------------------------------
# This script/function analysis the SPOT-VGT LAI changes
# First Date: 2018-xx-xx
# Revised Date: 2019-12-05, 2020-01-21
# -----------------------------------------------------------------

#fun.dyn.lai <- function ( target=c("WP"), Ref.wind=3.0, Del.wind=0, Ref.rain=200, Del.rain=0, 
#			  it1=36*2+1, it2=36*20,track.dir="/TRACK_DATA_TY_125D/",wrk.dir="/lfs/work/ychen/LAI_STUDY/",
#			  mon.str=1, mon.end=12) { 


# set working files directory  
wrk.dir <- c("/lfs/work/ychen/LAI_STUDY/")
track.dir <- c("/TRACK_DATA_TS_123D/") 
tc.pot.txt <-  substr( track.dir, start=16, stop=19)

target = c("WP")
#
#it1=90
#it2=90
#Del.wind = 0.0 
#Del.rain = 0.0
#set a reference wind speed  for maximum wind speed for a reference area
Ref.wind = 5.0
Ref.rain = 350
#ref.wind = 4.0 
#ref.rain = 0
#target=c("TW"); it1=36*12;it2=36*18; mon.str=3;mon.end=11

#2001 to 2018
target=c("WP"); it1=36*2+1;it2=36*19; mon.str=1;mon.end=12
#target=c("WP"); it1=36*19+1;it2=36*20; mon.str=3;mon.end=11
#target=c("TEST")

#test run
#target=c("WP"); it1=36*10+1;it2=36*11; mon.str=3;mon.end=11
#

#it1=36*12+1;it2=36*14; mon.str=3;mon.end=11
#test cases
#
#it1=275 #20060820
#it2=275
#it1=528 #20130831 Sulik
#it2=528
#it1=205 #20040910
#it2=205
it1=382 #20090810 Morakot
it2=382
neve=1
#neve=floor((it2-it1+1)/36)*(mon.end-mon.str+1)*3+1

#plot accumulative map 
#pdf(file="acc.plot.pdf",width=8, height=8)
#par(lwd=0.5, mar=c(5,5,1.5,1) ) #,mpg=c(2,1,0)
#layout(matrix(data=c(1,2), nrow=1, ncol=2, byrow=TRUE))
library("fields")
# The land cover maps were aggregated to 1km from ESA 300m land cover product from 1992 to 2015
# The leaf area index maps were the popernious 1km land product from 1999 to 2018 
# The typhoon track data was compiled from Join Typhoon Warnning Center (JTWC) from 1945 to 2017 
# load the R data after several pre-pocessing

# 1. TC track occurence from 1948 to 2017 (70 years) the domain should match the LAI map over the West Pacific region (5040 * 5040)
#    load("TC_return_frq_1948to2017.rda")
#    variable name: land.frq"
#    use a 2 year return frequency to create the TC mask over WP Basin 
#    set.frq=0.5
#    land.frq[ land.frq/70 <  set.frq]  <- NA
#    land.frq[ land.frq/70 >= set.frq]  <- 1

# set the TC mask radius
#radius <- c("100km")

# 2. ESA land cover map from 2015, the domain shold be also matched the LAI map over WP region
load(paste(wrk.dir,"/LANDCOVER_DATA/2001.esa.landcover.rda",sep=""))
# variable name: esa.lc  
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30 
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- NA 
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 20 ) ] <- 1
esa.lc[ (esa.lc >= 30)  & (esa.lc <= 100) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3 

# create the LC mask
lc.mask <- esa.lc
# agricuture mask 
lc.agr.mask <- lc.mask
lc.agr.mask[lc.mask!=1] <- NA
lc.agr.mask[lc.mask==1] <- 1
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- NA
lc.for.mask[lc.mask==2] <- 1
# other mask
lc.oth.mask <- lc.mask
lc.oth.mask[lc.mask!=3] <- NA
lc.oth.mask[lc.mask==3] <- 1


# We commemt out this, because we have a new approach to define 
# TC affected (tc_af1 & tc_af2) area and Reference area (tc_ref) 

# 3. JTWC TC track occurance with 100km radius mapping 
#load("TC_return_frq_1948to2017.rda")
#create the analyis area mask
#land.frq[(land.frq <  60) & (land.frq >  10) ] <- -1
#land.frq[land.frq < 60 ] <-  NA
#land.frq[land.frq >= 60 ] <-  1


go.raster <- TRUE
if (go.raster) {
# 4. load coastal boundary file
library("plyr")
library("fields")
library("rgdal")
library("maptools")
library("raster")
#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/data1/home/ychen/R/ne_10m_coastline/ne_10m_coastline.shp")
#convert arrary to raster obj for ploting
#xymn <- c(105.0089, 150, 5, 49.99107)
#tc.pot.raster <- raster( x=t(tc.pot.mask), 
#			  xmn=xymn[1], xmx=xymn[2], ymn=xymn[3], ymx=xymn[4],
#			  crs=CRS("+proj=longlat +datum=WGS84"))
}

# plot the mask 
#library(fields)
#image.plot( (tc.mask*lc.mask)[,5040:1], main="The Study Area" ) 

# create a table summaries the result of condictional LAI change based on wind and rainfall scales
#  table.comb <- data.frame()

# get LAI file list for loading  LAI map data this could be a loop of time period 
  source("read_ncdf4.R")
# get file list from the folder 
  lai.fnames <- list.files(path=paste(wrk.dir,"/LAI_DATA/",sep=""), pattern="*.nc") 

  
  wrk.yr.old = 0
  ieve = 0

  # 
   nx=5040
   ny=5040 

   nav.lat <- array(NA,dim=c(nx,ny)) 
   nav.lon <- array(NA,dim=c(nx,ny)) 
# zonal avage variables 
#   gain.zonal.for <- array(NA,dim=c(neve,ny))
#   gray.zonal.for <- array(NA,dim=c(neve,ny))
#   loss.zonal.for <- array(NA,dim=c(neve,ny))
#   gain.zonal.for.pot <- array(NA,dim=c(neve,ny))
#   gray.zonal.for.pot <- array(NA,dim=c(neve,ny))
#   loss.zonal.for.pot <- array(NA,dim=c(neve,ny))
 
#   gain.zonal.agr <- array(NA,dim=c(neve,ny))
#   loss.zonal.agr <- array(NA,dim=c(neve,ny)) 
   #nav.lat  <- array(NA, dim=c(ny)) 
   
   lai.1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[1], sep=""), var_st=2 )
   for (ix in 1:nx) {
      for (iy in 1:ny) {
         nav.lat[ix,iy] <- lai.1$lat[iy]
         nav.lon[ix,iy] <- lai.1$lon[ix] 
      }
   }

for ( it in it1:it2) { 
   
  # get obstime information 
   # lai.p0$time <- format(as.Date(lai.p0$time,origin="1970-01-01"),format="%Y%m%d")
   # print( paste("working date:", lai.p0$time, sep="") )
   wrk.date <- substr( lai.fnames[it], start=11, stop=18)
   wrk.date.bef <- substr( lai.fnames[it-1], start=11, stop=20)
   wrk.date.crt <- substr( lai.fnames[it], start=11, stop=20)
   wrk.yr  <- as.integer(substr(wrk.date, start=1, stop=4) )
   wrk.mon <- as.integer(substr(wrk.date, start=5, stop=6) ) 
   print( paste("working date:", wrk.date, sep="") )

   #load hourly track data
    hourly.track <- read.csv(file=paste("/lfs/work/ychen/LAI_STUDY/TRACK_DATA/hourly_data/wpb_",wrk.yr,"_hourly.txt",sep=""))


       # 3. JTWC TC track mask with radius mapping 
       # 
       # Define potential TC area as affected area 
       # This part is curical that we update the TC affected and reference area every 10 days
       #
       # update the tc track data 
       # load( paste(wrk.yr,"_TC_mask_",radius,".rda",sep="") )
       # tc.mask <- track.mask # loading from loading TC_mask.rda
       # load tc occurance map from 1999-2017 tc track data analysis as background reference if there was no tc activities   
       #load( paste(wrk.dir,"/TRACK_DATA/","TC_return_frq_1948to2017.rda",sep="") )
       #tc.bg <- land.frq/(2017-1948+1)
       # set a five year return frequence as back ground 
       #tc.bg[tc.bg >= 0.2 ] <- 1
       #tc.bg[tc.bg < 0.2 ]  <- NA
       #track.dir <- c("/TRACK_DATA/") 
       # load mask files for TC affect and reference area
       #tc.af1 <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_af1_",wrk.date,"00.nc",sep=""), var_st=1)
       #tc.af2 <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_af2_",wrk.date,"00.nc",sep=""), var_st=1)
       tc.pot <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_pot_",wrk.date,"00.nc",sep=""), var_st=1)
#
       tc.aff.mask <- tc.pot$tc_pot
       tc.aff.mask[tc.aff.mask==0] <- NA 
       
       if(  (length(which(!is.na(tc.aff.mask*lc.for.mask)))  < 5000)  ) {

       # go to the next step for the it loop because of s,all sample size pixels below 20   
         print( "Go to next it/time step!") 
         next
       # exit it loop
       } 
       # find tc.event id within a 10 days mask
       tc.eve <- levels(factor(tc.aff.mask))  
       print(paste("There is(are) ", length(tc.eve), " TC event(s) within 10 days!",sep=""))

       #load reference mask 
       #tc.ref <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_ref_",wrk.date,"00.nc",sep=""), var_st=1)
       #tc.ref.mask <- tc.ref$tc_ref  
       #tc.ref.mask[tc.ref.mask==0] <- NA
      
       # load maximum windspeed & accumulative rainfall data over 10 days from ERA5 Land hourly reanalysis  
       tc.mws <- fun_read_nc(arg1=paste(wrk.dir,"/WIND_DATA/","max_wind_", wrk.date,"00.nc",sep=""), var_st=1)   
       tc.acf <- fun_read_nc(arg1=paste(wrk.dir,"/PRECIP_DATA/","acc_rainf_", wrk.date,"00.nc",sep=""), var_st=1)   
     

       #tc.ref.mask[tc.acf$acc_rainf > ref.rain] <- NA  
       #tc.aff.mask[tc.acf$acc_rainf < ref.rain] <- NA

            #set a threshold pi*25^2 km2 of tc.aff.mask for analyising 
       #if ( length(which(tc.aff.mask != "NA")) < 3.1416*25**2.) {
         # set affected area on land  to NA
       #  tc.aff.mask[ which(tc.aff.mask != "NA")] <- NA 
         # use tc background map as refence map 
       #  tc.ref.mask <- tc.bg 
       #}            
      
   # 4. load LAI map data this could be in a loop of time period 
   # summer condition
   # if ( (wrk.mon >= mon.str) & (wrk.mon <= mon.end)) {  
      switch.yr <- FALSE 
      if ( (wrk.yr - wrk.yr.old) > 0 ) switch.yr <- TRUE 
      wrk.yr.old <- wrk.yr 
      print(paste("switch year:",switch.yr,sep="")) 
      
      ieve=ieve+1
      if (switch.yr) {
       #write table.comb for reference 
      # st.date <- substr( lai.fnames[it1], start=11, stop=18)
      # ed.date <- substr( lai.fnames[it], start=11, stop=18)
       #  home.dir <- c("/data1/home/ychen/")
      # save(table.comb, file = paste(target,"_",st.date,"to",ed.date,"_table.comb","w_",Ref.wind,"_p_",Ref.rain ,"_",tc.pot.txt,".rda",sep="") )

       # load the LAI dataset  
       # only update the newest year and shift the old one 
       # lai.n2 <- lai.n1
       # lai.n1 <- lai.p0
       # lai.p0 <- lai.p1
       # lai.p1 <- lai.p2
       # lai.p2 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it+2], sep="") )	      
       print(paste("update new year!!"))
       # lai.p0 <- lai.p1
       print("previous/current 10 days LAI load!")
       lai.p0 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it], sep=""), var_st=2 )
       print("next 10days  LAI load!")
       lai.p1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it+1], sep=""), var_st=2 )
    #  lai.p2 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it+2], sep="") )
    #  lai.n1 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it-1], sep="") )
    #  lai.n2 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it-2], sep="") )
      
       # update the land cover data 
       load( paste(wrk.dir,"/LANDCOVER_DATA/",wrk.yr,".esa.landcover.rda",sep="") )
       # set water(210) to NA
       esa.lc[ (esa.lc == 210 ) ] <- NA 
       esa.lc[ (esa.lc >= 10 ) & (esa.lc <= 20 ) ] <- 1
       esa.lc[ (esa.lc >= 30)  & (esa.lc <= 100) ] <- 2
      # esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
       esa.lc[ (esa.lc > 2) ]  <- 3
       # set water to NA
       esa.lc[ (esa.lc == 210 ) ] <- NA 
       # create the LC mask
       lc.mask <- esa.lc
       # update agricuture mask
       lc.agr.mask <- lc.mask
       lc.agr.mask[lc.mask!=1] <- NA
       lc.agr.mask[lc.mask==1] <- 1
       # update forest mask
       lc.for.mask <- lc.mask
       lc.for.mask[lc.mask!=2] <- NA
       lc.for.mask[lc.mask==2] <- 1
       # other mask
       lc.oth.mask <- lc.mask
       lc.oth.mask[lc.mask!=3] <- NA
       lc.oth.mask[lc.mask==3] <- 1

      } else{
       print("previous/current 10 days LAI load!")
       lai.p0 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it], sep=""), var_st=2)
       print("next 10 days LAI load!")  
       lai.p1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it+1], sep=""), var_st=2)
      }	 #end if      

   # find pixels for dek_m1 
   dek_m0 <- array(NA,dim=c(nx,ny)) 
   dek_m1 <- array(NA,dim=c(nx,ny)) 
   #dek_m2 <- array(NA,dim=c(nx,ny)) 
   #dek_m3 <- array(NA,dim=c(nx,ny)) 
   dek_p1 <- array(NA,dim=c(nx,ny)) 
   #dek_p2 <- array(NA,dim=c(nx,ny)) 
   #dek_p3 <- array(NA,dim=c(nx,ny)) 
   lai_be <- array(NA,dim=c(nx,ny))
   lai_af <- array(NA,dim=c(nx,ny))
   
# it.tc = 330
obs_cont = 6
cut.off = 0.1
# LAI Change analysis
# loop for spatil domain
# test for taiwan region 

if (target=="WP") {
nx1=500; nx2=5000 
ny1=1; ny2=4500 
}

if (target=="TW") {
nx1=1000; nx2=2000  
ny1=1800; ny2=3300 
}

if (target=="TEST") {
nx1=1000; nx2=1900  
ny1=2050; ny2=3550 
}


# test for Philippien region
if (target=="PH") {
nx1=  1600; nx2=2000
ny1=  3600; ny2=4100
}

# test for south of Japan region
if (target=="JP") {
nx1=  2700; nx2=3100
ny1=  1700; ny2=2100
}

# test region for Hainan island
if (target=="HN") {
nx1=200;  nx2=1200
ny1=2600; ny2=3600
}

   for ( iy in ny1:ny2) {
#   for ( iy in 2500:2500) {
	
	for ( ix in nx1:nx2) {

            # add if condiction only for potentail TC affect area
            #if ( (!is.na(lc.for.mask[ix,iy])) ) {       	

   #for ( ix in 1:nx) {
   #	for ( iy in 1:ny) {
              # get index for time
	      # test use 200
              # dek_m1
                be_cont = 60
                af_cont = 40
                if ( is.na(lai.p0$LENGTH_BEFORE[ix,iy])==FALSE &  is.na(lai.p0$QFLAG[ix,iy])==FALSE ) {
                if ( (lai.p0$LENGTH_BEFORE[ix,iy] <= be_cont) &
                     (lai.p0$LENGTH_AFTER[ix,iy] <=  af_cont) &
                     (lai.p0$QFLAG[ix,iy] == 0)  )  {
                        dek_m0[ix,iy] <- 1
                }else {
                        dek_m0[ix,iy] <- NA
                }
                }
                # calculate the lai before TC date

                lai_be[ix, iy] <- dek_m0[ix,iy]*lai.p0$LAI[ix,iy]
              # dek_p1
                be_cont = 40
                af_cont = 60
                if ( is.na(lai.p1$LENGTH_BEFORE[ix,iy])==FALSE &  is.na(lai.p1$QFLAG[ix,iy])==FALSE) {
                if ( (lai.p1$LENGTH_BEFORE[ix,iy] <= be_cont) &
                     (lai.p1$LENGTH_AFTER[ix,iy] <=  af_cont) &
                     (lai.p1$QFLAG[ix,iy] == 0) )  {
                        dek_p1[ix,iy] <- 1
                }else {
                        dek_p1[ix,iy] <- NA
                }
                }
                # calculate the lai after TC date  
		lai_af[ix, iy] <-  dek_p1[ix,iy]*lai.p1$LAI[ix,iy] 

          
     #}#end of if for land mask only 
    }#end of ix 

   # add zonal averge table for the analyis  
   # set cut.off 
   #  del.lai.zonal <- (lai_af[,iy] - lai_be[,iy]) * lc.for.mask[,iy] * tc.aff.mask[,iy] 
   #  tot.zonal.for <- length(which( del.lai.zonal != "NA"))     
   #  if (tot.zonal.for >= 20) {  
   #       # print(paste("total zonal pixcels:", tot.zonal.for,sep=""))
   #      gain.zonal.for[ieve,iy] <- length(which(del.lai.zonal >= cut.off)  ) /  tot.zonal.for *100.
   #       # print(paste("total gain pixcels:", gain.zonal.for[ieve,iy],sep=""))
   #      gray.zonal.for[ieve,iy] <- length(which( (del.lai.zonal < cut.off) & ( del.lai.zonal > -1*cut.off  ) )) / tot.zonal.for * 100.
   #       # print(paste("total gray pixcels:", gray.zonal.for[ieve,iy],sep=""))
   #      loss.zonal.for[ieve,iy] <- length(which(del.lai.zonal <= -1*cut.off)) / tot.zonal.for * 100.
   #       # print(paste("total loss pixcels:", loss.zonal.for[ieve,iy],sep=""))
   #  }

   #  del.lai.zonal.pot <- (lai_af[,iy] - lai_be[,iy]) * lc.for.mask[,iy] * (tc.ref.mask[,iy])  
   #  tot.zonal.for.pot <- length(which( del.lai.zonal.pot != "NA"))     
   #  if (tot.zonal.for.pot >= 20) {  
   #       # print(paste("total zonal pixcels:", tot.zonal.for,sep=""))
   #      gain.zonal.for.pot[ieve,iy] <- length(which(del.lai.zonal.pot >= cut.off)  ) /  tot.zonal.for.pot *100.
   #       # print(paste("total gain pixcels:", gain.zonal.for[ieve,iy],sep=""))
   #      gray.zonal.for.pot[ieve,iy] <- length(which( (del.lai.zonal.pot < cut.off) & ( del.lai.zonal.pot > -1*cut.off  ) )) / tot.zonal.for.pot * 100.
   #       # print(paste("total gray pixcels:", gray.zonal.for[ieve,iy],sep=""))
   #      loss.zonal.for.pot[ieve,iy] <- length(which(del.lai.zonal.pot <= -1*cut.off)) / tot.zonal.for.pot * 100.
   #       # print(paste("total loss pixcels:", loss.zonal.for[ieve,iy],sep=""))
   #  }


   # add information along the latitude zonal band 
   # summary the informatio to the table 
#   table.tmp <- data.frame( date.time = wrk.date,
#                            date.yr = wrk.yr,
#                            date.mon = wrk.mon,
#                            nav.lat = round(nav.lat[iy]),
#                            zonal.gain.for = gain.zonal.for[ieve,iy],
#   			    zonal.gray.for = gray.zonal.for[ieve,iy], 
#			    zonal.loss.for = loss.zonal.for[ieve,iy],
#			    zonal.gain.for.pot = gain.zonal.for.pot[ieve,iy],
#			    zonal.gray.for.pot = gray.zonal.for.pot[ieve,iy],
#			    zonal.loss.for.pot = loss.zonal.for.pot[ieve,iy] )
  # combine the table
#  table.comb <- rbind( table.comb, table.tmp ) 

  } #end of of iy (along  latitude) LAI analyis 
  # TC event loop Start 
  for ( ieve in 1: length(tc.eve)) {
 #  for ( ieve in 2:2 ) {
   # initialized the event based masks   
   #eve.tc.aff.mask <- array(NA, dim=c(5040,5040))
   eve.tc.ref.mask <- array(NA, dim=c(5040,5040))
   eve.tc.pot.mask <- array(NA, dim=c(5040,5040))
 
   #assign event based tc masks value 
   #eve.tc.aff.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
   eve.tc.ref.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
   eve.tc.pot.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
   
   #setp 1: select the event within 10 days & VMAX>65 for Typhoon
   event.tmp <-  subset(hourly.track, ((hourly.track$YYYYMMDDHH >= as.numeric(wrk.date.bef)) &
                                        (hourly.track$YYYYMMDDHH <  as.numeric(wrk.date.crt)) & (hourly.track$CY==tc.eve[ieve])) )
   event.tmp <- subset(event.tmp, (as.numeric(as.character(event.tmp$VMAX)) >= 35) )

    #mask out for having a 10_days max_wind  over the ref.wind in reference area
    tmp.arr <- array(NA,dim=c(nx,ny))
    tmp.arr[ ((tc.mws$max_wind <= Ref.wind)&(tc.acf$acc_rainf <= Ref.rain))] <- 1
    eve.tc.ref.mask <- (eve.tc.ref.mask )*tmp.arr

    eve.ref.pix.tot <- length(which(!is.na(eve.tc.ref.mask*lc.for.mask)))
    print(paste("total TC reference pixels:", eve.ref.pix.tot,sep=""))

    #------------------------
    #if( (eve.aff.pix.tot < 1000)  ) {
    # go to the next step for the it loop because of s,all sample size pixels below 20
    # print( "Go to next it/time step!")
    #next
    #}
    #------------------------

    # plot the results
    ld_go <- TRUE
     if(ld_go) {
     #dev.new()
     #get tC track to event base 
       
    # 
#     bitmap(file=paste("./plots_3d/LAI_delta_mask_",target,"_",wrk.date,".png",sep=""),type="png16m",
#              width = 8, height = 8, units = "in", res= 600 )
     layout(matrix(data=seq(1,4,1),nrow=2, ncol=2, byrow=TRUE))
     ylim <- c(10, 42)
     xlim <- c(105,140)
     #ylim <- c(5, 45)
     #xlim <- c(105,150)
     # wind speed & rainfall
     mws.raster <- raster( x=t(tc.mws$max_wind),
                             xmn=lai.p0$lon[1],  xmx=lai.p0$lon[nx],
                             ymn=lai.p0$lat[ny], ymx=lai.p0$lat[1],
                             crs=CRS("+proj=longlat +datum=WGS84"))
     acf.raster <- raster( x=t(tc.acf$acc_rainf),
                             xmn=lai.p0$lon[1],  xmx=lai.p0$lon[nx],
                             ymn=lai.p0$lat[ny], ymx=lai.p0$lat[1],
                             crs=CRS("+proj=longlat +datum=WGS84"))
   # par(xpd=TRUE)
     my.col    <- colorRampPalette(c("blue","cyan","green","yellow","orange","red","brown"))(20)
     plot(mws.raster, xlim=xlim, ylim=ylim, main=paste("10 days maximum wind-speed",sep=""),col=my.col,   
           legend.args=list(text="Ws.(m/s)", side=3, font=2, line=0.2, cex=0.8),legend.width=3,
                xlab= "Longitude", ylab="Latitude", cex.lab=1, cex.main=1 )
     points( y=as.numeric(as.character(event.tmp$LatNS)),
            x=as.numeric(as.character(event.tmp$LonEW)),
            col="black",cex=0.2 )
     plot(coastlines, add=T)
     # rainfall
     plot(acf.raster, xlim=xlim, ylim=ylim,main=paste("10 days accumulative rainfall",sep=""), col=my.col,
          legend.args=list(text="Rf.(mm)", side=3, font=2, line=0.2, cex=0.8),legend.width=3,
               xlab= "", ylab="", cex.lab=1, cex.main=1 )
     points( y=as.numeric(as.character(event.tmp$LatNS)),
            x=as.numeric(as.character(event.tmp$LonEW)),
            col="black",cex=0.2 )
     plot(coastlines, add=T)
     #  current lai plot 
     my.col    <- colorRampPalette(c("gray85","gray45","orange","yellow","green","forestgreen","purple"))(8)
     my.breaks <- seq(0, 7, length.out = 8)
     lai.raster <- raster( x=t(lai.p0$LAI),
                             xmn=lai.p0$lon[1],  xmx=lai.p0$lon[nx],
                             ymn=lai.p0$lat[ny], ymx=lai.p0$lat[1],
                             crs=CRS("+proj=longlat +datum=WGS84"))
     plot(lai.raster, ylim=ylim,xlim=xlim, legend=TRUE, col=my.col,breaks=my.breaks,
          legend.args=list(text="LAI", side=3, font=2, line=0.2, cex=0.8),legend.width=3,
          xlab="Longitude", ylab="Latitude", cex.lab=1.0, cex.main=1.0,
          main=paste("Current LAI: ",wrk.date,sep=""))
     points( y=as.numeric(as.character(event.tmp$LatNS)),
            x=as.numeric(as.character(event.tmp$LonEW)),
            col="black",cex=0.2 )
     plot(coastlines, add=T)
     # del1ta lai
     my.col    <- colorRampPalette(c("darkblue","blue","lightblue","white","pink","red","brown"))(6)
     my.breaks <- c( seq(-1.5, 1.5, length.out = 7)) 
     # lai_del in TC potentail area  mask 
     del_lai <- (lai_af - lai_be ) * lc.for.mask * eve.tc.pot.mask
     lai.raster.del <- raster( x=t(del_lai),
                             xmn=lai.p0$lon[1],  xmx=lai.p0$lon[nx],
                             ymn=lai.p0$lat[ny], ymx=lai.p0$lat[1],
                             crs=CRS("+proj=longlat +datum=WGS84"))
     plot(lai.raster.del, ylim=ylim,xlim=xlim, main="Potential affected forest area",zlim=c(-1.5,1.5), 
         legend.args=list(text="Del LAI", side=3, font=2, line=0.2, cex=0.8),legend.width=3,
         col=my.col,breaks=my.breaks,
              xlab= "", ylab="",cex.lab=1, cex.main=1 )
    #
     points( y=as.numeric(as.character(event.tmp$LatNS)),
            x=as.numeric(as.character(event.tmp$LonEW)),
            col="black",cex=0.2 )
     plot(coastlines, add=T)
  #   par(xpd=NA)
    # text(x=xlim[1]-5, y=lim[2]+5, paste("Obs. date:", wrk.date,sep=""),cex=1.5)
  
    # close png plot 
 #    dev.off()
   } #if ld_go 




 } #TC event loop  

  #release the memory 
 gc()
} #end of time loop 

# save the array 
#save( file=paste(target,"_gain.zonal.for_",radius,".rda",sep=""), gain.zonal.for )
#save( file=paste(target,"_loss.zonal.for_",radius,".rda",sep=""), loss.zonal.for )

#save( file=paste(target,"_gain.zonal.for_",radius,"_pot.rda",sep=""), gain.zonal.for.pot )
#save( file=paste(target,"_loss.zonal.for_",radius,"_pot.rda",sep=""), loss.zonal.for.pot )

#go.acc <- FALSE
#if  (go.acc) {
#
#   bitmap(file=paste("LAI_ACC_PLOT_",target,".png",sep=""),type="png16m",
#          width = 8, height = 8, units = "in", res= 600 )
#   layout(matrix(data=seq(1,2,1),nrow=1, ncol=2, byrow=TRUE))
	
# plot the time seriese data of acc.lai.drop
#plot( cumsum(table.comb$pix.loss.for), col="red",
#      yaxt="s", xaxt="n", type="l",
#      xlab="Observation date",
#      ylab="Acc. LAI drop/increase (km2, 10-days)", main = "Forest Area (affected)" )
#lines( cumsum(table.comb$pix.gain.for), col="blue")

# plot the time seriese data of acc.lai.drop
#plot( cumsum(table.comb$pix.gain.for.pot), col="blue",
#      yaxt="s", xaxt="n", type="l",
#      xlab="Observation date",
#      ylab="Acc. LAI drop/increase (km2, 10-days)", main= "Forest Area (unaffected)" )
#lines( cumsum(table.comb$pix.loss.for.pot), col="red")
#dev.off()
#}
#

#} # end of fun_dyn_lai
