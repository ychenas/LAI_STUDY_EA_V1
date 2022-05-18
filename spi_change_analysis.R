# -----------------------------------------------------------------
# This script/function analysis the SPOT-VGT LAI changes
# First Date:  -2018-
# Revised Date: 2019-12-05, 2020-01-21, 2020-04-23, 2020-04-28 
# Revised Date: 2020-02-21 add Effect size of precipitation 
# Revised Date: 2020-03-31 add offset days information
# Revised Date: 2020-05-22 bug fixed observation numbers issue
# Revised Date: 2020-05-26 bug fixed effect size of ET
# Revised Date: 2020-06-29 add quality check flags  
#          qc1: for checking mean difference between affect and reference before the TC event 
#          (remove) qc2: for checking mean difference caused by the TC event  
#          qc2: the error prppogation term as funtion of sample number and mean value
# Revised Date: 2020-07-03 add TC event intensity/maximum wind speed
# Revised Date: 2020-10-20 revised the simulation doamin to East Asia (0 to 60N, 90E,150E) 
# Revised Date: 2021-04-03 refine the quality check by incroplating the error propogation from samples. 
# Revised Date: 2021-07-24 add SPEI information for each event based analysis 
# -----------------------------------------------------------------

fun.dyn.lai <- function ( Ref.wind=3.0, Ref.rain=100, win.size=40, offdays=30, pdays=60, track.dir="TRACK_DATA_TS_125D")
{

# ====================================
# data for analysis from 1999 to 2018 
# ====================================

target=c("EA");
Del.wind=0; Del.rain=0; 

spei.fname = "spei01"

#all runs
it1=5; it2=36*20

#track.dir="/TRACK_DATA_2D/";
wrk.dir="/lfs/home/ychen/LAI_STUDY_EAsia/";
mon.str=1; mon.end=12;


#pdays: post TC-event days  
#offdays: offset days LAI observation  
test_ld <- FALSE
#test_ld <- TRUE


if (test_ld) {
  # set working files directory  
  wrk.dir <- c("/lfs/home/ychen/LAI_STUDY_EAsia/")
  track.dir <- c("/TRACK_DATA_2D/") 
  pdays=60
  win.size=60
  offdays=00
  #it1=90
  #it2=90
  Del.wind = 0.0 
  Del.rain = 0.0
  #set a reference wind speed  for maximum wind speed for a reference area
  #Sulik it=528
  Ref.wind = 0
  Ref.rain = 200
  #target=c("TEST"); it1=22; it2=22; mon.str=1;mon.end=12   #zero wind speed
  
  #target ="TEST"
  #target=c("TW"); it1=528; it2=528; mon.str=1;mon.end=12  #Sulik
  target=c("TW"); it1=382; it2=382; mon.str=1;mon.end=12   #Morakot
  # target=c("TW"); it1=274; it2=274; mon.str=1;mon.end=12
  #target=c("TW"); it1=275; it2=275; mon.str=1;mon.end=12
  #2018-06-20
  #target=c("EA"); it1=701; it2=701; mon.str=1;mon.end=12
}

Ref.wind <- as.integer(Ref.wind)
Ref.rain <- as.integer(Ref.rain)
win.size <- as.integer(win.size)
offdays  <- as.integer(offdays)
pdays    <- as.integer(pdays)

# pdays: days after tc events 
# pit: post itstep for 10days
#pdays = 60 
pit = as.integer(pdays/10)
  
# adjust for the final year 
#it2 = it2-pit 
#get output tc.track information and off-day
pos.txt <-  paste(substr(track.dir, start=12, stop=14),"_pos_",formatC(pdays, width=3, flag=0),sep="")
off.txt <-  paste("off_", formatC(offdays, width = 3,flag = 0),sep="")
#
#2001 to 2018
#target=c("WP"); it1=36*2+1;it2=36*20; mon.str=3;mon.end=11
#target=c("WP"); it1=36*19+1;it2=36*20; mon.str=3;mon.end=11
#target=c("TEST")

# get track data 
track.data <- data.frame()
for (iyr in 1999:2018) {
   tmp.data <- read.csv(file=paste(wrk.dir,"/TRACK_DATA/hourly_data/wpb_",iyr,"_hourly.txt",sep=""))
   track.data <- rbind(tmp.data,track.data)
}

#test run
#target=c("WP"); it1=36*10+1;it2=36*11; mon.str=3;mon.end=11
#it1=36*12+1;it2=36*14; mon.str=3;mon.end=11
#test cases
#it1=275 #20060820
#it2=275
#it1=528 #20130831 Sulik
#it2=528
#it1=205 #20040910
#it2=205
#it1=382 #20090810 Morakot
#it2=382
#neve=1
neve=floor((it2-it1+1)/36)*(mon.end-mon.str+1)*3+1

#get file list from the folder 
lai.fnames <- list.files(path=paste(wrk.dir,"/LAI_DATA/",sep=""), pattern="*.nc") 
#logic varibale ld.tc
ld.tc <- array(FALSE, dim=c(length(lai.fnames))) 
# assign the logic variable for with or without TC events every 10 days!
for (i in 1:length(lai.fnames) ) {
 # get lai date yyyydd and find the TC event   
 lai.stp <- substr(lai.fnames[i],start=11,stop=16)
 tmp.stp <- track.data$YYYYMMDDHH[which(substr(track.data$YYYYMMDDHH,start=1,stop=6) == lai.stp )][1]
 #assign logic to the vara
 if ( is.na(tmp.stp) == TRUE ) {
   ld.tc[i] <- FALSE
   print(paste("Without TC event:", substr(lai.fnames[i],start=11,stop=18),sep="") )
 }else{
   ld.tc[i] <- TRUE 
   #print(paste("Find TC event:", substr(lai.fnames[i],start=11,stop=18),sep="") )
 }
}

#plot accumulative map 
#pdf(file="acc.plot.pdf",width=8, height=8)
#par(lwd=0.5, mar=c(5,5,1.5,1) ) #,mpg=c(2,1,0)
#layout(matrix(data=c(1,2), nrow=1, ncol=2, byrow=TRUE))
#library("fields")
# The land cover maps were aggregated to 1km from ESA 300m land cover product from 1992 to 2015
# The leaf area index maps were the popernious 1km land product from 1999 to 2018 
# The typhoon track data was compiled from Join Typhoon Warnning Center (JTWC) from 1945 to 2017 
# load the R data after several pre-pocessing


# 1. TC track occurence from 1948 to 2017 (70 years) the domain should match the LAI map over the West Pacific region (5040 * 5040)
#    replaced by new definition     
#    load("TC_return_frq_1948to2017.rda")
#    variable name: land.frq"
#    use a 2 year return frequency to create the TC mask over WP Basin 
#    set.frq=0.5
#    land.frq[ land.frq/70 <  set.frq]  <- NA
#    land.frq[ land.frq/70 >= set.frq]  <- 1


# 2. ESA land cover map at 1km, the domain should be matched the LAI map over study area

load(paste(wrk.dir,"/LANDCOVER_DATA/1999.esa.landcover.east.asia.rda",sep=""))
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

go.raster <- FALSE
if (go.raster) {
# load coastal boundary file
library("plyr")
library("fields")
library("rgdal")
library("maptools")
library("raster")
#load the coastlines data by readOGR function from sp package
coastlines <- readOGR("/lfs/home/ychen/GIS//lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")
#convert arrary to raster obj for ploting
#xymn <- c(105.0089, 150, 5, 49.99107)
xymn <- c(90.00, 150, 0, 60)
tc.pot.raster <- raster( x=t(tc.pot.mask), 
			  xmn=xymn[1], xmx=xymn[2], ymn=xymn[3], ymx=xymn[4],
			  crs=CRS("+proj=longlat +datum=WGS84"))
}
# plot the mask 
#library(fields)
#image.plot( (tc.mask*lc.mask)[,5040:1], main="The Study Area" ) 
# create a table summaries the result of condictional LAI change based on wind and rainfall scales
  table.comb <- data.frame()
# get LAI file list for loading  LAI map data this could be a loop of time period 
  source("read_ncdf4.R")
  wrk.yr.old = 0
  ieve = 0
  # 
   nx=6722
   ny=6722 
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
   #   nav.lat  <- array(NA, dim=c(ny)) 
   
   lai.1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[1], sep=""), var_st=2 )
   for (ix in 1:nx) {
      for (iy in 1:ny) {
         nav.lat[ix,iy] <- lai.1$lat[iy]
         nav.lon[ix,iy] <- lai.1$lon[ix] 
      }
   }

  for ( it in it1:it2) { 
 
  # 
   # get obstime information 
   # lai.p0$time <- format(as.Date(lai.p0$time,origin="1970-01-01"),format="%Y%m%d")
   # print( paste("working date:", lai.p0$time, sep="") )
   wrk.date <- substr( lai.fnames[it], start=11, stop=18)
   wrk.yr  <- as.integer(substr(wrk.date, start=1, stop=4) )
   wrk.mon <- as.integer(substr(wrk.date, start=5, stop=6) ) 
   print( paste("working date:", wrk.date, sep="") )
 
   # skip the time step without TC activities 
   if (ld.tc[it] == FALSE){
      print("No TC events, so move to next step")  
      next
   } 
  
 # 3. JTWC TC track mask with radius mapping 
       # 
       # Define potential TC area as affected area 

       #
       # update the tc track data 
       # load( paste(wrk.yr,"_TC_mask_",radius,".rda",sep="") )
       # tc.mask <- track.mask # loading from loading TC_mask.rda
       # load tc occurance map from 1999-2018 tc track data analysis as background reference if there was no tc activities   
       load( "tc.occ.1999to2018.rda" )
       tc.occ <- arr.for.plot
       
       # load mask files for TC affect and reference area
       #tc.af1 <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_af1_",wrk.date,"00.nc",sep=""), var_st=1)
       #tc.af2 <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_af2_",wrk.date,"00.nc",sep=""), var_st=1)
#       tc.ref <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"tc_mask_ref_",wrk.date,"00.nc",sep=""), var_st=1)
       tc.pot <- fun_read_nc(arg1=paste(wrk.dir,track.dir,"/","tc_mask_pot_",wrk.date,"00.nc",sep=""), var_st=1)
#
       tc.aff.mask <- tc.pot$tc_pot
       tc.aff.mask[tc.aff.mask==0] <- NA 
       
       if(  (length(which(!is.na(tc.aff.mask*lc.for.mask)))  < 1000)  ) {
       # go to the next step for the it loop because of s,all sample size pixels below 20   
         print( "Go to next it/time step!") 
         next
       # exit it loop
       } 
       # find tc.event id within a 10 days mask
       tc.eve <- levels(factor(tc.aff.mask))  
       print(paste("There is(are) ", length(tc.eve), " TC event(s) within 10 days!",sep=""))
       
       # allocate the temporary array 
      # tmp.ace <- array(0,dim=c(nx,ny))
       tmp.acf <- array(0,dim=c(nx,ny))
       tmp.asw <- array(0,dim=c(nx,ny)) 
       # load maximum windspeed & accumulative rainfall & accumulative evapotranspiration from ERA5 Land hourly reanalysis  
       tc.mws <- fun_read_nc(arg1=paste(wrk.dir,"/WIND_DATA/","max_wind_", wrk.date,"00.nc",sep=""), var_st=1)   
       tc.acf <- fun_read_nc(arg1=paste(wrk.dir,"/PRECIP_DATA/","acc_rainf_", wrk.date,"00.nc",sep=""), var_st=1)   
       
       #before the TC event  set to 10-days      
       for (i in 1:1) {
       wrk.date.b1 <- substr(lai.fnames[it-i],   start=11, stop=18)
       #ace.b1 <- fun_read_nc(arg1=paste(wrk.dir,"/ET_DATA/","acc_et_", wrk.date.b1,"00.nc",sep=""), var_st=1)
       acf.b1 <- fun_read_nc(arg1=paste(wrk.dir,"/PRECIP_DATA/","acc_rainf_", wrk.date.b1,"00.nc",sep=""), var_st=1)
       asw.b1 <- fun_read_nc(arg1=paste(wrk.dir,"/SW_DATA/","soil_water_", wrk.date.b1,"00.nc",sep=""), var_st=1)   
       #asw.b1 <- fun_read_nc(arg1=paste(wrk.dir,"/SW_OBS_DATA/SUBSET/","soil_water_", wrk.date,"00.nc",sep=""), var_st=1)       

       #tmp.ace <- ace.b1$acc_et + tmp.ace
       tmp.acf <- acf.b1$acc_rainf + tmp.acf
       tmp.asw <- asw.b1$soil_water + tmp.asw
       }

       #calculate the average value over the pit events        
       #ace.b1$acc_et <- tmp.ace / as.numeric(1)   
       acf.b1$acc_rainf <- tmp.acf / as.numeric(1)
       asw.b1$soil_water <- tmp.asw / as.numeric(1)
 
       # reallocate the temporary array 
       #tmp.ace <- array(0,dim=c(nx,ny))
       tmp.acf <- array(0,dim=c(nx,ny))

       #after the TC event
       ld_go <- T
       if (ld_go) {
       for (i in 1:1) {
         wrk.date.p1 <- substr(lai.fnames[it+i],   start=11, stop=18)
         #ace.p1 <- fun_read_nc(arg1=paste(wrk.dir,"/ET_DATA/","acc_et_", wrk.date.p1,"00.nc",sep=""), var_st=1)
         acf.p1 <- fun_read_nc(arg1=paste(wrk.dir,"/PRECIP_DATA/","acc_rainf_", wrk.date.p1,"00.nc",sep=""), var_st=1)
         #tmp.ace <- ace.p1$acc_et + tmp.ace  
         tmp.acf <- acf.p1$acc_rainf + tmp.acf 
       }
       #calculate the average value over the pit events        
       #ace.p1$acc_et <- tmp.ace / as.numeric(pit)   
       acf.p1$acc_rainf <- tmp.acf / as.numeric(1)
       
  
       #clean the tmp array 
       rm(list= c("tmp.acf"))
       } #end ld_go if

# 4. load Monthly SPEI maps over the East Asia
       #wrk.date.spei <- paste(wrk.yr,"_",formatC(wrk.mon,width=2,flag=0),sep="") 
       wrk.date.spei <- paste(substr(lai.fnames[it + pit], start=11, stop=14),"_",substr(lai.fnames[it + pit], start=15, stop=16), sep="") 
       spei.data <- fun_read_nc(arg1=paste(wrk.dir,"/SPEI_DATA/",spei.fname,"_", wrk.date.spei,"_EA.nc",sep=""), var_st=1)
       print(paste("finish reading monthly SPEI.", wrk.date.spei,sep=""))

# 5. load LAI map data this could be in a loop of time period 
   
       # summer condition
       # if ( (wrk.mon >= mon.str) & (wrk.mon <= mon.end)) {  
       switch.yr <- FALSE 
       if ( (wrk.yr - wrk.yr.old) > 0 ) switch.yr <- TRUE 
       wrk.yr.old <- wrk.yr 
       print(paste("switch year:",switch.yr,sep="")) 
      
       ieve=ieve+1
       if (switch.yr) {
       #write table.comb for reference 
       st.date <- substr( lai.fnames[it1], start=11, stop=18)
       ed.date <- substr( lai.fnames[it], start=11, stop=18)
       # load the LAI dataset  
       # only update the newest year and shift the old one 
       print(paste("update new year!!"))
       # lai.p0 <- lai.p1
       print("current LAI load!")
       lai.b1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it-1], sep=""), var_st=2 )
       print("next n* 10days  LAI load!")
       lai.p1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it + pit], sep=""), var_st=2 )
    #  lai.p2 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it+2], sep="") )
    #  lai.n1 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it-1], sep="") )
    #  lai.n2 <- fun_read_nc( arg1=paste("LAI_SAMPLE/", lai.fnames[it-2], sep="") )
      
       # update the land cover data 
       load( paste(wrk.dir,"/LANDCOVER_DATA/",wrk.yr,".esa.landcover.east.asia.rda",sep="") )
       # set water(210) to NA
       esa.lc[ (esa.lc == 210 ) ] <- NA 
       esa.lc[ (esa.lc >= 10 ) & (esa.lc <= 40 ) ] <- 1
       esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
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
       print("current LAI load!")
       lai.b1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it-1], sep=""), var_st=2)
       print("next n*10 days LAI load!")  
       lai.p1 <- fun_read_nc( arg1=paste(wrk.dir,"/LAI_DATA/", lai.fnames[it + pit], sep=""), var_st=2)
      }	 #end if switch.yr      

   #  allocate the arries 
   dek_b1 <- array(NA,dim=c(nx,ny)) 
   dek_p1 <- array(NA,dim=c(nx,ny)) 
   #
   lai_be <- array(NA,dim=c(nx,ny))
   lai_af <- array(NA,dim=c(nx,ny))
   #
   #ace_be <- array(NA,dim=c(nx,ny))
   #ace_af <- array(NA,dim=c(nx,ny)) 
   #
   acf_be <- array(NA,dim=c(nx,ny))
   #acf_af <- array(NA,dim=c(nx,ny)) 
   #
   asw_be <- array(NA,dim=c(nx,ny))
 
   #Observation offset between two observations
   lai_off <- array(NA, dim=c(nx,ny)) 
   
   # SPEI data
   spei_aff <- array(NA, dim=c(nx,ny)) 
   spei_ref <- array(NA, dim=c(nx,ny))

# it.tc = 330
obs_cont = 6
cut.off = 0.0
# LAI Change analysis
# loop for spatil domain
# test for taiwan region 

if (target=="EA") {
nx1=1000; nx2=6500 
ny1=1500; ny2=6000 
}

if (target=="TW") {
nx1=3000; nx2=4000  
ny1=3500; ny2=4500 
}

if (target=="TEST") {
nx1=1000; nx2=4000  
ny1=2000; ny2=5000 
}


# test for Philippien region
#if (target=="PH") {
#nx1=  1600; nx2=2000
#ny1=  3600; ny2=4100
#}

# test for south of Japan region
#if (target=="JP") {
#nx1=  2700; nx2=3100
#ny1=  1700; ny2=2100
#}

# test region for Hainan island
#if (target=="HN") {
#nx1=200;  nx2=1200
#ny1=2600; ny2=3600
#}

   for ( iy in ny1:ny2) {
   #   for ( iy in 2500:2500) {
	
	for ( ix in nx1:nx2) {

            # add if condiction only for potentail TC affect area 
            # (this is important to save computational time for the analyis 
            if ( (!is.na(tc.aff.mask[ix,iy])) ) {       	

              # calculate the offset between two observations 
                if  (is.na(lai.b1$LENGTH_BEFORE[ix,iy])==FALSE & is.na(lai.b1$LENGTH_BEFORE[ix,iy])==FALSE ) {
                      b1.day <- seq( (0. - lai.b1$LENGTH_BEFORE[ix,iy]), lai.b1$LENGTH_AFTER[ix,iy], 1)      }
                else{ b1.day <- NA }
                # 
                if  (is.na(lai.p1$LENGTH_BEFORE[ix,iy])==FALSE & is.na(lai.p1$LENGTH_BEFORE[ix,iy])==FALSE ) {
                      p1.day <- seq( (0. - lai.p1$LENGTH_BEFORE[ix,iy]), lai.p1$LENGTH_AFTER[ix,iy], 1) + pdays}
                else{ p1.day <- NA }
                # case for overlap 
                     off.day <- -1* length (intersect ( b1.day, p1.day) ) 
                # case for offset 
                #if ( off.day == 0) { off.day <- ( ((min(p1.day)+max(p1.day))/2)   - ((min(p0.day)+max(p0.day))/2) )  }   
                if ( off.day == 0) { off.day <- (min(p1.day) - max(b1.day)) }              
 
                lai_off[ix,iy] <- off.day   

                # set window  
                #be_cont = 50
                #af_cont = 40
                be_cont = win.size+10
                af_cont = win.size
                
                if ( is.na(lai.b1$LENGTH_BEFORE[ix,iy])==FALSE &  is.na(lai.b1$QFLAG[ix,iy])==FALSE  & is.na(lai_off[ix,iy])==FALSE ) {
                if ( (lai.b1$LENGTH_BEFORE[ix,iy] <= be_cont) & (lai.b1$LENGTH_AFTER[ix,iy] <=  af_cont) & 
                     (lai_off[ix,iy] >= offdays  & lai_off[ix,iy] <= pdays+offdays) )   
                      {
                        dek_b1[ix,iy] <- 1
                }else {
                        dek_b1[ix,iy] <- NA
                }
                }
                # calculate the lai before TC date
                lai_be[ix,iy] <- dek_b1[ix,iy]*lai.b1$LAI[ix,iy]
                # calculate the accumulative ET over previous 10days           
                #ace_be[ix,iy] <- dek_b1[ix,iy]*ace.b1$acc_et[ix,iy]
                # calculate the accumulative Precipitation over previous 10days           
                acf_be[ix,iy] <- dek_b1[ix,iy]*acf.b1$acc_rainf[ix,iy]
                # calculate the average soil moisture            
                asw_be[ix,iy] <- dek_b1[ix,iy]*asw.b1$soil_water[ix,iy]
                
                # dek_p1
                #be_cont = 40
                #af_cont = 50
                be_cont = win.size
                af_cont = win.size+10
               
                if ( is.na(lai.p1$LENGTH_BEFORE[ix,iy])==FALSE &  is.na(lai.p1$QFLAG[ix,iy])==FALSE & is.na(lai_off[ix,iy])==FALSE ) {
                if ( (lai.p1$LENGTH_BEFORE[ix,iy] <= be_cont)  &  (lai.p1$LENGTH_AFTER[ix,iy]  <=  af_cont) & 
                     (lai_off[ix,iy] >= offdays  & lai_off[ix,iy] <= pdays+offdays) ) 
                      {
                        dek_p1[ix,iy] <- 1
                }else {
                        dek_p1[ix,iy] <- NA
                }
                }
                # calculate the lai after TC date  
		lai_af[ix,iy] <- dek_p1[ix,iy]*lai.p1$LAI[ix,iy] 
                # calculate the accumulative rainfall for current step 
                #ace_af[ix,iy] <- dek_p1[ix,iy]*ace.p1$acc_et[ix,iy]
                # calculate the accumulative Precipitation over previous 10days           
                #acf_af[ix,iy] <- dek_p1[ix,iy]*acf.p1$acc_rainf[ix,iy]
            
             }#end of if for land mask only 
    
        }#end of ix 
  } #end of of iy (along  latitude) LAI analyis 



  # TC event loop Start 
  for ( ieve in 1: length(tc.eve)) {
   # for ( ieve in 1:1 ) {
   # flags for the quality check ...
   qc1.score <- NA
   qc1.flag <- TRUE
   qc2.score <- NA 
   qc2.flag <- TRUE

   # initialized the event based masks   
   eve.tc.aff.mask <- array(0, dim=c(6722,6722))
   eve.tc.ref.mask <- array(0, dim=c(6722,6722))
   #assign event based tc masks value 
   eve.tc.aff.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
   eve.tc.ref.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
   #apply the mask from offset condiiction
   eve.tc.aff.mask <- eve.tc.aff.mask*dek_p1*dek_b1
   eve.tc.ref.mask <- eve.tc.ref.mask*dek_p1*dek_b1
    #mask out for having a 10_days max_wind  over the ref.wind in reference area
   # define the median wind and rainfall as reference
    # Ref.wind.1 <- mean(tc.mws$max_wind*tc.pot$tc_pot, na.rm=T)
    Ref.rain.1 <- as.numeric(Ref.rain) 
    Ref.wind.1 <- as.numeric(Ref.wind)  
    
    #assign event based tc masks value 
    eve.tc.aff.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0
    eve.tc.ref.mask[ tc.pot$tc_pot == as.integer(tc.eve[ieve]) ] <- 1.0

    #apply the mask from offset condition
    eve.tc.aff.mask <- eve.tc.aff.mask*dek_p1*dek_b1
    eve.tc.ref.mask <- eve.tc.ref.mask*dek_p1*dek_b1

    #mask out for having a 10_days max_wind  over the ref.wind in reference area
    #only rainfall condition
    if(Ref.wind.1 == 0) {
      print("Search affected/reference area, rainfall only mode") 
      #search/define reference area
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ (tc.acf$acc_rainf <=  Ref.rain.1)] <- 1
      eve.tc.ref.mask <- (eve.tc.ref.mask )*tmp.arr
      #search/define affected  area
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ (tc.acf$acc_rainf >  Ref.rain.1)] <- 1
      eve.tc.aff.mask <- (eve.tc.aff.mask )*tmp.arr
    }
    # only wind-speed condition
    if(Ref.rain.1 == 0) {
      print("Search affected/reference area, wind-speed only mode") 
      #search/define reference area
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ (tc.mws$max_wind <= Ref.wind.1)] <- 1
      eve.tc.ref.mask <- (eve.tc.ref.mask )*tmp.arr
      #search/define affected  area
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ (tc.mws$max_wind > Ref.wind.1)] <- 1
      eve.tc.aff.mask <- (eve.tc.aff.mask )*tmp.arr
    }
    # Combinded condition
    if ( (Ref.rain.1 > 0) & (Ref.wind.1 > 0) ) {  
      #search/define reference area 
      print("Search affected/reference area, combinded wind-speed and rainfall mode") 
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ ((tc.mws$max_wind <= Ref.wind.1)&(tc.acf$acc_rainf <= Ref.rain.1))] <- 1
      eve.tc.ref.mask <- (eve.tc.ref.mask )*tmp.arr
      #search for affected area
      tmp.arr <- array(NA,dim=c(nx,ny)) 
      tmp.arr[ ((tc.mws$max_wind > Ref.wind.1)|(tc.acf$acc_rainf > Ref.rain.1))] <- 1
      eve.tc.aff.mask <- ( eve.tc.aff.mask )*tmp.arr
    
    }

    print(paste("define ref.wind.eve:",Ref.wind.1,", define ref.rainf.eve:",Ref.rain.1,sep="")) 
 
    eve.aff.pix.tot <- length(which(!is.na(eve.tc.aff.mask*lc.for.mask)))
    print(paste("total TC affected pixels:", eve.aff.pix.tot,sep=""))
    
    eve.ref.pix.tot <- length(which(!is.na(eve.tc.ref.mask*lc.for.mask)))
    print(paste("total TC reference pixels:", eve.ref.pix.tot,sep=""))
       
    #------------------------
    if( (eve.aff.pix.tot < 20) | (eve.ref.pix.tot <= 20) ) {
    # go to the next step for the it loop because of s,all sample size pixels below 20   
     
     print( "sample < 100(aff) / 100(ref),  Go to next it/time step!") 
    next
    }
    # exit it loop
    
    #------------------------
 
   # get mean/max/min latitude and latitude for affected area within 10 days 
   eve.lat.med.aff <- median(eve.tc.aff.mask*lc.for.mask*nav.lat,na.rm=T) 
   eve.lat.max.aff <- max(eve.tc.aff.mask*lc.for.mask*nav.lat,na.rm=T)
   eve.lat.min.aff <- min(eve.tc.aff.mask*lc.for.mask*nav.lat,na.rm=T)
   #
   eve.lon.med.aff <- median(eve.tc.aff.mask*lc.for.mask*nav.lon,na.rm=T) 
   eve.lon.max.aff <- max(eve.tc.aff.mask*lc.for.mask*nav.lon,na.rm=T) 
   eve.lon.min.aff <- min(eve.tc.aff.mask*lc.for.mask*nav.lon,na.rm=T) 
   
   # get tc occurence information for TC affected area 
   eve.tc.occ.aff  <- median(eve.tc.aff.mask*lc.for.mask*tc.occ, na.rm=T) 

   # calculate pixels for each event
   eve.for.pix.aff <- length( which(!is.na(eve.tc.aff.mask*lc.for.mask)) ) 
   eve.agr.pix.aff <- length( which(!is.na(eve.tc.aff.mask*lc.agr.mask)) )
   eve.oth.pix.aff <- length( which(!is.na(eve.tc.aff.mask*lc.oth.mask)) )

   # calculate median of maxium wind speed over various land type in the affected area 
   eve.for.mws.aff <- median( eve.tc.aff.mask*lc.for.mask*tc.mws$max_wind, na.rm=T ) 
   eve.agr.mws.aff <- median( eve.tc.aff.mask*lc.agr.mask*tc.mws$max_wind, na.rm=T )
   eve.oth.mws.aff <- median( eve.tc.aff.mask*lc.oth.mask*tc.mws$max_wind, na.rm=T )

   # calculate median of accumulative rainfall over various land type in the affected area 
   eve.for.acf.aff <- median( eve.tc.aff.mask*lc.for.mask*tc.acf$acc_rainf, na.rm=T ) 
   eve.agr.acf.aff <- median( eve.tc.aff.mask*lc.agr.mask*tc.acf$acc_rainf, na.rm=T )
   eve.oth.acf.aff <- median( eve.tc.aff.mask*lc.oth.mask*tc.acf$acc_rainf, na.rm=T )

   # calculate median of accumulative rainfall over various land type in the affected area before the TC event 
   bef.for.acf.aff <- median( eve.tc.aff.mask*lc.for.mask*acf.b1$acc_rainf, na.rm=T ) 
   bef.agr.acf.aff <- median( eve.tc.aff.mask*lc.agr.mask*acf.b1$acc_rainf, na.rm=T )
   bef.oth.acf.aff <- median( eve.tc.aff.mask*lc.oth.mask*acf.b1$acc_rainf, na.rm=T )

   # calculate median of offset over various land type in the affected area 
   eve.for.off.aff <- median( eve.tc.aff.mask*lc.for.mask*lai_off, na.rm=T ) 
   eve.agr.off.aff <- median( eve.tc.aff.mask*lc.agr.mask*lai_off, na.rm=T )
   eve.oth.off.aff <- median( eve.tc.aff.mask*lc.oth.mask*lai_off, na.rm=T )

   # get mean/max/min latitude and latitude for refected area within 10 days 
   eve.lat.med.ref <- median( eve.tc.ref.mask*lc.for.mask*nav.lat,na.rm=T) 
   eve.lat.max.ref <- max( eve.tc.ref.mask*lc.for.mask*nav.lat,na.rm=T)
   eve.lat.min.ref <- min( eve.tc.ref.mask*lc.for.mask*nav.lat,na.rm=T)
   #
   eve.lon.med.ref <- median( eve.tc.ref.mask*lc.for.mask*nav.lon,na.rm=T) 
   eve.lon.max.ref <- max( eve.tc.ref.mask*lc.for.mask*nav.lon,na.rm=T) 
   eve.lon.min.ref <- min( eve.tc.ref.mask*lc.for.mask*nav.lon,na.rm=T) 

   # calculate pixels for each event
   eve.for.pix.ref <- length( which(!is.na( eve.tc.ref.mask*lc.for.mask)) ) 
   eve.agr.pix.ref <- length( which(!is.na( eve.tc.ref.mask*lc.agr.mask)) )
   eve.oth.pix.ref <- length( which(!is.na( eve.tc.ref.mask*lc.oth.mask)) )
 
   # calculate median of maxium wind speed over various land type in the affected area 
   eve.for.mws.ref <- median( eve.tc.ref.mask*lc.for.mask*tc.mws$max_wind, na.rm=T ) 
   eve.agr.mws.ref <- median( eve.tc.ref.mask*lc.agr.mask*tc.mws$max_wind, na.rm=T )
   eve.oth.mws.ref <- median( eve.tc.ref.mask*lc.oth.mask*tc.mws$max_wind, na.rm=T )

   # calculate median of accumulative rainfall over various land type in the affected area 
   eve.for.acf.ref <- median( eve.tc.ref.mask*lc.for.mask*tc.acf$acc_rainf, na.rm=T ) 
   eve.agr.acf.ref <- median( eve.tc.ref.mask*lc.agr.mask*tc.acf$acc_rainf, na.rm=T )
   eve.oth.acf.ref <- median( eve.tc.ref.mask*lc.oth.mask*tc.acf$acc_rainf, na.rm=T )

   # calculate median of offset days over various land type in the affected area 
   eve.for.off.ref <- median( eve.tc.ref.mask*lc.for.mask*lai_off, na.rm=T ) 
   eve.agr.off.ref <- median( eve.tc.ref.mask*lc.agr.mask*lai_off, na.rm=T )
   eve.oth.off.ref <- median( eve.tc.ref.mask*lc.oth.mask*lai_off, na.rm=T )
  
   # set a threshold to distinquish the difference retrived from satellite 
   lai_del_abs <- abs(lai_af - lai_be) 

   # remove the noise below the cut.off
   lai_af[lai_del_abs <= cut.off ] <- NA
   lai_be[lai_del_abs <= cut.off ] <- NA 

   #ace_af[lai_del_abs <= cut.off ] <- NA
   #ace_be[lai_del_abs <= cut.off ] <- NA 

  # acf_af[lai_del_abs <= cut.off ] <- NA
   acf_be[lai_del_abs <= cut.off ] <- NA 

   asw_be[lai_del_abs <= cut.off ] <- NA

  # add the check of the mean difference of LAI before the TC event in the reference and the affect area 
   lai.dif.bef <- abs(mean(lai_be*eve.tc.aff.mask*lc.for.mask,na.rm=T) - mean(lai_be*eve.tc.ref.mask*lc.for.mask,na.rm=T))
   print( paste("Absolute mean difference in LAI before the TC event between the reference and the affected area:", lai.dif.bef,sep="") )  
   if ( lai.dif.bef >= 0.25 ) {
       print( paste("Absolute mean difference over than 0.2:", lai.dif.bef," Skip the event !", sep="") )  
       qc1.flag <- FALSE
   }
   qc1.score <- lai.dif.bef
   
   # calculate effect size  
   #for forest area
   aa1 <- mean(lai_af*eve.tc.aff.mask*lc.for.mask,na.rm=T) - mean(lai_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)
     bef.for.lai.aff <- mean(lai_be  * lc.for.mask * eve.tc.aff.mask, na.rm=T) 
     aft.for.lai.aff <- mean(lai_af  * lc.for.mask * eve.tc.aff.mask, na.rm=T)
     print( paste("Forests LAI before in affected area:",bef.for.lai.aff) ) 
     print( paste("Forests LAI after in affected area:", aft.for.lai.aff) ) 
   bb1 <- sqrt( (sd(lai_af*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2. + sd(lai_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2.)/2.)
         print(paste("LAI STD pool: ", bb1,sep=""))
   std.for.lai.aff <- bb1 

   #for referance area    
   aa2 <- mean(lai_af*eve.tc.ref.mask*lc.for.mask,na.rm=T) - mean(lai_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)
     bef.for.lai.ref <- mean(lai_be  * lc.for.mask * eve.tc.ref.mask, na.rm=T) 
     aft.for.lai.ref <- mean(lai_af  * lc.for.mask * eve.tc.ref.mask, na.rm=T)
     print(paste("Forest LAI before in reference area:",bef.for.lai.ref) )
     print(paste("Forest LAI after in reference area:",  aft.for.lai.ref) )
   bb2 <- sqrt( (sd(lai_af*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2. + sd(lai_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2.)/2.)
         print(paste("LAI STD pool: ", bb2,sep=""))
   std.for.lai.ref <- bb2

   bb <- (bb1+bb2)/2.
   if (is.na(bb) != TRUE) {
      eff.size.for <- (aa1-aa2)/bb
  
      #calculation of the error propogation term (aa) 
      tot.n <- eve.aff.pix.tot + eve.ref.pix.tot
      delta.x <- 0.25   
      delta.y <-  sqrt(tot.n*(0.18)**2.)/tot.n
      big.x <- (aa1-aa2)    # LAI difference 
      big.y <- (bb1+bb2)/2. # LAI std
      aa <- abs(eff.size.for)  * sqrt( (delta.x/big.x)**2. + (delta.y/big.y)**2.)  

    }else{
      eff.size.for <- NA
      aa <- NA 
  } 
  
 
   print(paste("Effect size of LAI in forest area:", eff.size.for,sep=""))

  

  #==================================
  # Considering the error propogation,
  # we will only analysis the event with mean difference over than 0.25 = sqrt(0.185^2 + 0.185^2)  
   if ( is.na(aa) != TRUE ) {
        if(  (abs(aa) >= abs(eff.size.for)) | (is.na(eff.size.for) == TRUE) ) {
        # exit the TC event loop  
          print( "QC2! didn't pass") 
          qc2.flag <- FALSE 
       }
          print( paste("TC event, the normalized difference:", aa, sep="") )
         qc2.score <- aa 
   } 
  #=================================


   # mark the mask to spei_2m array 
   # calculate SPEI for the affect and reference area  

   #for forest area
   eve.for.spei.aff  <- mean( spei.data$SPEI*eve.tc.aff.mask*lc.for.mask,na.rm=T) 
     print( paste("Forests Mean SPEI in affected area:", eve.for.spei.aff) ) 
  
   #for referance forest area    
   eve.for.spei.ref  <- mean( spei.data$SPEI*eve.tc.ref.mask*lc.for.mask,na.rm=T) 
     print( paste("Forests Mean SPEI in reference area:", eve.for.spei.ref) )


   # calculate effect size for evapotranspiration  
   #for forest area
  # aa1 <- mean(ace_af*eve.tc.aff.mask*lc.for.mask,na.rm=T) - mean(ace_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)
  #   bef.for.ace.aff <- mean(ace_be  * lc.for.mask * eve.tc.aff.mask, na.rm=T) 
  #   aft.for.ace.aff <- mean(ace_af  * lc.for.mask * eve.tc.aff.mask, na.rm=T)
  #   print( paste("Forests ET before in affected area:",bef.for.ace.aff) ) 
  #   print( paste("Forests ET after in affected area:", aft.for.ace.aff) ) 
  # bb1 <- sqrt( (sd(ace_af*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2. + sd(ace_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2.)/2.)
  #       print(paste("ET STD pool, affected forest: ", bb1,sep=""))
  # std.for.ace.aff <- bb1 
   #for referance area    
  # aa2 <- mean(ace_af*eve.tc.ref.mask*lc.for.mask,na.rm=T) - mean(ace_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)
  #   bef.for.ace.ref <- mean(ace_be  * lc.for.mask * eve.tc.ref.mask, na.rm=T) 
  #   aft.for.ace.ref <- mean(ace_af  * lc.for.mask * eve.tc.ref.mask, na.rm=T)
  #   print(paste("Forests ET before in reference area:",bef.for.ace.ref) )
  #   print(paste("Forests ET after in reference area:",  aft.for.ace.ref) )
  # bb2 <- sqrt( (sd(ace_af*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2. + sd(ace_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2.)/2.)
  #       print(paste("ET STD pool, reference forest: ", bb2,sep=""))
  # std.for.ace.ref <- bb2
   #
  # aa <- (aa1-aa2)
  # bb <- (bb1+bb2)/2.
  # if (is.na(bb) != TRUE) {
  #    eff.size.ace.for <- aa/bb
  # }else{
  #    eff.size.ace.for <- NA
  # } 
  # print(paste("Effect size of ET in forest area:", eff.size.ace.for,sep=""))

ld_go <- F
if (ld_go) {
   # calculate effect size of Precipitation 
   #for forest area
   aa1 <- mean(acf_af*eve.tc.aff.mask*lc.for.mask,na.rm=T) - mean(acf_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)
     bef.for.acf.aff <- mean(acf_be  * lc.for.mask * eve.tc.aff.mask, na.rm=T) 
     aft.for.acf.aff <- mean(acf_af  * lc.for.mask * eve.tc.aff.mask, na.rm=T)
     print( paste("Forests Precip before in affected area:",bef.for.acf.aff) ) 
     print( paste("Forests Precip after in affected area:", aft.for.acf.aff) ) 
   bb1 <- sqrt( (sd(acf_af*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2. + sd(acf_be*eve.tc.aff.mask*lc.for.mask,na.rm=T)**2.)/2.)
         print(paste("Precip STD pool, affected forest: ", bb1,sep=""))
   std.for.acf.aff <- bb1 
   #for referance area    
   aa2 <- mean(acf_af*eve.tc.ref.mask*lc.for.mask,na.rm=T) - mean(acf_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)
     bef.for.acf.ref <- mean(acf_be  * lc.for.mask * eve.tc.ref.mask, na.rm=T) 
     aft.for.acf.ref <- mean(acf_af  * lc.for.mask * eve.tc.ref.mask, na.rm=T)
     print(paste("Forests Precip before in reference area:",bef.for.acf.ref) )
     print(paste("Forests Precip after in reference area:",  aft.for.acf.ref) )
   bb2 <- sqrt( (sd(acf_af*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2. + sd(acf_be*eve.tc.ref.mask*lc.for.mask,na.rm=T)**2.)/2.)
         print(paste("Precip STD pool, reference forest: ", bb2,sep=""))
   std.for.acf.ref <- bb2
   #
   aa <- (aa1-aa2)
   bb <- (bb1+bb2)/2.
   if (is.na(bb) != TRUE) {
      eff.size.acf.for <- aa/bb
   }else{
      eff.size.acf.for <- NA
   } 
   print(paste("Effect size of Precip in forest area:", eff.size.acf.for,sep=""))
}

   #calculate effect size for agricuture fields!  
   #agriculture affected area 
   aa1 <- mean(lai_af*eve.tc.aff.mask*lc.agr.mask,na.rm=T) - mean(lai_be*eve.tc.aff.mask*lc.agr.mask,na.rm=T)
     bef.agr.lai.aff <- mean(lai_be  * lc.agr.mask * eve.tc.aff.mask, na.rm=T)
     aft.agr.lai.aff <- mean(lai_af  * lc.agr.mask * eve.tc.aff.mask, na.rm=T)
     print(paste("Agr LAI before in affected area:", bef.agr.lai.aff))
     print(paste("Agr LAI after in affected area:",  aft.agr.lai.aff))
   bb1 <- sqrt( (sd(lai_af*eve.tc.aff.mask*lc.agr.mask,na.rm=T)**2. + sd(lai_be*eve.tc.aff.mask*lc.agr.mask,na.rm=T)**2.)/2.)
         print(paste("LAI STD pool:", bb1,sep=""))
   std.agr.lai.aff <- bb1

   #agriculture reference area   
   aa2 <- mean(lai_af*eve.tc.ref.mask*lc.agr.mask,na.rm=T) - mean(lai_be*eve.tc.ref.mask*lc.agr.mask,na.rm=T)
     bef.agr.lai.ref <- mean(lai_be  * lc.agr.mask * eve.tc.ref.mask, na.rm=T)
     aft.agr.lai.ref <- mean(lai_af  * lc.agr.mask * eve.tc.ref.mask, na.rm=T)
     print(paste("Agr LAI before in reference area:", bef.agr.lai.ref) )
     print(paste("Agr LAI after in reference area:",  aft.agr.lai.ref) )
   bb2 <- sqrt( (sd(lai_af*eve.tc.ref.mask*lc.agr.mask,na.rm=T)**2. + sd(lai_be*eve.tc.ref.mask*lc.agr.mask,na.rm=T)**2.)/2.)
     print(paste("LAI STD pool:", bb2,sep=""))
   std.agr.lai.ref <- bb2
 
   aa <- (aa1-aa2)
   bb <- (bb1+bb2)/2.
   
   if (is.na(bb) != TRUE) {
      eff.size.agr <- aa/bb
   }else{
      eff.size.agr <- NA
   } 
         print(paste("Effect size in  agriculture area:", eff.size.agr,sep=""))


  # calculate the mean soil water befor the TC event over 30 days
  
   bef.for.asw.aff <- mean( eve.tc.aff.mask*lc.for.mask*asw.b1$soil_water, na.rm=T )
  

  # search for the maximum windspeec for the event 
   eve.table  <-  subset (track.data, ((substr(track.data$YYYYMMDDHH,start=1,stop=4) == as.character(wrk.yr)) &
                                          (track.data$CY == tc.eve[ieve]) ))
   eve.intensity <- max(eve.table$VMAX)  
   

  # event base table 
  if ( is.na(eff.size.for) != TRUE ) {
     table.tmp <- data.frame( date.time = wrk.date,
                              date.yr = wrk.yr,
                              date.mon = wrk.mon,
                              eve.id = tc.eve[ieve],
                              eve.intensity = eve.intensity, 
                              eve.lat.med.aff = eve.lat.med.aff,
                              eve.lon.med.aff = eve.lon.med.aff,
                              eve.lat.max.aff = eve.lat.max.aff,
                              eve.lat.min.aff = eve.lat.min.aff,
                              eve.lon.max.aff = eve.lon.max.aff,
                              eve.lon.min.aff = eve.lon.min.aff,
                              eve.tc.occ.aff  = eve.tc.occ.aff,
                              eve.for.pix.aff = eve.for.pix.aff,
                              eve.agr.pix.aff = eve.agr.pix.aff,
                              eve.oth.pix.aff = eve.oth.pix.aff,
                              eve.for.mws.aff = eve.for.mws.aff,
                              eve.agr.mws.aff = eve.agr.mws.aff,
                              eve.oth.mws.aff = eve.oth.mws.aff,
                              eve.for.acf.aff = eve.for.acf.aff,
                              eve.agr.acf.aff = eve.agr.acf.aff,
                              eve.oth.acf.aff = eve.oth.acf.aff,
                              bef.for.acf.aff = bef.for.acf.aff,
                              bef.agr.acf.aff = bef.agr.acf.aff,
                              bef.oth.acf.aff = bef.oth.acf.aff,
                              bef.for.lai.aff = bef.for.lai.aff,
                              aft.for.lai.aff = aft.for.lai.aff,
                              std.for.lai.aff = std.for.lai.aff, 
                              bef.agr.lai.aff = bef.agr.lai.aff,
                              aft.agr.lai.aff = aft.agr.lai.aff,
                              std.agr.std.aff = std.agr.lai.aff,                           
                              bef.for.asw.aff = bef.for.asw.aff,
                              eve.lon.med.ref = eve.lon.med.ref,
                              eve.lat.med.ref = eve.lat.med.ref,
                              eve.lat.max.ref = eve.lat.max.ref,
                              eve.lat.min.ref = eve.lat.min.ref,
                              eve.lon.max.ref = eve.lon.max.ref,
                              eve.lon.min.ref = eve.lon.min.ref,
                              eve.for.pix.ref = eve.for.pix.ref,
                              eve.agr.pix.ref = eve.agr.pix.ref,
                              eve.oth.pix.ref = eve.oth.pix.ref,
                              eve.for.mws.ref = eve.for.mws.ref,
                              eve.agr.mws.ref = eve.agr.mws.ref,
                              eve.oth.mws.ref = eve.oth.mws.ref,
                              eve.for.acf.ref = eve.for.acf.ref,
                              eve.agr.acf.ref = eve.agr.acf.ref,
                              eve.oth.acf.ref = eve.oth.acf.ref,
                              eve.for.off.ref = eve.for.off.ref,
                              eve.agr.off.ref = eve.agr.off.ref,
                              eve.oth.off.ref = eve.oth.off.ref,
                              eve.for.spei.aff = eve.for.spei.aff,
                              eve.for.spei.ref = eve.for.spei.ref,
                              bef.for.lai.ref = bef.for.lai.ref,
                              aft.for.lai.ref = aft.for.lai.ref,
                              std.for.lai.ref = std.for.lai.ref, 
                              bef.agr.lai.ref = bef.agr.lai.ref,
                              aft.agr.lai.ref = aft.agr.lai.ref,
                              std.agr.std.ref = std.agr.lai.ref,                           
                              eff.size.for = eff.size.for,
                              eff.size.agr = eff.size.agr,
                              tot.n = tot.n,
                              qc1.score = qc1.score,
                              qc2.score = qc2.score,
                              qc1.flag = qc1.flag,
                              qc2.flag = qc2.flag  )
     print(table.tmp)
     # combine the table
     table.comb <- rbind( table.comb, table.tmp ) 

     print( paste( "Finish analysisng TC event: ",tc.eve[ieve],"; date:",wrk.date ,sep="") )
  
    # plot the results
     ld_go <- FALSE
     if(ld_go) {
     #dev.new()
 
     my.col    <- colorRampPalette(c("darkblue","blue","blue","lightblue","white","pink","red","red","brown"))(20)
     my.breaks <- seq(-1.25, 1.25, length.out = 21)
     # 
     # bitmap(file=paste(track.dir,"LAI_delta_mask_",target,"_",wrk.date,".png",sep=""),type="png16m",
     #        width = 8, height = 8, units = "in", res= 600 )
     layout(matrix(data=seq(1,6,1),nrow=3, ncol=2, byrow=TRUE))

     # mask plot 
     image.plot(eve.tc.aff.mask[nx1:nx2,ny2:ny1], main=paste("TC event:",tc.eve[ieve],"TC affected  mask",sep=""))
     image.plot(eve.tc.ref.mask[nx1:nx2,ny2:ny1], main=paste("TC event:",tc.eve[ieve],"TC reference mask",sep=""))
     # wind and rainfall plot 
     image.plot(tc.mws$max_wind[nx1:nx2,ny2:ny1]*lc.for.mask[nx1:nx2,ny2:ny1], main="10 days maximum wind fields (m/s)")
     image.plot(tc.acf$acc_rainf[nx1:nx2,ny2:ny1]*lc.for.mask[nx1:nx2,ny2:ny1], main="10 days accumulative rainfall (mm)")
     #
     ## forest area
     del_lai <- (lai_af - lai_be ) * lc.for.mask * eve.tc.aff.mask
     # gain.pix.for <- sum( del_lai[which(del_lai >= cut.off)], na.rm=T )  
     #loss.pix.for <- -1*sum( del_lai[which(del_lai <= -1*cut.off)], na.rm=T ) 
     # check plots
     image.plot(del_lai[nx1:nx2,ny2:ny1], main="TC affected forest area LAI change", zlim=c(-1.25,1.25) ,  col=my.col,breaks=my.breaks)
     # unaffected forest area
   
     # unaffected forest area
     del_lai <- (lai_af - lai_be ) * lc.for.mask* eve.tc.ref.mask
     #del_lai <- (lc.for.mask*tc.pot.mask)*(lai_af-lai_be)  
     #gain.pix.for.pot <- sum( del_lai[which(del_lai >= cut.off)], na.rm=T )  
     #loss.pix.for.pot <- -1*sum( del_lai[which(del_lai <= -1*cut.off)], na.rm=T ) 
     # check plots 
     image.plot(del_lai[nx1:nx2,ny2:ny1], main="TC reference forest area LAI change",zlim=c(-1.25,1.25),  col=my.col,breaks=my.breaks,
              xlab= paste("Obs. date:", wrk.date,sep="") )

     par(xpd=NA)
    #  text(x=0.5, y=1.2, paste("Obs. date:", wrk.date,sep=""))
  
    # close png plot 
    #dev.off()
   } #if ld_go 




  } else {
     print("No enough samples for  eff.size calculation!") 
  }

 } #TC event loop  

  ld_go <- FALSE
  if (ld_go) { 
  layout(matrix(data=seq(1,6,1),nrow=2, ncol=3, byrow=TRUE))
  plot( table.comb$nav.lat ~ table.comb$zonal.gain.for,xlim=c(0,100),col="green")
  plot( table.comb$nav.lat ~ table.comb$zonal.gray.for,xlim=c(0,100),col="gray")
  plot( table.comb$nav.lat ~ table.comb$zonal.loss.for,xlim=c(0,100),col="red")

  plot( table.comb$nav.lat ~ table.comb$zonal.gain.for.pot,xlim=c(0,100),col="green")
  plot( table.comb$nav.lat ~ table.comb$zonal.gray.for.pot,xlim=c(0,100),col="gray")
  plot( table.comb$nav.lat ~ table.comb$zonal.loss.for.pot,xlim=c(0,100),col="red")

  }
# calculate the accumulative LAI 
# set cut.off value
# agriculture area
#  del_lai <- (lai_af - lai_be ) * lc.agr.mask * tc.aff.mask
#  gain.pix.agr <- sum( del_lai[which(del_lai >= cut.off)], na.rm=T )  
#  loss.pix.agr <- -1*sum( del_lai[which(del_lai <= -1*cut.off)], na.rm=T ) 

#} #end of if for wrk.mon in summer conditions  

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

# set lat as factor for analysis
#  table.comb$lat <- as.factor(table.comb$nav.lat) 
# save the table for analyis
  st.date <- substr( lai.fnames[it1], start=11, stop=18)
  ed.date <- substr( lai.fnames[it2], start=11, stop=18)
#  home.dir <- c("/data1/home/ychen/")
  #save(table.comb, file = paste(target,"_",st.date,"to",ed.date,"_table.comb.rda",sep="") )
  save(table.comb, file = paste("./",spei.fname,"/",target,"_",st.date,"to",ed.date,"_table.comb",
       "w_",Ref.wind,"_p_",Ref.rain ,"_",pos.txt,"_",off.txt,"_cut_",cut.off,"_tc.occ.rda",sep="") )

  return(table.comb)

} # end of fun_dyn_lai

#========== Set the script to Auto RUN=========================== 
# Start convert the file by fun_compress with specific argunment
# If you want to submit the file by queue to obelix, you need to apply the 
# following lines, which allowed this R-script can be called under the shell script
# with the arggunmets sending by specific batch jobs  
#
args<-commandArgs(TRUE)
print(args)
fun.dyn.lai(args[1], args[2], args[3], args[4], args[5], args[6]) 
#
#========== End ================================================



