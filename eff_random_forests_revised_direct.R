

#####################################
# Quality Check function 
#####################################    
#QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.2, area.set=0) {
 
      #print( paste("Total events:", length(wrk.table$eve.for.pix.aff),sep=""))
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      all.eve <- length(wrk.table$eve.for.pix.aff)
      print( paste("Total events:", all.eve,sep=""))
 
      #QC1 mean difference  on LAI #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set) )
      tmp.eve <- length(wrk.table$eve.for.pix.aff)
      flag1.eve <- all.eve - length(wrk.table$eve.for.pix.aff)
      print( paste("Flag-1 events (remove Big LAI difference):", flag1.eve ,sep=""))
      
      wrk.table <- subset(wrk.table, ( (abs(wrk.table$qc1.score)<=qc1.set) &  (abs(wrk.table$eff.size.for) >= abs(wrk.table$qc2.score)))) 
      flag2.eve <- length(wrk.table$eve.for.pix.aff) 
      print( paste("Flag-2 events (Noise > difference & small LAI difference):", flag2.eve,sep=""))


      return(wrk.table)
      } 
####################################

#####################################
# Quality Check function (revised version by using 10% of difference between reference and TC-affected area 
#####################################    
#QC1 mean difference on LAI 10% of ref LAI
#QC2 mean change on LAI caused by TC  10% of ref LAI  
#fun_table.qc <- function(wrk.table, qc1.set=0.1, qc2.set=0.1, area.set=0) {
 
      #print( paste("Total events:", length(wrk.table$eve.for.pix.aff),sep=""))
      #subset the datset with QC and QA score
      #Avariable area
#      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
#      all.eve <- length(wrk.table$eve.for.pix.aff)
#      print( paste("Total events:", all.eve,sep=""))
 
      #QC1 mean difference on LAI 
      #QC2 mean change on LAI  
      #qc1.set <- 0.1
      #qc2.set <- 0.1
#      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score) <= qc1.set) )
#      tmp.eve <- length(wrk.table$eve.for.pix.aff)
#      flag1.eve <- all.eve - length(wrk.table$eve.for.pix.aff)
#      print( paste("Flag-1 events (remove Big LAI difference):", flag1.eve ,sep=""))
      
#      lai.c.aff <- wrk.table$bef.for.lai.aff  -  wrk.table$aft.for.lai.aff  
#      lai.c.ref <- wrk.table$bef.for.lai.ref  -  wrk.table$aft.for.lai.ref 
      
#      qc2.score.eq3 <- abs(  (lai.c.aff - lai.c.ref) / wrk.table$bef.for.lai.ref )     
     
      #add information of qc2 into the table  
#      wrk.table$qc2.score.eq3 <- qc2.score.eq3  
      
 
#      wrk.table <- subset(wrk.table, ( (abs(wrk.table$qc1.score)<=qc1.set) &  (abs(wrk.table$eff.size.for) <= qc2.set )) )
#      flag2.eve <- length(wrk.table$eve.for.pix.aff) 
#      print( paste("Flag-2 events (Noise > difference & small LAI difference):", flag2.eve,sep=""))


#      return(wrk.table)
#      } 
####################################





#set qc/qa parameters
qc1.set <- 0.5
qc2.set <- 0.1
area.set <- 0



#data.set <- c("./Rda_60/")
data.set <- c("./spei02/BG_R1/")

#data.set <- c("./spei02/BG_R1_QC/")


#get file-list 
in.files <- list.files(path=data.set, pattern="EA_*")
   
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   load(paste(data.set,in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
  #add case information to the table
  table.comb$run.case <- sub("_pos.*","",sub(".*table.comb","",in.files[id])) 
 
  # combine all table
   table.all <- rbind(table.comb,table.all) 
 
   #scatter.smooth(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30 ),
   #               col="lightgray", span=0.5, degree=1, 
   #               lpars=list(lwd=3,col="gray"), ylim=c(-0.5,0.5), xlab="",ylab="Effect size",  xaxt = "n")
}

# remove the NA cases
table.all <- table.all[complete.cases(table.all), ]
#table.all <- subset(table.all, table.all$eve.for.mws.aff>1.0)

# subset datataset (this is super critical!! we develop this QC flags)
table.all <- fun_table.qc(wrk.table=table.all, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)
 

library("dplyr")
library("randomForest")
#library("caret")


# read the enso index table 
enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)
pj.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ERA5_DATA/PJ_INDEX/PJ_index.csv",sep=",",header=T)  

# assign the positve/negtive value of eff_size as group B/A
table.all$group <- NA
table.all$group[table.all$eff.size.for > 0.] <- "Growth"
table.all$group[table.all$eff.size.for < 0.] <- "Damage"
table.all$group <- as.factor(table.all$group)

#assign the  index to the TC event 
table.all$oni.enso <- NA
table.all$pj.index <- NA

#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$oni.enso[it] <- enso.table$Ocean_Nino_Index[which(enso.table$Year == table.all$date.yr[it]  
                                                         & enso.table$Mon == table.all$date.mon[it] ) ] 
   table.all$pj.index[it] <- pj.table$PJ_index[which(pj.table$year == table.all$date.yr[it] 
                                               & pj.table$mon == table.all$date.mon[it] ) ] 
}


#define affect area 
eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex


#table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))
table.all$date.mon.day <-  as.numeric(table.all$date.mon+ table.all$date.yr) 





#mat1 <-matrix(c(1,1,2,3),2,2)
#layout(mat1, widths=c(2,2),height=c(2,2))

#hist.breaks <- c(seq(-5,5.0,length=100))
#hist.col <- ifelse( hist.breaks < 0, "orange2", "forestgreen") 

#hist(table.all$eff.size.for, col=hist.col
#     ,breaks=hist.breaks, main="All events", xlim=c(-1,1) )
#hist(table.all$eff.size.for[ table.all$enso >= 0.5 ], col= hist.col,
#     ,breaks=hist.breaks, main="El Nino (positive phase events)",xlim=c(-1,1))
#hist(table.all$eff.size.for[ table.all$enso <= -0.5 ], col= hist.col,
#     breaks=hist.breaks, main="La Nina (negtive phase events)",xlim=c(-1,1))


table.all <- subset( table.all, table.all$eve.for.mws.aff >1)  


eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100

# define the change of droughts index  
table.all$del.spei.for <- as.numeric(table.all$eve.for.spei.aff - table.all$eve.for.spei.ref)


ld.go <- T

if (ld.go) {


# prepare table for analysis

#table.ana <- select(table.all , group, eve.for.acf.aff, eve.lat.med.aff, eve.intensity, bef.for.asw.aff, date.mon.day, oni.enso)




# subset table
imp.table <- NA 

for (i in 1:10) {

#
#)

# "w_0_p_100_4D"  "w_0_p_60_2D"   "w_0_p_80_3D"   "w_10_p_0_3D"   "w_10_p_80_3D"  "w_12_p_0_4D"   "w_12_p_100_4D" "w_8_p_0_2D"    "w_8_p_60_2D"
#2D definition
if (i==1) table.ana <- subset(table.all , table.all$run.case=="w_0_p_60_2D")
if (i==2) table.ana <- subset(table.all , table.all$run.case=="w_8_p_0_2D")
if (i==3) table.ana <- subset(table.all , table.all$run.case=="w_8_p_60_2D")

#3D definition
if (i==4) table.ana <- subset(table.all , table.all$run.case=="w_0_p_80_3D")
if (i==5) table.ana <- subset(table.all , table.all$run.case=="w_10_p_0_3D")
if (i==6) table.ana <- subset(table.all , table.all$run.case=="w_10_p_80_3D")

#4D definition
if (i==7) table.ana <- subset(table.all , table.all$run.case=="w_0_p_100_4D")
if (i==8) table.ana <- subset(table.all , table.all$run.case=="w_12_p_0_4D")
if (i==9) table.ana <- subset(table.all , table.all$run.case=="w_12_p_100_4D")

#234D definition
if (i==10) table.ana <- table.all 

# prepare table for analysis
table.ana <- select(table.ana , group, 
                                eve.for.mws.aff,eve.for.acf.aff, eve.lat.med.aff,  eve.intensity, eff.cex, eve.tc.occ.aff,
                                date.mon.day, bef.for.acf.aff,  bef.for.lai.aff, pj.index, eve.for.spei.ref, del.spei.for)


# Set random seed to make results reproducible:
set.seed(17*i)


# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(table.ana))
# Mote Carlo sampling with n/2 of total data_set_size
#data_set_size <- as.integer(data_set_size/2)
#data_set_size = 131 

# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(table.ana), size = data_set_size)


# Assign the data to the correct sets
training <- table.ana[indexes,]


#Apply Classification with Random Forest:
# In the randomForest Package this would mean that you make a classification (type = "classification"), ntree = 2000 means that 2000 prediction trees are grown

rf_class = randomForest(group ~ ., data=training, ntree=666, mtry=20 , maxnodes=20,
                                    norm.votes=TRUE,  importance=TRUE,
                                    proximity=TRUE, type="classification", na.action=na.omit)
#
# Set your workingdirectory for the output here
#pdf(file='all_trees.pdf', paper = "a4r", width = 10) #Name of the PDF File. In this case it will be a PDF with 2000 pages
#for(k in 100:110){
  # the getTree function is part of the randomForest package; it extracts the structure of a tree from a randomForest object - I hope the biomod2 package has the same functionality
#  tree <- getTree(rf_class, k=k, labelVar = TRUE)
#  str(tree)
#  print(sum(tree[,1] +tree[,2]))
#  d <- to.dendrogram(tree)
#  str(d)
#  plot(d,center=TRUE,edgePar=list(t.cex=.55,p.col=NA,p.lty=0), yaxt = "n", 
#  main = paste("Tree Numbe:",k," Accuracy of classification: ",round(1-rf_class$err.rate[k],digits = 2)*100,"%"))
  #dev.new()
#}
#dev.off()

#dev.new()
#plot importance 
#varImpPlot(rf_class)
#export the importance from the model 

tmp.imp.table <- as.data.frame(importance(rf_class)) 
tmp.imp.table$vars <- row.names(tmp.imp.table)

print(tmp.imp.table)
# combine table 
imp.table <- rbind(tmp.imp.table, imp.table, make.row.names=FALSE)

}


# rm NA cases
imp.table <- imp.table[complete.cases(imp.table), ]
# set vars as factor 
imp.table$vars <- as.factor(as.character(imp.table$vars))


#if (data.set == "./Rda_60/") {
# adjust the orders 

# eve.for.acf.aff eve.for.spei.aff         date.mon  eve.lat.med.aff  bef.for.lai.aff  bef.for.acf.aff  eve.for.mws.aff          eff.cex    eve.intensity         pj.index
#      21.6946124       -1.0996020        7.9609186        1.7772625        4.9584312        3.3043613        2.0168017        0.7440358        3.4494686       10.6475595

imp.table$vars <- factor(imp.table$vars , 
                   levels=c("eve.for.acf.aff","pj.index","date.mon","bef.for.lai.aff","eve.intensity",
                            "bef.for.acf.aff","eve.for.mws.aff","eve.lat.med.aff","eff.cex", "eve.for.spei.aff"))
#
var.names <- rev(c("Accumulated rainfall(b)","Pacific Japan index","Month of landfall", "Prior leaf area(c)", "Cyclone intensity(a)", 
                   "Prior accumulated rainfall(c)","Maximum wind speed(b)", "Latitude of landfall(a)", "Affected area","Prior standarized precipitation evapotranspiration index(c)")) 
#
var.group <- rev(c("TCC","PC","PC","PC","TCC",
                    "PC","TCC","TCC","TCC","PC"))
#} else if(data.set == "./Rda_50") {


#imp.table$vars <- factor(imp.table$vars , 
#                   levels=c("eve.lat.med.aff","eve.for.acf.aff","eff.cex",
#                            "bef.for.acf.aff","bef.for.lai.aff","date.mon.day","eve.for.mws.aff",
#                            "eve.intensity","bef.for.asw.aff","oni.enso") )
#
#var.names <- c("Latitude(a)","AccRain(b)","AffArea",
#               "PriAccRain(c)","PriLAI(c)","Month","MaxWind(b)",
#               "IntTC(a)","PriSoilMoi(c)","ONI")

#var.group <- c("TCC","TCC","TCC",
#               "PC","PC","PC","TCC",
#               "TCC","PC","PC")
#} 

# get a new order of factor based on median of data
f.name <- imp.table$var
f.data <- imp.table$MeanDecreaseAccuracy
tmp.table <- data.frame(f.name, f.data)

new.order.table <- with(tmp.table, reorder(f.name, f.data, median, na.rm=T))

#boxplot(tmp.table$f.data ~ new.order.table)


go.pdf = F 
if(go.pdf) {pdf(file="Fig3_updated_PJ.pdf",width=12, height=9) } else{ print("go x-windows") }


#adjust label size 
par(cex.axis=0.5)

boundaries <- boxplot( tmp.table$f.data ~ new.order.table, ylim=c(-65,40),axes=F,main="",
              yaxt="n", medcol="black", horizontal=TRUE, outline=FALSE, cex.axis=1.5, cex.lab=1.5,cex.main=2.0,boxlwd=1.5,
              col=ifelse( var.group == "TCC", "lightgray", "white"))
#set axis texts 
# "Mean Decrease Accuracy"

axis(side=3, label=F, lwd.tick=0,at=c(-65,45),lwd=2) # axis, ticks
abline(v=-8,lwd=1.0)
abline(h=0.1,lwd=2.0)
axis(side=1, at=c(0,10,20,30,40), labels=c("0","10","20","30","40"), cex.axis=1.5)
par(xpd=T)
text(x=15,y=-1.0, label="Decrease in accuracy(%)",cex=1.5, font=2.0)
par(xpd=F)
#add legend to the plot 
#legend("topright",inset=.05, cex=1.5,title="",,box.lwd = 0,box.col = "white",bg = "white", 
#   c("TC Characteristic","Prior Condition"), fill=c("forestgreen","lightgray"), horiz=FALSE)

# Add Text on top of the box
nbGroup <- nlevels(imp.table$vars)
text(pos=4,# text to the right of reference coordinate   
  y=c(1:nbGroup), 
#  x=boundaries$stats[nrow(boundaries$stats),] + 7, 
  x= -65,
   paste(var.names[1:nbGroup],sep=""),cex=2.0)

#boxplot( tmp.table$f.data ~ new.order.table)

#import the package
#library("randomForest")
#library("reprtree")

# Perform training:
#rf_class = randomForest(group ~ ., data=training, ntree=100, mtry=2, importance=TRUE)
# Plot the tree 
#reprtree:::plot.getTree(rf_class)


#test <- training
#library(ggplot2)
#pGroup <- predict(rf_class,test,'vote')[,2]
#plotData <- lapply(names(test[,2:5]), function(x){
#  out <- data.frame(
#    var = x,
#    type = c(rep('Actual',nrow(test)),rep('Predicted',nrow(test))),
#    value = c(test[,x],test[,x]),
#    group = c(as.numeric(test$group)-1,pGroup)
#    )
#  out$value <- out$value-min(out$value) #Normalize to [0,1]
#  out$value <- out$value/max(out$value)
#  out
#})
#plotData <- do.call(rbind,plotData)
##qplot(value, group, data=plotData, facets = type ~ var, geom='smooth', span = 0.5)


} #ld_go

if(go.pdf) {dev.off()}



#dev.new()

# performing regrssion trees analysis based on the random forest analysis (assign the importance variables) 
#  
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees

library(rsample)     # data splitting 
library(dplyr)       # data wrangling


table.ana$group <- as.character(table.ana$group)
table.ana$group <- as.factor(table.ana$group)

# build the regression tree model
#tc.m1 <- rpart(formula = group ~.,  data = table.ana  )     
#rpart.plot (tc.m1)





