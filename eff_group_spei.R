## FUNCTIONS

#####################################
# Quality Check funtion 
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
 
      #QC1 mean difference on LAI #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      flag1.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set) )
      tmp.eve <- length(flag1.table$eve.for.pix.aff)
      flag1.eve <- all.eve - length(flag1.table$eve.for.pix.aff)
      print( paste("Flag-1 events (remove Big LAI difference):", flag1.eve ,sep=""))
      
      flag2.table <- subset(wrk.table, ((abs(wrk.table$qc2.score)<qc2.set) & (abs(wrk.table$qc1.score)<qc1.set) ))
      flag2.eve <- tmp.eve - length(flag2.table$eve.for.pix.aff) 
      print( paste("Flag-2 events (Neutral events):", flag2.eve,sep=""))
      fctor=1.25
      wrk.table <- subset(wrk.table, ( (abs(wrk.table$qc1.score)<=qc1.set) &  (abs(wrk.table$eff.size.for)*fctor >= abs(wrk.table$qc2.score)))) 
  
      return(list(qcqa=wrk.table,flag1=flag1.table,flag2=flag2.table))
      } 
####################################

#set qc/qa parameters
qc1.set <- 0.525
qc2.set <- 0.2
area.set <- 0

#rda.dir <- c("./Rda_60_tc.occ/")

rda.dir <- c("./spei02/")
spei.fname <- substr(rda.dir, start=3,stop=8)
#

# open the device for the figure output
#pdf(file="SPEI_plot.pdf", width=8,height=16)


library("dplyr")
library("randomForest")
library("caret")
#library("e1071")
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("tidyverse")

#get file-list 
#in.files <- list.files(path="./Rda_60/", pattern="EA_*")
in.files <- list.files(path=rda.dir, pattern="*_[2]D") 
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   print(in.files[id])
   load(paste(rda.dir,in.files[id],sep=""))
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
raw.table <- table.all
table.all <- subset(table.all, table.all$eve.for.mws.aff>0.0)
# subset datataset (this is super critical!! we develop this QC flags)
# combine qcqc table and flag2 table (neutral)  
table.tmp <- fun_table.qc(wrk.table=table.all, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)
# assign the positve/negtive value of eff_size as group B/A
table.tmp[['qcqa']]$group <- NA
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for >  0.18] <- "Positive"
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for < -0.18] <- "Negative"
table.tmp[['qcqa']]$group[abs(table.tmp[['qcqa']]$eff.size.for) <= 0.18] <- "Neutral"
table.tmp[['flag2']]$group <- "Neutral"
table.all <- rbind(table.tmp[['qcqa']],table.tmp[['flag2']])  
#table.all <- table.tmp[['qcqa']]
table.all$group <- as.factor(table.all$group)
# read the Oceanic Nino Index table 
oni.enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)
#assign the  index to the TC event 
table.all$oni.enso <- NA
#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$oni.enso[it] <- oni.enso.table$Ocean_Nino_Index[which(oni.enso.table$Year == table.all$date.yr[it]  
                                                         & oni.enso.table$Mon == table.all$date.mon[it] ) ] 
}
# convert the unit to the radius 
eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex
# add date stamp
table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))
# create table for spei analysis
#neg.pos.table <- subset(table.all, table.all$group!="Neutral" &  abs(table.all$eff.size.for) > 0.45)
neg.pos.table <- subset(table.all, abs(table.all$eff.size.for) > 0.3)


#calculate the median value
med_spei_df <- neg.pos.table %>%
  group_by(group) %>%
  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))

#spei.table <- table.all
p1 <- ggplot(neg.pos.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
p1 <- p1 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group), size=0.8, linetype="dashed")
p1 <- p1 + geom_density(alpha=0.7) 
p1 <- p1 + scale_color_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p1 <- p1 + scale_fill_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p1 <- p1 + theme(legend.position="top") 
p1 <- p1 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
p1 <- p1 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
p1 <- p1 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
p1 <- p1 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
p1 <- p1 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
p1 <- p1 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
p1 <- p1 + theme(legend.text = element_text(size = rel(1.0) ) ) 
p1 <- p1 + theme(legend.title = element_text(size = rel(1.0) ) ) 
p1 <- p1 + labs(x = expression(delta~"SPEI, 2D"), y = "Density", tag="A")
p1 <- p1 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
p1 <- p1 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
p1 <- p1 + guides(fill=guide_legend(title=""), color=guide_legend(title=""))


#  create table for spei analysis
#neu.table <- subset(table.all, table.all$group=="Neutral"  &  abs(table.all$eff.size.for) <= 0.3)
#calculate the median value
#med_spei_df <- neu.table %>%
#  group_by(group) %>%
#  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))
#p2 <- ggplot(neu.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
#p2 <- p2 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group), size=0.8, linetype="dashed")
#p2 <- p2 + geom_density(alpha=0.7) +  xlim(-1.5,1.5) + ylim(0,1.8)
#p2 <- p2 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
#p2 <- p2 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
#p2 <- p2 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
#p2 <- p2 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
#p2 <- p2 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
#p2 <- p2 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
#p2 <- p2 + theme(legend.text = element_text(size = rel(1.0) ) ) 
#p2 <- p2 + theme(legend.title = element_text(size = rel(1.0), face="bold" ) ) 
#p2 <- p2 + labs(x = expression(delta~"SPEI, 2D"), y = "", tag="A")
#p2 <- p2 + scale_color_discrete(name = "Effect Size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
#p2 <- p2 + scale_fill_discrete(name = "Effect size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
#p2 <- p2 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p2 <- p2 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p2 <- p2 + theme(legend.position="top") 
#p2 <- p2 + guides(fill=guide_legend(title="Effect Size Group:"), color=guide_legend(title="Effect Size Group:"))

#



#get file-list 
#in.files <- list.files(path="./Rda_60/", pattern="EA_*")
#in.files <- list.files(path=rda.dir, pattern="*_[2]D") 
in.files <- list.files(path=rda.dir, pattern="*_[3]D") 
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   print(in.files[id])
   load(paste(rda.dir,in.files[id],sep=""))
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
raw.table <- table.all
table.all <- subset(table.all, table.all$eve.for.mws.aff>0.0)
# subset datataset (this is super critical!! we develop this QC flags)
# combine qcqc table and flag2 table (neutral)  
table.tmp <- fun_table.qc(wrk.table=table.all, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)
# assign the positve/negtive value of eff_size as group B/A
table.tmp[['qcqa']]$group <- NA
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for >  0.18] <- "Positive"
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for < -0.18] <- "Negative"
table.tmp[['qcqa']]$group[abs(table.tmp[['qcqa']]$eff.size.for) <= 0.18] <- "Neutral"
table.tmp[['flag2']]$group <- "Neutral"
table.all <- rbind(table.tmp[['qcqa']],table.tmp[['flag2']])  
#table.all <- table.tmp[['qcqa']]
table.all$group <- as.factor(table.all$group)
# read the Oceanic Nino Index table 
oni.enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)
#assign the  index to the TC event 
table.all$oni.enso <- NA
#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$oni.enso[it] <- oni.enso.table$Ocean_Nino_Index[which(oni.enso.table$Year == table.all$date.yr[it]  
                                                         & oni.enso.table$Mon == table.all$date.mon[it] ) ] 
}
# convert the unit to the radius 
eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex
# add date stamp
table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))
# create table for spei analysis
#neg.pos.table <- subset(table.all, table.all$group!="Neutral" &  abs(table.all$eff.size.for) > 0.45)
neg.pos.table <- subset(table.all, abs(table.all$eff.size.for) > 0.3)




#calculate the median value
med_spei_df <- neg.pos.table %>%
  group_by(group) %>%
  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))

p3 <- ggplot(neg.pos.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
p3 <- p3 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group) , size=0.8, linetype="dashed")
p3 <- p3 + geom_density(alpha=0.7) 
p3 <- p3 + scale_color_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p3 <- p3 + scale_fill_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p3 <- p3 + theme(legend.position="none") 
p3 <- p3 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
p3 <- p3 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
p3 <- p3 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
p3 <- p3 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
p3 <- p3 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
p3 <- p3 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
p3 <- p3 + theme(legend.text = element_text(size = rel(1.0) ) ) 
p3 <- p3 + theme(legend.title = element_text(size = rel(1.0) ) ) 
p3 <- p3 + labs(x = expression(delta~"SPEI, 3D"), y = "Density", tag="B")
p3 <- p3 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
p3 <- p3 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))

#  create table for spei analysis
#neu.table <- subset(table.all, table.all$group=="Neutral"  &  abs(table.all$eff.size.for) <= 0.3)
#calculate the median value
#med_spei_df <- neu.table %>%
#  group_by(group) %>%
#  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))
#p4 <- ggplot(neu.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
#p4 <- p4 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group),size=0.8, linetype="dashed" )
#p4 <- p4 + geom_density(alpha=0.7) +  xlim(-1.5,1.5) + ylim(0,1.8)
#p4 <- p4 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
#p4 <- p4 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
#p4 <- p4 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
#p4 <- p4 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
#p4 <- p4 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
#p4 <- p4 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
#p4 <- p4 + theme(legend.text = element_text(size = rel(1.0) ) ) 
#p4 <- p4 + theme(legend.title = element_text(size = rel(1.0), face="bold" ) ) 
#p4 <- p4 + labs(x = expression(delta~"SPEI, 3D"), y = "Probablity Density", tag="C")
##p4 <- p4 + scale_color_discrete(name = "Effect Size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
##p4 <- p4 + scale_fill_discrete(name = "Effect size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
#p4 <- p4 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p4 <- p4 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p4 <- p4 + theme(legend.position="none") 
#p4 <- p4 + guides(fill=guide_legend(title="Effect Size Group:"), color=guide_legend(title="Effect Size Group:"))

#


#get file-list 
#in.files <- list.files(path=rda.dir, pattern="EA_*")
#in.files <- list.files(path=rda.dir, pattern="*_[4]D") 
in.files <- c("EA_19990220to20181231_table.combw_0_p_100_4D_pos_060_off_000_cut_0_tc.occ.rda") 
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   print(in.files[id])
   load(paste(rda.dir,in.files[id],sep=""))
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
raw.table <- table.all
table.all <- subset(table.all, table.all$eve.for.mws.aff>0.0)
# subset datataset (this is super critical!! we develop this QC flags)
# combine qcqc table and flag2 table (neutral)  
table.tmp <- fun_table.qc(wrk.table=table.all, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)
# assign the positve/negtive value of eff_size as group B/A
table.tmp[['qcqa']]$group <- NA
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for >  0.18] <- "Positive"
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for < -0.18] <- "Negative"
table.tmp[['qcqa']]$group[abs(table.tmp[['qcqa']]$eff.size.for) <= 0.18] <- "Neutral"
table.tmp[['flag2']]$group <- "Neutral"
table.all <- rbind(table.tmp[['qcqa']],table.tmp[['flag2']])  
#table.all <- table.tmp[['qcqa']]
table.all$group <- as.factor(table.all$group)
# read the Oceanic Nino Index table 
oni.enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)
#assign the  index to the TC event 
table.all$oni.enso <- NA
#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$oni.enso[it] <- oni.enso.table$Ocean_Nino_Index[which(oni.enso.table$Year == table.all$date.yr[it]  
                                                         & oni.enso.table$Mon == table.all$date.mon[it] ) ] 
}
# convert the unit to the radius 
eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex
# add date stamp
table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))
# create table for spei analysis
#neg.pos.table <- subset(table.all, table.all$group!="Neutral" &  abs(table.all$eff.size.for) > 0.2)

#neg.pos.table <- table.all
neg.pos.table <- subset(table.all, abs(table.all$eff.size.for) > 0.2)



#calculate the median value
med_spei_df <- neg.pos.table %>%
  group_by(group) %>%
  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))

p5 <- ggplot(neg.pos.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
p5 <- p5 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group), size=0.8, linetype="dashed")
p5 <- p5 + geom_density(alpha=0.7) 
p5 <- p5 + scale_color_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p5 <- p5 + scale_fill_discrete(name = "") + xlim(-1.5,1.5) + ylim(0,1.8) 
p5 <- p5 + theme(legend.position="none") 
p5 <- p5 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
p5 <- p5 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
p5 <- p5 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
p5 <- p5 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
p5 <- p5 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
p5 <- p5 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
p5 <- p5 + theme(legend.text = element_text(size = rel(1.0) ) ) 
p5 <- p5 + theme(legend.title = element_text(size = rel(1.0) ) ) 
p5 <- p5 + labs(x = expression( delta ~"SPEI, 4D"), y = "Density", tag="C" )   
p5 <- p5 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
p5 <- p5 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#
# create table for spei analysis
#neu.table <- subset(table.all, table.all$group=="Neutral"  &  abs(table.all$eff.size.for) <= 0.3)
#
#med_spei_df <- neu.table %>%
#  group_by(group) %>%
#  summarize(median=median( eve.for.spei.aff - eve.for.spei.ref ))
#p6 <- ggplot(neu.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) 
#p6 <- p6 + geom_vline(data = med_spei_df, aes(xintercept = median, color = group), size=0.8, linetype="dashed")
#p6 <- p6 + geom_density(alpha=0.7)+ xlim(-1.5,1.5) + ylim(0,1.8)
#p6 <- p6 + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18)) 
#p6 <- p6 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90)) 
#p6 <- p6 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
#p6 <- p6 + theme(axis.text = element_text(size = rel(1.0)  )  ) 
#p6 <- p6 + theme(axis.text.x = element_text(size = rel(1.0) ) )  
#p6 <- p6 + theme(axis.text.y = element_text(size = rel(1.0) ) ) 
#p6 <- p6 + theme(legend.text = element_text(size = rel(1.0) ) ) 
#p6 <- p6 + theme(legend.title = element_text(size = rel(1.0), face="bold" ) ) 
#p6 <- p6 + labs(x = expression(delta ~"SPEI, 4D"), y = "", tag="E" )
#p6 <- p6 + scale_color_discrete(name = "Effect Size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
#p6 <- p6 + scale_fill_discrete(name = "Effect size Group:") + xlim(-1.5,1.5) + ylim(0,1.8) 
#p6 <- p6 + scale_fill_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p6 <- p6 + scale_color_manual( breaks=c("Negative","Neutral","Positive"), values=c("orange2","lightblue","forestgreen"))
#p6 <- p6 + theme(legend.position="none")  

# show the multiple plot to the pdf page 
grid.arrange(p1, p3, p5, nrow=3, ncol=1 )  # Make a list of plots
pdf.fname <- paste("NEW_SPEI_Change_",spei.fname,".pdf",sep="")
ggsave( grid.arrange( p1,  p3 ,p5, nrow=3, ncol=1),
        filename=pdf.fname,  
        width=5, height=8) 
#dev.off()
#p2<-ggplot(spei.table, aes(x=eve.for.spei.ref+eve.for.spei.aff, fill=group)) + geom_density(alpha=0.4)
#p2


# output the statistic of different runs with differnet sinarious

#get runs
runs <- levels( as.factor(table.all$run.case) ) 


tot<-array(0,dim=length(runs))
neg<-array(0,dim=length(runs))
neu<-array(0,dim=length(runs))
pos<-array(0,dim=length(runs)) 

for ( irun in 1:length(runs) ) {  
  table.sub <- subset(table.all, table.all$run.case == as.character(runs[irun]))  
  
  tot[irun] <- length(table.sub$run.case)
  neg[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Negative"] ) 
  neu[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Neutral"] ) 
  pos[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Positive"] ) 

  print( paste("Run:", runs[irun] ))
  print( paste("Totoal n:",  tot[irun] ))  
  print( paste("Negative:",  neg[irun] ))  
  print( paste("Neutral:",   neu[irun] ))  
  print( paste("Postive:",   pos[irun] ))  
} 


print(paste("Mean Total-N:", format(mean(tot), digits=0),sep="")    )
print(paste("Std Total-N:",  format(sd(tot), digits=0)  ,sep="")    )

print(paste("Mean Negative-N:", format(mean(neg),digits=0),sep="")  )
print(paste("Std Nagative-N:",  format(sd(neg),digits=0)  ,sep="")  )

print(paste("Mean Negative-N:", format(mean(neu),digits=0) ,sep="") )
print(paste("Std Nagative-N:",  format(sd(neu), digits=0)  ,sep="") )

print(paste("Mean Negative-N:", format(mean(pos),digits=0),sep="")  )
print(paste("Std Nagative-N:",  format(sd(pos),digits=0), sep="")   )


