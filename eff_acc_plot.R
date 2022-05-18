


par(pty='s')  # force the plot to be square before we start



# load Rda files
#table.dir <- c("/lfs/home/ychen/scripts/R/Rscripts/LAI_STUDY_EA/Rda_60_tc.occ/")
table.dir <- c("/lfs/home/ychen/scripts/R/Rscripts/LAI_STUDY_EA/spei02/old_rda/")





fun_table.qc <- function(wrk.table, qc1.set=0.6, qc2.set=0.2, area.set=0) {

      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.mws.aff)> 0)  )
      #QC1 mean difference  on LAI
      #QC2 measurement uncertainty on LAI
      #qc1.set <- 0.5
      #qc2.set <-  ES > delta(ES)
      fact <- 1.25     
      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set &  ( abs(wrk.table$eff.size.for)*fact  >= abs(wrk.table$qc2.score))  ))
      return(wrk.table)
      }



# load different runs 

wrk.rda <-c("EA_19990220to20181231_table.combw_8_p_0_2D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_10_p_0_3D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_12_p_0_4D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_0_p_60_2D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_0_p_80_3D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_0_p_100_4D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_8_p_60_2D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_10_p_80_3D_pos_060_off_000_cut_0_tc.occ.rda",
"EA_19990220to20181231_table.combw_12_p_100_4D_pos_060_off_000_cut_0_tc.occ.rda")

run.lty <- c(1,1,1,1,1,1,1,1,1)
run.pch <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
run.txt <- c("wind(1a)","wind(1b)","wind(1c)", "rainfall(2a)","rainfall(2b)","rainfall(2c)", "combine(3a)","combine(3b)","combine(3c)") 
#run.col <- c("black","black","black","blue","blue","blue", "green","green","green")
run.col <- c(rep("darkgray",9))
run.lwd <- c(rep(2.5,9))



pdf("FigS2.pdf",height=16, width=16)

par(mfrow=c(1,1), mar=c(5,5,3,1))
 
for (irun in 1:8) {
## load the data 
load(paste(table.dir,"/",wrk.rda[irun], sep=""))

#convert knot to m/s 
table.comb$eve.intensity <- as.numeric(table.comb$eve.intensity)*0.51
all.table <- table.comb


qcqa.table <- fun_table.qc( table.comb, qc1.set=0.5, qc2.set=0.2, area.set=0) 
pos.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)> 0.2)  )
neg.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)< -0.2)  )
neu.qcqa.table <- subset(qcqa.table, (abs(qcqa.table$eff.size.for)<= 0.2)  )


if (irun == 1 ) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     xlim=c(15,85), 
     xlab="Tropical cyclones maximum intensity (m/s)",
     ylab="Cumulative proportion",main="",
     col="black",lwd=4,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
#title("Negative effect size TC intensity distributions (1999-2018)", line = 2.5, cex.main=2.0)
}
#add category line

abline(v=64*0.51,col="black",lty="dashed") 
abline(v=83*0.51,col="black",lty="dashed")
abline(v=96*0.51,col="black",lty="dashed")
abline(v=113*0.51,col="black",lty="dashed")
abline(v=136*0.51, col="black",lty="dashed")
axis(3, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("C-1","C-2","C-3","C-4","C-5"), las=1, cex.axis=1.5)

plot(ecdf(neg.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="orange", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)

#if (irun==9) {
#plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
#     col="black",lwd=1,add=T)
#}


par(xpd=TRUE)
text(x=5,y=1.1, label="",cex=2, font=2)  
par(xpd=FALSE) 



} # end of loop 

#legend('topleft', 
#       legend=run.txt,bg="white",  # text in the legend
#       col=c("black"),
#       lty=run.lty,
#       pch=run.pch,
#       lwd=c(1,1,1,1,1,1,1,1,1),
#       cex=1.5)  # specify the point type to be a square


for (irun in 1:8) {
## load the data 
load(paste(table.dir,"/",wrk.rda[irun], sep=""))
#convert knot to m/s 
table.comb$eve.intensity <- as.numeric(table.comb$eve.intensity)*0.51
all.table <- table.comb


all.table <- table.comb
qcqa.table <- fun_table.qc( table.comb, qc1.set=0.5, qc2.set=0.2, area.set=0) 
pos.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)> 0.2)  )
neg.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)< -0.2)  )
neu.qcqa.table <- subset(qcqa.table, (abs(qcqa.table$eff.size.for)<= 0.2)  )


if (irun == 1 ) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     xlim=c(15,85), 
     xlab="Tropical cyclones maximum intensity (m/s)",
     ylab="Cumulative proportion",main="",
     col="black",lwd=4,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,add=T)
#title("Positive effect size TC intensity distributions (1999-2018)", line = 2.5,cex.main=2.0)
}
#add category line
abline(v=64*0.51,col="black",lty="dashed") 
abline(v=83*0.51,col="black",lty="dashed")
abline(v=96*0.51,col="black",lty="dashed")
abline(v=113*0.51,col="black",lty="dashed")
abline(v=136*0.51, col="black",lty="dashed")
axis(3, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("C-1","C-2","C-3","C-4","C-5"), las=1, cex.axis=1.5)


plot(ecdf(pos.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="forestgreen", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)

#if (irun==9) {
#plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
#     col="black",lwd=1,add=T)
#}


table.comb$eve.intensity <- as.numeric(table.comb$eve.intensity)*0.51
all.table <- table.comb


par(xpd=TRUE)
text(x=5,y=1.1, label="",cex=2, font=2)  
par(xpd=FALSE)
 
} # end of loop 

for (irun in 1:9) {
## load the data 
load(paste(table.dir,"/",wrk.rda[irun], sep=""))
#convert knot to m/s 
table.comb$eve.intensity <- as.numeric(table.comb$eve.intensity)*0.51
all.table <- table.comb


all.table <- table.comb
qcqa.table <- fun_table.qc( table.comb, qc1.set=0.5, qc2.set=0.2, area.set=0) 
pos.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)> 0.2)  )
neg.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)< -0.2)  )
neu.qcqa.table <- subset(qcqa.table, (abs(qcqa.table$eff.size.for)<= 0.2)  )


if (irun == 1 ) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     xlim=c(15,85), 
     xlab="Tropical cyclones maximum intensity (m/s)",
     ylab="Cumulative proportion",main="",
     col="black",lwd=4,cex.lab=2, cex.axis=2, cex.main=1.5, cex.sub=2,add=T)
#title("Neutral effect size TC intensity distributions (1999-2018)", line = 2.5,cex.main=2.0)
}
#add category line
abline(v=64*0.51,col="black",lty="dashed") 
abline(v=83*0.51,col="black",lty="dashed")
abline(v=96*0.51,col="black",lty="dashed")
abline(v=113*0.51,col="black",lty="dashed")
abline(v=136*0.51, col="black",lty="dashed")
axis(3, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("C-1","C-2","C-3","C-4","C-5"), las=1, cex.axis=1.5)


plot(ecdf(neu.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="deepskyblue1", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)


if (irun==9) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     col="black",lwd=5,add=T)
}


par(xpd=TRUE)
text(x=5,y=1.1, label="",cex=2.0, font=2)  
par(xpd=FALSE)



} # end of loop 




dev.off()

