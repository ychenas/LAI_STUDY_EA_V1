


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

c_pecent <- c(23, 15, 10, 12, 18, 22 )
pdf("FigS2_3x3.pdf",height=16, width=16)

par(mfrow=c(3,3), mar=c(5,5,3,1))
for (irun in 1:9) {
## load the data 
load(paste(table.dir,"/",wrk.rda[irun], sep=""))

#convert knot to m/s 
table.comb$eve.intensity <- as.numeric(table.comb$eve.intensity)*0.51
all.table <- table.comb


qcqa.table <- fun_table.qc( table.comb, qc1.set=0.5, qc2.set=0.2, area.set=0) 
pos.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)> 0.2)  )
neg.qcqa.table <- subset(qcqa.table, ((qcqa.table$eff.size.for)< -0.2)  )
neu.qcqa.table <- subset(qcqa.table, (abs(qcqa.table$eff.size.for)<= 0.2)  )


#f ( (irun == 1)|(irun==4)|(irun==7) ) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     xlim=c(15,85), 
     xlab="Tropical cyclones maximum intensity (m/s)",
     ylab="Cumulative proportion",main="",
     col="white",lwd=0,cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
#title("Negative effect size TC intensity distributions (1999-2018)", line = 2.5, cex.main=2.0)
#}
#add legend 
leg.title<- c("2D runs(wind only)", "3D runs(wind only)", "4D runs(wind only)",
              "2D runs(rain only)", "3D runs(rain only)", "4D runs(rain only)",
              "2D runs(wind/rain)", "3D runs(wind/rain)", "4D runs(wind/rain)")
#if (irun <=if (irun <=if (irun <=3 ) {
     legend(title=leg.title[irun],60, 0.3, c("Positive", "Neutral", "Negative"), col = c("forestgreen","cyan","orange"),
            text.col = "black", lty = c(1, 1, 1),  merge = TRUE, bg = "gray90",cex=1.5,lwd=c(2,2,2))
#}


xx <- c(64*0.51,83*0.51,96*0.51,113*0.51,136*0.51)
yy <- c(35/100,47/100,55/100,66/100,87/100)

off.y <- c(0.03)
abline(v=64*0.51,col="gray",lty="dashed") 
abline(v=83*0.51,col="gray",lty="dashed")
abline(v=96*0.51,col="gray",lty="dashed")
abline(v=113*0.51,col="gray",lty="dashed")
abline(v=136*0.51, col="gray",lty="dashed")

abline(h=35/100,col="gray",lty="dashed") 
abline(h=47/100,col="gray",lty="dashed")
abline(h=55/100,col="gray",lty="dashed")
abline(h=66/100,col="gray",lty="dashed")
abline(h=87/100, col="gray",lty="dashed")
text(x=18,y=yy[1]-off.y,"35%",cex=2)
text(x=18,y=yy[2]-off.y,"47%",cex=2)
text(x=18,y=yy[3]-off.y,"55%",cex=2)
text(x=18,y=yy[4]-off.y,"66%",cex=2)
text(x=18,y=yy[5]-off.y,"87%",cex=2)






axis(3, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("C1","C2","C3","C4","C5"), las=1, cex.axis=1.5)

plot(ecdf(neg.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="orange", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)
plot(ecdf(neu.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="cyan", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)
plot(ecdf(pos.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="forestgreen", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)



#if (irun==9) {
plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     col="black",lwd=2,add=T)
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



dev.off()

