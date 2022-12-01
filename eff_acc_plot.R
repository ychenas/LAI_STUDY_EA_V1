


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
      fact <- 1.2   
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
run.txt <- c("wind(1a)","wind(1b)","wind(1c)", "rainfall(2a)","rainfall(2b)","rainfall(2c)", "combined(3a)","combined(3b)","combined(3c)") 
#run.col <- c("black","black","black","blue","blue","blue", "green","green","green")
run.col <- c(rep("darkgray",9))
run.lwd <- c(rep(2.5,9))

c_pecent <- c(23, 15, 10, 12, 18, 22 )
pdf("FigS2_3x3_cyclone_distribution.pdf",height=16, width=16)

par(mfrow=c(3,3), mar=c(5,5,3,1))
for (irun in 1:9) {

print(paste("current run:",run.txt[irun]))
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
     xlab="Tropical cyclones maximum wind speeds (m/s)",
     ylab="",main="",
     col="white",lwd=0,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
     yaxt="n",xaxt="n")
#title("Negative effect size TC intensity distributions (1999-2018)", line = 2.5, cex.main=2.0)
#}
#add legend 
leg.title<- c("2D (wind only)", "3D (wind only)", "4D (wind only)",
              "2D (rain only)", "3D (rain only)", "4D (rain only)",
              "2D (combined)", "3D (combined)", "4D runs(combined)")
#if (irun <=if (irun <=if (irun <=3 ) {
#}

#add climatology for C1-C5
axis(1, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("32","42","49","58","69"), las=1, cex.axis=1.5)

axis(2, at=c(0.31,0.45,0.55,0.66,0.87),
 labels=c("31%","45%","55%","66%","87%"), las=1, cex.axis=1.5)


#add C1-C5 text
axis(3, at=c(64*.51,83*.51,96*.51,113*.51,136*.51),
 labels=c("C1","C2","C3","C4","C5"), las=1, cex.axis=1.5)

plot(ecdf(neg.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="orange", lwd=2, lty=run.lty[irun],pch=run.pch[irun],add=T)
plot(ecdf(neu.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="cyan", lwd=2, lty=run.lty[irun],pch=run.pch[irun],add=T)
plot(ecdf(pos.qcqa.table[,"eve.intensity"]),do.points=T,verticals=TRUE,cex=2.0,
     col="forestgreen", lwd=4, lty=run.lty[irun],pch=run.pch[irun],add=T)



# add accumulative
plot( ecdf(qcqa.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     col="black",lwd=2,add=T)

plot( ecdf(all.table[,"eve.intensity"]), do.points=F,verticals=TRUE,
     col="gray",lwd=2,add=T)

my_ecdf <- ecdf(sort(qcqa.table[,"eve.intensity"])) 
tmp.table <- data.frame(x=sort(qcqa.table[,"eve.intensity"]), pdf=my_ecdf(sort(qcqa.table[,"eve.intensity"])))

#remove duplicated values
tmp.table[!duplicated(tmp.table),]


xx <- c(64*0.51,83*0.51,96*0.51,113*0.51,136*0.51)
#yy for each cdf distribution
#yy <- c(35/100,47/100,55/100,66/100,87/100)



#find C1-C5 acc-cdf
# x1---  x0---   x2---
# y1---  y0---   y2---

for ( i in 1:5) {
x0 <- xx[i]
tmp.1=0

# find x1, y1 for just below ix 
 for (j in 1: length(tmp.table$x)) {
   if ( ( tmp.1== 0) & (tmp.table$x[j] > x0) ) {
      tmp.1=1
      y2=tmp.table$pdf[j]
      x2=tmp.table$x[j]
      y1=tmp.table$pdf[j-1]
      x1=tmp.table$x[j-1]
   
    y0=((x0-x1)*(y2-y1))/(x2-x1)  + y1
   } 
 }

print( paste("x1:",formatC(x1,digits=2,format="f"), "x0:",x0,"x2:",formatC(x2,digits=2,format="f"),sep=","))
print( paste("y1:",formatC(y1,digits=2,format="f"), "y0:",formatC(y0,digits=2,format="f"),"y2:",formatC(y2,digits=2,format="f"),sep=","))



#add horizonal lines for x0
abline(h=y0,col="gray",lty="dashed") 
#add vertical line for y0
off.y <- c(0.03)
abline(v=x0, col="gray",lty="dashed")
text(x=18,y=y0+off.y, paste(formatC(y0*100,width=2,format="d"),"%",sep="") ,cex=1.5)



} 

#print(tmp.table)
#
#}

legend(title=leg.title[irun],62, 0.35, c("Positive", "Neutral", "Negative","All events","Climatology"), col = c("forestgreen","cyan","orange","black","gray"),
            text.col = "black", lty = c(1, 1, 1,1,1),  merge = TRUE, bg = "white",cex=1.5,lwd=c(1,1,1,1,1))

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

