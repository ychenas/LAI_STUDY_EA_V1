# multiple histogram plot 

x_del=0.1
#load the first case 
load("./Rda/WP_20010110to20181231_table.combw_3_p_200_122D.rda")
  #get density estimate
  dens=density(table.comb$eff.size.for)
  # Plot y-values scaled by number of observations against x values
  plot(x=dens$x-x_del,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       type="l", col="red",xlim=c(-1,1),ylim=c(0,45),
       xlab="TC effect size",ylab="Count estimate")

#load a case  
load("./Rda/WP_20010110to20181231_table.combw_3_p_200_123D.rda")
  #get density estimate
  dens=density(table.comb$eff.size.for)
  #plot the line 
  lines(x=dens$x-x_del,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       , col="red", lty=2)

#load a case  
load("./Rda/WP_20010110to20180110_table.combw_3_p_200_125D.rda")
  #get density estimate
  dens=density(table.comb$eff.size.for)
  #plot the line 
  lines(x=dens$x-x_del,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       , col="red", lty=3)

#load a case  
load("./Rda/WP_20010110to20180110_table.combw_3_p_200_127D.rda")
  #get density estimate
  dens=density(table.comb$eff.size.for)
  #plot the line 
  lines(x=dens$x-x_del,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       , col="red", lty=4)

#load a case  
load("./Rda/WP_20010110to20170620_table.combw_3_p_200_129D.rda")
  #get density estimate
  dens=density(table.comb$eff.size.for)
  #plot the line 
  lines(x=dens$x-x_del,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       , col="red", lty=5)

#add legend 


legend("topright", inset=.02,
       title="Potential TC Area:", 
       legend=c("2-D", "3-D", "5-D", "7-D", "9-D"),
       col=c("red"), 
       lty=1:5, cex=1.0)

