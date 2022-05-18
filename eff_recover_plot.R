# multiple histogram plot 

#sub_dir  <- c("./Rda_4050_Obs6/")
#run_date <- c("WP_20010110to20181231_table.comb")

#run_case <- c("w_3_p_350_TS_123D","w_3_p_350_TS_125D",
#              "w_3_p_350_TS_127D","w_3_p_350_TS_129D",
#              "w_2_p_350_TS_125D","w_3_p_350_TS_125D",
#              "w_5_p_350_TS_125D","w_7_p_350_TS_125D",
#              "w_5_p_40_TS_125D", "w_5_p_80_TS_125D",
#              "w_5_p_100_TS_125D","w_5_p_200_TS_125D",
#              "w_3_p_100_TS_125D" )
#legend.txt <- c("3m/s,p350-3D", "3m/s,p350-5D", "3m/s,p350-7D", "3m/s,p350-9D", 
#               "2m/s,p350-5D", "3m/s,p350-5D", "5m/s,p350-5D", "7m/s,p350-5D",
#                "5m/s,p40-5D", "5m/s,p80-5D","5m/s,p100-5D","5m/s,p200-5D","3m/s,p100-5D"),

sub_dir  <- c("./Rda_pos_off_2/")
run_date <- c("WP_20020110to20171231_table.comb")

run_case <- c( 
"w_3_p_100_TS_125Dpos_050_off_015",
"w_3_p_100_TS_125Dpos_070_off_015",
"w_3_p_100_TS_125Dpos_080_off_030",
"w_3_p_100_TS_125Dpos_100_off_030",
"w_3_p_100_TS_125Dpos_120_off_030",
"w_3_p_100_TS_125Dpos_150_off_030")

ncase = 6
nlag  = 1
# run_col <- c(rep("red",4),rep("cyan",4),rep("blue",4))
# run_lty <- c(seq(1,4,1),seq(1,4,1),seq(1,4,1))
 
# run_col <- c("red","cyan")
# run_lty <- c(1,2)


 all.table.comb <- data.frame()
 
  for (icase in 1:ncase) {
  #load the data 
  in_file <- paste(sub_dir,run_date,run_case[icase],".rda",sep="")
  print(in_file)
  load(in_file)
  
  #assign table to the dynamic text variable
  assign( value=table.comb, x=substr(run_case[icase], start=18, stop=24))  
  table.comb$eve.for.off.ass <- as.integer( substr(run_case[icase], start=22, stop=24))  
  all.table.comb <- rbind(table.comb, all.table.comb)
 
  } 
  # remove the small size events 
 
  #all.table.comb <- subset(all.table.comb, all.table.comb$eve.for.pix.aff >= 5000 )  

  # order by date and offset days 
  all.table.comb <- all.table.comb[order( all.table.comb$date.time,  all.table.comb$eve.id,all.table.comb$eve.for.off.aff),] 
 
 ylim=c(-2,2)
 xlim=c(0,200)

 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Days after the TC event (fast recovery)",ylab="Effect size (ES)")

# get TC event dati
  events <- as.character(unique(all.table.comb$date.time))
  line_col <- colorRampPalette(c("red","orange","yellow","green","blue","purple"))(length(events))
   for (ieve in 1:length(events)) {
#    for (ieve in 1:100) {
    eve.table <- subset(all.table.comb, as.character(all.table.comb$date.time) == events[ieve]) 
    # sub-event winthin 10days
      subevents <- as.character(unique(eve.table$eve.id))
      
      for (isubeve in 1: length(subevents) ) { 
          eve.table <- subset(eve.table, as.character(eve.table$eve.id) == subevents[isubeve])
          x_var <- eve.table$eve.for.off.aff
          y_var <- eve.table$eff.size.for
          if( (length(x_var)>=6) & (length(y_var)>=6) ) {   
             # add regression line
              fitline <- lm(as.vector(y_var[nlag:ncase]) ~ as.vector(x_var[nlag:ncase]), na.action=na.exclude )   
              if( abs(fitline$coefficients[2])< 0.025 & ( sign(y_var[1])!=sign(fitline$coefficients[2]))  ) { 
               points(x=x_var, y=y_var, col="black",bg=line_col[ieve],lwd=1.0,pch=21)
               abline(fitline,col=line_col[ieve],lty="dashed",add=T,lwd=1.5)
              } 
          }   
          #r_2   <- formatC(summary(fitline)$r.squared,width=4,digits=2,format="f") 
          #slope <- formatC(summary(fitline)$coefficients[2],width=4,digits=2,format="f") 
      }
   } 
 #add zero line 
 abline(h=0, col="lightgray",lty="dashed")
#open another window for plotting 
dev.new()

 ylim=c(-2,2)
 xlim=c(0,200)

 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Days after the TC event (slow recovery)",ylab="Effect size (ES)")

# get TC event dati
  events <- as.character(unique(all.table.comb$date.time))
  line_col <- colorRampPalette(c("red","orange","yellow","green","blue","purple"))(length(events))
   for (ieve in 1:length(events)) {
#    for (ieve in 1:100) {
    eve.table <- subset(all.table.comb, as.character(all.table.comb$date.time) == events[ieve]) 
    # sub-event winthin 10days
      subevents <- as.character(unique(eve.table$eve.id))
      
      for (isubeve in 1: length(subevents) ) { 
          eve.table <- subset(eve.table, as.character(eve.table$eve.id) == subevents[isubeve])
          x_var <- eve.table$eve.for.off.aff
          y_var <- eve.table$eff.size.for
          if( (length(x_var)>=6) & (length(y_var)>=6) ) {   
             # add regression line
              fitline <- lm(as.vector(y_var[nlag:ncase]) ~ as.vector(x_var[nlag:ncase]), na.action=na.exclude )   
              if( abs(fitline$coefficients[2])< 0.025 & ( sign(y_var[1])==sign(fitline$coefficients[2]))  ) { 
               points(x=x_var,y=y_var,col="black",bg=line_col[ieve],lwd=1.0,pch=21)
               abline(fitline,col=line_col[ieve],lty="dashed",add=T,lwd=1.5)
              } 
          }   
          #r_2   <- formatC(summary(fitline)$r.squared,width=4,digits=2,format="f") 
          #slope <- formatC(summary(fitline)$coefficients[2],width=4,digits=2,format="f") 
      }
   } 
 #add zero line 
 abline(h=0, col="lightgray",lty="dashed")

#open another window for plotting 
dev.new()

 xlim=c(-2,2)
 ylim=c(-0.025,0.025)

 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Effect size (ES)",ylab="Regression Slope (Beta)")

# get TC event dati
  events <- as.character(unique(all.table.comb$date.time))
  line_col <- colorRampPalette(c("red","orange","yellow","green","blue","purple"))(length(events))
   for (ieve in 1:length(events)) {
#    for (ieve in 1:100) {
    eve.table <- subset(all.table.comb, as.character(all.table.comb$date.time) == events[ieve]) 
    # sub-event winthin 10days
      subevents <- as.character(unique(eve.table$eve.id))
      
      for (isubeve in 1: length(subevents) ) { 
          eve.table <- subset(eve.table, as.character(eve.table$eve.id) == subevents[isubeve])
          x_var <- eve.table$eve.for.off.aff
          y_var <- eve.table$eff.size.for
          if( (length(x_var)>=6) & (length(y_var)>=6) ) {   
             # add regression line
              fitline <- lm(as.vector(y_var[nlag:ncase]) ~ as.vector(x_var[nlag:ncase]), na.action=na.exclude )   
          #    if( abs(fitline$coefficients[2])< 0.025 & ( sign(y_var[1])==sign(fitline$coefficients[2]))  ) { 
               points(x=y_var[1],y=fitline$coefficients[2],col="black",bg=line_col[ieve],lwd=1.0,pch=21)
          #     abline(fitline,col=line_col[ieve],lty="dashed",add=T,lwd=1.5)
          #    } 
          }   
          #r_2   <- formatC(summary(fitline)$r.squared,width=4,digits=2,format="f") 
          #slope <- formatC(summary(fitline)$coefficients[2],width=4,digits=2,format="f") 
      }
   } 
 #add zero line 
 abline(h=0, col="lightgray",lty="dashed")
 abline(v=0, col="lightgray",lty="dashed")
 
#open window for plotting
dev.new()


run_col <- colorRampPalette(c("red", "yellow"))(ncase)
 run_lty <- c(1)
 
 x_del=0.05
 #load the first case 
 ylim=c(0,50)
 xlim=c(-1.3,1.3)
   
 plot(x=0,  y=0,
      type="l", col="red",xlim=xlim,ylim=ylim,
      xlab="Effect size (ES)",ylab="Count estimate")
 #add zero line 
 abline(v=0, col="lightgray",lty="dashed")
 
  for (icase in 1:ncase) {
  #load data
  in_file <- paste(sub_dir,run_date,run_case[icase],".rda",sep="")
  print(in_file)
  load(in_file)
  #get density estimate
  dens=density(table.comb$eff.size.for)
  # Plot y-values scaled by number of observations against x values
  lines(x=dens$x-x_del,cex=1.5,
       y=length(table.comb$eff.size.for)*dens$y*dens$bw,
       type="l", col=run_col[icase], lty=run_lty, lwd=2.5)
  }

  #add legend 
#  grid(col="gray",lwd=1.5)
#  legend("topright", inset=.05, bg="white",
#       title="Search Area", 
 #      legend= legend.txt,
#       legend= run_case[1:ncase],
#       col=run_col, 
#       lty=run_lty,cex=1.0)

