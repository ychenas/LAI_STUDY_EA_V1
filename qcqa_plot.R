

# load table.comb table
load("EA_19990220to20181231_table.combw_8_p_60_2D_pos_050_off_000_cut_0.rda")
run1 <- table.comb

load("EA_19990220to20181231_table.combw_0_p_60_2D_pos_050_off_000_cut_0.rda")
run2 <- table.comb

load("EA_19990220to20181231_table.combw_8_p_0_2D_pos_050_off_000_cut_0.rda")
run3 <- table.comb

#load("EA_19990220to20181231_table.combw_0_p_80_3D_pos_050_off_000_cut_0.rda")
#run4 <- table.comb

#load("EA_19990220to20181231_table.combw_10_p_0_3D_pos_050_off_000_cut_0.rda")
#run5 <- table.comb

load("EA_19990220to20181231_table.combw_10_p_80_3D_pos_050_off_000_cut_0.rda")
run6 <- table.comb

#load("EA_19990220to20181231_table.combw_0_p_100_4D_pos_050_off_000_cut_0.rda")
#run7 <- table.comb

#load("EA_19990220to20181231_table.combw_12_p_0_4D_pos_050_off_000_cut_0.rda")
#run8 <- table.comb

#load("EA_19990220to20181231_table.combw_12_p_100_4D_pos_050_off_000_cut_0.rda")
#run9 <- table.comb


# combine table 
run_2D <- rbind( run1, run2, run3, run6)

#insert a qulity level color code

#      qc1 : T             F
#  qc2 
#   T        TT(green)     FT(yellow)
#   F        TF(yellow)    FF(green)
run_2D[, "qcqa.color"] <- "green"

n.total <- length(run_2D)

#for ( i in 1:n.total) {
#
#   if ( run_2D$qc1.score[i] <= 0.25 ) {
#        run_2D$qcqa.color[i] <- "yellow"
#   } 

#   if ( abs(run_2D$qc2.score[i]) >= abs(run_2D$eff.size.for[i]) ) {
#        run_2D$qcqa.color[i] <- "yellow"
#   } 

#}   



plot(abs(run_2D$eff.size.for) ~ abs(run_2D$qc2.score), xlim=c(0,4.5), ylim=c(0,4.5), pch=21, col="black", cex=1.0,bg=run_2D$qcqa.color)
abline(a=0,b=1, col="blue")
abline(a=0,b=0.5, col="blue", lty="dashed")
dev.new()
 
plot(abs(run_2D$eff.size.for) ~ abs(run_2D$qc1.score), xlim=c(0,4.5), ylim=c(0,4.5), pch=21, col="black", cex=1.0,bg=run_2D$qcqa.color)
abline(h=0,v=1,col="blue") 
abline(h=0,v=.5,col="blue",lty="dashed")



