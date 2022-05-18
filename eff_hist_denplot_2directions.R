# --------------------------------------------------
# multiple histogram plot and scattering plot 
# First Date: 2020-06-01
# Author: Yi-Ying Chen
# Email:  yiyingchen@gate.sinica.edu.tw
# --------------------------------------------------

go.pdf <- F
sub_dir  <- c("./")
# sub_dir  <- c("./")
  run_date <- c("WP_20010110to20171231_table.comb")
#  run_case <- c("w_4_p_150_S_125Dpos_050_off_000_cut_0")

#  run_case <- c("w_3_p_80_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_3_p_160_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_3_p_240_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_3_p_320_S_125Dpos_050_off_000_cut_0")
  
#  run_case <- c("w_3_p_80_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_5_p_80_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_7_p_80_S_125Dpos_050_off_000_cut_0")
#  run_case <- c("w_9_p_80_S_125Dpos_050_off_000_cut_0")

#  run_case <- c("w_6_p_100_S_122Dpos_050_off_000_cut_0")
#  run_case <- c("w_6_p_100_S_123Dpos_050_off_000_cut_0")
#  run_case <- c("w_6_p_100_S_125Dpos_050_off_000_cut_0")
  run_case <- c("w_6_p_100_S_127Dpos_050_off_000_cut_0") 
 
  plot_title <- paste(substr(run_case, start=1, stop=10),
                    substr(run_case, start=15,stop=16),sep="") 
  run_col <- c("red","blue")
  run_lty <- c(1,1)
  #5 load the data 
  in_file <- paste(sub_dir,run_date,run_case,".rda",sep="")
  print(in_file)
  load(in_file)

if (go.pdf) { pdf(file=paste("Scater_2D_Hist_ET_LAI_",plot_title,".pdf",sep=""),width=8, height=8) }

# 2 by 2 layout setting for combine plot
#par(oma = c(2, 0.5, 0.1, 0.1))
#layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1, 1), c(1, 3), TRUE)
#op = par(no.readonly = TRUE)


# subset datataset
 table.comb <- subset(table.comb, ((table.comb$eve.for.pix.aff)>3000 & as.integer(as.character(table.comb$date.time)) <= 20180101  )  )
 
 table.comb <- subset(table.comb, (abs(table.comb$qc1.score)<=0.25 &
                                   abs(table.comb$qc2.score)>=0.25 ))
 eff_cex <-  sqrt(table.comb$eve.for.pix.aff/3.1416)/100
 
 table.comb$eff.cex <- eff_cex


# 
 table.pos <- subset(table.comb, table.comb$eff.size.for > 0)
 table.neg <- subset(table.comb, table.comb$eff.size.for < 0)

table.pos <- table.pos[, c("eff.size.for","eve.intensity","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex")]
table.neg <- table.neg[, c("eff.size.for","eve.intensity","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex")]


go_pca <- TRUE

if (go_pca) {
# funtion for pair analysis
upper.panel <- function(x, y){
                points(x,y, cex=table.comb$eff.cex)
              }
# pair analysis
pairs(~eff.size.for+eve.for.mws.aff+eve.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+date.mon+eve.intensity+eve.lat.med.aff,upper.panel=upper.panel, 
lower.panel=NULL, data=table.comb, main="Pair Correlation")


# pca analysis 
pca <- prcomp(formula = ~ eve.for.mws.aff+eve.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+date.mon+eve.intensity+eve.lat.med.aff,
              data = table.comb,
              scale = TRUE)
#show pca 
pca
dev.new()

# pair analysis
pairs(~eff.size.for+eve.for.mws.aff+eve.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+date.mon+eve.intensity+eve.lat.med.aff,upper.panel=upper.panel,
lower.panel=NULL, data=table.comb, main="Pair Correlation")

#plot pca 
plot(pca,         
     type="line", # 
     main="Eigenvalue Plot of PCA") # 

#show variance =1.0 
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

vars <- (pca$sdev)^2  # 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars

# 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值
props <- vars / sum(vars)    
props

cumulative.props <- cumsum(props)  # 累加前n個元素的值
cumulative.props

#show top three PCs 
cumulative.props[3]

# 累積解釋比例圖
plot(cumulative.props)

# pca$rotation 
top3_pca.data <- pca$x[, 1:3]
top3_pca.data 

# 特徵向量(原變數的線性組合)
top3.pca.eigenvector <- pca$rotation[, 1:3]
top3.pca.eigenvector
pca$rotation

first.pca <- top3.pca.eigenvector[, 1]   #  第一主成份
second.pca <- top3.pca.eigenvector[, 2]  #  第二主成份
third.pca <- top3.pca.eigenvector[, 3]   #  第三主成份

dev.new()
# 第一主成份：由小到大排序原變數的係數
first.pca[order(first.pca, decreasing=FALSE)]  
##         RBI         H2B         H1B          HR         H3B          BB 
## -0.64629912 -0.51441491 -0.40991503 -0.34336124  0.01853759  0.05866272 
##          SB 
##  0.16722000
# 使用dotchart，繪製主成份負荷圖
dotchart(first.pca[order(first.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC1",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="red")         


dev.new()
# 第二主成份：由小到大排序原變數的係數
second.pca[order(second.pca, decreasing=FALSE)]  
##        H3B        H1B         SB        H2B        RBI         BB 
## -0.5595940 -0.4681242 -0.2741655 -0.2004156  0.1016251  0.2203673 
##         HR 
##  0.5417488
# 使用dotchart，繪製主成份負荷圖
dotchart(second.pca[order(second.pca, decreasing=FALSE)] ,  # 排序後的係數
         main="Loading Plot for PC2",                       # 主標題
         xlab="Variable Loadings",                          # x軸的標題
         col="blue")                                        # 顏色

dev.new()
# 第三主成份：由小到大排序原變數的係數
third.pca[order(third.pca, decreasing=FALSE)]  
##          BB          SB         RBI         H3B          HR         H2B 
## -0.78218985 -0.52853255 -0.25396353 -0.19427151 -0.03416307  0.01669591 
##         H1B 
##  0.07174689
# 使用dotchart，繪製主成份負荷圖
dotchart(third.pca[order(third.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC3",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="purple")                                     # 顏色



#assign the PCs to the table.comb


table.comb$PC1<-top3_pca.data[,1]
table.comb$PC2<-top3_pca.data[,2]
table.comb$PC3<-top3_pca.data[,3]

#do the pair analysis again

# pair analysis
pairs(~eff.size.for+PC1+PC2+PC3,upper.panel=upper.panel, 
lower.panel=NULL, data=table.comb, main="Pair Correlation")

}


# plot distribution plot of effect size with lai and latitude 

##### OPTION 1: hexbin from package 'hexbin' #######
#library(hexbin)

#x <- rnorm(mean=1.5, 5000)
#y <- rnorm(mean=1.6, 5000)
#df <- data.frame(x,y)
# Color housekeeping
#library(RColorBrewer)
#rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
#r <- rf(32)

# Create hexbin object and plot
#h <- hexbin(df)
#plot(h)
#plot(h, colramp=rf)




# select key variables 
#table.pos <- table.pos[, c("eve.intensity","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex")]
ld_go <-TRUE
if(ld_go) {
#dev.new()
es.lim =1.5 

#my.col    <- colorRampPalette(c("darkblue","blue","cyan","lightgray","yellow","pink","red","brown","brown"))(25)
eff.col    <- colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(16)
eff.breaks <- c( seq(-es.lim, es.lim, length.out = 15))
eff.breaks <- round(eff.breaks,2)

#assign the color based on effect size
table.comb$col <-  colorRampPalette(c("darkblue","blue","cyan","white","pink","red","brown"))(15)[as.numeric(cut(table.comb$eff.size.for ,breaks = eff.breaks))]

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=7, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(x=c(0,10), y=c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='')
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

par(oma = c(0.5, 0.5, 0.1, 0.1))
layout(matrix(c( 1, 2), 1, 2, byrow = TRUE), c(1, 1, 1), c(1, 3), TRUE)
op = par(no.readonly = TRUE)


# plot#1 scatter plot ES as funtion of LAI vs LATITUDE 
par(mar = c(4, 4, 1, 1))
  x_pos= table.comb$bef.for.lai.aff
  y_pos= table.comb$eve.lat.med.aff
  xy_bg= table.comb$col  
  ylim=c(10,48)
  xlim=c(2,6.5)
  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 21, bg = xy_bg, col = "black", lwd = 2.0, cex=eff_cex,
       ylab="",xlab=expression(paste(LAI[bef])), main="", )
  title(ylab=expression(paste(Latitude[aff])), line=2.0, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
  #axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
  #col.axis="black", las=2, cex.axis=1.0, tck=-.01)
  abline(v=0,h=0,lty="dashed")


#dev.new()

# plot#2 scatter plot ES as funtion of Soil Water and Months 
par(mar = c(4, 4, 1, 1))
  x_pos= table.comb$bef.for.asw.aff 
  y_pos= table.comb$eve.for.mws.aff
  xy_bg= table.comb$col  
  ylim=c(2, 12)
  xlim=c(0.25,0.40)
  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 21, bg = xy_bg, col = "black", lwd = 2.0, cex=eff_cex,
       ylab="",xlab=expression(paste(SoilWater[bef])), main="", )
  title(ylab=expression(paste(MaxWind[aff])), line=2.0, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
  #axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
  #col.axis="black", las=2, cex.axis=1.0, tck=-.01)
  abline(v=0,h=0,lty="dashed")







#dev.new()
# plot 3 legend on the top-right 
#par(mar = c(2, 5, 2, 2), xaxt = "n",cex=1.5)
#    color.bar(lut=eff.col,min=-es.lim,title="Effect size scale")
#

# plot#2 ES distribution
#par(mar = c(0, 4, 3, 0))
#  x_del=0.0
#  plot(x=0,  y=0,cex.main=1.0, cex.lab=1.0, cex.main=1.0, cex.sub=1.0,
#      type="l", col="red",xlim=xlim,ylim=c(0,30),
#      xlab="",ylab="Frequency",yaxt="n",xaxt="n")
  #get density estimation
#  table.comb.pos <- subset(table.comb, table.comb$eff.size.ace.for > 0.05 )  
#  table.comb.neg <- subset(table.comb, table.comb$eff.size.ace.for < -0.05) 
  #get density estimation
#  dens.pos = density(table.comb.pos$eff.size.for)
#  dens.neg = density(table.comb.neg$eff.size.for) 
#  dens.all = density(table.comb$eff.size.for) 


# positive ES_ET cases
#  lines(x=dens.pos$x-x_del,cex=1.5,
#        y=length(table.comb.pos$eff.size.for)*dens.pos$y*dens.pos$bw,
#        type="l", col="orange", lty="solid", lwd=4)
# negtive ES_ET cases
#  lines(x=dens.neg$x-x_del, cex=1.5,
#        y=length(table.comb.neg$eff.size.for)*dens.neg$y*dens.neg$bw,
#        type="l", col="forestgreen", lty="solid", lwd=4)
# all cases
#  dens.all = density(table.comb$eff.size.for)
#  lines(x=dens.all$x-x_del,cex=1.5,
#        y=length(table.comb$eff.size.for)*dens.all$y*dens.all$bw,
#        type="l", col="black", lwd=4)
#
#  abline(v=0, lty="dashed")
#  axis(2, at=seq(10,30,10),labels=seq(10,30,10), #pos=2.7,
#  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
# add zero line 
# abline(v=0, col="lightgray",lty="dashed")
#par(mar = c(4, 0, 0, 0))
  #load the data 
#  table.comb.pos <- subset(table.comb, table.comb$eff.size.for > 0.1 )  
#  table.comb.neg <- subset(table.comb, table.comb$eff.size.for < -0.1) 
  #get density estimate
#  dens.pos = density(table.comb.pos$eff.size.ace.for)
#  dens.neg = density(table.comb.neg$eff.size.ace.for) 

#plot(x=0, y=0, ylim=xlim, xlim=c(0,40), 
#     cex.main=1.0, cex.lab=1.0, cex.main=1.0, cex.sub=1.0, 
#     col="white",xlab="Frequency",  main="", yaxt = "n" )
#  abline(h=0, lty="dashed") 
# positive ES_LAI
#  lines(y=dens.pos$x-x_del,cex=1.5,
#        x=length(table.comb.pos$eff.size.for)*dens.pos$y*dens.pos$bw,
#        type="l", col=run_col[1], lty="solid", lwd=4)
## negtive ES_LAI
#  lines(y=dens.neg$x-x_del, cex=1.5,
#        x=length(table.comb.neg$eff.size.for)*dens.neg$y*dens.neg$bw,
#        type="l", col=run_col[2], lty="solid", lwd=4)  
# plot #4 
#par(mar = c(0, 0, 0, 0))
#par(cex.axis=1, cex.lab=2, cex.main=2.5, cex.sub=2)

#   plot(type="n",axes=FALSE,x=0)
#   legend(title=plot_title,"bottomleft",
#          legend=c(expression(paste(ES[LAI]," for ",ES[ET]," > 0")),
#            expression(paste(ES[LAI]," for ",ES[ET]," < 0")),
#            expression(paste(ES[ET]," for ",ES[LAI]," > 0")),
#            expression(paste(ES[ET]," for ",ES[LAI]," < 0"))),
#          col=c("orange","forestgreen","red", "blue"), 
#          lty=1,lwd=3, box.lty=0,cex=0.9) 
#
# open new device

#dev.new()

# 2 by 2 layout setting for combine plot
#par(oma = c(2, 0.5, 0.1, 0.1))
#layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1, 1), c(1, 3), TRUE)
#op = par(no.readonly = TRUE)


# subset datataset
# table.comb <- subset(table.comb, table.comb$eve.for.pix.aff >= 5000)
#  eff_cex <-  sqrt(table.comb$eve.for.pix.aff/3.1416)/100

# plot#1 scatter plot ES_ET vs ES_Rainf 
#par(mar = c(4, 4, 0, 0))
#  x_pos= table.comb$eff.size.acf.for
#  y_pos= table.comb$eff.size.ace.for 
#  ylim=c(-5,5)
#  xlim=c(-5,5)
#  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
#       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
#       pch = 21, bg = "lightgray", col = ifelse(  (table.comb$eff.size.for > 0.1)&(table.comb$eff.size.ace.for > 0.1)    ,"orange","black"),
#       lwd = 2.0, cex=eff_cex, ylab="",xlab=expression(paste(ES[Rainf])), main="", yaxt = "n" )
#  title(ylab=expression(paste(ES[ET])), line=2.5, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
#  axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
#  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
#  abline(v=0,h=0,lty="dashed")

#dev.new()
# plot#2 scatter plot ES_LAI vs ES_Rainf 
#par(mar = c(4, 4, 0, 0))
#  x_pos= table.comb$eff.size.acf.for
#  y_pos= table.comb$eff.size.for 
#  ylim=c(-5,5)
#  xlim=c(-5,5)
#  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
#       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
#       pch = 21, bg = "lightgray", col = ifelse(  (table.comb$eff.size.for > 0.1)&(table.comb$eff.size.ace.for > 0.1)    ,"orange","black"),
#       lwd = 2.0, cex=eff_cex, ylab="",xlab=expression(paste(ES[Rainf])), main="", yaxt = "n" )
#  title(ylab=expression(paste(ES[LAI])), line=2.5, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
#  axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
#  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
#  abline(v=0,h=0,lty="dashed")


#if (go.pdf) { dev.off() ; print("please check pdf file!" ) }

}

# select key variables 
#table.pos <- table.pos[, c("eve.intensity","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex")]
ld_go <-FALSE
if(ld_go) {
dev.new()

# plot#1 scatter plot ES_LAI vs ES_ET 
par(mar = c(4, 4, 0, 0))
  x_pos= table.comb$eff.size.for
  y_pos= table.comb$eff.size.ace.for 
  ylim=c(-5,3)
  xlim=c(-5,3)
  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 21, bg = "lightgray", col = "black", lwd = 2.0, cex=eff_cex,
       ylab="",xlab=expression(paste(ES[LAI])), main="", yaxt = "n" )
  title(ylab=expression(paste(ES[ET])), line=2.5, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
  axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
  abline(v=0,h=0,lty="dashed")

# plot#2 ES distribution
par(mar = c(0, 4, 3, 0))
  x_del=0.0
  plot(x=0,  y=0,cex.main=1.0, cex.lab=1.0, cex.main=1.0, cex.sub=1.0,
      type="l", col="red",xlim=xlim,ylim=c(0,30),
      xlab="",ylab="Frequency",yaxt="n",xaxt="n")
  #get density estimation
  table.comb.pos <- subset(table.comb, table.comb$eff.size.ace.for > 0.05 )  
  table.comb.neg <- subset(table.comb, table.comb$eff.size.ace.for < -0.05) 
  #get density estimation
  dens.pos = density(table.comb.pos$eff.size.for)
  dens.neg = density(table.comb.neg$eff.size.for) 
  dens.all = density(table.comb$eff.size.for) 


# positive ES_ET cases
  lines(x=dens.pos$x-x_del,cex=1.5,
        y=length(table.comb.pos$eff.size.for)*dens.pos$y*dens.pos$bw,
        type="l", col="orange", lty="solid", lwd=4)
# negtive ES_ET cases
  lines(x=dens.neg$x-x_del, cex=1.5,
        y=length(table.comb.neg$eff.size.for)*dens.neg$y*dens.neg$bw,
        type="l", col="forestgreen", lty="solid", lwd=4)
# all cases
  dens.all = density(table.comb$eff.size.for)
  lines(x=dens.all$x-x_del,cex=1.5,
        y=length(table.comb$eff.size.for)*dens.all$y*dens.all$bw,
        type="l", col="black", lwd=4)

  abline(v=0, lty="dashed")
  axis(2, at=seq(10,30,10),labels=seq(10,30,10), #pos=2.7,
  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
# add zero line 
# abline(v=0, col="lightgray",lty="dashed")
par(mar = c(4, 0, 0, 0))
  #load the data 
  table.comb.pos <- subset(table.comb, table.comb$eff.size.for > 0.1 )  
  table.comb.neg <- subset(table.comb, table.comb$eff.size.for < -0.1) 
  #get density estimate
  dens.pos = density(table.comb.pos$eff.size.ace.for)
  dens.neg = density(table.comb.neg$eff.size.ace.for) 

plot(x=0, y=0, ylim=xlim, xlim=c(0,40), 
     cex.main=1.0, cex.lab=1.0, cex.main=1.0, cex.sub=1.0, 
     col="white",xlab="Frequency",  main="", yaxt = "n" )
  abline(h=0, lty="dashed") 
# positive ES_LAI
  lines(y=dens.pos$x-x_del,cex=1.5,
        x=length(table.comb.pos$eff.size.for)*dens.pos$y*dens.pos$bw,
        type="l", col=run_col[1], lty="solid", lwd=4)
# negtive ES_LAI
  lines(y=dens.neg$x-x_del, cex=1.5,
        x=length(table.comb.neg$eff.size.for)*dens.neg$y*dens.neg$bw,
        type="l", col=run_col[2], lty="solid", lwd=4)  
# plot #4 
par(mar = c(0, 0, 0, 0))
#par(cex.axis=1, cex.lab=2, cex.main=2.5, cex.sub=2)

   plot(type="n",axes=FALSE,x=0)
   legend(title=plot_title,"bottomleft",
          legend=c(expression(paste(ES[LAI]," for ",ES[ET]," > 0")),
            expression(paste(ES[LAI]," for ",ES[ET]," < 0")),
            expression(paste(ES[ET]," for ",ES[LAI]," > 0")),
            expression(paste(ES[ET]," for ",ES[LAI]," < 0"))),
          col=c("orange","forestgreen","red", "blue"), 
          lty=1,lwd=3, box.lty=0,cex=0.9) 

# open new device

dev.new()

# 2 by 2 layout setting for combine plot
#par(oma = c(2, 0.5, 0.1, 0.1))
#layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 1, 1), c(1, 3), TRUE)
#op = par(no.readonly = TRUE)


# subset datataset
# table.comb <- subset(table.comb, table.comb$eve.for.pix.aff >= 5000)
  eff_cex <-  sqrt(table.comb$eve.for.pix.aff/3.1416)/100

# plot#1 scatter plot ES_ET vs ES_Rainf 
#par(mar = c(4, 4, 0, 0))
  x_pos= table.comb$eff.size.acf.for
  y_pos= table.comb$eff.size.ace.for 
  ylim=c(-5,5)
  xlim=c(-5,5)
  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 21, bg = "lightgray", col = ifelse(  (table.comb$eff.size.for > 0.1)&(table.comb$eff.size.ace.for > 0.1)    ,"orange","black"),
       lwd = 2.0, cex=eff_cex, ylab="",xlab=expression(paste(ES[Rainf])), main="", yaxt = "n" )
  title(ylab=expression(paste(ES[ET])), line=2.5, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
  axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
  abline(v=0,h=0,lty="dashed")

dev.new()
# plot#2 scatter plot ES_LAI vs ES_Rainf 
#par(mar = c(4, 4, 0, 0))
  x_pos= table.comb$eff.size.acf.for
  y_pos= table.comb$eff.size.for 
  ylim=c(-5,5)
  xlim=c(-5,5)
  plot(x=x_pos, y=y_pos , xlim=xlim,ylim=ylim,
       cex.main=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5,
       pch = 21, bg = "lightgray", col = ifelse(  (table.comb$eff.size.for > 0.1)&(table.comb$eff.size.ace.for > 0.1)    ,"orange","black"),
       lwd = 2.0, cex=eff_cex, ylab="",xlab=expression(paste(ES[Rainf])), main="", yaxt = "n" )
  title(ylab=expression(paste(ES[LAI])), line=2.5, cex.lab=1.5) #, family="Calibri Light")
  #y axis on the right
  axis(2, at=seq(-10,4),labels=seq(-10,4), #pos=2.7,
  col.axis="black", las=2, cex.axis=1.0, tck=-.01)
  abline(v=0,h=0,lty="dashed")


if (go.pdf) { dev.off() ; print("please check pdf file!" ) }

}
