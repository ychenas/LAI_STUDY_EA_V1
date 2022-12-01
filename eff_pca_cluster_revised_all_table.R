## FUNCTIONS

#################################################################
#**************************
#return the rules of a tree
#**************************
getConds<-function(tree){
  #store all conditions into a list
  conds<-list()
  #start by the terminal nodes and find previous conditions
  id.leafs<-which(tree$status==-1)
  j<-0
  for(i in id.leafs){
    j<-j+1
    prevConds<-prevCond(tree,i)
    conds[[j]]<-prevConds$cond
    while(prevConds$id>1){
      prevConds<-prevCond(tree,prevConds$id)
      conds[[j]]<-paste(conds[[j]]," & ",prevConds$cond)
      if(prevConds$id==1){
        conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
        break()
      }
    }
    
  }
  
  return(conds)
}

#**************************
#find the previous conditions in the tree
#**************************
prevCond<-function(tree,i){
  if(i %in% tree$right_daughter){
    id<-which(tree$right_daughter==i)
    cond<-paste(tree$split_var[id],">",tree$split_point[id])
  }
  if(i %in% tree$left_daughter){
    id<-which(tree$left_daughter==i)
    cond<-paste(tree$split_var[id],"<",tree$split_point[id])
  }
  
  return(list(cond=cond,id=id))
}

#remove spaces in a word
collapse<-function(x){
  x<-sub(" ","_",x)
  
  return(x)
}


#####################################################################################
to.dendrogram <- function(dfrep,rownum=1,height.increment=0.1){
  
  if(dfrep[rownum,'status'] == -1){
    rval <- list()
    
    attr(rval,"members") <- 1
    attr(rval,"height") <- 0.0
    attr(rval,"label") <- dfrep[rownum,'prediction']
    attr(rval,"leaf") <- TRUE
    
  }else{##note the change "to.dendrogram" and not "to.dendogram"
    left <- to.dendrogram(dfrep,dfrep[rownum,'left daughter'],height.increment)
    right <- to.dendrogram(dfrep,dfrep[rownum,'right daughter'],height.increment)
    rval <- list(left,right)
    
    attr(rval,"members") <- attr(left,"members") + attr(right,"members")
    attr(rval,"height") <- max(attr(left,"height"),attr(right,"height")) + height.increment
    attr(rval,"leaf") <- FALSE
    attr(rval,"edgetext") <- paste(dfrep[rownum,'split var'],"\n<",round(dfrep[rownum,'split point'], digits = 2),"=>", sep = " ")
  }
  
  class(rval) <- "dendrogram"
  
  return(rval)
}



#####################################
# Quality Check funtion 
#####################################    
#QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
fun_table.qc <- function(wrk.table, qc1.set, qc2.set, area.set=0) {
 
      #print( paste("Total events:", length(wrk.table$eve.for.pix.aff),sep=""))
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      all.eve <- length(wrk.table$eve.for.pix.aff)
      print( paste("Total events:", all.eve,sep=""))
 
      #QC1 mean difference on LAI #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.1

      qc1.score <-  abs((wrk.table$bef.for.lai.aff/wrk.table$bef.for.lai.ref)-1) 
      print(formatC(qc1.score,digits=2,format="g")) 
      
      flag1.table <- subset(wrk.table, qc1.score   <= qc1.set )


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
qc1.set <- 0.5
qc2.set <- 0.2
area.set <- 0

#rda.dir <- c("./Rda_60_tc.occ/")

#rda.dir <- c("./spei02/")
rda.dir <- c("./spei02/BG_R1/")
#


#get file-list 
in.files <- list.files(path=rda.dir, pattern="EA_*")
#in.files <- list.files(path=rda.dir, pattern="*_[2]D") 
#in.files <- list.files(path=rda.dir, pattern="*_[4]D") 
   
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
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for >  0.10] <- "Positive"
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for < -0.10] <- "Negative"
table.tmp[['qcqa']]$group[abs(table.tmp[['qcqa']]$eff.size.for) <= 0.10] <- "Neutral"

table.tmp[['flag2']]$group <- "Neutral"

table.all <- rbind(table.tmp[['qcqa']],table.tmp[['flag2']])  
#table.all <- table.tmp[['qcqa']]


table.all$group <- as.factor(table.all$group)


library("dplyr")
library("randomForest")
#library("caret")
#library("e1071")
library("ggplot2")
library("ggpubr")

library("gridExtra")

# read the enso index table 
enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)
pj.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ERA5_DATA/PJ_INDEX/PJ_index.csv",sep=",",header=T)  

# assign the positve/negtive value of eff_size as group B/A
#table.all$group <- NA
#table.all$group[table.all$eff.size.for >= 0.5 ] <- "Growth"
#table.all$group[abs(table.all$eff.size.for) < 0.5] <- "Neutral"
#table.all$group[table.all$eff.size.for <= 0.5] <- "Damage"
#table.all$group <- as.factor(table.all$group)

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


mat1 <-matrix(c(1,1,2,3),2,2)
layout(mat1, widths=c(2,2),height=c(2,2))



hist.breaks <- c(seq(-5,5.0,length=100))
hist.col <- ifelse( hist.breaks < 0, "orange2", "forestgreen") 
 

#hist(table.all$eff.size.for, col=hist.col
#     ,breaks=hist.breaks, main="All events", xlim=c(-1,1) )
#hist(table.all$eff.size.for[ table.all$oni.enso >= 0.5 ], col= hist.col,
#     ,breaks=hist.breaks, main="El Nino (positive phase events)",xlim=c(-1,1))
#hist(table.all$eff.size.for[ table.all$oni.enso <= -0.5 ], col= hist.col,
#     breaks=hist.breaks, main="La Nina (negtive phase events)",xlim=c(-1,1))



eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex

table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))


# define the change of droughts imdex  
table.all$del.spei.for <- as.numeric(table.all$eve.for.spei.aff - table.all$eve.for.spei.ref)


# prepare table for pca analysis
table.pca <- table.all[, c("eff.size.for",
                           "eve.for.spei.aff","date.mon.day","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","eve.lat.med.aff","eff.cex","pj.index")]

go_pca <- FALSE

if (go_pca) {
# funtion for pair analysis
upper.panel <- function(x, y){
                points(x,y, cex=table.comb$eff.cex)
              }
# pair analysis
#pairs(   ~eff.size.for+eve.for.mws.aff+eve.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+eve.lat.med.aff+enso,upper.panel=upper.panel,
#lower.panel=NULL, data=table.all, main="Pair Correlation")


# pca analysis
pca <- prcomp(formula = ~ eve.for.spei.aff+date.mon.day+eve.for.mws.aff+eve.for.acf.aff+bef.for.acf.aff+bef.for.lai.aff+eve.lat.med.aff+eff.cex+pj.index,
              data = table.pca,
              scale = TRUE)
#show pca
pca
dev.new()
layout(matrix(c(1,2), 1, 2,byrow = TRUE))
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
cumulative.props[6]

# 累積解釋比例圖
plot(cumulative.props)

# pca$rotation
pca.data <- pca$x[, 1:6]
pca.data

# 特徵向量(原變數的線性組合)
pca.eigenvector <- pca$rotation[, 1:6]
pca.eigenvector
pca$rotation

first.pca  <- pca.eigenvector[, 1]   #  
second.pca <- pca.eigenvector[, 2]   #  
third.pca  <- pca.eigenvector[, 3]   #  
fourth.pca <- pca.eigenvector[, 4]   # 
fiveth.pca <- pca.eigenvector[, 5]   #  
sixth.pca  <- pca.eigenvector[, 6]   #  


dev.new()
layout(matrix(c(1,2,3,4,5,6), 3, 2,byrow = TRUE))
#par(mfrow=c(2,3))
#
# 第一主成份：由小到大排序原變數的係數
first.pca[order(first.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(first.pca[order(first.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC1",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="red")
abline(v=0, col="gray",lty="dashed")

# 第二主成份：由小到大排序原變數的係數
second.pca[order(second.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(second.pca[order(second.pca, decreasing=FALSE)] ,  # 排序後的係數
         main="Loading Plot for PC2",                       # 主標題
         xlab="Variable Loadings",                          # x軸的標題
         col="blue")                                        # 顏色
abline(v=0, col="gray",lty="dashed")

# 第三主成份：由小到大排序原變數的係數
third.pca[order(third.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(third.pca[order(third.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC3",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="purple")                                     # 顏色
abline(v=0, col="gray",lty="dashed")

# 第四主成份：由小到大排序原變數的係數
fourth.pca[order(fourth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(fourth.pca[order(fourth.pca, decreasing=FALSE)] , # 排序後的係數
         main="Loading Plot for PC4",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="forestgreen")                                # 顏色
abline(v=0, col="gray",lty="dashed")

# 第五主成份：由小到大排序原變數的係數
fiveth.pca[order(fiveth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(fiveth.pca[order(fiveth.pca, decreasing=FALSE)] , # 排序後的係數
         main="Loading Plot for PC5",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="orange")                                     # 顏色
abline(v=0, col="gray",lty="dashed")

# 第六主成份：由小到大排序原變數的係數
sixth.pca[order(sixth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(sixth.pca[order(fiveth.pca, decreasing=FALSE)] ,  # 排序後的係數
         main="Loading Plot for PC6",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="brown")                                      # 顏色
abline(v=0, col="gray",lty="dashed")

#assign the PCs to the table.comb
table.all$pc1<-pca.data[,1]
table.all$pc2<-pca.data[,2]
table.all$pc3<-pca.data[,3]
table.all$pc4<-pca.data[,4]
table.all$pc5<-pca.data[,5]
table.all$pc6<-pca.data[,6]

#do the pair analysis again
dev.new()
# pair analysis
pairs(~eff.size.for+pc1+pc2+pc3+pc4+pc5+pc6,upper.panel=upper.panel,
,lower.panel=NULL, data=table.all, main="Pair Correlation")

}

# go cluster analysis
go.cluster <- FALSE
if (go.cluster) {
library(cluster)
data <- table.all[,c("pc1","pc2","pc3","pc4","pc5","pc6")] 
# pam = Partitioning Around Medoids
kmeans.cluster <- kmeans(data, centers=3)
#pam.cluster <- pam(data, k=3)
# comparison 
print(table(kmeans.cluster$cluster, table.all$group) )
}


# go for factor analysis
go.factor <- TRUE

if (go.factor) {

library(psych)


# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors,
# with varimax rotation

#prepare table for factor analysis
#table.fta <- table.all[, c("eff.size.for",
#                         "date.mon.day","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","eve.lat.med.aff","eff.cex","oni.enso")]

# group the data for Tropical Cyclone Characteristics
table.tcc <- table.all[, c("eve.for.mws.aff","eve.for.acf.aff","eve.lat.med.aff","eff.cex","eve.intensity")]


# group the data for Pre-event Characteristics
table.pc <- table.all[, c("eve.for.spei.ref","date.mon","bef.for.acf.aff","bef.for.lai.aff","pj.index")]

#SELECTED potentail variables from tcc & pc groups
table.pot <- table.all[, c("eve.for.mws.aff","eve.for.acf.aff","eve.lat.med.aff","eve.intensity","eff.cex", 
                           "date.mon","bef.for.acf.aff","bef.for.lai.aff","pj.index","eve.for.spei.ref","del.spei.for","eff.size.for")]


fta.pc <- principal(table.pc, nfactors=3, rotate="varimax")
print(fta.pc)

fta.tcc <- principal(table.tcc, nfactors=3, rotate="varimax")
print(fta.tcc)
#fta.tcc$scores[,"RC1"]

fta.pot <- principal(table.pot, nfactors=4, rotate="varimax")
print(fta.pot)


#fta.tcc <- principal(table.tcc, cor=TRUE)
#print(fta.tcc, digits=2, cutoff=.3, sort=TRUE)

#fta.pc <- principal(table.pc, cor=TRUE)
#print(fta.pc, digits=2, cutoff=.3, sort=TRUE)

} 

table.all$tcc.sc1 <-  fta.tcc$scores[,1]
table.all$tcc.sc2 <-  fta.tcc$scores[,2] 
table.all$tcc.sc3 <-  fta.tcc$scores[,3]

table.all$pc.sc1 <-  fta.pc$scores[,1]
table.all$pc.sc2 <-  fta.pc$scores[,2] 
table.all$pc.sc3 <-  fta.pc$scores[,3]


table.all$pot.sc1 <-  fta.pot$scores[,1]
table.all$pot.sc2 <-  fta.pot$scores[,2] 
table.all$pot.sc3 <-  fta.pot$scores[,3]
table.all$pot.sc4 <-  fta.pot$scores[,4]




#go classification of gression tree
go.tree <- TRUE

if (go.tree) {
dev.new()

# Classification Tree with rpart
library(rpart)
library(rpart.plot)
go.pdf = F
if(go.pdf) {pdf(file="Fig4_updated.pdf",width=12,height=9) } else{ print("go x-windows") }
#eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex","oni.enso"
# Six most important variables are 
# eve.for.acf.aff, date.mon.day, eve.lat.med.aff, bef.for.lai.aff, bef.for.acf.aff, oni.enso

#parameter seting 
#rpart.control(cp=0.0025, minsplit = 8, maxdepth = 5)

#fit.0 <- rpart(group ~ eve.for.acf.aff+
#                       bef.for.lai.aff+bef.for.asw.aff+
#                       eve.lat.med.aff+eve.intensity+
#                       oni.enso, method="class", data=table.all)

#fit.0 <- rpart(group ~  eve.for.acf.aff + bef.for.lai.aff +  pj.index +  
#                        eve.for.spei.ref , 
#                        method="class", data=table.all, control=rpart.control(cp=0, maxdepth=4,minsplit=100))

fit.0 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
                        method="class", data=table.all, control=rpart.control(cp=0, maxdepth=4,minsplit=100))


#2D
#fit.2d_w0_p60 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_0_p_60_2D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.2d_w0_p60, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )

#dev.new()

#fit.2d_w8_p0 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_8_p_0_2D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.2d_w8_p0, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )

#dev.new()
#fit.2d_w8_p60 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_8_p_60_2D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.2d_w8_p0, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )
#dev.new()

#3D
#fit.3d_w0_p80 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_0_p_80_3D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.3d_w0_p80, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )

#dev.new()

#fit.3d_w10_p0 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_10_p_0_3D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.3d_w10_p0, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )

#dev.new()
#fit.3d_w10_p80 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_10_p_80_3D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.3d_w10_p80, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )
#dev.new()

#4D
#fit.4d_w0_p100 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_0_p_100_4D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.4d_w0_p100, box.palette=list("white","white","white") )

#dev.new()

#fit.4d_w12_p0 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_12_p_0_4D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.4d_w12_p0, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",snip=F,box.palette=list("white","white","white") )

#dev.new()
#fit.4d_w12_p100 <- rpart(group ~ tcc.sc1   + tcc.sc2  +  tcc.sc3  + pc.sc1 + pc.sc2 + pc.sc3 , 
#                        method="class", data=subset(table.all, table.all$run.case=="w_12_p_100_4D"), 
#                        control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#rpart.plot(fit.4d_w12_p100, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No", snip=F,box.palette=list("white","white","white") )
#dev.new()

fit.pot <- rpart(group ~ pot.sc1   + pot.sc2  +  pot.sc3 + pot.sc4 , 
                        method="class", data=table.all, control=rpart.control(cp=0, maxdepth=4,minsplit=100))
#dev.new()
#rpart.plot(fit.pot, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No", snip=F,box.palette=list("white","white","white") )


fit.1 = prune(fit.pot, cp = 0.0035)
#dev.new()
#rpart.plot(fit.1, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No", snip=F,box.palette=list("white","white","white") )


#pruning the modeled tree05fit.1 = prune(fit.0, cp = 0.1)
#fit.1 = prune(fit.0, cp = 0.0025)
#fit.1 = prune(fit.0, cp = 0.008)
#fit.1 = prune(fit.0, cp = 0.01)


split.fun <- function(x, labs, digits, varlen, faclen)
{
    # replace variable names in the labels
    labs  <- sub("eve.for.acf.aff", "Accumulated rainfall(mm)", labs)
    labs  <- sub("pj.index",        "Pacific Japan index", labs)
    labs  <- sub("date.mon",        "Month", labs)
    labs  <- sub("eff.cex",         "Affected area (Mha)", labs)
    labs  <- sub("bef.for.lai.aff", "Prior leaf area(m^2/m^2)", labs)
    labs  <- sub("eve.intensity",   "Cyclone intensity(m/s)", labs)
    labs  <- sub("eve.for.spei.aff","Prior drought state(mm/mm)", labs)

    labs # return the modified labels
}

#rpart.plot(fit.0, uniform=TRUE, legend.cex = 0.5, legend.x=-1, legend.y=-1,tweak=1.25, yesno=2,yes.text="Yes", no.text="No",, split.fun=split.fun, snip=F,box.palette=list("white","white","white") )



#dev.off()

#fit.snip <- snip.rpart(fit.0, toss = 61) 

#rpart.plot(fit.snip)
#rpart.plot(fit.1, uniform=TRUE, box.palette=list("orange","forestgreen"))


# Grouping base on the pot.sc1 and pot.sc2 

table.gp1 <- subset( table.all, (table.all$pot.sc2 < 0.81 & table.all$pot.sc1 < 0.8))  
#table.gp1$gp <- c("gp1")
 
table.gp2 <- subset( table.all, (table.all$pot.sc2 < 0.81 & table.all$pot.sc1 >= 0.8))
#table.gp2$gp <- c("gp2")
 
table.gp3 <- subset( table.all, table.all$pot.sc2 >= 0.81 ) 
#table.gp3$gp <- c("gp3")
 

# c("eve.for.mws.aff","eve.for.acf.aff","eve.lat.med.aff","eve.intensity","eff.cex",

# "date.mon","bef.for.acf.aff","bef.for.lai.aff","pj.index","eve.for.spei.ref","del.spei.for","eff.size.for")]



# show some statistics 
print("Group1 Statistics:")
print(paste("Latitude of landfall :", median(table.gp1$eve.lat.med.aff), " +/- ", sd(table.gp1$eve.lat.med.aff) ) ) 
print(paste("Affected area :",        3.1416*((median(table.gp1$eff.cex)*100)**2.) , " +/- ", 3.1416*((sd(table.gp1$eff.cex)*100)**2.)))
print(paste("TC rainfall :",          median(table.gp1$eve.for.acf.aff), " +/- ", sd(table.gp1$eve.for.acf.aff) ))
print(paste("Maximum windspeed :",    median(table.gp1$eve.for.mws.aff), " +/- ", sd(table.gp1$eve.for.mws.aff) ))
print(paste("Intensity gust wind :",  median(table.gp1$eve.intensity)/3.6, " +/- ", sd(table.gp1$eve.intensity)/3.6 ))
print(paste("Pacific Japan Index  :", median(table.gp1$pj.index), " +/- ", sd(table.gp1$pj.index) ))
print(paste("Prior rainfall :",       median(table.gp1$bef.for.acf.aff), " +/- ", sd(table.gp1$bef.for.acf) ))
print(paste("Month of landfall :",    median(table.gp1$date.mon), " +/- ", sd(table.gp1$date.mon) ))
print(paste("Prior LAI :",      median(table.gp1$bef.for.lai.aff), " +/- ", sd(table.gp1$bef.for.lai.aff,na.rm=T) ))
print(paste("SPEI :",           median(table.gp1$eve.for.spei.ref), " +/- ", sd(table.gp1$eve.for.spei.ref) ))
print(paste("Change of SPEI :", median(table.gp1$del.spei.for), " +/- ", sd(table.gp1$del.spei.for) ))

print("Group2 Statistics:")
print(paste("Latitude of landfall :", median(table.gp2$eve.lat.med.aff), " +/- ", sd(table.gp2$eve.lat.med.aff) )  )
print(paste("Affected area :",        3.1416*((median(table.gp2$eff.cex)*100)**2.) , " +/- ", 3.1416*((sd(table.gp2$eff.cex)*100)**2.)))
print(paste("TC rainfall :",          median(table.gp2$eve.for.acf.aff), " +/- ", sd(table.gp2$eve.for.acf.aff) ))
print(paste("Maximum windspeed :",    median(table.gp2$eve.for.mws.aff), " +/- ", sd(table.gp2$eve.for.mws.aff) ))
print(paste("Intensity gust wind :",  median(table.gp2$eve.intensity)/3.6, " +/- ", sd(table.gp2$eve.intensity)/3.6 ))
print(paste("Pacific Japan Index  :", median(table.gp2$pj.index), " +/- ", sd(table.gp2$pj.index) ))
print(paste("Prior rainfall :",       median(table.gp2$bef.for.acf.aff), " +/- ", sd(table.gp2$bef.for.acf) ))
print(paste("Month of landfall :",    median(table.gp2$date.mon), " +/- ", sd(table.gp2$date.mon) ))
print(paste("Prior LAI :",      median(table.gp2$bef.for.lai.aff), " +/- ", sd(table.gp2$bef.for.lai.aff,na.rm=T) ))
print(paste("SPEI :",           median(table.gp2$eve.for.spei.ref), " +/- ", sd(table.gp2$eve.for.spei.ref) ))
print(paste("Change of SPEI :", median(table.gp2$del.spei.for), " +/- ", sd(table.gp2$del.spei.for) ))

print("Group3 Statistics:")
print(paste("Latitude of landfall :", median(table.gp3$eve.lat.med.aff), " +/- ", sd(table.gp3$eve.lat.med.aff) )  )
print(paste("Affected area :",        3.1416*((median(table.gp3$eff.cex)*100)**2.) , " +/- ", 3.1416*((sd(table.gp3$eff.cex)*100)**2.)))
print(paste("TC rainfall :",          median(table.gp3$eve.for.acf.aff), " +/- ", sd(table.gp3$eve.for.acf.aff) ))
print(paste("Maximum windspeed :",    median(table.gp3$eve.for.mws.aff), " +/- ", sd(table.gp3$eve.for.mws.aff) ))
print(paste("Intensity gust wind :",  median(table.gp3$eve.intensity)/3.6, " +/- ", sd(table.gp3$eve.intensity)/3.6 ))
print(paste("Pacific Japan Index  :", median(table.gp3$pj.index), " +/- ", sd(table.gp3$pj.index) ))
print(paste("Prior rainfall :",       median(table.gp3$bef.for.acf.aff), " +/- ", sd(table.gp3$bef.for.acf) ))
print(paste("Month of landfall :",    median(table.gp3$date.mon), " +/- ", sd(table.gp3$date.mon) ))
print(paste("Prior LAI :",      median(table.gp3$bef.for.lai.aff), " +/- ", sd(table.gp3$bef.for.lai.aff,na.rm=T) ))
print(paste("SPEI :",           median(table.gp3$eve.for.spei.ref), " +/- ", sd(table.gp3$eve.for.spei.ref) ))
print(paste("Change of SPEI :", median(table.gp3$del.spei.for), " +/- ", sd(table.gp3$del.spei.for) ))

print("ANOVA TESTS:")

table.gp1$gp="gp1"
table.gp2$gp="gp2"
table.gp3$gp="gp3"

table.anova<- rbind(table.gp1, table.gp2, table.gp3) 

print(" Latitude  ANOVA test p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(eve.lat.med.aff ~ gp, data=table.anova)
print(summary(aov.test))

#print("Latitude of landfall group1 vs group2 + group3") 

#x <- rep(table.gp1$eve.lat.med.aff,5) 
#y1 <- rep(table.gp2$eve.lat.med.aff,10)
#y2 <- rep(table.gp3$eve.lat.med.aff,10)
 
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y1),length(y2))
#length(x) <- n                      
#length(y1) <- n
#length(y2) <- n
#df = cbind.data.frame(x, y1, y2)

#one.way <- aov(x ~ y1+y2, data=df)
#print(summary(one.way))



#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

print(" Affect Area ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test <- aov(eff.cex ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))
#x<- table.gp1$eff.cex 
#y2<- table.gp2$eff.cex
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))


#print(" Affect Area group1 vs group3") 
#x<- table.gp1$eff.cex 
#y<- table.gp3$eff.cex
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))
print(" TC rainfall ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(eve.for.acf.aff ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))

#print(" TC rainfall group1 vs group2") 
#x<- table.gp1$eve.for.acf.aff 
#y<- table.gp2$eve.for.acf.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)



#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" TC rainfall group1 vs group3") 
#x<- table.gp1$eve.for.acf.aff 
#y<- table.gp3$eve.for.acf.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))


print(" maxwind ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(eve.for.mws.aff ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))

#print(" Maxwind group1 vs group2") 
#x<- table.gp1$eve.for.mws.aff 
#y<- table.gp2$eve.for.mws.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" Maxwind group1 vs group3") 
#x<- table.gp1$eve.for.mws.aff 
#y<- table.gp3$eve.for.mws.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#rint(summary(one.way))

print(" gust ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<- aov(eve.intensity ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))

#
#print(" Gust wind group1 vs group2") 
#x<- table.gp1$eve.intensity 
#y<- table.gp2$eve.intensity
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)

#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" Gust wind group1 vs group3") 
#x<- table.gp1$eve.intensity 
#y<- table.gp3$eve.intensity
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))
print(" PJI ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(pj.index ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))

#
#print(" PJI group1 vs group2") 
#x<- table.gp1$pj.index 
#y<- table.gp2$pj.index
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" PJI group1 vs group3") 
#x<- table.gp1$pj.index 
#y<- table.gp3$pj.index
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

print(" Prior rianfall ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(bef.for.acf.aff ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))

#
#print(" Proior rainfall group1 vs group2") 
#x<- table.gp1$bef.for.acf.aff 
#y<- table.gp2$bef.for.acf.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" Proior rainfall group1 vs group3") 
#x<- table.gp1$bef.for.acf.aff 
#y<- table.gp3$bef.for.acf.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))
print(" Month  ANOVA  p<0.05 for the sinifigance of difference in groups") 
aov.test<- aov(date.mon ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))
#
#print(" month group1 vs group2") 
#x<- table.gp1$date.mon 
#y<- table.gp2$date.mon
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)
#
#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" month group1 vs group3") 
#x<- table.gp1$date.mon 
#y<- table.gp3$date.mon
#res <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)



#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

print(" Prior LAI ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(bef.for.lai.aff ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))
#


#
#print(" Prior LAI group1 vs group2") 
#x<- table.gp1$bef.for.lai.aff 
#y<- table.gp2$bef.for.lai.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" Prior LAI group1 vs group3") 
#x<- table.gp1$bef.for.lai.aff 
#y<- table.gp3$bef.for.lai.aff
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))
print("SPEI ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(eve.for.spei.ref ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))
#


#print(" SPEI group1 vs group2") 
#x<- table.gp1$eve.for.spei.ref 
#y<- table.gp2$eve.for.spei.ref
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)



#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print(" SPEI group1 vs group3") 
#x<- table.gp1$eve.for.spei.ref 
#y<- table.gp3$eve.for.spei.ref
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))
print("Delta SPEI ANOVA p<0.05 for the sinifigance of difference in groups") 
aov.test<-aov(del.spei.for ~ gp, data=table.anova)
print(summary(aov.test))
print(TukeyHSD(aov.test,conf.level=0.95))
#



#print("Chnage of  SPEI group1 vs group2") 
#x<- table.gp1$del.spei.for
#y<- table.gp2$del.spei.for
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)

#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))

#print("Change of  SPEI group1 vs group3") 
#x<- table.gp1$del.spei.for
#y<- table.gp3$del.spei.for
#res<-t.test(x, y, alternative = "two.sided", var.equal = FALSE)
#print(res)


#n <- max(length(x), length(y))
#length(x) <- n                      
#length(y) <- n
#df = cbind.data.frame(x, y)
#one.way <- aov(x ~ y, data=df)
#print(summary(one.way))



# grow tree
#fit <- rpart(group ~ pc1 + pc2 + pc3 + pc4,
#   method="class", data=table.all)
# plot tree
#rpart.plot(fit.0, uniform=TRUE,
#         main="Classification Tree for Effect Size")
#text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postscript plot of tree
#post(fit, file = "rgression_tree.ps",
#   title = "Classification Tree for Effect Size")

pdf(file="Fig2_Rainfall_SPEI.pdf", width=8,height=8)


par(mfrow = c(2, 2),     # 2x2 layout
    oma = c(2, 2, 1, 1), # two rows of text at the outer left and bottom margin one row at top and right
    mar = c(2, 1, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA)            # allow content to protrude into outer margin (and beyond)

#set color palette
my.color<- colorRampPalette(c("brown","yellow","gray","gray","green","forestgreen"))(16)
#  my.color <- viridis(n=12, alpha=1, direction=1, option="H") # D for viridus  H for turbo/rainbow
#  my.color[1] <- "gray"
#  my.breaks<- round(seq(0, 6.0, length.out = 13), digits=1)
  my.breaks<- c(-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0)
  my.breaks.txt<- c("-2.0","-1.5","-1.0","0.0","0.5","1.0","1.5","2.0","2.5")

table.gp1 <- arrange(table.gp1,factor(group, levels = c("Neutral", "Negative", "Positive")))
table.gp2 <- arrange(table.gp2,factor(group, levels = c("Neutral", "Negative", "Positive")))
table.gp3 <- arrange(table.gp3,factor(group, levels = c("Neutral", "Negative", "Positive")))

#olor selection
colors <- c("gray", # gray
            "#FDAE61", # Orange
            "#66BD63") # Darker green
rbPal <- colorRampPalette(c("gray","#66BD63"))

col_neu <- colors[1]
col_neg <- colors[2]
col_pos <- colors[3]

cex_set=2.0

# Scatter plot

table.plot <- table.anova
#This adds a column of color values
# based on the y values
table.plot$col  <- rbPal(20)[as.numeric(cut(c(0.0,1.0),breaks = 20))]
# Reorder the factor levels
#reordered_groups <- factor(table.plot$group , levels = c("Neutral","Negative","Positive"))
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Positive"  )

plot( mean(table.plot$bef.for.acf.aff) ~ mean(table.plot$eve.for.spei.ref),
     pch = c(21), cex= cex_set, xlim=c(-.5,.5), ylim=c(0,200), bg=NA,
     xlab="",ylab="Rainfall (mm)")
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5

print("all positve cases")    
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )


for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_pos)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_pos)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_pos)
}
points( y=y2, x=x2,  pch = 21, cex=cex_set, bg=col_pos)

#
table.plot <- table.anova
# Reorder the factor levels
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Neutral" )

points(y= mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref),
     pch = c(21), cex= cex_set,  ylim=c(0,200),   bg=NA)

x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5

print("all neutral cases")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )


for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neu)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neu)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neu)
}
points( y=y2, x=x2,  pch = c(21), cex= cex_set, bg=col_neu)


#dev.new()
table.plot <- table.anova
# Reorder the factor levels
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Negative"  )

#table.plot <- subset(table.plot, as.character(table.plot$group)=="Negative")
#table.plot <- subset(table.plot, as.character(table.plot$group)=="Neutral")
points( y=mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref),
     pch = c(21), cex= cex_set,  ylim=c(0,200), bg=NA)
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5

print("all negative cases")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neg)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neg)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neg)
}
points( y=y2, x=x2,  pch = c(21), cex= cex_set, bg=col_neg)

# Legend
legend("topright",
       legend = c("Prior condition", "Post period","Neutral event pathway","Positive event pathway","Negative event pathway"),
       pch = c(1,21,NA,NA,NA),cex=c(cex_set,cex_set),lty=c(NA,NA,"dashed","dashed","dashed"),
       col = c("black","black","gray","#66BD63","#FDAE61"),pt.bg = c("NA","gray"))
#col="#FDAE61")
# col="#66BD63")
#       col = colors[factor(levels(reordered_groups))])




table.plot <- table.anova
#This adds a column of color values
# based on the y values
table.plot$col  <- rbPal(20)[as.numeric(cut(c(0.0,1.0),breaks = 20))]
# Reorder the factor levels
#reordered_groups <- factor(table.plot$group , levels = c("Neutral","Negative","Positive"))
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Neutral" & as.character(table.plot$gp) == "gp1" )

plot(y= mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref),
     pch = c(0), cex=cex_set, xlim=c(-1.0,1.0), ylim=c(0,200), col=col_neu,
     xlab="",ylab="")

x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
print("neutral group1") 
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neu)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neu)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neu)
}
points( y=y2, x=x2,  pch = c(15), cex=cex_set, col=col_neu)


#overlape group 2
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Neutral" & as.character(table.plot$gp) == "gp2" )
points(y=mean(table.plot$bef.for.acf.aff),x= mean(table.plot$eve.for.spei.ref), 
pch = c(1), cex= cex_set, col=col_neu)

x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
print("neutral group2") 
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neu)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neu)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neu)
}
points( y=y2, x=x2,  pch = c(16), cex= cex_set, col=col_neu)

#overlape group 3
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Neutral" & as.character(table.plot$gp) == "gp3" )
points(y=mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref),
 pch = c(2), cex= cex_set, col=col_neu)
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
 print("neutral group3")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neu)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neu)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neu)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neu)
}
points( y=y2, x=x2,  pch = c(17), cex= cex_set, col=col_neu)


# Legend
legend("topright",
       legend = c("Group-1", "Group-2", "Group-3"),
       pch = c(0,1,2),  col=col_neu)


table.plot <- table.anova
#This adds a column of color values
# based on the y values
table.plot$col  <- rbPal(20)[as.numeric(cut(c(0.0,1.0),breaks = 20))]
# Reorder the factor levels
#reordered_groups <- factor(table.plot$group , levels = c("Neutral","Negative","Positive"))
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Positive" & as.character(table.plot$gp) == "gp1" )

plot( mean(table.plot$bef.for.acf.aff) ~ mean(table.plot$eve.for.spei.ref),
      pch = c(0), cex= cex_set, xlim=c(-1.0,1.0), ylim=c(0,200), col=col_pos,
     xlab="Standardized Precipitation and Evapotranspiration Index", ylab="Rainfall (mm)" )
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
   print("postive group1")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_pos)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_pos)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_pos)
}
points( y=y2, x=x2,  pch = c(15), cex= cex_set, col=col_pos)
#    col = colors[reordered_groups])

#overlape group 2
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Positive" & as.character(table.plot$gp) == "gp2" )
points(mean(table.plot$bef.for.acf.aff) ~ mean(table.plot$eve.for.spei.ref), 
pch = c(1), cex=cex_set, col=col_pos)
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
   print("positive group2")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_pos)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_pos)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_pos)
}
points( y=y2, x=x2,  pch = c(16), cex= cex_set, col=col_pos)
#    col = colors[reordered_groups])


#overlape group 3
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Positive" & as.character(table.plot$gp) == "gp3" )
points(mean(table.plot$bef.for.acf.aff) ~ mean(table.plot$eve.for.spei.ref),
 pch = c(2), cex= cex_set, col=col_pos)
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
   print("positive group3")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_pos)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_pos)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_pos)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_pos)
}
points( y=y2, x=x2,  pch = c(17), cex= cex_set, col=col_pos)


# Legend
legend("topright",
       legend = c("Group-1", "Group-2", "Group-3"),
       pch = c(0,1,2), col=col_pos)

#dev.new()


#dev.new()
table.plot <- table.anova
#This adds a column of color values
# based on the y values
table.plot$col  <- rbPal(20)[as.numeric(cut(c(0.0,1.0),breaks = 20))]
# Reorder the factor levels
#reordered_groups <- factor(table.plot$group , levels = c("Neutral","Negative","Positive"))
reordered_groups <- factor(table.plot$gp , levels = c("gp1","gp2","gp3"))
table.plot <- subset(table.plot, as.character(table.plot$group)=="Negative" & as.character(table.plot$gp) == "gp1" )

plot( y=mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref),
     pch = c(0), cex= cex_set, xlim=c(-1.0,1.0), ylim=c(0,200), col=col_neg,
     xlab="Standardized Precipitation and Evapotranspiration Index", ylab="")
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
   print("negative group1")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neg)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neg)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neg)
}
points( y=y2, x=x2,  pch = c(15), cex= cex_set, col=col_neg)
#    col = colors[reordered_groups])


#overlape group 2
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Negative" & as.character(table.plot$gp) == "gp2" )

points(y=mean(table.plot$bef.for.acf.aff), x= mean(table.plot$eve.for.spei.ref), 
pch = c(1), cex= cex_set, col=col_neg)

x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
   print("negative group2")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neg)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neg)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neg)
}
points( y=y2, x=x2,  pch = c(16), cex=cex_set, col=col_neg)

#overlape group 3
table.plot <- table.anova
table.plot <- subset(table.plot, as.character(table.plot$group)=="Negative" & as.character(table.plot$gp) == "gp3" )
points(y=mean(table.plot$bef.for.acf.aff),x= median(table.plot$eve.for.spei.ref),
 pch = c(2), cex= cex_set, col="#FDAE61")
x1<- mean(table.plot$eve.for.spei.ref)
x2<- mean(table.plot$eve.for.spei.ref) + mean(table.plot$del.spei.for)*mean(table.plot$eff.size.for)
del_x1 <- sd(table.plot$eve.for.spei.ref)*.5
del_x2 <- sd(table.plot$del.spei.for*table.plot$eff.size.for)*0.5

y1<- mean(table.plot$bef.for.acf.aff)
y2<- mean(table.plot$eve.for.acf.aff) + mean(table.plot$bef.for.acf.aff)
del_y1 <- sd(table.plot$bef.for.acf.aff)*.5
del_y2 <- sd(table.plot$eve.for.acf.aff)*.5
    print("negative group3")
print(paste("prior SPEI mean:", x1, "prior SPEI std:", del_x1) )
print(paste("post  SPEI mean:", x2, "post  SPEI std:", del_x2) )
print(paste("prior rainfall mean:", y1, "prior rainfall std:", del_y1) )
print(paste("post  rainfall mean:", y2, "post  rainfall std:", del_y2) )

for (i in 1:length(x1)) {
 lines(x=c(x1[i],x2[i]), y=c(y1[i],y2[i]),type="l",lty="dashed",col=col_neg)
 #prior uncetainty 
 lines(x=c(x1[i] - del_x1 , x1[i] + del_x1), y=c(y1[i],y1[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x1[i],x1[i]), y=c(y1[i] - del_y1, y1[i] + del_y1),type="l",lty="solid",col= col_neg)
 #poest period uncertainty
 lines(x=c(x2[i] - del_x2 , x2[i] + del_x2), y=c(y2[i],y2[i]),type="l",lty="solid",col= col_neg)
 lines(x=c(x2[i],x2[i]), y=c(y2[i] - del_y2, y2[i] + del_y2),type="l",lty="solid",col= col_neg)
}
points( y=y2, x=x2,  pch = c(17), cex= cex_set, col=col_neg)

# Legend
legend("topright",
       legend = c("Group-1", "Group-2", "Group-3"),
       pch = c(0,1,2),  col=col_neg)

dev.off()




# subset the table for different group

table.pos <- subset(table.all, table.all$group == "Positive") 
table.neu <- subset(table.all, table.all$group == "Neutral") 
table.neg <- subset(table.all, table.all$group == "Negative") 


# get the statistic from boxplot funtion
#  a matrix, each column contains the extreme of the lower
#          whisker[1], the lower hinge[2], the median[3], the upper hinge[4] and the
#          extreme of the upper whisker[5] for one group.

#stat.pos.aff <- boxplot(aft.for.lai.aff ~ date.mon, data=table.pos, ylim=c(1,6)) 
#stat.pos.ref <- boxplot(aft.for.lai.ref ~ date.mon, data=table.pos, ylim=c(1,6)) 

#stat.neu.aff <- boxplot(aft.for.lai.aff ~ date.mon, data=table.neu, ylim=c(1,6)) 
#stat.neu.ref <- boxplot(aft.for.lai.ref ~ date.mon, data=table.neu, ylim=c(1,6)) 

#stat.neg.aff <- boxplot(aft.for.lai.aff ~ date.mon, data=table.neg, ylim=c(1,6)) 
#stat.neg.ref <- boxplot(aft.for.lai.ref ~ date.mon, data=table.neg, ylim=c(1,6)) 



# plot the box plot
library(ggplot2)
library(tidyverse)
#theme_set(theme_bw(16))
#table.all %>% 
 
ld_go <- F
if(ld_go) { 
lai.plot <- ggplot(table.all ,aes(x=date.mon,y=aft.for.lai.ref)) + 
          geom_boxplot(aes(group=date.mon),  outlier.size=0) + 
          geom_smooth(data=table.neu, col="lightblue" ,  aes(x=date.mon,y=aft.for.lai.aff), lwd=2 ) +
          geom_smooth(data=table.neg, col="orange",      aes(x=date.mon,y=aft.for.lai.aff), lwd=2 ) + 
          geom_smooth(data=table.pos ,col="forestgreen", aes(x=date.mon,y=aft.for.lai.aff), lwd=2 ) +
          theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20),
                axis.title.y = element_text(size = rel(1.5), angle = 90),
                axis.title.x = element_text(size = rel(1.5), angle = 00), 
                axis.text = element_text(size = rel(1.2)  ),
                axis.text.x = element_text(size = rel(1.2) ),
                axis.text.y = element_text(size = rel(1.2) ),  
                legend.text = element_text(size = rel(1.2) ),
                legend.title = element_text(size = rel(1.2), face="bold" ), legend.position="top"  ) +
         #labs(x ="Month", y = "Dynamics of LAI after TC disturbances", tag="")+
         scale_x_continuous(name="Month", breaks=seq(0,12,1), labels=seq(0,12,1), limits=c(0.5,12.5) )+
         scale_y_continuous(name="Vegetation state (LAI) after TC landfall", breaks=seq(0,6,1), labels=seq(0,6,1), limits=c(0,6.2) )

acf.plot <- ggplot(table.all ,aes(x=date.mon,y= eve.for.acf.aff)) + 
          geom_boxplot(aes(group=date.mon),  outlier.size=0) + 
          geom_smooth(data=table.neu, col="lightblue" ,  aes(x=date.mon,y= eve.for.acf.aff), lwd=2 ) +
          geom_smooth(data=table.neg, col="orange",      aes(x=date.mon,y= eve.for.acf.aff), lwd=2 ) + 
          geom_smooth(data=table.pos ,col="forestgreen", aes(x=date.mon,y= eve.for.acf.aff), lwd=2 ) +
          theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20),
                axis.title.y = element_text(size = rel(1.5), angle = 90),
                axis.title.x = element_text(size = rel(1.5), angle = 00), 
                axis.text = element_text(size = rel(1.2)  ),
                axis.text.x = element_text(size = rel(1.2) ),
                axis.text.y = element_text(size = rel(1.2) ),  
                legend.text = element_text(size = rel(1.2) ),
                legend.title = element_text(size = rel(1.2), face="bold" ), legend.position="top"  ) +
         #labs(x ="Month", y = "Dynamics of LAI after TC disturbances", tag="")+
         scale_x_continuous(name="Month", breaks=seq(0,12,1), labels=seq(0,12,1), limits=c(0.5,12.5) )+
         scale_y_continuous(name="Area average TC rainfall (mm)", breaks=seq(0,150,20), labels=seq(0,150,20), limits=c(0,150) )


mws.plot <- ggplot(table.all ,aes(x=date.mon,y= eve.for.mws.aff)) + 
          geom_boxplot(aes(group=date.mon),  outlier.size=0) + 
          geom_smooth(data=table.neu, col="lightblue" ,  aes(x=date.mon,y= eve.for.mws.aff), lwd=2 ) +
          geom_smooth(data=table.neg, col="orange",      aes(x=date.mon,y= eve.for.mws.aff), lwd=2 ) + 
          geom_smooth(data=table.pos ,col="forestgreen", aes(x=date.mon,y= eve.for.mws.aff), lwd=2 ) +
          theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20),
                axis.title.y = element_text(size = rel(1.5), angle = 90),
                axis.title.x = element_text(size = rel(1.5), angle = 00), 
                axis.text = element_text(size = rel(1.2)  ),
                axis.text.x = element_text(size = rel(1.2) ),
                axis.text.y = element_text(size = rel(1.2) ),  
                legend.text = element_text(size = rel(1.2) ),
                legend.title = element_text(size = rel(1.2), face="bold" ), legend.position="top"  ) +
         #labs(x ="Month", y = "Dynamics of LAI after TC disturbances", tag="")+
         scale_x_continuous(name="Month", breaks=seq(0,12,1), labels=seq(0,12,1), limits=c(0.5,12.5) )+
         scale_y_continuous(name="Area average maximum wind speed (m/s)", breaks=seq(5,25,5), labels=seq(5,25,5), limits=c(5,25) )



  #          mapping = aes(x = grp, y = average, group=1))
#ggsave("boxplot_with_line_connecting_mean_values_ggplot2_R.png")


# analysis the LAI signal 

# Difference of spei plot
 
#dev.new()

#table.all$group   
  
# Use semi-transparent fill
#p <- ggplot(df, aes(x=weight, fill=sex)) +
#            geom_density(alpha=0.4)
#p
# Add mean lines
#p+geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
#             linetype="dashed")
spei.table <- subset(table.all, table.all$group!="Neutral" &  abs(table.all$eff.size.for) > 0.45)
#spei.table <- table.all
#dev.new()
pdf(file="SPEI_plot.pdf", width=8,height=8)
p1 <- ggplot(spei.table, aes(x= eve.for.spei.aff - eve.for.spei.ref , fill=group)) + geom_density(alpha=0.4) + xlim(-1.5,1.5) 
#p1 <- p1 + theme_bw()
p1 <- p1 + ggtitle("Standard Precpipitation - Evaporation Index (SPEI) Analysis") + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=18))
p1 <- p1 + labs(x = "SPEI Change (Affected Area - Reference Area)", y = " Density ")
p1 <- p1 + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))
p1 <- p1 + theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) 
p1 <- p1 + theme(axis.text = element_text(size = rel(1.0)  )  )
p1 <- p1 + theme(axis.text.x = element_text(size = rel(1.0) ) ) 
p1 <- p1 + theme(axis.text.y = element_text(size = rel(1.0) ) )
p1 <- p1 + theme(legend.text = element_text(size = rel(1.0) ) )
p1 <- p1 + theme(legend.title = element_text(size = rel(1.0) ) )
#
#p1
grid.arrange(p1, p1, p1, p1, nrow=2, ncol=2)  # Make a list of plots
#ggsave( plots, 
#        filename="SPEI_plot.pdf", 
#         width=8, height=8) 
dev.off()
#p2<-ggplot(spei.table, aes(x=eve.for.spei.ref+eve.for.spei.aff, fill=group)) + geom_density(alpha=0.4)
#p2


} #end ld_go

# add a note of  drought-busting & flooding events in the end of table
 
spei_aft  <-  table.all$eve.for.spei.aff
spei_pri  <-  table.all$eve.for.spei.aff.pri


#table.all$drought.buster <- ifelse( ((spei_pri <= -0.455555)  ((spei_aft - spei_pri) > 0.2) )   , "TRUE", "FALSE") 
table.all$drought.buster <- ifelse( ((spei_pri) >= -5.455555  )   , "TRUE", "FALSE") 
table.all$flooding <- ifelse( (spei_aft >= 0.5 ), "TRUE", "FALSE") 
table.all$increase_spei  <- ifelse( ( spei_aft >= spei_pri), "TRUE", "FALSE") 





if(go.pdf) {dev.off()}
#printcp(fit) # display the results
#dev.new()
#plotcp(fit.0) # visualize cross-validation results
#summary(fit) # detailed summary of splits


# output the table as csv.file
write.table(table.all, file="all.table.csv", sep=",",row.names=F, na = "NA")


} 




ld_go <- F

if (ld_go) {

# output the statistic of different runs with differnet sinarious

#get runs
runs <- levels( as.factor(table.all$run.case) ) 


tot<-array(0,dim=length(runs))
neg<-array(0,dim=length(runs))
neu<-array(0,dim=length(runs))
pos<-array(0,dim=length(runs)) 
bst<-array(NA,dim=length(runs)) 
fld<-array(NA,dim=length(runs))

neg.bst <- array(NA,dim=length(runs))
neu.bst <- array(NA,dim=length(runs))
pos.bst <- array(NA,dim=length(runs))

neg.fld <- array(NA,dim=length(runs))
neu.fld <- array(NA,dim=length(runs))
pos.fld <- array(NA,dim=length(runs))


 
for ( irun in 1:length(runs) ) {  
  table.sub <- subset(table.all, table.all$run.case == as.character(runs[irun]))  
  #print(table.sub) 
  tot[irun] <- length(table.sub$run.case)
  neg[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Negative"] ) 
  neu[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Neutral"] ) 
  pos[irun] <- length(  table.sub$run.case[as.character(table.sub$group) == "Positive"] ) 
  bst[irun] <- length(  table.sub$run.case[as.character(table.sub$drought.buster) == "TRUE"] ) 
  fld[irun] <- length(  table.sub$run.case[as.character(table.sub$flooding) == "TRUE"] ) 

  neg.bst[irun] <- length( table.sub$run.case[ (as.character(table.sub$drought.buster) == "TRUE") & (as.character(table.sub$group) == "Negative") ] )
  neu.bst[irun] <- length( table.sub$run.case[ (as.character(table.sub$drought.buster) == "TRUE") & (as.character(table.sub$group) == "Neutral")  ] )
  pos.bst[irun] <- length( table.sub$run.case[ (as.character(table.sub$drought.buster) == "TRUE") & (as.character(table.sub$group) == "Positive")  ]  )

  neg.fld[irun] <- length( table.sub$run.case[ (as.character(table.sub$flooding) == "TRUE") & (as.character(table.sub$group) == "Negative") ] )
  neu.fld[irun] <- length( table.sub$run.case[ (as.character(table.sub$flooding) == "TRUE") & (as.character(table.sub$group) == "Neutral")  ] )
  pos.fld[irun] <- length( table.sub$run.case[ (as.character(table.sub$flooding) == "TRUE") & (as.character(table.sub$group) == "Positive")  ]  )


  print( paste("Run:", runs[irun] ))
  print( paste("Totoal n:",  tot[irun] ))  
  print( paste("Negative:",  neg[irun] ))  
#  print( paste("Negative_Flooding:",        neg.fld[irun] ))  
#  print( paste("Negative_Drought-busting:", neg.bst[irun] ))  

  print( paste("Neutral:",   neu[irun] ))  
#  print( paste("Neutral_Flooding:",         neu.fld[irun] ))  
#  print( paste("Neutral_Drought-busting:",  neu.bst[irun] ))  

  print( paste("Postive:",   pos[irun] ))  
#  print( paste("Postive_Flooding:",         pos.fld[irun] ))  
#  print( paste("Postive_Drought-busting:",  pos.bst[irun] ))  

#  print( paste("Buster event:", bst[irun] ))
#  print( paste("Flood event:", fld[irun]  ))



} 


print(paste("Mean Total-N:", format(mean(tot), digits=0),sep="")    )
print(paste("Std Total-N:",  format(sd(tot), digits=0)  ,sep="")    )

print(paste("Mean Negative-N:", format(mean(neg),digits=0),sep="")  )
print(paste("Std Nagative-N:",  format(sd(neg),digits=0)  ,sep="")  )

print(paste("Mean Negative-N:", format(mean(neu),digits=0) ,sep="") )
print(paste("Std Nagative-N:",  format(sd(neu), digits=0)  ,sep="") )

print(paste("Mean Negative-N:", format(mean(pos),digits=0),sep="")  )
print(paste("Std Nagative-N:",  format(sd(pos),digits=0), sep="")   )

print(paste("Total Events:",    format(sum(tot),digits=0),sep="")  )
print(paste("Buster Evnets:",   format(sum(bst),digits=0), sep="")   )
print(paste("Floods Evnets:",   format(sum(fld),digits=0), sep="")   )




# add bar plot to show the pri and post SPEI and ES_LAI

# sorting the table all by prior SPEI value

#sort by mpg (ascending) and cyl (descending)
#newdata <- mtcars[order(mpg, -cyl),]
table.all <- table.all[order(-table.all$eve.for.spei.aff.pri),]
# subset to QC/QA pass eve
table.qcqa <- subset(table.all, (table.all$drought.buster==TRUE))
#table.qcqa <- subset(table.all, ( table.all$flooding == TRUE ) )
#table.qcqa <- subset(table.all, ( table.all$drought.buster== TRUE )&(table.all$qc1.flag==TRUE | table.all$qc2.flag==TRUE) )
#table.qcqa <- subset(table.all, ( table.all$flooding == FALSE  & table.all$drought.buster== FALSE) )
#



table.qcqa$col <- as.character(table.qcqa$group)
for (i in 1:length(table.qcqa$col)) {
   if ( (table.qcqa$col[i]) == "Neutral") table.qcqa$col[i] <- "cadetblue1"  
   if ( (table.qcqa$col[i]) == "Positive") table.qcqa$col[i] <- "forestgreen"  
   if ( (table.qcqa$col[i]) == "Negative") table.qcqa$col[i] <- "orange"  
}

# plot the figure 
dev.off()
dev.new()

csa <- function(x1,y1,x2,y2,arr.size,first.col,second.col, ...) {
    cols <- colorRampPalette( c(first.col,second.col) )(250)
    x <- approx(c(0,1),c(x1,x2), xout=seq(0,1,length.out=251))$y
    y <- approx(c(0,1),c(y1,y2), xout=seq(0,1,length.out=251))$y
    
    arrows(x[250],y[250],x[251],y[251], col=cols[250],length=arr.size, ...)
    segments(x[-251],y[-251],x[-1],y[-1],col=cols, ...)

}


color.scale.arrow <- Vectorize(csa, c('x1','y1','x2','y2') )


neve = length(table.qcqa$col)
#plot(0,0, xlim=c(-2.,2.), ylim=c(-10*neve,neve*20),col="white",xlab="SPEI", ylab="", yaxt="n")
tmp=neve*20
for (i in neve:1) {
  #add  retangular bar
   if (table.qcqa$eve.for.spei.aff.pri[i] > table.qcqa$eve.for.spei.aff[i]) {
      x_left = table.qcqa$eve.for.spei.aff[i] 
      x_right = table.qcqa$eve.for.spei.aff.pri[i] 
   }else {
      x_left = table.qcqa$eve.for.spei.aff.pri[i]
      x_right = table.qcqa$eve.for.spei.aff[i]
   }
  y_bottom = tmp - (abs(table.qcqa$eff.size.for[i])*10) - 25.0
  y_top = tmp 
  
  #rect(xleft =  x_left,  ybottom = y_bottom,
  #     xright = x_right, ytop = y_top ,
  #       col = table.qcqa$col[i] , border = NA)

#  color.scale.arrow(x1=table.qcqa$eve.for.spei.aff.pri[i],y1=y_bottom,
#                    x2=table.qcqa$eve.for.spei.aff[i],    y2=y_bottom, 
#                    table.qcqa$col[i], table.qcqa$col[i], lwd=(table.qcqa$tot.n[i])/400000,
#                    arr.size=0.05)
 
#print(paste("x_left:",x_left, "x_right:",x_right, "ES:",abs(table.qcqa$eff.size.for[i])))

 tmp = y_bottom
} 



# add 



par(mar = c(5, 6, 1, 2) , mgp=c(3,1,0))
 
par(oma = c(1, 1, 0, 0))

#go.pdf = F
#if(go.pdf) {pdf(file="FigS4.pdf",width=8,height=6) } else{ print("go x-windows") }


#boxobj <- boxplot( table.all$eff.size.for ~ round(table.all$eve.tc.occ.aff,digits=0), outline=FALSE, 
#                       ylim=c(-0.8,1.5),cex.lab=1.5, cex.axis=1.2, 
#                       names = c("(0,0.5]","(0.5,1.5]","(1.5,2.5]","(2.5,3.5]","(3.5,4.5]"," > 4.5"),
#                       ylab="Effect size (unitless)", xlab=expression("Return frequency (" * yr^-1 * ")"))

pos.table <- subset(table.all, table.all$eff.size.for>0 & table.all$eve.tc.occ.aff < 4.5)
occ.pos <- pos.table$eve.tc.occ.aff
eff.pos <- pos.table$eff.size.for

neg.table <- subset(table.all, table.all$eff.size.for<0 & table.all$eve.tc.occ.aff < 4.5 ) 
occ.neg <- neg.table$eve.tc.occ.aff
eff.neg <- neg.table$eff.size.for

## gg scatter plot 
library(ggpubr)
library(gridExtra)

pdf(file="FigS4_pos_neg.pdf",width=8,height=6)


ytxtleft = c("(a) (b)")


p1 <- ggscatter(data=pos.table, y="eff.size.for", x="eve.tc.occ.aff", add="loess", conf.int=TRUE, 
           ylab="Positve effect size (unitless)", xlab=expression("Return frequency (" * yr^-1 * ")"))

p2 <- ggscatter(data=neg.table, y="eff.size.for", x="eve.tc.occ.aff", add="loess", conf.int=TRUE, 
           ylab="Negative effect size (unitless)", xlab=expression("Return frequency (" * yr^-1 * ")"))

grid.arrange(p1, p2, ncol=1, nrow=2, left=ytxtleft, heights=c(1.2, 1))           

dev.off()




#smooth lines 
#fit.3.pos <- lm(eff.pos ~ poly(occ.pos,3)) 
#smooth lines
#fit.3.neg <- lm(eff.neg ~ poly(occ.neg,3)) 

#new.occ.pos <- seq(min(occ.pos), max(occ.pos), length.out=length(eff.pos))
#new.occ.neg <- seq(min(occ.neg), max(occ.neg), length.out=length(eff.neg))

#new.eff.pos <- predict.lm(fit.3.pos, newdata=(list(x=new.occ.pos)))
#new.eff.neg <- predict.lm(fit.3.neg, newdata=(list(x=new.occ.neg)))



#plot(table.all$eff.size.for ~ table.all$eve.tc.occ.aff,
#     ylim=c(-2.,2.),cex.lab=1.5, cex.axis=1.2, xlim=c(0,4.5),
#     ylab="Effect size (unitless)", xlab=expression("Return frequency (" * yr^-1 * ")"))

#lines(new.occ.pos, new.eff.pos, col="blue", lty=4)

#lines(new.occ.neg, new.eff.neg, col="red", lty=4) 

#abline(h=0, col="gray",lty="dashed")

# Add sample size on top
#nbGroup <- nlevels(as.factor(round(table.all$eve.tc.occ.aff,digits=0)))

#gp.txt <- c("a","a","b","c","a","a")

#text( 
#  x=c(1:nbGroup), 
#  y=boxobj$stats[nrow(boxobj$stats),] + 0.1, 
#  paste( gp.txt ,sep=""), font=2,cex=1.5,  
#)


novaD<-aov( table.all$eff.size.for ~ as.factor(round(table.all$eve.tc.occ.aff,digits=0)))



if(go.pdf) {dev.off()}



} # end ld_go 

