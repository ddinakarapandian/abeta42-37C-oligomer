# AFP Functions

# #Set heatmap3 color bar parameters
# breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
# barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red")
# 
# breakBarColorsCor=c(-200,seq(-1, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
# barColorsCor = colorpanel(length(breakBarColorsCor)-1, "black", "white", "orange")

# colorblind-friendly palette for figures
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
palette(cbPalette)


# load AFP's favorite packages ###############
AFP_load <- function(){
  if (!require("pacman")) install.packages("pacman")
  if (!require("BiocManager")) install.packages("BiocManager")
  # Load contributed packages with pacman
  pacman::p_load(pacman, party, rio,gplots,readxl,matrixStats,stringi,heatmap3,ClassDiscovery,Hmisc,RColorBrewer,stats,ropls,rstatix,ggpubr,ggforce)
  library("reshape")
  library("tidyverse")
}


# rotate.opls ##################

#Used to rotate scores and loadings by specified degrees, then perform axis swaps and/or flips. Output is list of scores, loadings, and rotation matrix.
rotate.opls <- function(plsOut,degrees,swap=F,flip1=F,flip2=F){
  #rotate by specified degrees
  theta=degrees*pi/180
  rotmat=rbind(c(cos(theta),-sin(theta)),c(sin(theta), cos(theta)))
  scores = plsOut@scoreMN %*% rotmat
  T1=scores[,1]
  T2=scores[,2]
  loadings = plsOut@loadingMN %*% rotmat
  P1=loadings[,1]
  P2=loadings[,2]
  
  #swap axis 1 and 2 if specified
  if(swap==TRUE){
    P1temp=P1; P2temp=P2
    P1=P2temp; P2=P1temp
    T1temp=T1; T2temp=T2
    T1=T2temp; T2=T1temp
  }
  
  #Flip axes 1 and/or 2 if specified
  if(flip1==TRUE){
    P1=-P1
    T1=-T1
  }
  if(flip2==TRUE){
    P2=-P2
    T2=-T2
  }
  out <- list()
  out$T1 <- T1
  out$T2 <- T2
  out$P1 <- as.matrix(P1)
  out$P2 <- as.matrix(P2)
  out$rotmat <- rotmat
  return(out)
}


# scores.plot ############################

#Used to create scores plot colored by categorical variable of interest, saved to working directory. No output.
scores.plot <- function(T1,T2,color,color.str="",analysis="PLSDA",plot.title=NULL){
  
  if (!require("tidyverse")) install.packages("tidyverse")
  if(is.null(plot.title)){plot.title=analysis}
  if(analysis=="PCA"){axis.label="PC"}else{axis.label="LV"}
  ScoresPlot=data.frame(T1,T2,as.factor(color))
  colnames(ScoresPlot)=c("T1","T2","color")
  plotOut <- ggplot(ScoresPlot, aes(x=T1, y=T2, color=as.factor(color))) + 
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    geom_point(size=6)+
    xlab(paste0("Scores on ",axis.label,"1"))+
    ylab(paste0("Scores on ",axis.label,"2"))+
    #stat_ellipse()+
    ggtitle(plot.title)+
    xlim(-1.1*max(abs(T1)),1.1*max(abs(T1)))+
    ylim(-1.1*max(abs(T2)),1.1*max(abs(T2)))+
    theme(panel.background = element_rect(fill = 'white',colour='black'), 
          text = element_text(size=20),
          panel.border=element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.ticks=element_line(color="black")
    )+
    labs(color = color.str)
  return(plotOut)
}


# scores.plot.ellipse ############################

#Used to create scores plot colored by categorical variable of interest, saved to working directory. No output.
scores.plot.ellipse <- function(T1,T2,color,color.str="",analysis="PLSDA",plot.title=NULL,level=0.95){
  
  if (!require("tidyverse")) install.packages("tidyverse")
  if(is.null(plot.title)){plot.title=analysis}
  if(analysis=="PCA"){axis.label="PC"}else{axis.label="LV"}
  ScoresPlot=data.frame(T1,T2,as.factor(color))
  colnames(ScoresPlot)=c("T1","T2","color")
  plotOut <- ggplot(ScoresPlot, aes(x=T1, y=T2, color=as.factor(color))) + 
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    geom_point(size=6)+
    xlab(paste0("Scores on ",axis.label,"1"))+
    ylab(paste0("Scores on ",axis.label,"2"))+
    stat_ellipse(level=level)+
    ggtitle(plot.title)+
    xlim(-1.1*max(abs(T1)),1.1*max(abs(T1)))+
    ylim(-1.1*max(abs(T2)),1.1*max(abs(T2)))+
    theme(panel.background = element_rect(fill = 'white',colour='black'), 
          text = element_text(size=20),
          panel.border=element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.ticks=element_line(color="black")
    )+
    labs(color = color.str)
  return(plotOut)
}


# scores.plot.gradient ##########################

#Used to create scores plot colored by quantitative variable of interest, saved to working directory. No output.
scores.plot.gradient <- function(T1,T2,color,color.str,analysis="PLSDA",plot.title=NULL){
  if (!require("tidyverse")) install.packages("tidyverse")
  if(is.null(plot.title)){plot.title=analysis}
  if(analysis=="PCA"){axis.label="PC"}else{axis.label="LV"}
  breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
  barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red2")
  sc = scale_colour_gradientn(colours = barColors, limits=c(min(color), max(color)))
  ScoresPlot=data.frame(T1,T2,color)
  colnames(ScoresPlot)=c("T1","T2","color")
  plotOut <-    ggplot(ScoresPlot, aes(x=T1, y=T2, color=color)) + 
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    geom_point(size=6)+
    sc+
    xlab(paste0("Scores on ",axis.label,"1"))+
    ylab(paste0("Scores on ",axis.label,"2"))+
    ggtitle(plot.title)+
    xlim(-1.1*max(abs(T1)),1.1*max(abs(T1)))+
    ylim(-1.1*max(abs(T2)),1.1*max(abs(T2)))+
    theme(panel.background = element_rect(fill = 'white',colour='black'), 
          text = element_text(size=20),
          panel.border=element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.ticks=element_line(color="black")
    )+
    labs(color = color.str)
  return(plotOut)
}


# loadings.chart #####################


#Used to create loadings bar plots sorted by value and scaled by diving by max value, saved to working directory. No output.
loadings.chart <- function(loadings,component_str,lim=1,cex=1){
  P1Sort=sort(t(loadings))/max(abs(loadings))
  indP1=sort(loadings, index.return=TRUE)$ix
  plotOut <- barplot(P1Sort,
          main = paste0("Signals in ",component_str),
          xlab = "",
          ylab = "",
          cex.names = cex,
          names.arg = rownames(loadings)[indP1],
          col = "lightblue",
          horiz = FALSE,
          las=2,
          ylim=c(-lim,lim))
  return(plotOut)
}

loadings.chart2 <- function(loadings,component_str,x_size=10){
  if (!require("ggpubr")) install.packages("ggpubr")
  library("ggpubr")
  P1Sort=sort(t(loadings))/max(abs(loadings))
  indP1=sort(loadings, index.return=TRUE)$ix
  
  df <- data.frame(P1Sort,rownames(loadings)[indP1])
  colnames(df) <- c("Loadings","Analytes")
  plotOut <- ggbarplot(df,x="Analytes",y="Loadings",fill="lightblue",xlab = "",ylab="",width = 0.8) +
    rotate_x_text() +
    ggtitle(paste0("Signals in ",component_str)) +
    theme(plot.title = element_text(hjust=0.5,face = "bold"),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = x_size))
  return(plotOut)
}



# loadings.chart.sd ############################

#Used to create loadings bar plots with sd error bars sorted by value and scaled by diving by max value, saved to working directory. No output.
loadings.chart.sd <- function(loadings,sd,component_str,lim=1){
  P1Sort=sort(t(loadings))/max(abs(loadings))
  sd=sd/max(abs(loadings))
  indP1=sort(loadings, index.return=TRUE)$ix
  barCenters <- barplot(P1Sort,
                        main = paste0("Signals in ",component_str),
                        xlab = "",
                        ylab = "",
                        names.arg = rownames(loadings)[indP1],
                        col = "lightblue",
                        horiz = FALSE,
                        las=2,
                        ylim=c(-lim,lim))
  segments(barCenters,P1Sort-sd[indP1],barCenters,P1Sort+sd[indP1])
  arrows(barCenters,P1Sort-sd[indP1],barCenters,P1Sort+sd[indP1],code=3,angle=90,length=0.05)
  return(barCenters)
}

# color.bar ########################

#Define function for making color bars (taken from http://www.colbyimaging.com/wiki/statistics/color-bars)
color.bar = function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


# opls.LOOCV ###########################

# Function opls.LOOCV
# executes function opls() for PLSR/PLSDA and conducts a leave-n out cross validation
# note: the term LOOCV (leave one out cross validation) is used broadly below, but changing value of n in input will leave-n out for all output

#Input: 
# x = x input for opls(), the numeric matrix of samples (rows) and analytes (columns)
# y = y input for opls(), either an array of regression variable values (PLSR) or an array of factors for group membership (PLSDA)
# n = 1, the number of randomized samples to be left out per iteration, must be an integer
# runs = 100, the number of iterations to conduct leave-n out cross validation, must be an integer

#Output:
# A list with the following objects:
# loadings = loadings in LV1 (first column) and LV2 (second column) from opls() of all samples (not LOOCV)
# scores = scores in LV1 (first column) and LV2 (second column) from opls() of all samples (not LOOCV)
# loadings_LOOCV = mean values of loadings in LV1 (first column) and LV2 (second column) from LOOCV iterations
# LV1_mat = matrix of LV1 output of each LOOCV iteration (rows), for use in calculating standard deviation after rotation
# LV2_mat = matrix of LV2 output of each LOOCV iteration (rows), for use in calculating standard deviation after rotation

opls.LOOCV <- function(x,y,n=1,runs=100){
  require("ropls"); require("matrixStats")
  oplsOut <- opls(x,y,predI=2)
  Load_LV1_LOOCV=matrix(0, nrow=runs, ncol=dim(x)[2])
  Load_LV2_LOOCV=matrix(0, nrow=runs, ncol=dim(x)[2])
  #iterations of opls
  for (i in 1:runs)
  {
    leftOut=sample(1:dim(x)[1],n,replace=FALSE) #creates randomized array of sample indices
    x_LOOCV=apply(x[-leftOut,],2,scale)
    indConstantColumns=colSds(x_LOOCV)<1e-10  #Remove constant columns and then run the PLS
    opls.out_LOOCV <- opls(x_LOOCV[,!indConstantColumns], y[-leftOut],predI=2,fig.pdfC="none",info.txtC="none")    #bioconductor PLS-DA   
    P1_LOOCV=matrix(0,ncol = 1, nrow=dim(x)[2])
    P2_LOOCV=matrix(0,ncol = 1, nrow=dim(x)[2])
    P1_LOOCV[!indConstantColumns]=opls.out_LOOCV@loadingMN[,1]
    P2_LOOCV[!indConstantColumns]=opls.out_LOOCV@loadingMN[,2]
    rownames(P1_LOOCV)=colnames(x)
    rownames(P2_LOOCV)=colnames(x)
    
    #test for axis flip from original loadings vector, assign to loading matrix
    Load_LV1_LOOCV[i,]=P1_LOOCV * sign(sum(P1_LOOCV * oplsOut@loadingMN[,1]))
    Load_LV2_LOOCV[i,]=P2_LOOCV * sign(sum(P2_LOOCV * oplsOut@loadingMN[,2]))
    
    # #test for axis swap from original loadings vectors
    # if(sum(Load_LV1_LOOCV[i,]*oplsOut@loadingMN[,1]) < sum(Load_LV2_LOOCV[i,]*oplsOut@loadingMN[,2])){
    #   temp1 <- Load_LV1_LOOCV[i,]
    #   temp2 <- Load_LV2_LOOCV[i,]
    #   Load_LV1_LOOCV[i,] <- temp2
    #   Load_LV2_LOOCV[i,] <- temp1
    # }
  }   #end LOOCV opls iterations
  
  #compute mean and standard deviation for stored loadings over all iterations
  P1_LOOCV=colMeans(Load_LV1_LOOCV)
  P2_LOOCV=colMeans(Load_LV2_LOOCV)
  
  listOut=list()
  listOut$loadings <- oplsOut@loadingMN
  listOut$scores <- oplsOut@scoreMN
  listOut$loadings_LOOCV <- base::cbind(P1_LOOCV,P2_LOOCV)
  listOut$LV1_mat <- Load_LV1_LOOCV
  listOut$LV2_mat <- Load_LV2_LOOCV
  
  return(listOut)
}


# rotate.opls.LOOCV ###########################

# Function rotate.opls.LOOCV
# rotates output from function opls.LOOCV by specified degrees and performs specified axis swap/flips

#Input: 
# opls.LOOCV_out = object output from opls.LOOCV()
# degrees = specified degrees of clockwise rotation
# swap = F, change to TRUE for axis swap (LV1 to LV2, LV2 to LV1)
# flip1 = F, change to TRUE to multiply LV1 by -1 (flip across its axis)
# flip2 = F, change to TRUE to multiply LV2 by -1 (flip across its axis)

#Output:
# A list with the following objects:
# T1 = scores in rotated LV1
# T2 = scores in rotated LV2
# P1 = loadings in rotated LV1
# P2 = loadings in rotated LV2
# P1_LOOCV = mean value of loadings in LV1, calculated across LOOCV iterations
# P2_LOOCV = mean value of loadings in LV2, calculated across LOOCV iterations
# sd_LV1_LOOCV = standard deviation of each loading in LV1, calculated across LOOCV iterations
# sd_LV2_LOOCV = standard deviation of each loading in LV2, calculated across LOOCV iterations
# rotmat = matrix of rotation by specified degrees

rotate.opls.LOOCV <- function(opls.LOOCV_out,degrees,swap=F,flip1=F,flip2=F){
  require("matrixStats")
  #rotate by specified degrees
  theta=degrees*pi/180
  rotmat=rbind(c(cos(theta),-sin(theta)),c(sin(theta), cos(theta)))
  scores = opls.LOOCV_out$scores %*% rotmat
  T1=scores[,1]
  T2=scores[,2]
  loadings = opls.LOOCV_out$loadings %*% rotmat
  P1=loadings[,1]
  P2=loadings[,2]
  
  #rotate matrix of loadings from LOOCV
  P1_mat <- opls.LOOCV_out$LV1_mat
  P2_mat <- opls.LOOCV_out$LV2_mat
  for(i in 1:dim(opls.LOOCV_out$LV1_mat)[1]) {
    load_mat <- base::cbind(opls.LOOCV_out$LV1_mat[i,],opls.LOOCV_out$LV2_mat[i,])
    rot_load_mat <- load_mat %*% rotmat
    P1_mat[i,] <- rot_load_mat[,1]
    P2_mat[i,] <- rot_load_mat[,2]
  }
  
  #swap axis 1 and 2 if specified
  if(swap==TRUE){
    P1temp=P1; P2temp=P2
    P1_mat_temp = P1_mat; P2_mat_temp = P2_mat
    P1=P2temp; P2=P1temp
    P1_mat = P2_mat_temp; P2_mat = P1_mat_temp
    T1temp=T1; T2temp=T2
    T1=T2temp; T2=T1temp
  }
  
  #Flip axes 1 and/or 2 if specified
  if(flip1==TRUE){
    P1=-P1
    P1_mat=-P1_mat
    T1=-T1
  }
  if(flip2==TRUE){
    P2=-P2
    P2_mat=-P2_mat
    T2=-T2
  }
  
  P1_LOOCV=colMeans(P1_mat)
  names(P1_LOOCV)=names(P1)
  P2_LOOCV=colMeans(P2_mat)
  names(P2_LOOCV)=names(P2)
  stdevLV1=colSds(P1_mat) 
  stdevLV2=colSds(P2_mat)
  
  
  out <- list()
  out$T1 <- T1
  out$T2 <- T2
  out$P1 <- as.matrix(P1)
  out$P2 <- as.matrix(P2)
  out$P1_LOOCV <- as.matrix(P1_LOOCV)
  out$P2_LOOCV <- as.matrix(P2_LOOCV)
  out$sd_LV1_LOOCV <- stdevLV1
  out$sd_LV2_LOOCV <- stdevLV2
  out$rotmat <- rotmat
  return(out)
}


# R corr adjust FDR by Levi ################

rcorr.adjustFDR<-function (x, type = c("pearson", "spearman"), use = c("complete.obs", 
                                                                       "pairwise.complete.obs")) 
{
  opt <- options(scipen = 5)
  on.exit(options(opt))
  type <- match.arg(type)
  use <- match.arg(use)
  x <- if (use == "complete.obs") 
    as.matrix(na.omit(x))
  else as.matrix(x)
  R <- rcorr(x, type = type)
  P <- P.unadj <- R$P
  p <- P[lower.tri(P)]
  adj.p <- p.adjust(p, method = "fdr")
  P[lower.tri(P)] <- adj.p
  P[upper.tri(P)] <- 0
  P <- P + t(P)
  #P <- ifelse(P < 1e-04, 0, P)
  #P <- format(round(P, 4))
  diag(P) <- NA
  #P[c(grep("0.0000", P), grep("^ 0$", P))] <- "<.0001"
  #P[grep("0.000$", P)] <- "<.001"
  #P.unadj <- ifelse(P.unadj < 1e-04, 0, P.unadj)
  #P.unadj <- format(round(P.unadj, 4))
  diag(P.unadj) <- NA
  #P.unadj[c(grep("0.0000$", P.unadj), grep("^ 0$", P.unadj))] <- "<.0001"
  #P.unadj[grep("0.000$", P.unadj)] <- "<.001"
  result <- list(R = R, P = P, P.unadj = P.unadj, type = type)
  class(result) <- "rcorr.adjustFDR"
  result
}



# Heatmap_Func by Brae ###################

Heatmap_Func <- function(heatMapdata, title, x_str = "", y_str = "", color1 = "blue", color2 = "white", color3 = "red2", verticalcluster = TRUE, horizontalcluster = TRUE,cex_r=1,cex_c=0.6,mar=c(5,2)) {
  if (!require("gplots")) install.packages("gplots")
  if (!require("heatmap3")) install.packages("heatmap3")
  breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
  barColors = colorpanel(length(breakBarColors)-1, color1, color2, color3)
  if (horizontalcluster & verticalcluster) {
    hc = hclust(dist(t(heatMapdata), method = "euclidean"), method = "average")
    hcv = hclust(dist(heatMapdata, method = "euclidean"), method = "average")
    dataHeatMap = data.frame(heatMapdata)
    ColLabels = colnames(heatMapdata)
    RowLabels = rownames(heatMapdata)
    
    output_hm <- heatmap3(heatMapdata,
                          col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=2.5),
                          Rowv=as.dendrogram(hcv),
                          Colv=as.dendrogram(hc),
                          scale="none", cexCol=cex_c,cexRow=cex_r, margins=mar,   #5,2
                          highlightCell=data.frame(rep(1:dim(heatMapdata)[1],each=(dim(heatMapdata)[2]+1)),rep(1:(dim(heatMapdata)[2]+1),times=dim(heatMapdata)[1]),'black',1),
                          labCol=ColLabels, labRow=RowLabels, main = title, xlab = x_str, ylab = y_str
    )
  } else if (horizontalcluster) {
    hc = hclust(dist(t(heatMapdata), method = "euclidean"), method = "average")
    dataHeatMap = data.frame(heatMapdata)
    ColLabels = colnames(heatMapdata)
    RowLabels = rownames(heatMapdata)
    
    output_hm <- heatmap3(heatMapdata,
                          col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=2.5),
                          Rowv=NA,
                          Colv=as.dendrogram(hc),
                          scale="none", cexCol=cex_c,cexRow=cex_r, margins=mar,   #5,2
                          highlightCell=data.frame(rep(1:dim(heatMapdata)[1],each=(dim(heatMapdata)[2]+1)),rep(1:(dim(heatMapdata)[2]+1),times=dim(heatMapdata)[1]),'black',1),
                          labCol=ColLabels, labRow=RowLabels, main = title, xlab = x_str, ylab = y_str
    )
  } else if (verticalcluster) {
    hcv = hclust(dist(heatMapdata, method = "euclidean"), method = "average")
    dataHeatMap = data.frame(heatMapdata)
    ColLabels = colnames(heatMapdata)
    RowLabels = rownames(heatMapdata)
    
    output_hm <- heatmap3(heatMapdata,
                          col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=2.5),
                          Rowv=as.dendrogram(hcv),
                          Colv=NA,
                          scale="none", cexCol=cex_c,cexRow=cex_r, margins=mar,   #5,2
                          highlightCell=data.frame(rep(1:dim(heatMapdata)[1],each=(dim(heatMapdata)[2]+1)),rep(1:(dim(heatMapdata)[2]+1),times=dim(heatMapdata)[1]),'black',1),
                          labCol=ColLabels, labRow=RowLabels, main = title, xlab = x_str, ylab = y_str
    )
  } else {
    dataHeatMap = data.frame(heatMapdata)
    ColLabels = colnames(heatMapdata)
    RowLabels = rownames(heatMapdata)
    
    output_hm <- heatmap3(heatMapdata,
                          col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=2.5),
                          Rowv=NA,
                          Colv=NA,
                          scale="none", cexCol=cex_c,cexRow=cex_r, margins=mar,   #5,2
                          highlightCell=data.frame(rep(1:dim(heatMapdata)[1],each=(dim(heatMapdata)[2]+1)),rep(1:(dim(heatMapdata)[2]+1),times=dim(heatMapdata)[1]),'black',1),
                          labCol=ColLabels, labRow=RowLabels, main = title, xlab = x_str, ylab = y_str
    )
  }
  return(output_hm)
}

Heatmap_SideCol <- function(hm_mat,SideCol,title="",x.str="",y.str="",color1="blue",color2="white",color3="red2",cex_r=1,cex_c=1,mar=c(5,2)){
  if (!require("gplots")) install.packages("gplots")
  if (!require("heatmap3")) install.packages("heatmap3")
  breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
  barColors = colorpanel(length(breakBarColors)-1, color1, color2, color3)
  hc = hclust(dist(t(hm_mat), method = "euclidean"), method = "average")
  ColLabels = colnames(hm_mat)
  RowLabels = rownames(hm_mat)
  
  output_hm <- heatmap3(hm_mat,
                        col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA),cex=2.5),
                        Rowv=NA,
                        Colv=as.dendrogram(hc),
                        scale="none", cexCol=cex_c,cexRow=cex_r, margins=mar,   #5,2
                        highlightCell=data.frame(rep(1:dim(hm_mat)[1],each=(dim(hm_mat)[2]+1)),rep(1:(dim(hm_mat)[2]+1),times=dim(hm_mat)[1]),'black',1),
                        labCol=ColLabels, labRow=RowLabels, main = title, xlab = x.str, ylab = y.str, RowSideColors = SideCol
  )
}


# Correlation Heatmap ######################
Heatmap_Cor <- function(hmData, title="", x.str = "", y.str = "", 
                        color1 = "black", color2 = "white", color3 = "orange", 
                        cex_r=1,cex_c=0.6,mar=c(12,13)){
  if (!require("gplots")) install.packages("gplots")
  if (!require("heatmap3")) install.packages("heatmap3")
  breakBarColors=c(-200,seq(-1, 1, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
  barColors = colorpanel(length(breakBarColors)-1, color1, color2, color3)
  
  corOut=cor(hmData)
  corOut[is.na(corOut)] = 0
  hCor=hclust(as.dist((1-corOut)/2))
  
  heatmap3(corOut,
           col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
           Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",margins=mar,
           cexCol=cex_c,cexRow=cex_r,
           main = title, xlab = x.str, ylab = y.str)
}

Heatmap_p <- function(hmData, title="", x.str = "", y.str = "", 
                        color1 = "red",color2= "green", color3 = "black", 
                        cex_r=1,cex_c=0.6,mar=c(12,13),
                        loglimit=-5){
  if (!require("gplots")) install.packages("gplots")
  if (!require("heatmap3")) install.packages("heatmap3")
  breakBarColors=c(0,10^seq(loglimit, 0, 0.05)) 
  barColors = colorpanel(length(breakBarColors)-1, color1, color2, color3)
  
  corOut_FDR= rcorr.adjustFDR(hmData)
  
  corOut=cor(hmData)
  corOut[is.na(corOut)] = 0
  hCor=hclust(as.dist((1-corOut)/2))
  
  heatmap3(corOut_FDR$P,
           col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
           Rowv=as.dendrogram(hCor), Colv=as.dendrogram(hCor),  scale="none",margins=mar,
           cexCol=cex_c,cexRow=cex_r,
           main = title, xlab = x.str, ylab = y.str)
}


# MeanSEM_function by Brae #####################


errorplot.meanSE <- function(x, y, colors=x, x.str = "", y.str = "", line_color = 'black', title = y.str) {
  require("ggpubr")
  dataPlot=data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggerrorplot(dataPlot,x="x_dot",y="y_dot",
                         combine=TRUE,add=c("mean","mean_se","dotplot"),
                         add.params=list(fill="colors_dot",size=1),color=line_color,
                         legend="none",error.plot = "errorbar")+    
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_line(color="gray",linetype=2,size=0.2),
          axis.line=element_line(color="black"),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16),
          plot.title=element_text(hjust=0.5))+
    xlab(x.str)+
    ylab(y.str)+
    ggtitle(paste(title))
  return(plotOut)
}

errorplot.meanSD.compare <- function(x, y, colors=x, x.str = "", y.str = "", line_color = 'black', title = y.str, comparisonMethod = "t.test",my_comparisons=list(c(1,2)),p.display="p.signif",paired=F) {
  require("ggpubr")
  dataPlot=data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggerrorplot(dataPlot,x="x_dot",y="y_dot",
                         combine=TRUE,add=c("mean","violin"),
                         add.params=list(fill="colors_dot",size=0.5),color=line_color,
                         legend="none",error.plot = "errorbar",desc_stat = "mean_sd")+    
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_line(color="gray",linetype=2,size=0.2),
          axis.line=element_line(color="black"),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16),
          plot.title=element_text(hjust=0.5))+
    xlab(x.str)+
    ylab(y.str)+
    stat_compare_means(comparisons=my_comparisons, method = comparisonMethod,label=p.display,paired=paired) +
    ggtitle(paste(title))
  return(plotOut)
}

errorplot.meanSE.anova <- function(x, y, colors=x, x.str = "", y.str = "", line_color = 'black', title = y.str,p.display="p.signif") {
  require("ggpubr")
  dataPlot=data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggerrorplot(dataPlot,x="x_dot",y="y_dot",
                         combine=TRUE,add=c("mean","mean_se","dotplot"),
                         add.params=list(fill="colors_dot",size=1),color=line_color,
                         legend="none",error.plot = "errorbar")+    
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_line(color="gray",linetype=2,size=0.2),
          axis.line=element_line(color="black"),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16),
          plot.title=element_text(hjust=0.5))+
    xlab(x.str)+
    ylab(y.str)+
    stat_compare_means(method = "anova",label=p.display) +
    ggtitle(paste(title))
  return(plotOut)
}

#Used to create dot plots with mean +/- SEM bars, saved to working directory. No output.
dotplot.mean.sem <- function(x,y,colors=x,x.str="",y.str="",my_comparisons=list(c(1,2)),title="",line_color="black"){
  if (!require("ggpubr")) install.packages("ggpubr")
  library("ggpubr")
  dataPlot=data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggline(dataPlot,x="x_dot",y="y_dot",combine=TRUE,add=c("mean_se","dotplot"),add.params=list(fill="colors_dot",size=1.2),color=line_color)+
    theme(panel.background = element_blank(),
          panel.grid.major=element_line(color='gray',linetype=2,size=0.2),
          plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          axis.line=element_line(color='black'),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16))+
    xlab(x.str)+
    ylab(y.str)+
    stat_compare_means(comparisons = my_comparisons,method="t.test") +
    ggtitle(title)
  return(plotOut)
}

dotplot.mean.sem.gradient <- function(x,y,colors=x,x.str,y.str,my_comparisons=list(c(1,2)),title=""){
  if (!require("ggpubr")) install.packages("ggpubr")
  library("ggpubr")
  dataPlot=data.frame(factor(x),y,colors)
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggline(dataPlot,x="x_dot",y="y_dot",color="colors_dot",combine=TRUE,
                    add=c("mean_se","dotplot"),position = position_dodge(),
                    size=2)+
    theme(panel.background = element_blank(),
          panel.grid.major=element_line(color='gray',linetype=2,size=0.2),
          plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          axis.line=element_line(color='black'),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16))+
    xlab(x.str)+
    ylab(y.str)+
    stat_compare_means(comparisons = my_comparisons,method="t.test") +
    gradient_color(c("blue","white","red"))
    ggtitle(title)
  return(plotOut)
}

#Used to create dot plots with mean +/- SEM bars, saved to working directory. No output.
dotplot.mean.sem.facet <- function(x,y,colors=x,facet,x.str="",y.str="",my_comparisons=list(c(1,2)),title=""){
  if (!require("ggpubr")) install.packages("ggpubr")
  library("ggpubr")
  dataPlot=data.frame(factor(x),y,factor(colors),facet)
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot","facet")
  plotOut <- ggline(dataPlot,x="x_dot",y="y_dot",combine=TRUE,add=c("mean_se","dotplot"),add.params=list(fill="colors_dot",size=1.2),color="black",palette=cbPalette,facet.by = "facet")+
    theme(panel.background = element_blank(),
          panel.grid.major=element_line(color='gray',linetype=2,size=0.2),
          plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.line=element_line(color='black'),
          text = element_text(size=20),
          axis.text.x = element_text(face = "bold", color = "black",
                                     size =10,angle=45,vjust=1,hjust=1),
          axis.text.y = element_text(face = "bold", color = "black",
                                     size = 16))+
    xlab(x.str)+
    ylab(y.str)+
    stat_compare_means(method = "anova") +
    ggtitle(title)
  return(plotOut)
}



# scatter_reg by Brae, edited AFP 1/24/21 ########



scatter_reg <- function(x,y,x.str="",y.str="",title="",point_color="black",line_color="red", 
                        showRSQ = TRUE, showPVAL = TRUE, 
                        line_size = 1.5, point_size = 4, line_format = "longdash", axisfontsize = 14, 
                        titlefontsize = 20,
                        text_position = "top right") {
  dataReg=data.frame(x,y)
  colnames(dataReg) <- c("x","y")
  linmodel=lm(y~x,dataReg)
  RSQ=summary(linmodel)$r.squared
  f <- summary(linmodel)$fstatistic
  pval <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(pval) <- NULL
  
  if (text_position == "top right"){p_x=0.8; p_y=1;r_x=0.8; r_y=0.9}
  if (text_position == "top left"){p_x=0.2; p_y=1;r_x=0.2; r_y=0.9}
  if (text_position == "bottom right"){p_x=0.8; p_y=0.15;r_x=0.8; r_y=0.05}
  if (text_position == "bottom left"){p_x=0.2; p_y=0.15;r_x=0.2; r_y=0.05}

  
  if(showRSQ == FALSE){if(r_y==0.05){p_y=r_y}; r_x = NaN}
  if(showPVAL== FALSE){if(p_y==1){r_y=p_y};p_x = NaN}
  
  output_plot <- ggplot(dataReg, aes(x=x, y=y)) +
    geom_smooth(method='lm',se=FALSE,color=line_color,linetype=line_format, size = line_size)+
    geom_point(size=point_size,color=point_color)+
    theme(panel.background = element_rect(fill = 'transparent'),
          text = element_text(size=20),
          axis.text.x = element_text(color = "black",
                                     size = axisfontsize),
          axis.text.y = element_text(color = "black",
                                     size = axisfontsize),
          plot.title = element_text(hjust = 0.5,face='bold',size=titlefontsize),
          panel.border=element_blank(),
          axis.line=element_line(color='black'),
          plot.margin=margin(10,20,10,10)
    )+
    xlab(x.str)+
    ylab(y.str)+
    annotate("text",size=6, x={r_x*max(x,na.rm=T)+(1-r_x)*min(x,na.rm=T)}, y={r_y*max(y,na.rm=T)+(1-r_y)*min(y,na.rm=T)}, label=bquote(R^2==.(formatC(RSQ, digits=2))))+
    annotate("text",size=6, x=p_x*max(x,na.rm=T)+(1-p_x)*min(x,na.rm=T), y=p_y*max(y,na.rm=T)+(1-p_y)*min(y,na.rm=T), label= bquote(p==.(formatC(pval, digits=2))))+
    ggtitle(title)
  
  return(output_plot)
}


# scatter_reg_group

scatter_reg_group <- function(x,y,grouping,x.str="",y.str="",title="",
                              colors=c("black","red"), 
                        showRSQ = TRUE, showPVAL = TRUE, 
                        line_size = 1.5, point_size = 4, 
                        line_format = "longdash", axisfontsize = 14, 
                        text_position = "top right") {
  dataReg=data.frame(x,y,factor(grouping))
  colnames(dataReg) <- c("x","y","grouping")
  
  # for(i in 1:length(unique(dataReg$grouping))){
  #   linmodel=lm(y~x,dataReg,which(grouping==unique(grouping)[i]))
  #   RSQ=summary(linmodel)$r.squared
  #   f <- summary(linmodel)$fstatistic
  #   pval <- pf(f[1],f[2],f[3],lower.tail=F)
  #   attributes(pval) <- NULL
  # }

  
  if (text_position == "top right"){p_x=0.8; p_y=1;r_x=0.8; r_y=0.9}
  if (text_position == "top left"){p_x=0.2; p_y=1;r_x=0.2; r_y=0.9}
  if (text_position == "bottom right"){p_x=0.8; p_y=0.15;r_x=0.8; r_y=0.05}
  if (text_position == "bottom left"){p_x=0.2; p_y=0.15;r_x=0.2; r_y=0.05}
  
  if(showRSQ == FALSE){if(r_y==0.05){p_y=r_y}; r_x = NaN}
  if(showPVAL== FALSE){if(p_y==1){r_y=p_y};p_x = NaN}
  
  output_plot <- ggplot(dataReg, aes(x=x, y=y, group=grouping, color=grouping)) +
    geom_smooth(method='lm',se=FALSE,size = line_size,linetype=line_format)+
    geom_point(size=point_size)+
    scale_color_manual(values=colors)+
    theme(panel.background = element_rect(fill = 'transparent'),
          text = element_text(size=20),
          axis.text.x = element_text(color = "black",
                                     size = axisfontsize),
          axis.text.y = element_text(color = "black",
                                     size = axisfontsize),
          plot.title = element_text(hjust = 0.5,face='bold'),
          panel.border=element_blank(),
          axis.line=element_line(color='black'),
          plot.margin=margin(10,20,10,10)
    )+
    xlab(x.str)+
    ylab(y.str)+
    annotate("text",size=6, x={r_x*max(x,na.rm=T)+(1-r_x)*min(x,na.rm=T)}, y={r_y*max(y,na.rm=T)+(1-r_y)*min(y,na.rm=T)}, label=bquote(R^2==.(formatC(RSQ, digits=2))))+
    annotate("text",size=6, x=p_x*max(x,na.rm=T)+(1-p_x)*min(x,na.rm=T), y=p_y*max(y,na.rm=T)+(1-p_y)*min(y,na.rm=T), label= bquote(p==.(formatC(pval, digits=2))))+
    ggtitle(title)
  
  return(output_plot)
}

