#----Preliminaries----
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
###if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

  if (!require("pacman")) install.packages("pacman")
  if (!require("BiocManager")) install.packages("BiocManager")
  pacman::p_load(readxl,tidyverse,ropls)
  
  library('fastcluster')
  library("gplots")
  library("dendextend")
  library(reshape)
  library(matrixStats)
  library("readxl")
  library("heatmap3")
  library("compare")
  library("grid")
  library("gridGraphics")
  library("ggplot2")
  library("scatr")
  library("ggpubr")
  library("RColorBrewer")
  library(boot)
  library(caret)
  library(hrbrthemes)
  library("ClassDiscovery")
  library("plsdof")
  library("gplots")
  library("lme4")
  library("foreign")
  library("corrplot")
  library("Hmisc")
  library("PerformanceAnalytics")
  library("dplyr")
  library("caret")
  library("ggpubr")
  library("ggpmisc")
  library("RColorBrewer")
  library("gplots")
  library("SAMtool")
  library("impute")
  library("latex2exp")
  library("base")
  library("RColorBrewer")
  library("tibble")
  source("WoodLabFunctions.R")


#source("summarySE.R")

#Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) #Outside numbers clip outliers. This is for zscoring.
# barColors = colorpanel(length(breakBarColors)-1, "#58539f", "white", "#b83732")
barColors = colorpanel(length(breakBarColors)-1, "#264385", "white", "#b83732")
barColors = colorpanel(length(breakBarColors)-1, "darkgrey", "white", "mediumorchid")

pdfPath=paste("Scale_Bar_ZScore.pdf", sep="")
pdf(pdfPath, width=2,height=5,pointsize = 14)
color.bar(barColors,-1.5,1.5,nticks=5)
dev.off()


#----Read data----
df_in = data.frame(read_excel("Monomer and Oligomer WT ABeta.xlsx", sheet="PV 1.6 Copy"))
firstprotein = "Eotaxin"
lastprotein = "TNFa"
Analytes <- colnames(df_in)[{which(colnames(df_in)==firstprotein)}:{which(colnames(df_in)==lastprotein)}]
# Analytes = colnames(df_in )[2:9]

##if (filter_cyt){df_in <- select(df_in,Condition,Sample,all_of(filtered_cyts))}

analyte_range=which(colnames(df_in) %in% Analytes)

#create array of average background for each analyte
df_b = df_in %>%
  dplyr::filter(Sample == "Background") %>%
  dplyr::select(all_of(analyte_range))
df_b_mean = colMeans(df_b)


#create data frame of sample data
df_samples <- df_in %>%
  dplyr::filter(Sample != "Background")

#subtract background from analytes
df_samples[,analyte_range] <- sweep(df_samples[,analyte_range],2,df_b_mean)

#set negative values to 0
df_samples[df_samples<0] <- 0
# 
# df_samples = df_in %>% filter(Sample != "Timed")
df_samples <- df_samples %>%
  filter(Sample %in% c(
    "Vehicle", "WT Monomer",
    "WT Oligomer" 
    ))

df_samples <- df_samples %>%
  mutate(Sample = factor(Sample, levels = c(
    "Vehicle", "WT Monomer",
    "WT Oligomer" 
  ))) %>%
  arrange(Sample)


# Get data for Run 3
# Load new data
df_in2 <- data.frame(read_excel("Monomer and Oligomer WT ABeta.xlsx", sheet="Tarunya Data"))

# Define first and last analytes
firstprotein2 <- "Eotaxin"
lastprotein2  <- "TNFa"
Analytes2 <- colnames(df_in2)[
  which(colnames(df_in2) == firstprotein2):which(colnames(df_in2) == lastprotein2)
]

# Get analyte column indices
analyte_range2 <- which(colnames(df_in2) %in% Analytes2)

# Calculate average background values
df_b2 <- df_in2 %>%
  dplyr::filter(Sample == "Background") %>%
  dplyr::select(all_of(analyte_range2))
df_b_mean2 <- colMeans(df_b2)

# Extract non-background samples
df_samples2 <- df_in2 %>%
  dplyr::filter(Sample != "Background")

# Subtract background
df_samples2[, analyte_range2] <- sweep(df_samples2[, analyte_range2], 2, df_b_mean2)

# Set negative values to zero
df_samples2[df_samples2 < 0] <- 0

# Filter and arrange selected samples
df_samples2 <- df_samples2 %>%
  filter(Sample %in% c(    "Vehicle", "WT Monomer",
                           "WT Oligomer" 
  )) %>%
  mutate(Sample = factor(Sample, levels = c(    "Vehicle", "WT Monomer",
                                                "WT Oligomer" 
  ))) %>%
  arrange(Sample)


# # Get CoStim data
# df_in3 <- data.frame(read_excel("Monomer and Oligomer WT ABeta.xlsx", sheet="PV 1.1"))
# 
# firstprotein3 <- "Eotaxin"
# lastprotein3  <- "TNFa"
# Analytes3 <- colnames(df_in3)[
#   which(colnames(df_in3) == firstprotein3):which(colnames(df_in3) == lastprotein3)
# ]
# analyte_range3 <- which(colnames(df_in3) %in% Analytes3)
# 
# # Calculate background mean
# df_b3 <- df_in3 %>%
#   dplyr::filter(Sample == "Background0") %>%
#   dplyr::select(all_of(analyte_range3))
# df_b_mean3 <- colMeans(df_b3)
# 
# # Sample data
# df_samples3 <- df_in3 %>%
#   dplyr::filter(Sample != "Background0")
# df_samples3[, analyte_range3] <- sweep(df_samples3[, analyte_range3], 2, df_b_mean3)
# df_samples3[df_samples3 < 0] <- 0
# 
# # Filter and format
# df_samples3 <- df_samples3 %>%
#   filter(Sample %in% c("Vehicle", "WT Monomer",
#                        "WT Oligomer" )) %>%
#   mutate(Sample = factor(Sample, levels = c(
#     "Vehicle", "WT Monomer",
#     "WT Oligomer" 
#   ))) %>%
#   arrange(Sample)
# 


# df_Z3 = df_samples3
# df_Z3[,analyte_range] = apply(df_samples3[,analyte_range],2,"scale")
# df_Z3[is.na(df_Z3)] = 0
# 
# df_Z3 = df_Z3 %>%
#   arrange(Sample)
# 
# df_Mat3 <- df_Z3[,analyte_range] %>%
#   as.matrix()
# 
# rownames(df_Mat3) = df_samples3$Sample



normalize_to_vehicle <- function(df_plate, analyte_range) {
  # Filter for Vehicle rows
  vehicle_rows <- df_plate %>% filter(Sample == "Vehicle")
  
  # Check if there are any Vehicle samples
  if (nrow(vehicle_rows) == 0) {
    warning("No Vehicle samples found in this plate.")
    return(df_plate)
  }
  
  # Compute mean across Vehicle samples per analyte
  vehicle_means <- colMeans(vehicle_rows[, analyte_range], na.rm = TRUE)
  
  # Replace 0 or NA in vehicle_means with 1 to avoid division by zero/NA
  vehicle_means[is.na(vehicle_means) | vehicle_means == 0] <- 1
  
  # Normalize
  df_normalized <- df_plate
  df_normalized[, analyte_range] <- sweep(
    df_plate[, analyte_range],
    2,
    vehicle_means,
    FUN = "/"
  )
  
  return(df_normalized)
}

# Normalize each plate separately
df_samples_norm1 <- normalize_to_vehicle(df_samples, analyte_range)
df_samples_norm2 <- normalize_to_vehicle(df_samples2, analyte_range2)
# df_samples_norm3 <- normalize_to_vehicle(df_samples3, analyte_range3)

# (Optional: ensure column alignment before merging)
df_samples_norm2 <- df_samples_norm2[, colnames(df_samples_norm1)]
# df_samples_norm3 <- df_samples_norm3[, colnames(df_samples_norm1)]

# Merge all normalized data
df_samplesc <- bind_rows(df_samples_norm1, df_samples_norm2)

# Arrange by desired sample order
df_samplesc <- df_samplesc %>%
  mutate(Sample = factor(Sample, levels = c("Vehicle", "WT Monomer", "WT Oligomer"))) %>%
  arrange(Sample)



df_Mat <- as.matrix(df_samplesc[, analyte_range])
df_Z = apply(df_Mat, 2, scale)
rownames(df_Z) = df_samplesc$Sample





# # 3. Batch correction
# Batch <- as.factor(df_samplesc$Run)
# df_Mat_batch_corrected <- removeBatchEffect(t(df_Mat), batch = Batch)
# df_Mat_batch_corrected <- t(df_Mat_batch_corrected)

# #z=score
# df_Z = apply(df_Mat_batch_corrected, 2, scale)
# rownames(df_Z) = df_samplesc$Sample


# Combine
# df_samplesc2 = rbind(df_samplesc, df_samples3)
# 
# df_Matboth = rbind(df_Z,df_Mat3)
# 
# #Remove IL-10
# df_samplesc2 = rbind(df_samplesc, df_samples3)
# 
# df_Matboth = rbind(df_Z,df_Mat3)
# 
# 
# df_Matboth_df <- as.data.frame(df_Matboth)  
# 
# 
# df_Matboth_df$Sample = factor(rownames(df_Matboth), levels = 
#                                 c("Vehicle","ASN","N27D","N27D CoStim", "Q15E","Q15E CoStim", "WT", "WT CoStim"))
# 
# df_Matboth_df <- df_Matboth_df %>%
#   arrange(Sample)
# 
# df_MatBoth_df2 = as.matrix(df_Matboth_df[,-9])
# 
# # df_Zboth$Sample = df_samplesc2$Sample
# 
# # df_Zboth <- as.data.frame(df_Zboth) %>%
# #   tibble::rownames_to_column(var = "Sample") %>%
# #   filter(Sample %in% c("Vehicle", "N27D", "N27D CoStim", "Q15E", "Q15E CoStim", "WT", "WT CoStim")) %>%
# #   mutate(Sample = factor(Sample, levels = c("Vehicle", "N27D", "N27D CoStim", "Q15E", "Q15E CoStim", "WT", "WT CoStim"))) %>%
# #   arrange(Sample) %>%
# #   column_to_rownames(var = "Sample")
# 
# 
# 
# # 
# df_samplesc2 = rbind(df_samplesc, df_samples3)
# df_samplesc2 <- df_samplesc2 %>%
#   filter(Sample %in% c("Vehicle","ASN","N27D","N27D CoStim", "Q15E","Q15E CoStim", "WT", "WT CoStim")) %>%
#   mutate(Sample = factor(Sample, levels = c("Vehicle","ASN","N27D","N27D CoStim", "Q15E","Q15E CoStim", "WT", "WT CoStim"))) %>%
#   arrange(Sample)
# # 
# # 
co=cor(df_Z, method="pearson")
hc=hcludf_Zhc=hclust(as.dist(1-co), method="complete")

ConditionColorBar <- case_when(
  df_samplesc$Sample == "Vehicle"     ~ brewer.pal(9, "Oranges")[7],
  df_samplesc$Sample == "WT Monomer"  ~ "#1E90FF",   # brighter non-palette blue
  df_samplesc$Sample == "WT Oligomer" ~ brewer.pal(9, "Blues")[9]
)



pdfPath=paste("./Heatmap.pdf", sep="")
pdf(pdfPath, width=20,height=14,pointsize = 20)

heatmap3(df_Z,
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, Colv=NA,  scale="none",
         labRow = NA, RowSideColors = ConditionColorBar, RowSideLabs = c("Conditions"),
         labCol = colnames(df_Z), main="", cexRow=1, cexCol = 1.2) 

dev.off()



## Making Heatmaps
pdfPath <- "./Legends_Conditions.pdf"
pdf(pdfPath, width = 8, height = 6, pointsize = 10)

# Create a blank plotting area
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', 
     xlim = c(0, 1), ylim = c(0, 1))

# Updated legend labels (must match Sample column)
legend_labels <- c(
  "Vehicle", 
  "WT Monomer", "WT Oligomer"
)

# Updated colors for the legend
legend_colors <- c(
  "Vehicle"     = brewer.pal(9, "Oranges")[7],
  "WT Monomer"  = "#1E90FF",
  "WT Oligomer" = brewer.pal(9, "Blues")[9]
)

# Draw the legend
legend("topleft", legend = legend_labels, 
       pch = 16, pt.cex = 3, cex = 1.5, bty = 'n',
       col = legend_colors[legend_labels])

# Add main title
mtext("Conditions", at = 0.2, cex = 2)

# Close the PDF device
dev.off()


#### easier way to make dotplots 
library(ggbeeswarm)
dataIn = df_samplesc
dataIn <- rbind(df_samples, df_samples2)
dataIn$Type <- ifelse(dataIn$Sample == "Vehicle", "Vehicle",
                      ifelse(grepl("Monomer", dataIn$Sample), "Monomer",
                             ifelse(grepl("Oligomer", dataIn$Sample), "Oligomer", NA)))


Analytes <- colnames(df_samplesc)[{which(colnames(df_samplesc)==firstprotein)}:{which(colnames(df_samplesc)==lastprotein)}]


for (i in Analytes) {
  dataBar <- data.frame("Group" = dataIn$Sample, "Values" = dataIn[[i]])
  
  # Corrected PDF output
  output <- paste0(i, " Normalized.pdf")
  pdf(file = output, width = 7, height = 8.5)  # Replaced `quartz()` with `pdf()`
  print({
    e <- ggplot(data = dataBar, aes(x = Group, y = Values, color = Group, shape = Group)) +
      geom_beeswarm(aes(fill = Group, shape = Group), size = 4, stroke = 1, cex = 2,
                    groupOnX = TRUE, dodge.width = 0) +
      stat_summary(geom = "crossbar", fun = mean,
                   position = position_dodge(0.7), colour = "darkblue",
                   linewidth = 0.4, width = 0.6, show.legend = FALSE) +
      stat_summary(geom = "errorbar", color = "black",
                   position = position_dodge(0.7), linewidth = 0.8, 
                   width = 0.4, fun.data = mean_se) +
      stat_compare_means(exact = FALSE, method = "wilcox",
                         comparisons = list(c(1,2), c(3,1), c(2,3)),
                         symnum.args = list(), 
                         p.adjust.method = "bonferroni", label = "p.value",
                         bracket.size = 0.85, tip.length = 0,
                         # label.y = c(0.95 * max(dataIn[[i]]), 
                         #             1.08 * max(dataIn[[i]]), 
                         #             1.2 * max(dataIn[[i]]), 
                         #             1.32 * max(dataIn[[i]])),
                         size = 10, face = "bold.italic", ref.group = ".all.") +
      scale_color_manual(values = c(
        "Vehicle"     = brewer.pal(9, "Oranges")[7],
        "WT Monomer"  = "#1E90FF",
        "WT Oligomer" = brewer.pal(9, "Blues")[9]
      )) +
      scale_fill_manual(values = c(
        "Vehicle"     = brewer.pal(9, "Oranges")[7],
        "WT Monomer"  = "#1E90FF",
        "WT Oligomer" = brewer.pal(9, "Blues")[9]
      )) +      scale_shape_manual(values = c(
        "Vehicle"         = 22,  # Square
        "WT Monomer"      = 24,
        "WT Oligomer"     = 21
      )) +
      scale_x_discrete(labels = scales::wrap_format(5)) +
      scale_y_continuous(limits = c(0, 1.6 * max(dataIn[[i]])), expand = c(0, 0)) +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_text(size = 30, color = "black", face = "bold.italic"),
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.line = element_line(colour = 'black', linewidth = 0.65, linetype = 'solid'),
        text = element_text(size = 30),
        axis.text.x = element_text(face = "bold", color = "black", vjust = 0.6, size = 35,angle = 45),
        axis.text.y = element_text(face = "bold", color = "black", size = 25),
        axis.title.y = element_text(size = 35, margin = margin(t = 0, r = 14, b = 0, l = 0)),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = "none", 
        legend.spacing.x = unit(5, "pt")
      ) +
      xlab("") +
      ylab(paste("Fluorescent Intensity [a.u.]")) +
      ggtitle(i)
  })
  
  dev.off()  # Close PDF device
}



#----Discriminant Partial least squares----

indConstantColumns=GMCM:::colSds(df_Z)==0   

theta=0

# if(ageVal==1){theta=0}
# if(ageVal==2){theta=45}
# if(ageVal==6){theta=0}

PCA=ropls::opls(x=df_Z[,!indConstantColumns], predI=2, crossvalI = 5)
PCAOut=rotate.opls(PCA, degrees=theta)

stdevLV1=colSds(PCAOut$P1) 
stdevLV2=colSds(PCAOut$P2)


scores.plot(T1=PCAOut$T1,T2=PCAOut$T2,color=df_samplesc$Sample)


# predictor = factor(df_samples$Abeta)
# 
# plsOut = opls.LOOCV(x=df_Mat, y = predictor, n = 1, runs =100)
# 
# 
# plsOut=opls.LOOCV(x=df_Mat[,!indConstantColumns],y=predictor,n=1,runs=100)
# 
# 
# plsOut=rotate.opls.LOOCV(plsOut,degrees=theta,swap=F,flip1=F,flip2=F)
# scores.plot(T1=plsOut$T1,T2=plsOut$T2,color=df_samples$Sample)
# 
# #loading chart PCA
# load1_LOOCV = loadings.chart.sd(loadings=plsOut$P1, 
#                                 sd = analyte_sd_df$SD,
#                                 component_str = "LV1",
#                                 lim = 2.5)
# 
# load1_LOOCV = loadings.chart.sd(loadings=plsOut$P1, 
#                                 sd = plsOut$sd_LV1_LOOCV,
#                                 component_str = "LV1",
#                                 lim = 2.5)

##-------Initialize the loadings out----------
Data_P1=matrix(0,ncol = 1, nrow=dim(df_Z)[2])
Data_P2=matrix(0,ncol = 1, nrow=dim(df_Z)[2])
Data_P1Std=matrix(0,ncol = 1, nrow=dim(df_Z)[2])
Data_P2Std=matrix(0,ncol = 1, nrow=dim(df_Z)[2])
rownames(Data_P1)=colnames(df_Z); rownames(Data_P2)=colnames(df_Z)
rownames(Data_P1Std)=colnames(df_Z); rownames(Data_P2Std)=colnames(df_Z)


T1=plsOut$T1
T2=plsOut$T2
Data_P1[!indConstantColumns]=plsOut$P1_LOOCV
Data_P2[!indConstantColumns]=plsOut$P2_LOOCV
Data_P1Std[!indConstantColumns]=as.matrix(plsOut$sd_LV1_LOOCV); rownames(Data_P1Std)=rownames(Data_P1)
Data_P2Std[!indConstantColumns]=plsOut$sd_LV2_LOOCV; rownames(Data_P2Std)=rownames(Data_P2)


#Bar plot the loadings and heatmaps Sorted the same way
# library(latex2exp)
# rownames(Data_P1)=colnames(df_Mat)
# indP1=sort(Data_P1, index.return=TRUE)$ix
# 
# Data_P1Sort=as.matrix(Data_P1[indP1,1])/max(abs(Data_P1))
# Data_P1SortStd=as.matrix(Data_P1Std[indP1,1])/max(abs(Data_P1))
# 
# pdf("./LV1Sort_PLSR_MonovsOligo.pdf", width=12,height=6, pointsize = 18)
# 
# BP=barplot(t(Data_P1Sort),
#            main = "",
#            xlab = "",
#            ylab = "",
#            names.arg = TeX(rownames(Data_P1Sort)),
#            col = "#A7C0DE",
#            horiz = FALSE,
#            las=2,
#            ylim=c(-1.5,0.7))
# arrows(BP,Data_P1Sort+Data_P1SortStd, BP, Data_P1Sort-Data_P1SortStd, angle=90, code=3, length=0.04,lwd=3)
# dev.off()


# for PCA
T1=PCAOut$T1
T2=PCAOut$T2
Data_P1[!indConstantColumns]=PCAOut$P1
Data_P2[!indConstantColumns]=PCAOut$P2

rownames(Data_P1)=colnames(df_Z)
indP1=sort(Data_P1, index.return=TRUE)$ix

Data_P1Sort=as.matrix(Data_P1[indP1,1])/max(abs(Data_P1))
Data_P1SortStd=as.matrix(Data_P1Std[indP1,1])/max(abs(Data_P1))

pdf("./PCA_PC1sort.pdf", width=12,height=6, pointsize = 18)

BP=barplot(t(Data_P1Sort),
           main = "",
           xlab = "",
           ylab = "",
           names.arg = TeX(rownames(Data_P1Sort)),
           col = "#A7C0DE",
           horiz = FALSE,
           las=2,
           ylim=c(-1.5,1.5))
# arrows(BP,Data_P1Sort+Data_P1SortStd, BP, Data_P1Sort-Data_P1SortStd, angle=90, code=3, length=0.04,lwd=3)
dev.off()



pdf("./PCA_PC1sort shading.pdf", width = 12, height = 6, pointsize = 18)

k <- nrow(t(Data_P1Sort))
dens <- rep(c(15, 25, 35, 20), length.out = k)
angs <- rep(c(45, 135, 0, 90), length.out = k)

BP <- barplot(
  t(Data_P1Sort),
  beside = TRUE,
  names.arg = TeX(rownames(Data_P1Sort)),
  col     = rep("black", k),
  density = dens,
  angle   = angs,
  border  = "black",
  las     = 2,
  ylim    = c(-1.5, 1.5)
)

dev.off()

# 
# #Scores Plot fancy version
# #PLSR
# blah3 <- data.frame(
#   T1 = T1,
#   T2 = T2,
#   Sample = factor(df_samplesc$Sample, 
#                   levels = c("WT Monomer", "WT Oligomer", "N27D Monomer", "N27D Oligomer", "Vehicle"))) 
# 
# 
# # Rename columns for clarity
# colnames(blah3) <- c("T1", "T2",  "Sample")
# 
# pdf("./PLSR_Scores_MonovsOligo.pdf", width=6,height=4, pointsize = 16)
# 
# ggplot(blah3, aes(x=T1, y=T2, color=Sample)) + 
#   geom_point(size=4, shape = 24, color = "black", stroke = 1,aes(fill = Sample))+
#   # geom_point(size=5, shape = 17,aes(color = Sample))+
#   xlab("Scores on LV1")+
#   ylab("Scores on LV2")+
#   ggtitle("")+
#   xlim(-6,6)+
#   ylim(-6,6)+
#   geom_vline(xintercept=0)+
#   geom_hline(yintercept=0)+
#   theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
#         text = element_text(size=20),
#         axis.text.x = element_text(face = "bold", color = "black", 
#                                    size = 18),
#         axis.text.y = element_text(face = "bold", color = "black", 
#                                    size = 16))+
#   scale_color_manual(values = c("#E08214", "#FDB863", "#B2ABD2", "#8073AC", "#542788")) +
#   scale_fill_manual(values = c("#E08214", "#FDB863", "#B2ABD2", "#8073AC", "#542788")) +
#   labs(color = "")
# 
# dev.off()


# scores fancier version 
#PCA
T1 = PCAOut$T1
T2 = PCAOut$T2
Data_P1 = PCAOut$P1
Data_P2 = PCAOut$P2

blah3 <- data.frame(
  T1 = T1,
  T2 = T2,
  Sample = factor(df_samplesc$Sample, levels = c(
    "Vehicle", "WT Monomer", "WT Oligomer"
  ))
)


# Rename columns for clarity
colnames(blah3) <- c("T1", "T2", "Sample")

blah3$Type <- ifelse(blah3$Sample == "Vehicle", "Vehicle",
                     ifelse(grepl("Monomer", blah3$Sample), "Monomer",
                            ifelse(grepl("Oligomer", blah3$Sample), "Oligomer", NA)))

colnames(blah3) <- c("T1", "T2", "Sample", "Type")
# Open PDF device
pdf("./PCA_Scores.pdf", width = 6, height = 4, pointsize = 16)

# Plot
ggplot(blah3, aes(x = T1, y = T2, fill = Sample, shape = Sample)) + 
  geom_point(size = 5, stroke = 1, aes(fill = Sample, shape = Sample)) +
  xlab("Scores on PC1") +
  ylab("Scores on PC2") +
  ggtitle("") +
  xlim(-5, 5) +
  ylim(-5, 5) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    text = element_text(size = 20),
    axis.text.x = element_text(face = "bold", color = "black", size = 18),
    axis.text.y = element_text(face = "bold", color = "black", size = 16)
  ) +
  scale_color_manual(values = c(
    "Vehicle"     = brewer.pal(9, "Oranges")[7],
    "WT Monomer"  = "#1E90FF",
    "WT Oligomer" = brewer.pal(9, "Blues")[9]
  )) +
  scale_fill_manual(values = c(
    "Vehicle"     = brewer.pal(9, "Oranges")[7],
    "WT Monomer"  = "#1E90FF",
    "WT Oligomer" = brewer.pal(9, "Blues")[9]
  )) +
  scale_shape_manual(values = c(
    "Vehicle"         = 22,  # Square
    "WT Monomer"      = 24,
    "WT Oligomer"     = 21
  )) +
  labs(color = "", shape = "", fill = "")
# Close PDF
dev.off()


# PCA scre dot plot
pdf("./PCA_PC1_BarPlot.pdf", width = 7, height = 7, pointsize = 16)
  ggplot(blah3, aes(x = Sample, y = T1)) +
  stat_summary(geom = "bar",aes(fill = Sample), fun = mean,position = position_dodge(0.7),
               linewidth = 0.4, width = 0.5,show.legend = FALSE, colour = "black",
  )+# plot mean(column)
  geom_beeswarm(size = 4, aes(shape = Sample), fill = "black", stroke = 1, cex = 2, groupOnX = TRUE, dodge.width = 0) +
  stat_summary(geom = "errorbar", color = "black",
               position = position_dodge(0.7), linewidth = 0.8,
               width = 0.4, fun.data = mean_se) +
  stat_compare_means(exact = FALSE, method = "wilcox",
                     comparisons = list(c(1,2), c(3,1), c(2,3)),
                     symnum.args = list(), 
                     p.adjust.method = "bonferroni", label = "p.value",
                     bracket.size = 0.85, tip.length = 0,
                     # label.y = c(0.95 * max(dataIn[[i]]), 
                     #             1.08 * max(dataIn[[i]]), 
                     #             1.2 * max(dataIn[[i]]), 
                     #             1.32 * max(dataIn[[i]])),
                     size = 9, face = "bold.italic", ref.group = ".all.") +
    scale_color_manual(values = c(
      "Vehicle"     = brewer.pal(9, "Oranges")[7],
      "WT Monomer"  = "#1E90FF",
      "WT Oligomer" = brewer.pal(9, "Blues")[9]
    )) +
    scale_fill_manual(values = c(
      "Vehicle"     = brewer.pal(9, "Oranges")[7],
      "WT Monomer"  = "#1E90FF",
      "WT Oligomer" = brewer.pal(9, "Blues")[9]
    )) +
    scale_shape_manual(values = c(
    "Vehicle"         = 22,  # Square
    "WT Monomer"      = 24,
    "WT Oligomer"     = 21
  )) +
  scale_x_discrete(labels = scales::wrap_format(5)) +
  # scale_y_continuous(limits = c(0, 10)) +
  ylab("PC1 Score") +
  xlab("") +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'white'),
    strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"),
    strip.background = element_rect(colour = "white", fill = "white"),
    axis.line = element_line(colour = 'black', linewidth = 0.65, linetype = 'solid'),
    text = element_text(size = 25),
    axis.text.x = element_text(face = "bold", color = "black", vjust = 0.6, size = 27),
    axis.text.y = element_text(face = "bold", color = "black", size = 25),
    axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 14, b = 0, l = 0)),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "none", 
    legend.spacing.x = unit(5, "pt")
  )
  
dev.off()



#PC2
# PCA scre dot plot
pdf("./PCA_PC2_BarPlot.pdf", width = 7, height = 7, pointsize = 16)

ggplot(blah3, aes(x = Sample, y = T2, color = Sample)) +
  stat_summary(geom = "bar",aes(fill = Sample), fun = mean,position = position_dodge(0.7),
               linewidth = 0.4, width = 0.5,show.legend = FALSE, colour = "black",
  )+# plot mean(column)
  geom_beeswarm(size = 4, aes(shape = Sample), fill = "black", stroke = 1, cex = 2, groupOnX = TRUE, dodge.width = 0) +
  stat_summary(geom = "errorbar", color = "black",
               position = position_dodge(0.7), linewidth = 0.8,
               width = 0.4, fun.data = mean_se) +
  stat_compare_means(exact = FALSE, method = "wilcox",
                     comparisons = list(c(1,2), c(3,1), c(2,3)),
                     symnum.args = list(), 
                     p.adjust.method = "bonferroni", label = "p.value",
                     bracket.size = 0.85, tip.length = 0,
                     # label.y = c(0.95 * max(dataIn[[i]]), 
                     #             1.08 * max(dataIn[[i]]), 
                     #             1.2 * max(dataIn[[i]]), 
                     #             1.32 * max(dataIn[[i]])),
                     size = 9, face = "bold.italic", ref.group = ".all.") +
  scale_color_manual(values = c(
    "Vehicle"     = brewer.pal(9, "Oranges")[7],
    "WT Monomer"  = "#1E90FF",
    "WT Oligomer" = brewer.pal(9, "Blues")[9]
  )) +
  scale_fill_manual(values = c(
    "Vehicle"     = brewer.pal(9, "Oranges")[7],
    "WT Monomer"  = "#1E90FF",
    "WT Oligomer" = brewer.pal(9, "Blues")[9]
  )) +
  scale_x_discrete(labels = scales::wrap_format(5)) +
  # scale_y_continuous(limits = c(0, 10)) +
  ylab("PC1 Score") +
  xlab("") +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'white'),
    strip.text.x = element_text(size = 20, color = "black", face = "bold.italic"),
    strip.background = element_rect(colour = "white", fill = "white"),
    axis.line = element_line(colour = 'black', linewidth = 0.65, linetype = 'solid'),
    text = element_text(size = 25),
    axis.text.x = element_text(face = "bold", color = "black", vjust = 0.6, size = 23, angle = 90),
    axis.text.y = element_text(face = "bold", color = "black", size = 25),
    axis.title.y = element_text(size = 30, margin = margin(t = 0, r = 14, b = 0, l = 0)),
    axis.ticks.length = unit(.25, "cm"),
    legend.position = "none", 
    legend.spacing.x = unit(5, "pt")
  )

dev.off()



# Sort heatmap and loarding chart the same way 
indP1=sort(Data_P1, index.return=TRUE)$ix
pdf("./Heatmap LV1Sort.pdf", width=20,height=16, pointsize = 18)

heatmap3(df_Z[,indP1],
         col=barColors, breaks=breakBarColors,legendfun=function()showLegend(legend=c(NA),col=c(NA)), 
         Rowv=NA, Colv=NA,  scale="none",
         labRow = NA, RowSideColors = ConditionColorBar, RowSideLabs = c("Conditions"),
         labCol = colnames(df_Z[,indP1]), main="", cexRow=1, cexCol = 1.2) 

dev.off()
