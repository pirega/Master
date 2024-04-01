#####################
# Código para análisis de expresión diferencial
#####################


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
#BiocManager::install("org.Hs.eg.db")

# 1. Carga de librerias

library (edgeR)
library(ggplot2)
library(org.Hs.eg.db)
library(statmod)
library(ggrepel)
library(UpSetR)
library(venn)
library(pheatmap)
library(matrixStats)
library(limma)
library(grid)

proyectos <- c("PRJNA558102", "PRJNA542148", "PRJNA794736")
for (i in proyectos){
  setwd(paste0("/media/scratch1/09MBIF/1_datos/3_anotacion/conteo/",i,"/"))
  
  #2. CARGA DE DATOS
  meta_data <- read.table("metadata.txt", header = TRUE, sep = "\t")
  names(meta_data)
  group <- factor(meta_data$categorizacion)
  enfermedad <- factor(meta_data$enfermedad)
  sex <- factor(meta_data$sex)
  project <- factor(meta_data$BioProject)
  data_counts <- readDGE(list.files(pattern = "SRR*"), group = group, labels = meta_data$Run)
  data_counts <-data_counts[,which(data_counts$samples$lib.size!="NA")]
  group <- group[which(data_counts$samples$lib.size!="NA")]   
  enfermedad <-    enfermedad[which(data_counts$samples$lib.size!="NA")]  
  project <- project[which(data_counts$samples$lib.size!="NA")]
  # 3. Adding information
  columns(org.Hs.eg.db)
  ann <- select(org.Hs.eg.db, keys = row.names(data_counts$counts), 
                columns = c("ENTREZID","SYMBOL", "GENENAME"), keytype = "SYMBOL")
  head (ann)
  table(ann$SYMBOL==row.names(data_counts$counts))
  
  # 4. Add the groups info
  data_counts$samples$group <- group
  
  #5. Add gene annotation
  data_counts$genes <- ann
  
  
  #6. Filtering to remove low counts
  keep <- rowSums(cpm(data_counts)>0) >=2
  table(keep)
  
  data_counts2 <- data_counts[keep,keep.lib.sizes=FALSE]
  
  #7. Qué proporción de genes se mantienen?
  dim(data_counts2)[1]/dim(data_counts)[1]*100
  
  # 8. Normalization
  data_counts2 <- calcNormFactors(data_counts2)
  data_counts2$samples
  
  #9. Multidimesional scaling plots
  #comparaciones <- c("disease","group")
  pch <- c(0,1,2,15,16,17)
  colors <- rep(c("darkgreen", "red", "blue"),2)
  pdf("MDS_plot.pdf")
  #par(mfrow=c(2,1))
  plotMDS(data_counts2,col=colors[group], pch=pch[group], main= "MDS plot")
  legend("topleft", legend=levels(group),pch=pch,col=colors, ncol=2, cex=0.8)
  plotMDS(data_counts2,col=colors[enfermedad], pch=pch[enfermedad], main= "MDS plot")
  legend("topleft", legend=levels(enfermedad),pch=pch,col=colors, ncol=2, cex=0.8)
  plotMDS(data_counts2,col=colors[sex], pch=pch[sex], main= "MDS plot")
  legend("topleft", legend=levels(sex),pch=pch,col=colors, ncol=2, cex=0.8)
  plotMDS(data_counts2,col=colors[project], pch=pch[project], main= "MDS plot")
  legend("topleft", legend=levels(project),pch=pch,col=colors, ncol=2, cex=0.8)
  dev.off()
  
  #10. Datos normalizados?
  normcounts <- cpm(data_counts2,normalized.lib.sizes = TRUE)
  
  
  
  #11.a Design matrix model    
  #Enfermedad  
  design <- model.matrix(~0 + enfermedad)
  colnames(design)<- levels(enfermedad)
  head(design)
  #12.a Estimate Dispersion
  data_counts2
  data_counts2 <- estimateDisp(data_counts2, design, robust = TRUE)
  data_counts2$common.dispersion
  pdf(paste0("enfermedad_BCV_plot.pdf"))
  plotBCV(data_counts2)
  dev.off()
  
  #13.a GLM Model Adjustment
  fit <- glmQLFit(data_counts2, design, robust = TRUE)
  pdf("enfermedad_QLDisp_plot.pdf")
  plotQLDisp(fit, ylab = "dispersion QL")
  dev.off()
  
  #14.a Prueba de significancia
  comparativa <- makeContrasts(NASH-NAFL, levels = design)
  res <- glmQLFTest(fit, contrast = comparativa)
  head(res$table)
  res_corrected <- topTags(res, n=Inf)
  head(res_corrected$table)
  res_sub <- res$table
  
  nrow(res_corrected$table[res_corrected$table$FDR <=0.05 & abs(res_corrected$table$logFC) >= 2, ])
  is.de <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 2)
  head(is.de)
  summary(is.de)
  pdf("MD_plot_enfermedad.pdf")
  plotMD(res, status = is.de, col=c("red", "blue"), legend="topright")
  dev.off()
  
  #15.b VolcanoPlot
  data <- res_corrected$table
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC >2 & data$FDR < 0.05]<- "UP"
  data$diffexpressed[data$logFC < -2 & data$FDR < 0.05]<- "DOWN"
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$SYMBOL[data$diffexpressed != "NO"]

  pdf("Volcano_plot_enfermedad.pdf")
  ggplot(data=data, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel()+
    scale_color_manual(values=c("blue","grey", "red"))+
    ggtitle("VolcanoPlot enfermedad")
  dev.off()
  data<- data[data$diffexpressed != "NO",]
  write.csv(data, file=paste0(i,"_DGE_enfermedad.csv"))
  
  #11.b Design matrix model    
  #Grupos del paper  
  design <- model.matrix(~0 + group)
  colnames(design)<- levels(group)
  head(design)
  #12.b Estimate Dispersion
  data_counts2
  data_counts2 <- estimateDisp(data_counts2, design, robust = TRUE)
  data_counts2$common.dispersion
  pdf(paste0("groups_BCV_plot.pdf"))
  plotBCV(data_counts2)
  dev.off()
  
  #13.b GLM Model Adjustment
  fit <- glmQLFit(data_counts2, design, robust = TRUE)
  pdf("Group_QLDisp_plot.pdf")
  plotQLDisp(fit, ylab = "dispersion QL")
  dev.off()
  
  #14.b Prueba de significancia
  #comparacion <- c("NASH_F0_F1-NASH_F2", "NASH_F2-NASH_F3", "NASH_F3-NASH_F4")
  comparativa <- makeContrasts(NASH_F2-NASH_F0_F1, levels = design)
  res <- glmQLFTest(fit, contrast = comparativa)
  head(res$table)
  res_corrected <- topTags(res, n=Inf)
  head(res_corrected$table)
  res_sub <- res$table
  
  nrow(res_corrected$table[res_corrected$table$FDR <=0.05 & abs(res_corrected$table$logFC) >= 2, ])
  is.de <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 2)
  head(is.de)
  summary(is.de)
  pdf("Grupos_MD_plot_NASHF0_1vsF2.pdf")
  plotMD(res, status = is.de, col=c("red", "blue"), legend="topright")
  dev.off()
  
  #15.b VolcanoPlot
  data <- res_corrected$table
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC >2 & data$FDR < 0.05]<- "UP"
  data$diffexpressed[data$logFC < -2 & data$FDR < 0.05]<- "DOWN"
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$SYMBOL[data$diffexpressed != "NO"]

  pdf("Volcano_plot_NASH_F2vsF0-1.pdf")
  ggplot(data=data, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel()+
    scale_color_manual(values=c("blue","grey", "red"))+
    ggtitle("Group_VolcanoPlot")
  dev.off()
  data<- data[data$diffexpressed != "NO",]
  write.csv(data, file=paste0(i,"_grupos_DGE_NASH_F2vsF0-1.csv"))
  
  #14.c Prueba de significancia
  #comparacion <- c("NASH_F0_F1-NASH_F2", "NASH_F2-NASH_F3", "NASH_F3-NASH_F4")
  comparativa <- makeContrasts(NASH_F3-NASH_F2, levels = design)
  res <- glmQLFTest(fit, contrast = comparativa)
  head(res$table)
  res_corrected <- topTags(res, n=Inf)
  head(res_corrected$table)
  res_sub <- res$table
  
  nrow(res_corrected$table[res_corrected$table$FDR <=0.05 & abs(res_corrected$table$logFC) >= 2, ])
  is.de <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 2)
  head(is.de)
  summary(is.de)
  pdf("Grupos_MD_plot_NASH_F3vsF2.pdf")
  plotMD(res, status = is.de, col=c("red", "blue"), legend="topright")
  dev.off()
  
  #15.c VolcanoPlot
  data <- res_corrected$table
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC >2 & data$FDR < 0.05]<- "UP"
  data$diffexpressed[data$logFC < -2 & data$FDR < 0.05]<- "DOWN"
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- data$SYMBOL[data$diffexpressed != "NO"]
  pdf("Volcano_plot_NASH_F3vsF2.pdf")
  ggplot(data=data, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel()+
    scale_color_manual(values=c("blue","grey", "red"))+
    ggtitle("Group_VolcanoPlot")
  dev.off()
  data<- data[data$diffexpressed != "NO",]
  write.csv(data, file=paste0(i,"_grupos_DGE_NASH_F3vsF2.csv"))
  
  
  #14.d Prueba de significancia
  #comparacion <- c("NASH_F0_F1-NASH_F2", "NASH_F2-NASH_F3", "NASH_F3-NASH_F4")
  comparativa <- makeContrasts(NASH_F4-NASH_F3, levels = design)
  res <- glmQLFTest(fit, contrast = comparativa)
  head(res$table)
  res_corrected <- topTags(res, n=Inf)
  head(res_corrected$table)
  res_sub <- res$table
  
  nrow(res_corrected$table[res_corrected$table$FDR <=0.05 & abs(res_corrected$table$logFC) >= 2, ])
  is.de <- decideTests(res, adjust.method = "BH", p.value = 0.05, lfc = 2)
  head(is.de)
  summary(is.de)
  pdf("Grupos_MD_plot_NASH_F4vsF3.pdf")
  plotMD(res, status = is.de, col=c("red", "blue"), legend="topright")
  dev.off()
  
  #15.d VolcanoPlot
  data <- res_corrected$table
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC >2 & data$FDR < 0.05]<- "UP"
  data$diffexpressed[data$logFC < -2 & data$FDR < 0.05]<- "DOWN"
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO"] <- rownames(data)[data$diffexpressed != "NO"]
  pdf("Volcano_plot_NASH_F4vsF3.pdf")
  ggplot(data=data, aes(x=logFC, y=-log10(FDR), col=diffexpressed, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel()+
    scale_color_manual(values=c("blue","grey", "red"))+
    ggtitle("Group_VolcanoPlot")
  dev.off()
  data<- data[data$diffexpressed != "NO",]
  write.csv(data, file=paste0(i,"_grupos_DGE_NASH_F4vsF3.csv"))
  
  save.image(paste0(i,".RData"))
}
