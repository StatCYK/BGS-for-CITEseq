library(Seurat)
library(ggplot2)
library(SeuratData)
library(cowplot)
library(dplyr)
library(leaps)
library(SeuratDisk)
library(tidyverse)
library(scales)
library(foreach)
library(doParallel)

clr_function <- function(x) {
  return(log(x = (1+x)/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
}


setwd("./pbmc_HIV/")


PBMC <- LoadH5Seurat("./pbmc_multimodal.h5seurat",assay = "ADT")

marker.names <- read.csv("./data/ADTmarkers.csv")

df <- data.frame(
  cell.type = as.factor(PBMC$celltype.l3)
)


PBMC <- NormalizeData(PBMC, normalization.method = 'CLR', margin = 2) %>% ScaleData()


SaveH5Seurat(
  PBMC,
  filename = "./PBMC_ADT.h5Seurat",
  overwrite = FALSE
)


for(cell.name in unique(PBMC$celltype.l3)){
  # the 10 markers find by WNN paper for cell markers are 
  if(length(which(marker.names$celltype.l3 == cell.name)) >0 ){
  cell.marker <- marker.names$protein[which(marker.names$celltype.l3 == cell.name)]
  for(i in 1:length(cell.marker)){
    cell.marker[i] <-  gsub("`","",cell.marker[i])
  }
  cell.marker.data <- PBMC.subset[cell.marker]
  PBMC.cell = PBMC[cell.marker]
  cell.marker2 <- c()
  for(name in cell.marker){
    cell.marker2 <- c(cell.marker2, gsub("-",".",name))
  }
  y <- (PBMC$celltype.l3==cell.name)+0
  ADT.data <- t(PBMC.cell@assays$ADT@scale.data)
  
  colnames(ADT.data) = rownames(PBMC.cell@assays$ADT@scale.data)
  dat <- data.frame(cbind(ADT.data,y))
  write.csv(dat,paste0("./data/",cell.name,"_160000.csv"))
  }
}
