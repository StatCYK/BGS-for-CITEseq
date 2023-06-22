# Name: seurat_compare.R
# Function: Under current problem setting, comparing BSS with Seurat & scanpy
# ------------------------------------------------------------------------------
# rm(list = ls())

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(dbscan)
library(leaps)
library(caret) # PCA
library(jackstraw)
# BiocManager::install("MAST")
# BiocManager::install("DESeq2")
setwd("../CBMC analysis/code")
set.seed(123)
t1 = Sys.time()
# ------- Functions ------------------------------------------------------------
clr_function = function(x) { # For normalizing ADT (x)
  return(log(x = (1+x)/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
}

r_sq <- function(predicted, actual) { # For computing R-squared
  RSS <- sum((actual - predicted)^2)
  TSS <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (RSS / TSS)
  return(R_squared)
}

a_r_sq<- function(predicted, actual, p) {
  RSS <- sum((actual - predicted)^2)
  TSS <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (RSS / TSS)
  n <- length(actual)
  adjusted_R_squared <- 1 - ((1 - R_squared) * (n - 1) / (n - p - 1))
  return(adjusted_R_squared)
}


# ------ read and preprocess ---------------------------------------------------
cbmc.rna <- as.sparse(read.csv(file = "./real data/CBMC/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",",
                               header = TRUE, row.names = 1))

cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "./real data/CBMC/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
                               header = TRUE, row.names = 1))

all.equal(colnames(cbmc.rna), colnames(cbmc.adt))

# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

DefaultAssay(cbmc) <- 'RNA'
cbmc <- NormalizeData(cbmc,scale.factor = 1e3) %>% FindVariableFeatures() %>% ScaleData() 

CD14 = cbmc["rna_CD14"]
y = matrix(CD14@assays$RNA@data)
hist(y)

DefaultAssay(cbmc) <- 'ADT'
VariableFeatures(cbmc) <- rownames(cbmc[["ADT"]])
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

y = matrix(CD14@assays$RNA@data)
DefaultAssay(cbmc) <- 'ADT' 
x = t(cbmc@assays$ADT@scale.data)

data_8000 = as.data.frame(cbind(y,x))
colnames(data_8000) <- c("y",colnames(x))

write.csv(data_8000,"~/Dropbox/Projects/quantum computing/QAS/DATA/CiteSeq/cbmc/data_8000.csv")


### Saving un-normalized data ###
unnormalized_ADT = t(as.matrix(cbmc@assays$ADT@counts))
data_8000_unnormalized = as.data.frame(cbind(y, unnormalized_ADT))
row.names(data_8000_unnormalized) = NULL
write.csv(data_8000_unnormalized,"../data/data_8000_unnormalized.csv", row.names = F)

