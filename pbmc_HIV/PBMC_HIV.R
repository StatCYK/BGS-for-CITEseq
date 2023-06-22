library(Seurat)
library(ggplot2)
library(SeuratData)
library(cowplot)
library(dplyr)
library(harmony)
library(leaps)
library(foreach)
library(doParallel)

setwd("~/Dropbox/Projects/quantum computing/QAS/DATA/pbmc_HIV/")

# Load in RNAsqc Data
Cite.rna <- Matrix::readMM('./GSM5008737_RNA_3P-matrix.mtx')
rna.names <- read.csv(file = "./GSM5008737_RNA_3P-features.tsv", sep = "\t",header = FALSE)$V1
rna.bar<-read.csv(file = "./GSM5008737_RNA_3P-barcodes.tsv.gz", sep = ",",header = FALSE)

# Load in the ADT UMI matrix
Cite.adt <- Matrix::readMM(file = "./GSM5008738_ADT_3P-matrix.mtx.gz")
adt.bar <- read.csv(file = "./GSM5008738_ADT_3P-barcodes.tsv.gz", sep = ",",header = FALSE)
adt.names <- read.csv(file = "./GSM5008738_ADT_3P-features.tsv", sep = "\t", header = FALSE)$V1

# Load in the HTO data
HTO.mat <- Matrix::readMM(file = "./GSM5008739_HTO_3P-matrix.mtx.gz")
HTO.info <-read.csv("./GSE164378_Antibody_HTO_info.csv",sep = ",",header = TRUE)
HTO.names <- read.csv(file = "./GSM5008739_HTO_3P-features.tsv.gz", sep = "\t", header = FALSE)$V1


Donor.list <- c()
date.list <- c()
batch.list <- c()
for (i in 1:24) {
  idx = which(HTO.info$X3..HTO == HTO.names[i])
  Donor.list[i] = HTO.info$Donor[idx]
  date.list[i] = HTO.info$time[idx]
  batch.list[i] = HTO.info$Batch[idx]
}

sample.ID = apply(HTO.mat, 2, which.max)
sample.donor = Donor.list[sample.ID]
sample.date = date.list[sample.ID]
sample.batch = batch.list[sample.ID]

row.names(Cite.rna) = rna.names
row.names(Cite.adt) = adt.names

colnames(Cite.rna) = colnames(Cite.adt) = adt.bar$V1

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
Cite.rna.rdc <- CollapseSpeciesExpressionMatrix(Cite.rna,ncontrols = 1000)

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(Cite.rna), colnames(Cite.adt))

# creates a Seurat object based on the scRNA-seq data
Cite <- CreateSeuratObject(counts = Cite.rna.rdc)

# We can see that by default, the Cite object contains an assay storing RNA measurement
Assays(Cite)

# create a new assay to store ADT information
adt_assay <- CreateAssayObject(counts = Cite.adt)

# add this assay to the previously created Seurat object
Cite[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(Cite)

DefaultAssay(Cite) <- 'RNA'
Cite <- NormalizeData(Cite) %>% FindVariableFeatures() %>% ScaleData()
Cite <- RunPCA(Cite,reduction.name = 'pca')
Cite <- RunUMAP(Cite, reduction='pca',reduction.name = "RNA_UMAP", dims = 1:30, reduction.key = "rna.UMAP_")
Cite <- FindNeighbors(Cite, reduction='pca')
Cite <- FindClusters(Cite, resolution = 0.8)

Cite$donor <- sample.donor
Cite$batch <- sample.batch

#### visualize batch effect ######
p1=DimPlot(Cite, reduction = 'RNA_UMAP',label = TRUE, group.by = "donor", 
        repel = TRUE, label.size = 3)+ggtitle(NULL)+xlab("RNA_UMAP1")+ylab("RNA_UMAP2")
p2=DimPlot(Cite, reduction = 'RNA_UMAP',label = TRUE, group.by = "batch", 
           repel = TRUE, label.size = 3)+ggtitle(NULL)+xlab("RNA_UMAP1")+ylab("RNA_UMAP2") 
p1+p2

### batch effect correction #######
DefaultAssay(Cite) <- 'RNA'
Cite <- RunHarmony(Cite,"donor", plot_convergence = TRUE)

#### downstream analysis after harmony correction ####
DefaultAssay(Cite) <- 'RNA'
Cite <- RunUMAP(Cite, reduction = "harmony", dims = 1:30,reduction.name = "umap") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

#### visualize after correct batch effect ######
p3=DimPlot(Cite, reduction = 'umap',label = TRUE, group.by = "donor", 
           repel = TRUE, label.size = 3)+ggtitle(NULL)+xlab("RNA_UMAP1")+ylab("RNA_UMAP2")
p4=DimPlot(Cite, reduction = 'umap',label = TRUE, group.by = "batch", 
           repel = TRUE, label.size = 3)+ggtitle(NULL)+xlab("RNA_UMAP1")+ylab("RNA_UMAP2") 

p5=DimPlot(Cite, reduction = 'umap',label = TRUE, group.by = "RNA_snn_res.1", 
           repel = TRUE, label.size = 3)+ggtitle(NULL)+xlab("RNA_UMAP1")+ylab("RNA_UMAP2") 

p3+p4+p5


DefaultAssay(Cite) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(Cite) <- rownames(Cite[["ADT"]])
Cite <- NormalizeData(Cite, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() #%>% RunPCA(reduction.name = 'apca')

clr_function <- function(x) {
  return(log(x = (1+x)/(exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE)/length(x = x)))))
}

ADT.data <- t(apply(Cite@assays$ADT@counts,2,clr_function))
colnames(ADT.data) = adt.names
adt.sd <- apply(ADT.data, 2, sd) # show the adt.sd
hist(adt.sd,40)

x = ADT.data[,order(adt.sd,decreasing = TRUE)[1:100]]
x = scale(x)
### HUMAN_CD14 as response ##

rna.CD14.idx = which(row.names(Cite.rna.rdc) == "CD14")

FeatureScatter(Cite, feature1 = "adt_CD14", feature2 = "rna_CD14")

DefaultAssay(Cite) <- 'RNA'
CD14 = Cite["rna_GZMB"]#Cite["rna_CD14"]
y = matrix(CD14@assays$RNA@scale.data)
write.csv(y,"./RNA_CD14.csv")
DefaultAssay(Cite) <- 'ADT'
#x = t(Cite@assays$ADT@scale.data)
write.csv(x,"./ADT228.csv")


fit = lm(y~x)
res = summary(fit)
res

models <- regsubsets(Fertility~., data = swiss, nvmax = 5)

sample.info = data.frame(date = sample.date,donor =sample.donor)
RNA.ADT.data = data.frame(cbind(y,x,sample.info))

fmla = as.formula(paste("y ~", paste(colnames(x), collapse = " + ")))
models <- regsubsets(y~., data = RNA.ADT.data, really.big = TRUE)
summary(models)
write.csv(RNA.ADT.data,"./RNA_CD14_ADT.csv")



gene.name <- "CDKN1C"
CD14 = Cite[gene.name]#Cite["rna_CD14"]
y = matrix(CD14@assays$RNA@scale.data)
#### summary changes of rna CD14 in each sample ####
rna.CD14.smp.lvl = c()
for (i in 1:24) {
  sample.i = which(sample.ID == i)
  y.sample = y[sample.i]
  rna.CD14.smp.lvl[i] = mean(y.sample[y.sample>0])#mean(y[sample.i]>0)
}
sample.CD14.data = data.frame(rna.CD14 = rna.CD14.smp.lvl, date = date.list, donor = Donor.list) #rna.CD14.smp.lvl

ggplot(sample.CD14.data, aes(x=date, y=rna.CD14) ) + ylab(gene.name)+
  #geom_boxplot(alpha = 0.2,width=0.5,aes(x=date, y=rna.CD14,color = date,fill = date))+
  geom_line(size=1.2,aes(group = donor,color = donor,linetype = donor))+
  geom_point(size = 2) +
  theme_bw()




x.sample = x[sample.donor == "P2",]
y.sample = y[sample.donor == "P2"]
x.sample = x.sample[y.sample>0,]
y.sample = y.sample[y.sample>0]

fit = lm(y~x)
res = summary(fit)

fit = lm(y.sample~x.sample)
res = summary(fit)
res
adt.test = res$coefficients[,4]

### subsampling to reduce compute cost ####
#x = read.csv("./ADT228.csv",row.names = 1)


y = read.csv("./RNA_CD14.csv",row.names = 1)
set.seed(1234)

sample.idx = sample(c(1:nrow(x)), size = 1000, replace = FALSE)
adt.test = c()
for (p in 1:ncol(x)) {
  adt.test[p] = cor.test(y[sample.idx],x[sample.idx,p],method = "spearman")$p.value
}


colnames(x)[order(abs(adt.test),decreasing = FALSE)[1:30]]
# 
adt.test.adj = p.adjust(adt.test)
colnames(x)[which(adt.test.adj<0.01)] 

