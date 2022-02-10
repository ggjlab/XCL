setwd("E:\\Xenopus\\Xenopus _DGE\\Adult")
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(reshape2)
message("read dge files")
aa<-list.files()
dgemeta<-colsplit(aa,"_",c(1,2,3,4))
dgemeta$batch<-paste0(dgemeta$`1`,"_",dgemeta$`2`)
dgemeta$tissue<-dgemeta$`1`
dgemeta$date<-dgemeta$`3`

qq=1
data <- read.table(aa[qq],header = T,row.names = 1)
colnames(data)<-paste0(dgemeta$batch[qq],".",colnames(data))

message("start rmbatch")
data_500less <- data[,colSums(data)<500]
data_500more<-data[,colSums(data)>=500]
more500<- data_500more
less<- data_500less
raw<-data
allumi<-data.frame(umi=colSums(less))
ssa<-allumi[with(allumi,order(umi,decreasing = F)),]
if (length(ssa)>500){
  ss<-rownames(allumi)[with(allumi,order(umi,decreasing = F))][1:500]
  less<-less[,ss]
} else{ss<-rownames(allumi)[with(allumi,order(umi,decreasing = F))][length(ssa)]}
less_data<-data.frame(gene=rowSums(less))
#table(less_data$gene>10)
usegene<-rownames(less_data)[less_data$gene>10]
more500<-more500[usegene,]
less<-less[usegene,]
raw<-raw[usegene,]
background <- data.frame(var=replicate(1,n = nrow(more500)),
                         cellnum_express =rowSums(more500>0),
                         rowMean_500more =rowMeans(more500),
                         row.names = rownames(more500)
                         ,rowMeans_all=rowMeans(raw)						
)
temp <- merge(background,data.frame(rowMean_less =rowMeans(less)),all.x=F, by="row.names")
background <- data.frame(temp[,-1],row.names = temp[,1])
for (m in rownames(background)){
  background[m,"var"] <- var(as.numeric(more500[m,]))
  background[m,"sd"] <- sqrt(background[m,"var"])
}
background <- background[with(background,order(-rowMean_less,-rowMean_500more,-cellnum_express, -sd)),]
background$multi<-background$rowMean_less*background$sd
background<-background[background$multi>=1  ,]
med<-median(background$rowMean_500more/background$rowMean_less)
background[,"batchValue"] <- background[,"rowMean_less"]*med
background$batchValue <- round(background$batchValue) 
background <- background[background$batchValue>0,]
dge_m<- data_500more
result_dge <- dge_m
for (i in rownames(background)) { result_dge[i,] <- result_dge[i,]-background[i,"batchValue"] }
result_dge[result_dge<0] <- 0   
seurat_object<-CreateSeuratObject(result_dge)
message("rmbatch is done")
message("doubletfinder start")
message("seurat process start")
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt:")
seurat_object <- subset(seurat_object, subset = percent.mt < 20)
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- FindNeighbors(seurat_object, dims = 1:50)
seurat_object <- FindClusters(seurat_object, resolution = 1)
#save(seurat_object, file = "seurat_object_rmbatch.rdata")
message("seurat process done")

message("DoubletFinder process start")
library(DoubletFinder)
doublettest <- paramSweep_v3(seurat_object, PCs = 1:50, sct = FALSE)
#head(sweep.res.list_kidney)
doublettest2 <- summarizeSweep(doublettest, GT = FALSE)
#head(sweep.res.list_kidney)
doublettest3 <- find.pK(doublettest2)
mpK<-as.numeric(as.vector(doublettest3$pK[which.max(doublettest3$BCmetric)]))
seurat_object$use<-Idents(seurat_object)
annotations <- seurat_object$use
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*length(seurat_object$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
final <- doubletFinder_v3(seurat_object, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#DimPlot(final,group.by=colnames(final@meta.data)[grep("DF",colnames(final@meta.data))] )
#FeaturePlot(final,features = "nCount_RNA")

cellname.use<-rownames(final@meta.data)[final@meta.data[9]=="Singlet"]
use<-subset(final,cells=cellname.use)
use$tissue<-dgemeta$tissue[qq]
use$batch<-dgemeta$batch[qq]
use$date<-dgemeta$date[qq]

mergeuse<-use

for (qq in 26:length(aa)) {
  data <- read.table(aa[qq],header = T,row.names = 1)
  colnames(data)<-paste0(dgemeta$batch[qq],".",colnames(data))
  message("start rmbatch")
  data_500less <- data[,colSums(data)<500]
  data_500more<-data[,colSums(data)>=500]
  more500<- data_500more
  less<- data_500less
  raw<-data
  allumi<-data.frame(umi=colSums(less))
  ssa<-allumi[with(allumi,order(umi,decreasing = F)),]
  if (length(ssa)>500){
    ss<-rownames(allumi)[with(allumi,order(umi,decreasing = F))][1:500]
    less<-less[,ss]
  } else{ss<-rownames(allumi)[with(allumi,order(umi,decreasing = F))][length(ssa)]}
  less_data<-data.frame(gene=rowSums(less))
  #table(less_data$gene>10)
  usegene<-rownames(less_data)[less_data$gene>10]
  more500<-more500[usegene,]
  less<-less[usegene,]
  raw<-raw[usegene,]
  background <- data.frame(var=replicate(1,n = nrow(more500)),
                           cellnum_express =rowSums(more500>0),
                           rowMean_500more =rowMeans(more500),
                           row.names = rownames(more500)
                           ,rowMeans_all=rowMeans(raw)						
  )
  temp <- merge(background,data.frame(rowMean_less =rowMeans(less)),all.x=F, by="row.names")
  background <- data.frame(temp[,-1],row.names = temp[,1])
  for (m in rownames(background)){
    background[m,"var"] <- var(as.numeric(more500[m,]))
    background[m,"sd"] <- sqrt(background[m,"var"])
  }
  background <- background[with(background,order(-rowMean_less,-rowMean_500more,-cellnum_express, -sd)),]
  background$multi<-background$rowMean_less*background$sd
  background<-background[background$multi>=1  ,]
  med<-median(background$rowMean_500more/background$rowMean_less)
  background[,"batchValue"] <- background[,"rowMean_less"]*med
  background$batchValue <- round(background$batchValue) 
  background <- background[background$batchValue>0,]
  dge_m<- data_500more
  result_dge <- dge_m
  for (i in rownames(background)) { result_dge[i,] <- result_dge[i,]-background[i,"batchValue"] }
  result_dge[result_dge<0] <- 0   
  seurat_object<-CreateSeuratObject(result_dge)
  message("rmbatch is done")
  message("doubletfinder start")
  message("seurat process start")
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt:")
  seurat_object <- subset(seurat_object, subset = percent.mt < 20)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2500)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- FindNeighbors(seurat_object, dims = 1:50)
  seurat_object <- FindClusters(seurat_object, resolution = 1)
  #save(seurat_object, file = "seurat_object_rmbatch.rdata")
  message("seurat process done")
  
  message("DoubletFinder process start")
  library(DoubletFinder)
  doublettest <- paramSweep_v3(seurat_object, PCs = 1:50, sct = FALSE)
  #head(sweep.res.list_kidney)
  doublettest2 <- summarizeSweep(doublettest, GT = FALSE)
  #head(sweep.res.list_kidney)
  doublettest3 <- find.pK(doublettest2)
  mpK<-as.numeric(as.vector(doublettest3$pK[which.max(doublettest3$BCmetric)]))
  seurat_object$use<-Idents(seurat_object)
  annotations <- seurat_object$use
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.05*length(seurat_object$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  final <- doubletFinder_v3(seurat_object, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  #DimPlot(final,group.by=colnames(final@meta.data)[grep("DF",colnames(final@meta.data))] )
  #FeaturePlot(final,features = "nCount_RNA")
  
  cellname.use<-rownames(final@meta.data)[final@meta.data[9]=="Singlet"]
  use<-subset(final,cells=cellname.use)
  use$tissue<-dgemeta$tissue[qq]
  use$batch<-dgemeta$batch[qq]
  use$date<-dgemeta$date[qq]
  
  mergeuse<-merge(mergeuse,use)
  
}
table(mergeuse$batch)
save(mergeuse,file = "Adultraw1-25.rdata")

mergeuse<-SCTransform(mergeuse)
data<-CreateSeuratObject(mergeuse@SCT$scale.data)
SaveH5Seurat(data, filename = "adult.h5Seurat")
Convert("adult.h5Seurat", dest = "h5ad")