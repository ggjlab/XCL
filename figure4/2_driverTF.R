setwd("/media/ggj/ggjlab2/LiaoYuan/")
library(reshape2)
library(base)
####process JSDR result
#NF48
JSDR <- read.table("./github_scripts/figure4/result/St48_TF_JSDR.out",header = T,sep="\t")
colnames(JSDR) <- c("Gene","Celltype","JSDR")
JSDR_1 <- dcast(JSDR,Gene~Celltype) 
JSDR_1[1:2,1:2]
rownames(JSDR_1) <- JSDR_1$Gene
JSDR_1 <- JSDR_1[,-1]
which(is.na(JSDR_1[1,]))
which(is.na(JSDR_1[,1]))
JSDR_1 <- scale(JSDR_1,center = TRUE)
JSDR_1[JSDR_1<0]=0
St48_Averge <- JSDR_1

#NF54
JSDR <- read.table("./github_scripts/figure4/result/St54_TF_JSDR.out",header = T,sep="\t")
colnames(JSDR) <- c("Gene","Celltype","JSDR")
JSDR_1 <- dcast(JSDR,Gene~Celltype) 
JSDR_1[1:2,1:2]
rownames(JSDR_1) <- JSDR_1$Gene
JSDR_1 <- JSDR_1[,-1]
which(is.na(JSDR_1[1,]))
which(is.na(JSDR_1[,1]))
JSDR_1 <- scale(JSDR_1,center = TRUE)
JSDR_1[JSDR_1<0]=0
St54_Averge <- JSDR_1

#NF59
JSDR <- read.table("./github_scripts/figure4/result/St59_TF_JSDR.out",header = T,sep="\t")
colnames(JSDR) <- c("Gene","Celltype","JSDR")
JSDR_1 <- dcast(JSDR,Gene~Celltype) 
JSDR_1[1:2,1:2]
rownames(JSDR_1) <- JSDR_1$Gene
JSDR_1 <- JSDR_1[,-1]
which(is.na(JSDR_1[1,]))
which(is.na(JSDR_1[,1]))
JSDR_1 <- scale(JSDR_1,center = TRUE)
JSDR_1[JSDR_1<0]=0
St59_Averge <- JSDR_1

#NF66
JSDR <- read.table("./github_scripts/figure4/result/St66_TF_JSDR.out",header = T,sep="\t")
colnames(JSDR) <- c("Gene","Celltype","JSDR")
JSDR_1 <- dcast(JSDR,Gene~Celltype) 
JSDR_1[1:2,1:2]
rownames(JSDR_1) <- JSDR_1$Gene
JSDR_1 <- JSDR_1[,-1]
which(is.na(JSDR_1[1,]))
which(is.na(JSDR_1[,1]))
JSDR_1 <- scale(JSDR_1,center = TRUE)
JSDR_1[JSDR_1<0]=0
St66_Averge <- JSDR_1


colnames(St48_Averge) <- paste(colnames(St48_Averge),"NF48",sep = "|")
colnames(St54_Averge) <- paste(colnames(St54_Averge),"NF54",sep = "|")
colnames(St59_Averge) <- paste(colnames(St59_Averge),"NF59",sep = "|")
colnames(St66_Averge) <- paste(colnames(St66_Averge),"NF66",sep = "|")

Gene <- intersect(rownames(St48_Averge), rownames(St54_Averge))
Gene <- intersect(Gene, rownames(St59_Averge))
Gene <- intersect(Gene, rownames(St66_Averge))
length(Gene)
St48_Averge <- St48_Averge[Gene,]
St54_Averge <- St54_Averge[Gene,]
St59_Averge <- St59_Averge[Gene,]
St66_Averge <- St66_Averge[Gene,]

NF48_Averge <- St48_Averge
NF54_Averge <- St54_Averge
NF59_Averge <- St59_Averge
NF66_Averge <- St66_Averge

####process annotation
stage <- c("NF48","NF54","NF59","NF66")
anno <- NULL
for(i in 1:4){
  temp <- readxl::read_xlsx("./github_scripts/figure4/celltype_annotation.xlsx",i)
  colnames(temp) <- c("cluster","celltype")
  temp <- temp[order(temp$cluster),]
  temp$stage <- stage[i]
  temp$type <- paste0(temp$cluster,"|",temp$stage)
  
  anno <- rbind(anno,temp)
}

anno$celltype2 <- colsplit(anno$celltype,"_",names = c("n1","n2"))$n1
anno$celltype2 <- gsub("naïve T cell|T helper cell|Proliferating T cell|Regulatory T cell","T cell",anno$celltype2)
anno$celltype2 <- gsub("Proliferating macrophage","Macrophage",anno$celltype2)
anno$celltype2 <- gsub("Proliferating epithelial cell","Epithelial cell",anno$celltype2)
anno$celltype2 <- gsub("Proliferating muscle satellite cell","Muscle satellite cell",anno$celltype2)
anno$celltype2 <- gsub("Excitatory neuron","Neuron",anno$celltype2)
anno$celltype2 <- gsub("Proliferating erythrocyte","Erythrocyte",anno$celltype2)
anno$celltype2 <- gsub("Distal chondrocyte","Chondrocyte",anno$celltype2)
anno$celltype2 <- gsub("Vascular endothelial cell","Endothelial cell",anno$celltype2)
unique(anno$celltype2)
anno$type <- paste0(anno$celltype2,"|",anno$stage)
celltype1 <- unique(anno$celltype2)

####driverTF
for(i in 1:length(celltype1)){
  anno.temp <- anno[anno$celltype2==celltype1[i],]
  stage.num <- unique(anno.temp$stage)
  stage_matrix <- c("NF48_Averge","NF54_Averge","NF59_Averge","NF66_Averge")
  if(length(stage.num)==4){
    ID_48 <- anno.temp[anno.temp$stage==stage.num[1],]$type
    ID_54 <- anno.temp[anno.temp$stage==stage.num[2],]$type
    ID_59 <- anno.temp[anno.temp$stage==stage.num[3],]$type
    ID_66 <- anno.temp[anno.temp$stage==stage.num[4],]$type
    
    Gene_set1 <-  as.matrix(St48_Averge[,unique(ID_48)] * St54_Averge[,unique(ID_54)])
    Gene_set2 <-  as.matrix(St54_Averge[,unique(ID_54)] * St59_Averge[,unique(ID_59)])
    Gene_set3 <-  as.matrix(St59_Averge[,unique(ID_59)] * St66_Averge[,unique(ID_66)])
    Gene_set <- cbind(Gene_set1,Gene_set2,Gene_set3)
    
    Gene_set_TF <- Gene_set
    Gene_set_TF <- as.data.frame(Gene_set_TF)
    A=Gene_set_TF[Gene_set_TF$V1 > 1,]
    B=Gene_set_TF[Gene_set_TF$V2 > 1,]
    C=Gene_set_TF[Gene_set_TF$V3 > 1,]
    tf_m=union(rownames(A),rownames(B))
    tf_m=union(tf_m,rownames(C))
    Gene_set_TF_USE <- Gene_set_TF[tf_m,]
    Gene_set_TF_USE <- Gene_set_TF_USE[order(Gene_set_TF_USE$V1,decreasing = T),]
    
    colnames(Gene_set_TF_USE) <- c("NF48_NF54","NF54_NF59","NF59_NF66")
    dir.create("./github_scripts/figure4/driver_TF_result",showWarnings = F)
    
    filename <- paste0("./github_scripts/figure4/driver_TF_result/",celltype1[i],"_driverTF_result.csv")
    write.csv(Gene_set_TF_USE,file = filename)
  }
  if(length(stage.num)==3){
    ID1 <- anno.temp[anno.temp$stage==stage.num[1],]$type
    ID2 <- anno.temp[anno.temp$stage==stage.num[2],]$type
    ID3 <- anno.temp[anno.temp$stage==stage.num[3],]$type
    
    num1 <- which(colsplit(stage_matrix,"_",names = c("n1","n2"))$n1==stage.num[1])
    num2 <- which(colsplit(stage_matrix,"_",names = c("n1","n2"))$n1==stage.num[2])
    num3 <- which(colsplit(stage_matrix,"_",names = c("n1","n2"))$n1==stage.num[3])
    Gene_set1 <-  as.matrix(get(stage_matrix[num1])[,unique(ID1)] * get(stage_matrix[num2])[,unique(ID2)])
    Gene_set2 <-  as.matrix(get(stage_matrix[num2])[,unique(ID2)] * get(stage_matrix[num3])[,unique(ID3)])
    Gene_set <- cbind(Gene_set1,Gene_set2)
    
    Gene_set_TF <- Gene_set
    Gene_set_TF <- as.data.frame(Gene_set_TF)
    A=Gene_set_TF[Gene_set_TF$V1 > 1,]
    B=Gene_set_TF[Gene_set_TF$V2 > 1,]
    
    tf_m=union(rownames(A),rownames(B))
    tf_m=union(tf_m,rownames(C))
    Gene_set_TF_USE <- Gene_set_TF[tf_m,]
    Gene_set_TF_USE <- Gene_set_TF_USE[order(Gene_set_TF_USE$V1,decreasing = T),]
    
    colnames(Gene_set_TF_USE) <- c(paste0(stage.num[1],"_",stage.num[2]),paste0(stage.num[2],"_",stage.num[3]))
    dir.create("./github_scripts/figure4/driver_TF_result",showWarnings = F)
    
    filename <- paste0("./github_scripts/figure4/driver_TF_result/",celltype1[i],"_driverTF_result.csv")
    write.csv(Gene_set_TF_USE,file = filename)
  }
  
}





