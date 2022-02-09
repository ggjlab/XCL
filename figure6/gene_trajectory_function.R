require(WGCNA)
require(flashClust)
require(Hmisc)
require(dplyr)
require(openxlsx)
require(ggplot2)
require(cowplot)


## avg for variation in module score as a function of time. Compares many samplings of scores between pre-infection and each time point.
TestModuleTemporalVariation <- function(tps = c("St54","St59","St66"), #which time points to run the function on
                                        meta.data, sample.size = 50, ntest = 1000, name.of.feature = "M"){ #name.of.feature refers to prefix for columns with module scores
  
  #apply to run over one module at a time
  mod.min.pvals <- apply(meta.data[,grep(name.of.feature, colnames(meta.data), value=TRUE)], 2, function(mod){
    #get the scores for this cluster separated by time points in a list
    scores.by.tp = lapply(tps, function(time) mod[meta.data$stage == time])
    
    #run the wilcox test ntest times with sample.size between each timepoint and pre-infection, report the average p-val
    wilcox.pval.by.tp = mapply(scores = scores.by.tp[-1], tp = tps[-1], SIMPLIFY=FALSE, function(scores,tp){
      test.pval = replicate(ntest, expr={
        wilcox.test(x=sample(scores.by.tp[[1]], sample.size),y=sample(scores, sample.size))$p.value})
      return(mean(test.pval))
    })
    #return the min p-value for all comparison within that cluster, FDR corrected
    min.pval = min(p.adjust(unlist(wilcox.pval.by.tp)))
    return(min.pval)
  })
  return(mod.min.pvals) #returns the vector of minimum average p-value across all time point comparisons for each module
}

##
GetCelltypeSeuratObject <- function(St48cluster,St54cluster,St59cluster,St66cluster){
  load("/media/ggj/ggjlab2/LiaoYuan/xiaotu/St48/St48.RData")
  #Idents(pbmc) <- pbmc$RNA_snn_res.0.6
  DimPlot(pbmc,label=T)
  st48 <- subset(pbmc,idents = St48cluster)
  
  load("/media/ggj/ggjlab2/LiaoYuan/xiaotu/St54/St54.RData")
  #Idents(pbmc) <- pbmc$RNA_snn_res.0.6
  DimPlot(pbmc,label=T)
  st54 <- subset(pbmc,idents = St54cluster)
  #st54 <- as.data.frame(GetAssayData(st54@assays$RNA,slot="counts"))
  
  #st59
  load("/media/ggj/ggjlab2/LiaoYuan/xiaotu/St59/St59.RData")
  DimPlot(pbmc,label=T)
  st59 <- subset(pbmc,idents = St59cluster)
  #st59 <- as.data.frame(GetAssayData(st59@assays$RNA,slot="counts"))
  
  #st66
  load("/media/ggj/ggjlab2/LiaoYuan/xiaotu/St66/St66.RData")
  DimPlot(pbmc,label=T)
  st66 <- subset(pbmc,idents = St66cluster)
  #st66 <- as.data.frame(GetAssayData(st66@assays$RNA,slot="counts"))
  
  data <- merge(st48,st54)
  data <- merge(data,st59)
  data <- merge(data,st66)
  data <- NormalizeData(data)
  Idents(data) <- data$stage
  return(data)
} 

stageDEGcaculate <- function(data=data,logfc.threshold = 0.1,min.pct = 0.1){
  St48_St54 <- FindMarkers(data,ident.1 = "St48",ident.2 = "St54",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St48_St54 <- St48_St54[St48_St54$p_val_adj<=0.05,]
  
  St48_St59 <- FindMarkers(data,ident.1 = "St48",ident.2 = "St59",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St48_St59 <- St48_St59[St48_St59$p_val_adj<=0.05,]
  
  St48_St66 <- FindMarkers(data,ident.1 = "St48",ident.2 = "St66",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St48_St66 <- St48_St66[St48_St66$p_val_adj<=0.05,]
  
  St54_St59 <- FindMarkers(data,ident.1 = "St54",ident.2 = "St59",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St54_St59 <- St54_St59[St54_St59$p_val_adj<=0.05,]
  
  St54_St66 <- FindMarkers(data,ident.1 = "St54",ident.2 = "St66",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St54_St66 <- St54_St66[St54_St66$p_val_adj<=0.05,]
  
  St59_St66 <- FindMarkers(data,ident.1 = "St59",ident.2 = "St66",logfc.threshold = logfc.threshold,min.pct = min.pct)
  St59_St66 <- St59_St66[St59_St66$p_val_adj<=0.05,]
  
  gene.use <- c(rownames(St48_St54),rownames(St48_St59),rownames(St48_St66),rownames(St54_St59),rownames(St54_St66),rownames(St59_St66))
  gene.use <- unique(gene.use)
  
  return(gene.use)
}

SoftClustering <- function(data=data,cluster_num=30,celltype="Stromal",
                           gene.use=gene.use,usemestimate=T,m=1.25){
  library(Mfuzz)
  avg <- AverageExpression(data)
  #avg <- as.data.frame(avg$RNA)
  avg <- log(avg$RNA+1)
  avg <- as.data.frame(avg)
  
  #gene.use <- c(rownames(St48_St54),rownames(St48_St59),rownames(St48_St66),rownames(St54_St59),rownames(St54_St66),rownames(St59_St66))
  gene.use <- unique(gene.use)
  
  avg <- avg[gene.use,]
  
  #mfuzz
  min.mem = 0
  mfuzz <- as.matrix(avg)
  mfuzz_class <- new('ExpressionSet',exprs=mfuzz)
  
  mfuzz_class <- standardise(mfuzz_class)
  avg_standard <- mfuzz_class@assayData$exprs
  #dim(avg_standard)
  
  set.seed(123)
  
  cluster_num <- cluster_num
  if(usemestimate==T){
    m=mestimate(mfuzz_class)
  }
  else{
    m=m
  }
  mfuzz_cluster <- mfuzz(mfuzz_class,c=cluster_num,m=m)
  
  ##grep anno
  cluster_size <- mfuzz_cluster$size
  names(cluster_size) <- 1:cluster_num
  cluster_size
  
  head(mfuzz_cluster$cluster)
  
  avg_cluster <- mfuzz_cluster$cluster
  avg_standard <- mfuzz_class@assayData$exprs
  avg_standard_cluster <- cbind(avg_standard[names(avg_cluster),],avg_cluster)
  
  avg_standard_cluster <- as.data.frame(avg_standard_cluster)
  cluster.use <- avg_standard_cluster
  cluster.use$gene <- rownames(cluster.use)
  cluster.use <- cluster.use[,c(5,6)]
  
  cl=mfuzz_cluster
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  
  memship <- t(memship)
  memship[1:5,1:5]
  dim(memship)
  memship <- as.data.frame(memship)
  memship$cluster <- rownames(memship)
  memship.wide <- melt(memship,id="cluster")  
  cluster.use$type <- paste0(cluster.use$avg_cluster,"_",cluster.use$gene)
  memship.wide$type <- paste0(memship.wide$cluster,"_",memship.wide$variable)
  cluster.use <- merge(cluster.use,memship.wide,by="type")
  cluster.use <- cluster.use[,c(3,4,6)]
  
  dir.create(paste0("/media/ggj/ggjlab2/LiaoYuan/figure6/",celltype),showWarnings = F)
  
  csvname <- paste0("/media/ggj/ggjlab2/LiaoYuan/figure6/",celltype,"/module.csv")
  write.csv(cluster.use,csvname)
  
  cluster <- list(cluster=cluster.use,exp_matrix=avg_standard_cluster,membership=memship.wide)
  return(cluster)
}

plotSoftClusteringResults <-function(data=cluster.use,cluster_num=30,celltype="Stromal",avg_standard_cluster=avg_standard_cluster,memship.wide=memship.wide){
  library(gridExtra)
  plot <- list()
  cluster <- 1:cluster_num
  for(i in 1:length(cluster)){
    temp <- avg_standard_cluster[avg_standard_cluster$avg_cluster==cluster[i],]
    temp <- temp[,-c(5)]
    temp <- t(temp)
    temp <- as.data.frame(temp)
    
    
    temp <- temp[,data[data$cluster==i,]$gene]
    temp$stage <- rownames(temp)
    use <- melt(temp,id="stage")
    cluster2 <- memship.wide[memship.wide$cluster==cluster[i],]
    use <- merge(use,cluster2,by='variable')
    use$stage <- substr(use$stage,1,4)
    p1 <-  ggplot(data=use, aes(x=stage, y=value.x,group=variable,color=value.y))+scale_color_gradientn(colors = c("#49FF00","#FFFF00","#FF0000"),limits=c(0,1))+
      geom_line()+theme_classic()+ylab("Expression changes")+labs(col="Membership value")+ggtitle(paste0("Cluster",cluster[i]))+NoLegend()
    
    plot[[i]] <- p1
  }
  
  plot[['nrow']] <- 6
  plot[['ncol']] <- 5
  
  dir.create(paste0("/media/ggj/ggjlab2/LiaoYuan/figure6/",celltype),showWarnings = F)
  pdfname <- paste0("/media/ggj/ggjlab2/LiaoYuan/figure6/",celltype,"/stage_all.pdf")
  pdf(pdfname,width = 20,height = 15)
  do.call('grid.arrange', plot)
  dev.off()
} 


TestModuleTimeVariation <- function(type="Module1",data=data,cluster_num=30,
                                    gene.clusters.names.final=gene.clusters.names.final){
  if(type=="Module1"){
    test_result <- matrix(NA,nrow = length(gene.clusters.names.final),ncol = 4)
    test_result <- as.data.frame(test_result)
    colnames(test_result) <- c("wilcox1","wilcox2","wilcox3","isSig")
    for(i in 1:length(gene.clusters.names.final)){
      temp <- data@meta.data[gene.clusters.names.final[i]]
      temp$stage <- colsplit(rownames(temp),"_",names=c("n1","n2"))$n1
      
      colnames(temp)[1] <- "Modulescore"
      test1 <- wilcox.test(temp[temp$stage=="St48",]$Modulescore,temp[temp$stage=="St54",]$Modulescore,
                           alternative = "less")
      
      test2 <- wilcox.test(temp[temp$stage=="St54",]$Modulescore,temp[temp$stage=="St59",]$Modulescore,
                           alternative = "less")
      
      test3 <- wilcox.test(temp[temp$stage=="St59",]$Modulescore,temp[temp$stage=="St66",]$Modulescore,
                           alternative = "greater")
      
      test_result[i,1] <- test1$p.value
      test_result[i,2] <- test2$p.value
      test_result[i,3] <- test3$p.value
      if(test1$p.value<1e-3&test2$p.value<1e-3&test3$p.value<1e-3){
        test_result[i,4] <- "True"
      }
      else{
        test_result[i,4] <- "False"
      }
    }
  }     
  if(type=="Module2"){
    test_result <- matrix(NA,nrow = length(gene.clusters.names.final),ncol = 4)
    test_result <- as.data.frame(test_result)
    colnames(test_result) <- c("wilcox1","wilcox2","wilcox3","isSig")
    for(i in 1:length(gene.clusters.names.final)){
      temp <- data@meta.data[gene.clusters.names.final[i]]
      temp$stage <- colsplit(rownames(temp),"_",names=c("n1","n2"))$n1
      
      colnames(temp)[1] <- "Modulescore"
      test1 <- wilcox.test(temp[temp$stage=="St48",]$Modulescore,temp[temp$stage=="St54",]$Modulescore,
                           alternative = "greater")
      
      test2 <- wilcox.test(temp[temp$stage=="St54",]$Modulescore,temp[temp$stage=="St59",]$Modulescore,
                           alternative = "greater")
      
      test3 <- wilcox.test(temp[temp$stage=="St59",]$Modulescore,temp[temp$stage=="St66",]$Modulescore,
                           alternative = "less")
      
      test_result[i,1] <- test1$p.value
      test_result[i,2] <- test2$p.value
      test_result[i,3] <- test3$p.value
      if(test1$p.value<1e-3&test2$p.value<1e-3&test3$p.value<1e-3){
        test_result[i,4] <- "True"
      }
      else{
        test_result[i,4] <- "False"
      }
    }
    
  }
  if(type=="Module3"){
    test_result <- matrix(NA,nrow = length(gene.clusters.names.final),ncol = 4)
    test_result <- as.data.frame(test_result)
    colnames(test_result) <- c("wilcox1","wilcox2","wilcox3","isSig")
    for(i in 1:length(gene.clusters.names.final)){
      temp <- data@meta.data[gene.clusters.names.final[i]]
      temp$stage <- colsplit(rownames(temp),"_",names=c("n1","n2"))$n1
      
      colnames(temp)[1] <- "Modulescore"
      test1 <- wilcox.test(temp[temp$stage=="St48",]$Modulescore,temp[temp$stage=="St54",]$Modulescore,
                           alternative = "greater")
      
      test2 <- wilcox.test(temp[temp$stage=="St54",]$Modulescore,temp[temp$stage=="St59",]$Modulescore,
                           alternative = "greater")
      
      test3 <- wilcox.test(temp[temp$stage=="St59",]$Modulescore,temp[temp$stage=="St66",]$Modulescore,
                           alternative = "greater")
      
      test_result[i,1] <- test1$p.value
      test_result[i,2] <- test2$p.value
      test_result[i,3] <- test3$p.value
      if(test1$p.value<1e-3&test2$p.value<1e-3&test3$p.value<1e-3){
        test_result[i,4] <- "True"
      }
      else{
        test_result[i,4] <- "False"
      }
    }
    
  }
  if(type=="Module4"){
    test_result <- matrix(NA,nrow = length(gene.clusters.names.final),ncol = 4)
    test_result <- as.data.frame(test_result)
    colnames(test_result) <- c("wilcox1","wilcox2","wilcox3","isSig")
    for(i in 1:length(gene.clusters.names.final)){
      temp <- data@meta.data[gene.clusters.names.final[i]]
      temp$stage <- colsplit(rownames(temp),"_",names=c("n1","n2"))$n1
      
      colnames(temp)[1] <- "Modulescore"
      test1 <- wilcox.test(temp[temp$stage=="St48",]$Modulescore,temp[temp$stage=="St54",]$Modulescore,
                           alternative = "less")
      
      test2 <- wilcox.test(temp[temp$stage=="St54",]$Modulescore,temp[temp$stage=="St59",]$Modulescore,
                           alternative = "less")
      
      test3 <- wilcox.test(temp[temp$stage=="St59",]$Modulescore,temp[temp$stage=="St66",]$Modulescore,
                           alternative = "less")
      
      test_result[i,1] <- test1$p.value
      test_result[i,2] <- test2$p.value
      test_result[i,3] <- test3$p.value
      if(test1$p.value<1e-3&test2$p.value<1e-3&test3$p.value<1e-3){
        test_result[i,4] <- "True"
      }
      else{
        test_result[i,4] <- "False"
      }
    }
    
  }
  module <- rownames(test_result)[which(test_result$isSig=="True")]
  
  return(module) 
}

GetModuleGene <- function(module=module,data=data,cluster_num=30,Sig.Var.Modules=Sig.Var.Modules){
  gene <- NULL
  for(i in 1:length(module)){
    gene.temp <- Sig.Var.Modules[[as.numeric(module[i])]]
    gene <- unique(c(gene,gene.temp))
  }
  return(gene)
}

GetOrthgene <- function(gene=gene,celltype="Stromal",type="module1"){
  orthgene <- read.table("./orth/human_to_frog.txt",sep="\t",header=T)
  gene <- gsub("\\.[L|S|g]","",gene)
  gene.use <- orthgene[orthgene$Xenbase.gene.symbol%in%gene,]
  savepath <- paste0("./figure6/",celltype,"/",type,"_orth.csv")
  write.csv(gene.use$Human.gene.symbol,savepath)
  return(gene.use$Human.gene.symbol)
}