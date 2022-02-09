setwd("/media/ggj/ggjlab2/LiaoYuan/")

library(Seurat)
library(Mfuzz)
library(reshape2)
library(gridExtra)
library(org.Xl.eg.db)    
library(ggplot2)   
library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(openxlsx)
source("./github_scripts/figure6/gene_trajectory_function.R")

# Neuron
# get neuron cell clusters from four stages
data <- GetCelltypeSeuratObject(St48cluster = c(5),
                                St54cluster = c(21),
                                St59cluster = c(29),
                                St66cluster = c(45))
gc()

# calculate differential expressed genes
# use 'FindMarkers' function in Seurat 
gene.use <- stageDEGcaculate(data,logfc.threshold = 0.1,min.pct = 0.1)

# We performed soft clustering for differential expressed genes
cluster <- SoftClustering(data,gene.use = gene.use,usemestimate = T,cluster_num = 30,celltype="Neuron")

# plot 
plotSoftClusteringResults(data=cluster$cluster,cluster_num=30,celltype="Neuron",
                          avg_standard_cluster = cluster$exp_matrix,memship.wide = cluster$membership)

Sig.Modules <- list()
for(i in 1:length(unique(cluster$cluster$cluster))){
  Sig.Modules[[i]] <- cluster$cluster[cluster$cluster$cluster==i,]$gene
}

# Score each cell for each module using Seurat's built in AddModuleScore
data = AddModuleScore(data, features = Sig.Modules, ctrl.size=5,name = "M")
gene.clusters.names <- colnames(data@meta.data)[grep("M",colnames(data@meta.data))]

# Look at how many cells are present at each stage to choose the appropriate sample size
table(data@meta.data$stage)

# Determine the smallest pval for wilcox.test for every timepoint compared to reference (pre-infection) for each module.
# tps is a userdefined list of time points, whose values are present in data@meta.data $TimePoint. 
# The first time point in tps is used as the reference time point to which all others are compared
clust.wilcox.tests.minp = TestModuleTemporalVariation(meta.data = data@meta.data, sample.size = 150, ntest=1000, tps = c("St48","St54","St59","St66"))
clust.wilcox.tests.minp #print the results

# keep any modules with a minimum pval < 1e-10
Sig.Var.Modules = Sig.Modules[clust.wilcox.tests.minp <= 1e-10]

# remove the original module scores, and add scores for only those modules that were variant in time
data@meta.data = data@meta.data[,-grep("M", colnames(data@meta.data), ignore.case = FALSE)]
data = AddModuleScore(data, features = Sig.Var.Modules, ctrl.size=5, name = "M")
gene.clusters.names.final = grep("M", colnames(data@meta.data), ignore.case = FALSE, value = TRUE)

# rename the modules with the appropraite naming scheme
names(Sig.Var.Modules) = gene.clusters.names.final


# plot box of modules
geneclusts.time.final = lapply(gene.clusters.names.final, function(c_name){
  p = ggplot(data@meta.data, aes_string(x="stage", y=c_name, fill="stage")) +
    geom_boxplot()  + theme(legend.position = "none")
  return(p)
})
plot_grid(plotlist = geneclusts.time.final)

module1 <- TestModuleTimeVariation(type="Module1",gene.clusters.names.final = gene.clusters.names.final,
                                   data=data)
module1
module2 <- TestModuleTimeVariation(type="Module2",gene.clusters.names.final = gene.clusters.names.final,
                                   data=data)
module2
module3 <- TestModuleTimeVariation(type="Module3",gene.clusters.names.final = gene.clusters.names.final,
                                   data=data)
module3
module4 <- TestModuleTimeVariation(type="Module4",gene.clusters.names.final = gene.clusters.names.final,
                                   data=data)
module4

######
module1.gene <- GetModuleGene(module = module1,Sig.Var.Modules = Sig.Var.Modules)
module2.gene <- GetModuleGene(module = module2,Sig.Var.Modules = Sig.Var.Modules)
module3.gene <- GetModuleGene(module = module3,Sig.Var.Modules = Sig.Var.Modules)
module4.gene <- GetModuleGene(module = module4,Sig.Var.Modules = Sig.Var.Modules)

######orthgene
module1.orthgene <- GetOrthgene(gene=module1.gene,celltype="Neuron",type="module1")
module2.orthgene <- GetOrthgene(gene=module2.gene,celltype="Neuron",type="module2")
module3.orthgene <- GetOrthgene(gene=module3.gene,celltype="Neuron",type="module3")
module4.orthgene <- GetOrthgene(gene=module4.gene,celltype="Neuron",type="module4")

###check
Sig.Modules_use <- list()
module.use <- c("module1.gene","module2.gene","module3.gene","module4.gene")
for(i in 1:4){
  Sig.Modules_use[[i]] <- get(module.use[i])
}
data = AddModuleScore(data, features = Sig.Modules_use, ctrl.size=5,name = "Module")
gene.clusters.names.check <- colnames(data@meta.data)[grep("Module",colnames(data@meta.data))]
geneclusts.time.final.check = lapply(gene.clusters.names.check, function(c_name){
  p = ggplot(data@meta.data, aes_string(x="stage", y=c_name, fill="stage")) +
    geom_boxplot()  + theme(legend.position = "none")
  return(p)
})
plot_grid(plotlist = geneclusts.time.final.check)


save.image("./figure6/Neuron/result.RData")
