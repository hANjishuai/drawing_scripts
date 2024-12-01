subcellanno <- function(newcluster){
  library(Seurat)
  library(ggplot2)
  if (!dir.exists("../output/04.subcell/03.annotation/")) {
    dir.create("../output/04.subcell/03.annotation/")
  }
  celldata <- load(paste0("../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata"))
  cellcluster <- read.csv(newcluster,header = TRUE)
  
  # level2
  #取出注释好的csv文件
  newcluster = cellcluster$manual_level2
  #对应注释时的清晰度
  names(newcluster) <- cellcluster$cluster
  #seurat对象同步调频
  seuratdata$seurat_clusters <- seuratdata$SCT_snn_res.0.2
  Idents(seuratdata) <- seuratdata$seurat_clusters
  #增加注释
  seuratdata <- RenameIdents(seuratdata, newcluster)
  seuratdata$manual_L2 <- Idents(seuratdata) %>% factor(.,levels=levels(seuratdata))
  cellplot0 <-  DimPlot (seuratdata, reduction = "umap", 
                         label = TRUE, pt.size = 0.5)
  ggsave(paste0("../output/04.subcell/03.annotation/2.cellplot.",
                "manual_L2",".pdf"),plot=cellplot0,width = 15,height =10)
  save(seuratdata,file = paste0("../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata"))
  
  if(!dir.exists("../output/04.subcell/04.pyscenic")){
    dir.create("../output/04.subcell/04.pyscenic")
  }#为后面pyscenic分析准备
  #1、矩阵
  pyscenic_counts=LayerData(seuratdata,assay = "RNA",layer = 'counts')
  write.csv(t(as.matrix(pyscenic_counts)),
            file = paste0("../output/04.subcell/04.pyscenic/","sce_exp.csv"))
  #2、meta数据
  cellInfo <- seuratdata@meta.data[,c("manual_L2","condition","nCount_RNA","nFeature_RNA")]
  colnames(cellInfo) <- c('CellType','condition','nGene' ,'nUMI')
  write.csv(cellInfo, 
            file = paste0("../output/04.subcell/04.pyscenic/","cellInfo.csv"))
  
}