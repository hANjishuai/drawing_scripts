cellanno <- function(projectname,species) {
    #注释细胞类型
    library(SingleR)
    library(celldex)
    library(ggplot2)
    library(Seurat)
    plan("multisession", workers =20) ####多线程运算
    options(future.globals.maxSize = 5000 * 1024^2)
    if (!dir.exists(paste0("../output/03.celltype"))) dir.create(paste0("../output/03.celltype"))
    if (species == "human") {
        #immgen<- ImmGenData()
        cellmarker <- celldex::HumanPrimaryCellAtlasData()
      }
    if (species == "mouse") {
        cellmarker<- celldex::ImmGenData()
      }
    clusterdata <- load(paste0("../output/01.QC/5.",projectname,".cluster.Rdata"))
    singler.results <- SingleR(test = GetAssayData(seuratdata,assay="RNA"),
                               clusters = Idents(seuratdata),
                               ref=cellmarker,assay.type.test = 1,
                               labels = cellmarker$label.main)
  
  #手动矫正
  seuratdata[["SingleR.cluster.labels"]] <- 
    singler.results$labels[match(seuratdata[[]][["seurat_clusters"]], rownames(singler.results))]
  cellplot <-  DimPlot (seuratdata, reduction = "umap", group.by="SingleR.cluster.labels",
                        label = TRUE, pt.size = 0.5) 
  #+ NoLegend()
  cellplot
  ggsave(paste0("../output/03.celltype/1.cellplot.pdf"),plot=cellplot,width = 15,height =15)
  #提取细胞类型
  cellcluster <- cbind(rownames(singler.results),singler.results$labels) 
  write.csv(cellcluster,paste0("../output/03.celltype/0.singleR_cellplot.csv"),quote = TRUE,row.names = TRUE)
  
  save(seuratdata,file=paste0("../output/01.QC/6.",projectname,".celltype.Rdata") )
}

cellanno.check <- function(projectname,newcluster) {
  library(Seurat)
  library(ggplot2)
  celldata <- load(paste0("../output/01.QC/6.",projectname,".celltype.Rdata"))
  #自己去修改0.singleR_cellplot.csv文件，添加manual_level1列后，再往下
  cellcluster <- read.csv(newcluster,header = TRUE)
  #newcluster<-"../output/03.celltype/0.singleR_cellplot.csv"
  ####定义新的亚群
  ####
  # level1
  newcluster = cellcluster$manual_level1
  names(newcluster) <- cellcluster$cluster
  Idents(seuratdata) <- seuratdata$seurat_clusters
  seuratdata <- RenameIdents(seuratdata, newcluster)
  seuratdata$manual_L1 <- Idents(seuratdata)
  cellplot0 <-  DimPlot (seuratdata, reduction = "umap", 
                         label = TRUE, pt.size = 0.5)
  ggsave(paste0("../output/03.celltype/2.cellplot.","manual_L1",".pdf"),plot=cellplot0,width = 15,height =10)
  if (F) {
    # level2
    newcluster = cellcluster$manual_level2
    names(newcluster) <- cellcluster$cluster
    Idents(seuratdata) <- seuratdata$seurat_clusters
    seuratdata <- RenameIdents(seuratdata, newcluster)
    seuratdata$manual_L2 <- Idents(seuratdata)
    cellplot0 <-  DimPlot (seuratdata, reduction = "umap", 
                           label = TRUE, pt.size = 0.5)
    ggsave(paste0("../output/03.celltype/2.cellplot.","manual_L2",".pdf"),plot=cellplot0,width = 15,height =10)
    
    # level3
    newcluster = cellcluster$manual_level3
    names(newcluster) <- cellcluster$cluster
    Idents(seuratdata) <- seuratdata$seurat_clusters
    seuratdata <- RenameIdents(seuratdata, newcluster)
    seuratdata$manual_L3 <- Idents(seuratdata)
    cellplot0 <-  DimPlot (seuratdata, reduction = "umap", 
                           label = TRUE, pt.size = 0.5)
    ggsave(paste0("../output/03.celltype/3.cellplot.","manual_L3",".pdf"),plot=cellplot0,width = 15,height =10)
    cellmate <- seuratdata@meta.data
    write.csv(cellmate,"../output/03.celltype/cellmate.csv",quote = TRUE,row.names = FALSE)
    
  }
  save(seuratdata,file=paste0("../output/01.QC/7.",projectname,".celltype.check.Rdata") )
}


