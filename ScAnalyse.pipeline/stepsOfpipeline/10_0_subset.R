subcell.qc <- function(species){
  library(Seurat)
  library(ggplot2)
  library(glmGamPoi)
  library(harmony)
  options(future.globals.maxSize = 50000 * 1024^2)  # 设置为 1000 MB
  
  if (!dir.exists(paste0("../output/04.subcell/"))) 
    dir.create(paste0("../output/04.subcell/"))
  rawdata <- load(paste0("../output/01.QC/subset_","B_cell.Rdata"))
  metaData <- subcells@meta.data
  #metaData$orig.ident <- gsub("_.*$","",rownames(metaData))
  subcells@meta.data <-metaData
  subcells <- LayerData(subcells,assay = "RNA",layer = "counts")
  #最好用这个layerdata函数提取，不然会失去细胞名和基因名
  seuratdata <- CreateSeuratObject(counts = subcells,
                     meta.data = metaData[,c(1,8,19)],
                     project = "Bcell" )
  colnames(seuratdata) 
  rownames(seuratdata)
  ####切换数据集
  DefaultAssay(seuratdata) <- "RNA"
  #####检查质量
  if (species=="human"){
    seuratdata[["percent.mt"]] <- PercentageFeatureSet(seuratdata, pattern = "^MT-")
    ####计算核糖体比例
    rb.genes <- rownames(seuratdata)[grep("^RP[SL]",rownames(seuratdata))]
    C <- seuratdata[["RNA"]]@layers$counts
    rownames(C) <- rownames(seuratdata)
    colnames(C) <- colnames(seuratdata)
    percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
    seuratdata <- AddMetaData(seuratdata, percent.ribo, col.name = "percent.ribo")
  }
  if (species=="mouse"){
    #线粒体基因
    #rownames(seuratdata)[grepl('^MT-',rownames(seuratdata),ignore.case = T)]
    #核糖体蛋白基因
    #rownames(seuratdata)[grepl('^Rp[sl]',rownames(seuratdata),ignore.case = T)]
    ##计算上述两种基因在所有细胞中的比例
    seuratdata[["percent.mt"]] <- PercentageFeatureSet(seuratdata, pattern = "^mt-")
    ####计算核糖体比例
    rb.genes <- rownames(seuratdata)[grep("^Rp[sl]",rownames(seuratdata))]
    C <- as.data.frame(seuratdata[["RNA"]]@layers$counts)
    rownames(C) <- rownames(seuratdata)
    colnames(C) <- colnames(seuratdata)
    percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
    seuratdata <- AddMetaData(seuratdata, percent.ribo, col.name = "percent.ribo")
  }
  if (!dir.exists(paste0("../output/04.subcell/01.QC/"))) 
    dir.create(paste0("../output/04.subcell/01.QC/"))
  p1 <- VlnPlot(seuratdata, alpha = 0.5,
                features = c("percent.ribo", "percent.mt",
                             "nFeature_RNA", "nCount_RNA"),
                group.by =  "condition",
                ncol = 4,pt.size = 0.01)
  p1
  ggsave(paste0("../output/04.subcell/01.QC/0.raw.vlnplot.pdf"),plot=p1,width = 10,height =6)
  plot1 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2
  p
  ggsave(paste0("../output/04.subcell/01.QC/1.farture.ncount.mt.scatter.pdf"),plot=p,width = 10,height =5)
  
  #### 过滤低质量细胞
  seuratdata.sub <- subset(seuratdata,
                           ####过滤测序深度基因数目>200（低质量细胞），
                           ####基因数目小于6000（doublets细胞）
                           ####线粒体比例小于10（破损细胞）
                           ####
                           subset = nFeature_RNA > 50 & nFeature_RNA < 6000 & percent.mt < 8
  )
  ####质控后数据
  p <- VlnPlot(seuratdata.sub, 
               features = c("percent.ribo", "percent.mt",
                            "nFeature_RNA", "nCount_RNA"), 
               group.by =  "condition",
               ncol = 4,pt.size = 0.01)
  p
  ggsave(paste0("../output/04.subcell/01.QC/1.QC.vlnplot.pdf"),plot=p,width = 20,height =6)
  
  plot1 <- FeatureScatter(seuratdata.sub, 
                          feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuratdata.sub,
                          feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2
  p
  ggsave(paste0("../output/04.subcell/01.QC/1.QC.farture.ncount.mt.scatter.pdf"),plot=p,width = 10,height =5)
  
  ####添加文件
  
  write("cellnumbers",file = "../output/04.subcell/01.QC/cellnums.txt")
  sink(file = "../output/04.subcell/01.QC/cellnums.txt",append = T);
  print("rawdata")
  print(table(Idents(seuratdata)));
  print("QCdata")
  print(table(Idents(seuratdata.sub)));
  sink()
  #### 标准化数据
  # store mitochondrial percentage in object meta data
  seuratdata <- PercentageFeatureSet(seuratdata, pattern = "^MT-", col.name = "percent.mt")
  
  # run sctransform
  seuratdata  <- SCTransform(seuratdata , vars.to.regress = "percent.mt", verbose = FALSE)
  
  # These are now standard steps in the Seurat workflow for visualization and clustering
  seuratdata  <- RunPCA(seuratdata , verbose = FALSE)
  save(seuratdata,file=paste0("../output/04.subcell/01.QC/3.",projectname,".QC.Rdata") )
  ElbowPlot(seuratdata)#15
  seuratdata <- RunHarmony(seuratdata,"orig.ident")
  seuratdata  <- RunUMAP(seuratdata , 
                         reduction = "harmony",
                         dims = 1:15,
                         seed.use = 10, 
                         verbose = FALSE)
  seuratdata  <- FindNeighbors(seuratdata , 
                               reduction = "harmony",
                               dims = 1:15, 
                               verbose = FALSE)
  seuratdata  <- FindClusters(seuratdata ,
                              resolution = seq(0,1,0.1),
                              random.seed = 10,
                              verbose = FALSE)
  DimPlot(seuratdata , label = TRUE)
  
  if (!dir.exists("../output/04.subcell/02.cluster/")) {
    dir.create("../output/04.subcell/02.cluster/")
  }
  
  save(seuratdata,file = "../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata")
}