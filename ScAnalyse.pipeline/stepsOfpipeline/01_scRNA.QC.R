scRNA.qc <- function(projectname,species){
  library(Seurat)
  library(ggplot2)
  library(glmGamPoi)
  options(future.globals.maxSize = 50000 * 1024^2)  # 设置为 1000 MB
    
  if (!dir.exists(paste0("../output/01.QC"))) dir.create(paste0("../output/01.QC"))
  rawdata <- load(paste0("../output/01.QC/2.",projectname,".raw.Rdata"))
  ####切换数据集
  DefaultAssay(seuratdata) <- "RNA"
  
  #####检查质量
  if (species=="human"){
    seuratdata[["percent.mt"]] <- PercentageFeatureSet(seuratdata, pattern = "^MT-")
    seuratdata[["percent.ribo"]] <- PercentageFeatureSet(seuratdata, pattern = "^RP[SL]")
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
  Idents(seuratdata) <- seuratdata$orig.ident
  
  p1 <- VlnPlot(seuratdata, alpha = 0.5,
                features = c("percent.ribo", "percent.mt",
                             "nFeature_RNA", "nCount_RNA"),
                group.by =  "condition",
                ncol = 4,pt.size = 0.00001)
  p1
  ggsave(paste0("../output/01.QC/0.raw.vlnplot.pdf"),plot=p1,width = 10,height =6)

  plot1 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2
  p
  ggsave(paste0("../output/01.QC/1.farture.ncount.mt.scatter.pdf"),plot=p,width = 10,height =5)
  
  #### 过滤低质量细胞
  seuratdata.sub <- subset(seuratdata,
                       ####过滤测序深度基因数目>200（低质量细胞），
                       ####基因数目小于6000（doublets细胞）
                       ####线粒体比例小于10（破损细胞）
                       ####
                       subset = nFeature_RNA > 50 & nFeature_RNA < 6000 & percent.mt < 12.5
                       )
  ####质控后数据
  p <- VlnPlot(seuratdata.sub, 
               features = c("percent.ribo", "percent.mt",
                            "nFeature_RNA", "nCount_RNA"), 
               group.by =  "condition",
               ncol = 4,pt.size = 0.01)
  p
  ggsave(paste0("../output/01.QC/1.QC.vlnplot.pdf"),plot=p,width = 20,height =6)
  
  plot1 <- FeatureScatter(seuratdata.sub, 
                          feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuratdata.sub,
                          feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p <- plot1 + plot2
  p
  ggsave(paste0("../output/01.QC/1.QC.farture.ncount.mt.scatter.pdf"),plot=p,width = 10,height =5)
  
  ####添加文件
  
  write("cellnumbers",file = "../output/01.QC/cellnums.txt")
  sink(file = "../output/01.QC/cellnums.txt",append = T);
  print("rawdata")
  print(table(Idents(seuratdata)));
  print("QCdata")
  print(table(Idents(seuratdata.sub)));
  sink()
  #### 标准化数据
if(T){
  cellnum <- dim(seuratdata)
  
  seuratdata.sub <- SCTransform(seuratdata.sub,
                                variable.features.n=cellnum[1])
  seuratdata.sub <- SCTransform(seuratdata.sub, 
                                vars.to.regress = "percent.mt",
                                verbose = FALSE)
  
  seuratdata=seuratdata.sub
 
  save(seuratdata,file=paste0("../output/01.QC/3.",projectname,".QC.Rdata") )
}
  #save(seuratdata.sub,file=paste0("../output/01.QC/3.",projectname,".QC_failtosct.Rdata") )
}

###数据太大了，内存不够（256G都不够！！），只能去scanpy了
if(F){
  library(Seurat)
  library(sceasy)
  library(reticulate)
  load("../output/01.QC/3.B1.QC_failtosct.Rdata")
  # 在R语言中加载python环境
  use_condaenv('Rstudio-server')
  loompy <- reticulate::import('loompy')
  # Seurat to AnnData
  seuratdata.sub[["RNA"]] <- as(seuratdata.sub[["RNA"]], "Assay")
  sceasy::convertFormat(seuratdata.sub, from="seurat", to="anndata",
                        outFile='../output/01.QC/3.B1.QC_failtosct_V5.h5ad')
  #原文链接：https://blog.csdn.net/zfyyzhys/article/details/142603822
  #ha5d格式数据转化成seruat对象
  library(sceasy)
  library(reticulate)
  library(Seurat)
  library(BiocParallel)
  register(MulticoreParam(workers = 10, progressbar = TRUE))
  # h5ad转为Seurat
  sceasy::convertFormat(obj = "../output/02.cluster/HC_DLE_SLE_SKIN.h5ad", 
                        from="anndata",to="seurat",
                        outFile = '3.B1.Scanpy_2_seurat.rds')
  #这种方法得到的数据是SeruatV4版本的，所以如果要用于SeruatV5的话还需要再转化一下。
  #还有细胞数很多的话sceasy就不好用了，这个时候可以用dior包。
  
}







