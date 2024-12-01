#首先我们跑一遍标准的Seurat_V4流程，获取单细胞转录组分析结果seurat对象

run_Seurat <- function(input_type = c("dataframe", "10X","h5"),#输入数据类型，这里只提供了3
                       file_path,#储存目标文件的文件夹，最好setwd到这个目录，文件名命名为样本名称
                       species=c('human','mouse'),#物种
                       percent.mt,#线粒体基因过滤比例
                       percent.rb,#核糖体基因过滤比例，如果不需要过滤，这个数值设置为100
                       min_Counts,#最低过滤counts
                       max_Counts,#最高过滤counts
                       min_feature,#最低过滤基因
                       max_feature,#最高低过滤基因
                       intergrated.method = c("CCA","Harmony"),#数据整合去批次方式
                       vars.to.regress,#需要回归的参数（ScaleData那里）
                       dims,#选择降维的PC数
                       resolution,#降维分辨率
                       reduction=c("UMAP","TSNE")#UMAP或者TSNE降维
){
  
  folders=list.files(file_path)
  
  if(input_type=="dataframe"){
    files = list.files()[1]
    files_sufix = substring(files,nchar(files)-2)
    
    
    if(files_sufix == "csv"){
      
      names = gsub(".csv",'',folders)
      matrxi_list = lapply(folders,function(x){ read.csv(x, header = T, row.names = 1)})
      names(matrxi_list) <- names
      
    }else{
      
      names = gsub(".txt",'',folders)
      matrxi_list = lapply(folders,function(x){ read.table(x,header = TRUE, row.names = 1) })
      names(matrxi_list) <- names
    }
    
    
    SeuratObj = list()
    
    for (i in 1:length(matrxi_list)) {
      
      Obj = CreateSeuratObject(counts = matrxi_list[[i]],  min.cells=3, min.features = 200)
      SeuratObj[[i]] <- Obj
      names(SeuratObj)[i] <- names[i]
    }
    
  }
  
  if(input_type=="10X"){
    
    
    SeuratObj = lapply(folders,function(x){ 
      CreateSeuratObject(counts = Read10X(x), project = x, min.cells = 3, min.features = 200)
    })
    
    
  }
  
  if(input_type=="h5"){
    
    names = gsub(".h5",'',folders)
    matrxi_list = lapply(folders,function(x){Read10X_h5(x)})
    names(matrxi_list) <- names
    
    SeuratObj = list()
    for (i in 1:length(matrxi_list)) {
      
      Obj = CreateSeuratObject(counts = matrxi_list[[i]], project = names(matrxi_list)[i], min.cells=3, min.features = 200)
      SeuratObj[[i]] <- Obj
      
    }
    
  }
  
  
  if(species=='human'){
    
    mito.pattern <- "^MT-"
    ribo.pattern <- "^RBS|RPL"
    hb.pattern <- "^HB[^(P)]"
    
    s.genes <- Seurat::cc.genes.updated.2019$s.genes
    g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
    
  }else{
    
    mito.pattern <- "^mt-"
    ribo.pattern <- "^Rbs|Rpl"
    hb.pattern <- "^Hb[^(p)]"
    s.genes <- c("Mcm4","Exo1","Slbp","Gmnn","Cdc45","Msh2","Mcm6","Rrm2","Pold3","Blm", "Ubr7","Mcm5","Clspn","Hells","Nasp","Rpa2","Rad51ap1","Tyms","Rrm1","Rfc2","Prim1",   
                 "Brip1","Usp1","Ung", "Pola1","Mcm2","Fen1","Tipin","Pcna","Cdca7","Uhrf1","Casp8ap2","Cdc6","Dscc1","Wdr76","E2f8","Dtl", "Ccne2","Atad2","Gins2","Chaf1b","Pcna-ps2")
    g2m.genes <- c("Nuf2","Psrc1","Ncapd2","Ccnb2","Smc4","Lbr", "Tacc3","Cenpa","Kif23","Cdca2","Anp32e","G2e3","Cdca3","Anln","Cenpe","Gas2l3","Tubb4b","Cenpf","Dlgap5","Hjurp","Cks1brt", "Gtse1","Bub1","Birc5",
                   "Ube2c","Rangap1","Hmmr","Ect2","Tpx2","Ckap5","Cbx5","Nek2","Ttk", "Cdca8","Nusap1","Ctcf","Cdc20","Cks2","Mki67","Tmpo","Ckap2l","Aurkb","Kif2c","Cdk1","Kif20b","Top2a","Aurka","Ckap2","Hmgb2","Cdc25c","Ndc80","Kif11" )
    
  }
  
  
  
  for(x in seq_along(SeuratObj)) {
    SeuratObj[[x]] <- CellCycleScoring(SeuratObj[[x]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    SeuratObj[[x]][["percent.mt"]] <- PercentageFeatureSet(SeuratObj[[x]], pattern = mito.pattern)
    SeuratObj[[x]][["percent.rb"]] <- PercentageFeatureSet(SeuratObj[[x]], pattern = ribo.pattern)
    SeuratObj[[x]][["percent.hb"]] <- PercentageFeatureSet(SeuratObj[[x]], pattern = hb.pattern)
  }
  
  preQC <- list()
  for(x in seq_along(SeuratObj)){
    
    p = VlnPlot(SeuratObj[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"), ncol = 5, group.by = "orig.ident", pt.size = 0)
    preQC[[x]] <- p
  }
  
  Pre_QC = Seurat::CombinePlots(preQC,ncol = 1)
  high = length(SeuratObj)
  ggsave(Pre_QC,file="Pre_QC.pdf",width=8,height=4*high)
  
  
  for(x in seq_along(SeuratObj)){
    
    a=min_Counts 
    b=max_Counts 
    c=min_feature 
    d=max_feature 
    e=percent.mt
    f=percent.rb
    
    SeuratObj[[x]] <- subset(SeuratObj[[x]], 
                             subset = nCount_RNA > a & 
                               nCount_RNA < b &
                               nFeature_RNA > c &
                               nFeature_RNA < d &
                               percent.mt < e &
                               percent.hb < 1 &
                               percent.rb <f)
  }
  
  
  postQC <- list()
  for(x in seq_along(SeuratObj)){
    
    p = VlnPlot(SeuratObj[[x]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hb"), ncol = 5, group.by = "orig.ident", pt.size = 0)
    postQC[[x]] <- p
  }
  
  Post_QC = Seurat::CombinePlots(postQC,ncol = 1)
  high = length(SeuratObj)
  ggsave(Post_QC,file="Post_QC.pdf",width=8,height=4*high)
  
  
  SeuratObj <- lapply(X = SeuratObj, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
  })
  
  
  
  if (intergrated.method=="CCA") {
    
    SeuratObj <- FindIntegrationAnchors(object.list = SeuratObj, dims = 1:50)
    SeuratObj_combine <- IntegrateData(anchorset = SeuratObj, dims = 1:50)
    
    saveRDS(SeuratObj_combine, file='./SeuratObj_cca_intergrate.rds')
    
    SeuratObj_combine <- ScaleData(SeuratObj_combine, vars.to.regress = vars.to.regress, verbose = T)
    SeuratObj_combine <- RunPCA(SeuratObj_combine, npcs = 100, verbose = FALSE)
    
    PCA_plot <- Seurat::ElbowPlot(SeuratObj_combine, ndims=100)
    ggsave(PCA_plot,file="PCA_plot.pdf",width=15,height=12)
    
    
    if(reduction=="UMAP"){
      
      SeuratObj_combine <- RunUMAP(SeuratObj_combine, reduction = "pca", dims = 1:dims)
      
    }else{
      
      SeuratObj_combine <- RunTSNE(SeuratObj_combine, reduction = "pca", dims = 1:dims)
      
    }
    
    
    SeuratObj_combine <- FindNeighbors(SeuratObj_combine, reduction = "pca", dims = 1:dims)
    SeuratObj_combine <- FindClusters(SeuratObj_combine,resolution = resolution)
    
    return(SeuratObj_combine)
    
    
  }
  
  if(intergrated.method == "Harmony"){
    
    SeuratObj_combine = merge(SeuratObj[[1]], y= SeuratObj[-1], add.cell.ids = names(SeuratObj))
    
    SeuratObj_combine <- ScaleData(SeuratObj_combine, vars.to.regress = vars.to.regress, verbose = T)
    SeuratObj_combine <- FindVariableFeatures(SeuratObj_combine, selection.method = "vst", nfeatures = 4000)
    SeuratObj_combine <- RunPCA(SeuratObj_combine, npcs = 100, verbose = FALSE)
    
    PCA_plot <- Seurat::ElbowPlot(SeuratObj_combine, ndims=100)
    ggsave(PCA_plot,file="PCA_plot.pdf",width=15,height=12)
    
    SeuratObj_combine <- RunHarmony(SeuratObj_combine, "orig.ident")
    
    if(reduction=="UMAP"){
      
      SeuratObj_combine <- RunUMAP(SeuratObj_combine, reduction = "harmony", dims = 1:dims)
      
    }else{
      
      SeuratObj_combine <- RunTSNE(SeuratObj_combine, reduction = "harmony", dims = 1:dims)
      
    }
    
    
    SeuratObj_combine <- FindNeighbors(SeuratObj_combine, reduction = "pca", dims = 1:dims)
    SeuratObj_combine <- FindClusters(SeuratObj_combine,resolution = resolution)
    
    return(SeuratObj_combine)
    
    
  }
  
  
}
