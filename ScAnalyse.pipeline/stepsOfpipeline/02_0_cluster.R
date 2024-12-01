

#InstallData("ifnb")
#projectname <- "GSE168944"

cluster_marker <- function(projectname) {
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(dplyr)
  #library(future)
  plan("multisession", workers =10) ####多线程运算
  options(future.globals.maxSize = 500000 * 1024^2)
  if (!dir.exists(paste0("../output/02.cluster"))) dir.create(paste0("../output/02.cluster"))
  
  #载入数据
  rawdata <- load(paste0("../output/01.QC/3.",projectname,".QC.Rdata"))
  DimPlot(seuratdata, reduction = "umap",label = TRUE,
          repel = FALSE, pt.size = 0.01,label.size = 3,
          seed = 10, group.by = "condition",raster = FALSE)

  seuratdata <- RunPCA(seuratdata, verbose = FALSE,seed.use = 10)##避免用irlba计算PCA
  # Examine and visualize PCA results a few different ways
  PCAnum <- length(seuratdata@reductions$pca)
  
  write("summary",file = "../output/02.cluster/summary.txt")
  sink(file = "../output/02.cluster/summary.txt",append = TRUE);
  print(c('cellnumber',dim(seuratdata)));
  print(seuratdata[["pca"]], dims = 1:PCAnum, nfeatures = 10);
  sink()
  ####可视化聚类
  p <- DimPlot(seuratdata, reduction = "pca",raster = FALSE)
  p
  ggsave(paste0("../output/02.cluster/01.PCA.dim.pdf"),
         plot=p,width = 15,height =13, scale = 0.5)
  pdf("../output/02.cluster/02.PCA.heatmap.pdf",width = 15,height = 5)
  for (i in 1:c(ceiling(PCAnum/5))) {
    if (5*i< PCAnum | 5*i ==PCAnum) { 
    DimHeatmap(seuratdata, dims = c(5*i-4):c(5*i), nfeatures = 20,ncol = 5,
               cells = 1000, balanced = TRUE)
    } else {
      DimHeatmap(seuratdata, dims = c(5*i-4):PCAnum, nfeatures = 20,ncol = 5,
                 cells = 1000, balanced = TRUE)
    }
  }
  dev.off()
  
  ####确定数据维度
  #主成分累积贡献大于80%
  #PC本身对方差贡献小于5%
  #两个连续PCs之间差异小于0.1%
  PC1 <- seuratdata@reductions$pca@feature.loadings
    # Determine percent of variation associated with each PC
  pct <- seuratdata [["pca"]]@stdev / sum( seuratdata [["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% 
  # and % variation associated with the PC as less than 5
  co1 <- which(cumu > 80 & pct < 5)[1]
  co1
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
              decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  co2
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  pcs
  
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
  
  
  # Elbow plot to visualize 
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  ggsave(paste0("../output/02.cluster/03.Determine_PCnums.pdf"),
         plot=p,width = 7,height =5,scale = 0.5)
  
  #### 其他方法
  #  seuratdata <- JackStraw(seuratdata,num.replicate = 100)
  #  seuratdata <- ScoreJackStraw(seuratdata,dims = 1:20)
  #  JackStrawPlot(seuratdata,dims = 1:20)
  #  p <- ElbowPlot(seuratdata,ndims = PCAnum)
  #  p$data
  
  best.num <- pcs ##观察上图获得最佳聚类数目
  
  ####细胞聚类  
  seuratdata <- RunHarmony(seuratdata,"orig.ident")
  seuratdata <- RunUMAP(seuratdata, 
                        reduction = "harmony",
                        dims = 1:best.num,
                        seed.use = 10)
  seuratdata <- FindNeighbors(seuratdata,
                              reduction = "harmony",
                              dims = 1:best.num)
  seuratdata <- FindClusters(seuratdata, 
                             resolution = seq(0.1,1.0,0.2),#清晰度越多，越占内存，谨慎选择
                             random.seed = 10)
  seuratdata$seurat_clusters <-seuratdata$SCT_snn_res.0.5
  cell_cluster <- seuratdata$seurat_clusters
  cell_ID <- names(seuratdata$seurat_clusters)
  cell_cluster <- data.frame(cell_ID=cell_ID,cell_cluster=cell_cluster)
  write.csv(cell_cluster, "../output/02.cluster/cell_cluster.csv",row.names = FALSE)
 
  if ("Cluster" %in% colnames(seuratdata@meta.data) ) {
    Idents(seuratdata) <- seuratdata$Cluster
  } else {
    Idents(seuratdata) <- seuratdata$seurat_clusters
  }
  embed_umap <- Embeddings(seuratdata,"umap") #提取umap图坐标
  write.csv(embed_umap,"../output/02.cluster/embed_umap.csv") 
  # Visualization
  pdf("../output/02.cluster/04.Visualization.cluster.pdf",width = 6,height = 6)
  DimPlot(seuratdata, reduction = "umap",label = TRUE,
          repel = FALSE, pt.size = 0.01,label.size = 3,
          seed = 10, group.by = "condition",raster = FALSE)
  DimPlot(seuratdata, reduction = "umap", label = TRUE,
          label.size = 3, 
          repel = FALSE,raster = FALSE)
  dev.off()
  
  p3 <- DimPlot(seuratdata, reduction = "umap", 
                split.by = "condition",ncol = 2,raster = FALSE)
  p3
  ggsave(paste0("../output/02.cluster/05.Visualization.cluster.by.ident.pdf"),
         plot=p3,width = 8,height =4)
  
  sink(file = "../output/02.cluster/summary.txt",append = TRUE);
  print("cluster_cells" );
  print(table(seuratdata$seurat_clusters));
  #print(table(seuratdata@active.ident))
  print("condition_cells" );
  print(table(seuratdata$condition));
  print(best.num)#12
  sink()
  #sorted
  #seuratdata$orig.ident %>% table() %>% names() 
  #%>% paste0(collapse = "','")#嫌弃麻烦，做得工具
  seuratdata$orig.ident <- factor(seuratdata$orig.ident,
                                  levels = c('HC1_B','HC2_B','HC3_B',
                                             'DLE1_B','DLE2_B','DLE3_B','DLE4_B',
                                             'SLE1_B','SLE2_B','SLE3_B','SLE4_B'))
  seuratdata$condition <- factor(seuratdata$condition,levels = c('HC','DLE','SLE'))
  save(seuratdata,file=paste0("../output/01.QC/4.",projectname,".inter.cluster.Rdata") )
  # Identify the 10 most highly variable genes
  top10_frist <- head(VariableFeatures(seuratdata), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seuratdata)
  plot2 <- LabelPoints(plot = plot1, points = top10_frist, repel = T)
  plot2
  ggsave(paste0("../output/02.cluster/06.top10.VF.pdf"),plot=plot2,width = 11,height =12,
         scale = 0.5)
 
  # find markers for every cluster compared to all remaining cells, report only the positive
  # ones
  
  seuratdata <- ScaleData(seuratdata, verbose = FALSE)
  seuratdata.markers <- FindAllMarkers(seuratdata, only.pos = TRUE,
                                       min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(seuratdata.markers,"../output/02.cluster/cluster.markers.csv",row.names = FALSE) 
 
  seuratdata.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  write.csv(top10,"../output/02.cluster/cluster.markers.top10.csv",row.names = FALSE) 
  
  p_top10 <- DoHeatmap(subset(seuratdata,downsample=100), 
                       features = top10$gene) + NoLegend()
  p_top10
  ggsave(paste0("../output/02.cluster/06.Marker.heatmap.pdf"),
         plot=p_top10,width = 6,height =12)
  
  seuratdata.markers <- read.csv("../output/02.cluster/cluster.markers.csv",header = TRUE)
  cluster_marker_row <- data.frame(cluster=0:(length(unique(seuratdata$seurat_clusters))-1),markers=NA)
  for (c in cluster_marker_row$cluster) {
    temp_c = seuratdata.markers$gene[seuratdata.markers$cluster == c]
    temp_c = paste0(temp_c,collapse = ",")
    cluster_marker_row$markers[which(cluster_marker_row$cluster== c)] = temp_c
  }  
  write.csv(cluster_marker_row,"../output/02.cluster/cluster.markers_row.csv",row.names = FALSE) 
  save(seuratdata,seuratdata.markers,file=paste0("../output/01.QC/5.",projectname,".cluster.Rdata") )
}


