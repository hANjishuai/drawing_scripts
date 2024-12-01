#01.差异分析####
find_trans_Profile <- function() {
  #A:Env_setting
  cell<-"Bcell"
  PutIn <- '../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata'
  AnnotationType <- "manual_L2"
  
  #B:所需包以及函数载入
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(Seurat)
  library(ggpubr)
  library(tidyverse)
  theme_set(theme_cowplot())
  
  library(ComplexHeatmap)
  library(circlize)
  library(paletteer)
  library(tidyr)
  library(cowplot)
  library(ggplotify)
  library(tibble)
  
  #C:正式pipline
  if (!dir.exists(paste0("../output/04.subcell/06.DEG/"))) dir.create(paste0("../output/04.subcell/06.DEG/"))
  if (!dir.exists(paste0("../output/04.subcell/06.DEG/",cell))) dir.create(paste0("../output/04.subcell/06.DEG/",cell))
  #load
  celldata=load(PutIn)
  #选定亚群注释集
  Idents(seuratdata) <- AnnotationType
  #数据实验选择“RNA”，正态和标准化，差异分析
  DefaultAssay(seuratdata) <- "RNA"
  all.genes<-rownames(seuratdata)
  subcell <- seuratdata %>% NormalizeData() %>% 
    ScaleData(.,features = all.genes)
  all.marker.RNA <-FindAllMarkers(subcell,assay = "RNA",only.pos = F,
                                  min.pct = 0.5,logfc.threshold = 0.5)
  write.csv(all.marker.RNA,
            file = paste0("../output/04.subcell/06.DEG/Bcell/",cell,"_marker.csv"))
  #筛选差异基因
  filter_marker.RNA <- all.marker.RNA %>% 
    mutate(.,status=ifelse(all.marker.RNA$avg_log2FC>0,"top","down")) %>%
    filter(p_val_adj<0.05) %>% group_by(cluster)%>%arrange(cluster,desc(avg_log2FC))
  write.csv(filter_marker.RNA,
            file = paste0("../output/04.subcell/06.DEG/Bcell/","all_marker.csv"))
  
  #筛选前10的基因,并画热图
  ##提取前10基因
  top_10<- all.marker.RNA %>% 
    mutate(.,status=ifelse(all.marker.RNA$avg_log2FC>0,"top","down")) %>%
    filter(p_val_adj<0.05) %>% group_by(cluster)%>%top_n(10,wt = avg_log2FC)
  gene_10 <- top_10[,c('gene','cluster')] 
  gene_10 <- gene_10[which(!duplicated(gene_10$gene)),]

  ##取出前10基因表达矩阵
  express <- LayerData(seuratdata,assay = 'SCT',layer = 'data') %>% data.frame()
  express <- subset(express,subset = rownames(express) %in% gene_10$gene) 
  express <- rownames_to_column(express,var = 'gene')
  express <- express[match(gene_10$gene,express$gene),]
  rownames(express) <- express$gene
  express <- express[,-1] %>% t()
  ##添加细胞注释信息
  meta <- seuratdata@meta.data
  express <- cbind(data.frame(meta$manual_L2),express)
  express2 <- aggregate(express[,-1],by = list(express$meta.manual_L2),mean)
  express2 <- column_to_rownames(express2,var = 'Group.1')-1.25
  matrix <- t(express2)
  matrix <- cbind(data.frame(gene_10$cluster),matrix)
  #热图
  ##设置颜色：
  col_fun <- colorRamp2(c(-1.5,0,1.5),c("#e8f6e5","#fca069","#d43627"))
  pdf("../output/04.subcell/06.DEG/Bcell/TOP10_Marker_gene.pdf",
      width = 2,
      height = 6.5)
  Heatmap(matrix[,-1],
          cluster_columns = F,
          cluster_rows = F,
          show_row_names = T,
          row_names_gp = gpar(fontface="bold",fontsize=3.5),
          col = col_fun,
          show_heatmap_legend = F,
          border = F,
          gap = unit(1, "mm"),
          rect_gp = gpar(col="white",lwd=1),
          column_names_gp = gpar(fontface="bold",fontsize=5.5),
          column_names_centered = F,
          column_names_rot = 90,
          left_annotation = c(),
          #split = matrix$gene_10.cluster,
          row_title_rot = 90,
          column_title_rot = 0,
          row_title_gp = gpar(fontsize=3,fontface="bold"))
  dev.off()
  }

