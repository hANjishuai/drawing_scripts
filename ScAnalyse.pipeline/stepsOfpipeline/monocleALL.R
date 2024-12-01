#version 2.0.0 2024/09/29
#history
# 1.更改seurat对象内部结构，更加适合seurat v5
# 2.增加format参数，确定随着差异时序基因变化的细胞类型
# 3.自动下载所需的包
# 4.不再自动画图
if (!require(AnnotationDbi)) {
  BiocManager::install("AnnotationDbi")
}
if (!require(monocle)) {
  BiocManager::install("monocle")
}
if (!require(ggsci)) {
  BiocManager::install("ggsci")
}
library(AnnotationDbi)
library(monocle)
library(ggsci)
if (!dir.exists("../output/04.subcell/08.monocle")) {
  dir.create("../output/04.subcell/08.monocle")}
monocleAll <- function(seu=NULL,format=NA){
  DefaultAssay(seu) <- "RNA"   #拟时序需要用rnaslot
  SeuratObject <- seu
  #提取转录组counts信息--表达矩阵
  #  expr_matrix <-as(as.matrix(SeuratObject@assays$RNA@counts),"sparseMatrix")
  expr_matrix <-as(as.matrix(SeuratObject@assays$RNA@layers$counts),"sparseMatrix")
  #提取表型信息--pdata
  p_data <- SeuratObject@meta.data
  #p_data$celltype <- SeuratObject@active.ident#整合每个细胞的细胞鉴定信息到p_data,若有,可跳
  #提取基因信息
  f_data <- data.frame(gene_short_name = row.names(SeuratObject),row.names = row.names(SeuratObject))
  #构建CDS对象
  pd <- new('AnnotatedDataFrame', data = p_data)
  fd <- new('AnnotatedDataFrame', data = f_data)
  cds <- newCellDataSet(expr_matrix,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
  #估计size factor和离散度
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  #过滤低表达的基因         
  cds <- detectGenes(cds,min_expr = 0.1)
  head(fData(cds))
  expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
  
  #dpFeature获取拟时序高变基因
  diff <- differentialGeneTest(cds[expressed_genes,],
                               fullModelFormulaStr = format,
                               cores=1)
  save()
  #从算法上来说,~后可以是pdata中的任意一列
  deg <- subset(diff,qval<0.01)
  deg <- deg[order(deg$qval,decreasing = F),]
  write.table(deg,file=paste0("train.monocle.DEG_",format(Sys.time(),"%H_%M_%S"),".xls"),col.name=T,row.names = F,sep = "\t",quote=F)  
  #得到轨迹构建基因并可视化                           
  ordergene <- rownames(deg)
  cds <- setOrderingFilter(cds,ordergene)
  #可通过 cds@featureData@data$use_for_ordering,table()后查询
  #可视化轨迹构建基因
  pdf(paste0("train.ordergenes",format(Sys.time(),"%H_%M_%S"),".pdf"))
  plot_ordering_genes(cds)
  dev.off()
  #降维
  print("正在降维，请耐心等待")
  cds <- reduceDimension(cds,
                         max_components = 2,
                         method = "DDRTree" )
  #拟时序排列细胞
  print("正在排列细胞，请耐心等待")
  cds <- orderCells(cds)
  
  #保存cds文件
  save(cds,file =paste0('cds_',format(Sys.time(),"%y_%m_%d"),".RData"))
  if(F){
    plot1 <- plot_cell_trajectory(cds,color_by =unlist(strsplit(id,"~"))[2],size = 1,show_backbone = T)+ scale_color_npg()
    ggsave(paste0("train.monocle.pseudotime_",format(Sys.time(),"%H_%M_%S"),".pdf")
           , plot = plot1, width = 7, height = 7)
    
    plot2 <- plot_cell_trajectory(cds,color_by = unlist(strsplit(id,"~"))[2])+ scale_color_npg()+ facet_wrap(unlist(strsplit(id,"~"))[2], nrow = 1)
    ggsave(paste0("train.monocle.pseudotime2_",format(Sys.time(),"%H_%M_%S"),".pdf")
           , plot = plot2, width = 7, height = 7)
    
    plot3 <- plot_cell_trajectory(cds,color_by = "State",size = 1,show_backbone = T)+ scale_color_npg()
    ggsave(paste0("train.monocle.pseudotime3_",format(Sys.time(),"%H_%M_%S"),".pdf")
           , plot = plot3, width = 7, height = 7)
    
    plot4 <- plot_cell_trajectory(cds,color_by = "State") + facet_wrap("~State", nrow = 1)+ scale_color_npg()
    ggsave(paste0("train.monocle.pseudotime4_",format(Sys.time(),"%H_%M_%S"),".pdf")
           , plot = plot4, width = 7, height = 7)
    
    df <- pData(cds)
    ## pData(cds)取出的是cds对象中cds@phenoData@data的内容 View(df)
    pdf("train.monocle.pseudotime7.pdf",width=7,height=7)
    ggplot(df,aes(Pseudotime, colour = NA, fill=NA)) + 
      #
      geom_density(bw=0.5,linewidth=1,alpha =0.5)+theme_classic2()
    dev.off()
  }
}