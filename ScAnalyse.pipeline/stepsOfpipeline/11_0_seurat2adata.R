#R中的操作，在单细胞seurat对象中提取数据，用来在python中构建adata
#首先在R中我们定义一个函数，提取需要的文件

#我们将seurat中的seurat单细胞对象转化为scanpy对象，能够在python中使用(pyscenic,cellphoneDB)
#这个操作和我们在RNA速率中的分析一样，现在R中使用一个函数将所有的文件分离出来，然后在python中构建

#2024-10-25 Fanghan_Jay Version::0.0.1
seurat_to_adata <- function(object,#seurat对象
                            Dimension=c('UMAP','TSNE'),#降维方式
                            path){#文件保存路径
  library(Seurat)
  library(Matrix)
  
  seurat_obj <- object
  seurat_obj$barcode <- Cells(seurat_obj)
  if(Dimension=='UMAP'){
    cell.embeddings<- seurat_obj@reductions$umap@cell.embeddings
    seurat_obj$UMAP_1 <- cell.embeddings[,1]
    seurat_obj$UMAP_2 <- cell.embeddings[,2]
  }else{
    
    cell.embeddings<- seurat_obj@reductions$tsne@cell.embeddings
    seurat_obj$TSNE_1 <- cell.embeddings[,1]
    seurat_obj$TSNE_2 <- cell.embeddings[,2]
  }
  #保存metadat
  write.csv(seurat_obj@meta.data, file=paste0(path,'metadata.csv'), quote=F, row.names=F)
  #保存matrix
  counts_matrix <- LayerData(seurat_obj,assay = 'RNA',layer = 'counts')
  writeMM(counts_matrix, file=paste0(path, 'counts.mtx'))
  #PCA降维信息(harmony过，就使用harmony降维数据)
  #write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(path,'pca.csv'), quote=F,row.names=F)
  write.csv(seurat_obj@reductions$harmony@cell.embeddings, file=paste0(path,'harmony.csv'), quote=F,row.names=F)
  #保存gene name
  write.table(data.frame('gene'=rownames(counts_matrix)),file=paste0(path,'gene_names.csv'),
              quote=F,row.names=F,col.names=F)
}

#测试集
if(F){
testassays <- load("../output/01HC-DLE-SLE_blood-result/04.subcell/02.cluster/1.CD79A_Bcell.cluster.Rdata")
seurat_to_adata(seuratdata,"UMAP","../output/01HC-DLE-SLE_blood-result/test/")
}

