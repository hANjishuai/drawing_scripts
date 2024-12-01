
sample_seuratdata <- function() {
    library(Seurat)
    loaddata=load(paste0("../output/01.QC/7.",projectname,".celltype.check.Rdata"))
    # 定义每个聚类希望下采样的数量
    n_sample_per_cluster <- 100  # 根据需要设置
    # 创建一个列表来存储下采样后的细胞
    sampled_cells <- unlist(lapply(unique(seuratdata$manual_L3), function(cluster) {
        cells <- WhichCells(seuratdata, idents = cluster)
        sample(cells, min(n_sample_per_cluster, length(cells)))  # 随机抽样
    }))
    # 根据下采样的细胞创建新对象
    downsampled_obj <- subset(seuratdata, cells = sampled_cells)
    save(downsampled_obj,file=paste0("../output/01.QC/8.",projectname,".celltype.sample.Rdata") )
    
    
    
}
