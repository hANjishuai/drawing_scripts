subcell <- function(projectname){
  library(Seurat)
  celldata <- load(paste0("../output/01.QC/7.",projectname,".celltype.check.Rdata"))
  Idents(seuratdata) <- seuratdata$manual_L1
  celltype <- unique(seuratdata$manual_L1)
  for (cell in celltype) {
    subcells <- subset(seuratdata, idents = cell)
    subcells@meta.data$seurat_clusters <- subcells@meta.data$manual_L1
    subcells@meta.data$seurat_clusters <- droplevels(subcells@meta.data$seurat_clusters,
                                                     exclude=setdiff(levels(subcells@meta.data$seurat_clusters),
                                                                     unique(subcells@meta.data$seurat_clusters)))
    save(subcells,file=paste0("../output/01.QC/subset_",cell,".Rdata") )
  }
}
