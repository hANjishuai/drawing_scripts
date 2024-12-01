
clustermaker_table <- function() {
  library(Seurat)
  library(dplyr)
  celldata <- load(paste0("../output/01.QC/7.",projectname,".celltype.check.Rdata"))
  plan("multisession", workers =24) ####多线程运算

  seuratdata.markers <- read.csv("../output/02.cluster/cluster.markers.csv")
  cellcluster <- read.csv("../output/03.celltype/0.singleR_cellplot.csv",header = TRUE)
  
  seuratdata.markers$cell <- NA
  for (i in 1:nrow(cellcluster)) {
    index <- seuratdata.markers$cluster %in% cellcluster$cluster[i-1]
    seuratdata.markers$cell[index] <- cellcluster$manual_level1[i-1]
  }
    # save(seuratdata,seuratdata.markers,file=paste0("../output/01.QC/",projectname,".celltype.check_marker.Rdata") )
  
  write.csv(seuratdata.markers,"../output/03.celltype/cell.markers.csv",row.names = FALSE) 
  
  cellnum <- length(cellcluster$manual_level1)
  cluster_marker_row <- data.frame(cluster=cellcluster$manual_level1,markers=NA)
  for (c in cluster_marker_row$cluster) {
    temp_c = seuratdata.markers$gene[seuratdata.markers$cell == c]
    temp_c = paste0(temp_c,collapse = ",")
    cluster_marker_row$markers[which(cluster_marker_row$cluster== c)] = temp_c
  }  
  write.csv(cluster_marker_row,"../output/03.celltype/cell.markers_row.csv",row.names = FALSE) 
  
  seuratdata.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  write.csv(top10,"../output/03.celltype/cell.markers.top10.csv",row.names = FALSE) 
  }


