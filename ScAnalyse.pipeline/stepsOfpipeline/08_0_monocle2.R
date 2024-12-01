
# rm(list = ls())
# gc()
#cell="Microglia"
monocle2_self <- function(cell,T0name){
    library(Seurat)
    library(monocle)
    options(future.globals.maxSize = 50000 * 1024^2)  # 设置为 1000 MB
    
if (!dir.exists(paste0("../output/07.TI"))) dir.create(paste0("../output/07.TI"))
if (!dir.exists(paste0("../output/07.TI/monocle2"))) dir.create(paste0("../output/07.TI/monocle2"))
if (!dir.exists(paste0("../output/07.TI/monocle2/",cell))) dir.create(paste0("../output/07.TI/monocle2/",cell))

savename <- paste0("../output/07.TI/monocle2/",cell,"/")
#Load Seurat Obj
celldata=load(paste0("../output/01.QC/subset_",cell,".Rdata"))

table(subcells$manual_L1)

if (ncol(subcells)>80000) {
  subcells <- subset(subcells,downsample=3000)
}
# 选择指定聚类的细胞子集
DefaultAssay(subcells) <- "RNA"
subcells <- NormalizeData(subcells)

#Transfer into cds
cds <- as.CellDataSet(subcells)
rm(subcells)
# Monocle2 process
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## ordering by marker gene per cluster
deg <- read.csv(paste0("../output/03.celltype/cell.markers.top10.csv"))
index <- grep(cell,deg$cell)
deg <- deg[index,]
sel.gene <- unique(deg$gene)
cds <- monocle::setOrderingFilter(cds, sel.gene)

## dimension reduciton
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

## ordering cells
cds <- monocle::orderCells(cds)

## ordering cells by assigning root nodes
cds$manual_L3 <- factor(cds$manual_L3,unique(cds$manual_L3))
GM_state <- function(cds){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds$manual_L3)[,T0name]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
if (!is.null(T0name)) {
    cds <- monocle::orderCells(cds, root_state =  GM_state(cds))
}

# Visualization
p1 <- monocle::plot_cell_trajectory(cds, color_by = "manual_L3")
p1
ggsave(plot = p1,
       filename = paste0("../output/07.TI/monocle2/",cell,"/trajetory.pdf"),
       width = 8,height =6)
# Visualization by tissue
p2 <- monocle::plot_cell_trajectory(cds, color_by = "manual_L3")  + facet_wrap(~manual_L3)
p2
ggsave(plot = p2,
       filename = paste0("../output/07.TI/monocle2/",cell,"/trajetory_tissue.pdf"),
       width = 8,height =6)

# Visualization by State
p3 <- monocle::plot_cell_trajectory(cds, color_by = "State")
p3
ggsave(plot = p3,
       filename = paste0("../output/07.TI/monocle2/",cell,"/trajetory_state.pdf"),
       width = 8,height =5)
p4 <- monocle::plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State)
p4
ggsave(plot = p4,
       filename = paste0("../output/07.TI/monocle2/",cell,"/trajetory_state_tisse.pdf"),
       width = 8,height =5)
# Visualization by Pseudotime
p5 <- monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
p5
ggsave(plot = p5,
       filename = paste0("../output/07.TI/monocle2/",cell,"/trajetory_Pseudotime.pdf"),
       width = 8,height =5)
}




