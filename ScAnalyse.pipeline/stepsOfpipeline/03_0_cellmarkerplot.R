manual_marker_plot <- function(projectname) {
  library(gdata)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(readxl)
  plan("multisession", workers =20) ####多线程运算
  options(future.globals.maxSize = 5000 * 1024^2)
  if (!dir.exists(paste0("../output/03.celltype"))) dir.create(paste0("../output/03.celltype"))
  #### load data ####
  celldata=load(paste0("../output/01.QC/5.",projectname,".cluster.Rdata") )
  # cellsignal cellmarker
  cellsignal <- read_xls("../db/immune_Cell_markers.xls",sheet = 1,
                         col_names = F)
  cs <- apply(cellsignal,1,function(x) {
    temp_x <- x[-1] %>% as.data.frame
    return(temp_x)
  })
  names(cs) <- cellsignal$...1
  if (F) {
  # microglia marker 
  microgliasignal <- read_xls("../db/neural_markers.xls",sheet = 1)
  
  ms <- apply(microgliasignal, 1, function (x) {
    temp_s <- x[-1]
    temp_info <- convert_human_to_mosue(gene_list =temp_s )
    return(temp_info$Symbol)
  })
  names(ms) <- microgliasignal$Cells
  }
  # merge data
  ms <- list()
  all_genes <- c(cs,ms)
  all_genes0 <- unlist(all_genes)
  all_genes0 <- all_genes0[all_genes0!=""]
  all_genes0 <- all_genes0[!is.na(all_genes0)]
  all_genes0 <- all_genes0[!duplicated(all_genes0)]
  p_all_genes <- DoHeatmap(subset(seuratdata,downsample=8000),
                           features = all_genes0)
  ggsave(paste0("../output/03.celltype/0.manual_all_marker_heatmap.png"),
         plot=p_all_genes,width = 8,height =15)
  #### plot heatmap ####
  scaled_data <- seuratdata@assays$SCT@scale.data %>% rownames
  pdf(paste0("../output/03.celltype/0.manual_single_marker_heatmap.pdf"),
      width = 8,height =8)
  for (l in 1:length(all_genes)) {
    temp_g <- all_genes[[l]]
    temp_g <- temp_g[temp_g!=""]
    temp_g <- temp_g[!is.na(temp_g)]
    temp_g <- temp_g[!duplicated(temp_g)]
    if (length(intersect(scaled_data ,temp_g))>0 ) {
      p_temp_genes <- DoHeatmap(subset(seuratdata,downsample=1000),
                                features = temp_g,assay = "SCT")+
                                labs(title = names(all_genes)[l])
      #ggsave(plot = p_temp_genes,filename = paste0("../output/03.celltype/0.",
      #                                             names(all_genes)[l],
      #                                             ".png"))
      plot(p_temp_genes)
    }
    
  }
  dev.off()
}

