

totalcell_statistics <- function(celltype,specie) {
  
  library(Seurat)
  library(ggplot2)
  library(ggthemes)
  library(reshape2)
  library(pheatmap)
  library(patchwork)
  library(gtools)
  if (!dir.exists(paste0("../output/03.celltype/statistics"))) dir.create(paste0("../output/03.celltype/statistics"))
  
  cellmate <- read.csv("../output/03.celltype/cellmate.csv",header = TRUE)
  cellmate$Tissue <- gsub("-","_",cellmate$Tissue)
  
  #### all tissues and cells
  p <- ggplot(cellmate,aes(x=Tissue,y=manual_L3))
  p1 <- p+geom_count(aes(color = after_stat(n)),alpha=.5,shape=16) +
    scale_color_gradient(low = "blue",
                         high = "red" )+theme_minimal()
  p1
  ggsave(paste0("../output/03.celltype/statistics/1.tissue_cell_manual_point.pdf"),
         plot = p1,width = 12,height = 9)

  #### sicellmatengle cell
  
  tissue_cell_m <- dcast(cellmate,Tissue~manual_L3 )
  tissue_cell_m <- tissue_cell_m[mixedorder(tissue_cell_m$Tissue,decreasing = TRUE),]
  rownames(tissue_cell_m) <- tissue_cell_m$Tissue
  tissue_cell_m1 <- tissue_cell_m[,-1]
  #tissue_cell_m1 <- log2(tissue_cell_m1+1)
  pdf("../output/03.celltype/statistics/1.tissue_cell_manual_heatmap.pdf",width = 8,height = 6)
  pheatmap::pheatmap(tissue_cell_m1,cellwidth = 12,cellheight = 10,
                     color = colorRampPalette(c("blue", "white", "red"))(50),
                     cluster_cols = FALSE,cluster_rows = FALSE)
  dev.off()
  
  # calculate cell rate
  cell_rate <- function(singletissue){
    rate_a <- lapply(singletissue,function (x) {
      x/sum(singletissue)*100
    })
    return(as.numeric(rate_a))
  }
  tissue_cell_rate <- apply(tissue_cell_m[,-1],1,cell_rate) %>% as.data.frame
  rownames(tissue_cell_rate) <- paste0(colnames(tissue_cell_m)[-1],"_rate")
  colnames(tissue_cell_rate) <- tissue_cell_m$Tissue
  tissue_cell_rate <- as.data.frame(tissue_cell_rate)
  tissue_cell_rate <- cbind(Tissue=rownames(tissue_cell_rate),tissue_cell_rate)
  #merge count and rate
  
  tissue_cell_m_r <- cbind(tissue_cell_m,t(tissue_cell_rate[,-1]))

  #### histogram plot rate
  rate_m <- melt(tissue_cell_rate)
  rate_m$G <- gsub("-.*","",rate_m$variable)
  rate_m$Tissue <- gsub("_rate","",rate_m$Tissue)
  p1 <- ggplot(rate_m,aes(x=variable,y=value,fill=G))+geom_bar(stat = "identity")+
    facet_wrap(~ Tissue,ncol=5)+theme_few()+
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,face = "bold"),
          text= element_text(size = 14,face = "bold"))
  p1
  ggsave(paste0("../output/03.celltype/statistics/3.bar_single.pdf"),
         plot = p1,width = 16,height = 12)
 
  #### histogram plot numbers
  num_m <- melt(tissue_cell_m)
  num_m$G <- gsub("-.*","",num_m$Tissue)
         
  p1 <- ggplot(num_m,aes(x=variable,y=value,fill=G))+geom_bar(stat = "identity")+
    facet_wrap(~ Tissue,ncol=3)+theme_few()+
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,face = "bold"),
          text= element_text(size = 14,face = "bold"))
  p1
  ggsave(paste0("../output/03.celltype/statistics/3.bar_single_numbers.pdf"),
         plot = p1,width = 20,height = 12)
}
