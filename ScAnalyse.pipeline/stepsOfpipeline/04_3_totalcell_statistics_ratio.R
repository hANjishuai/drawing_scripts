

totalcell_statistics_ratio <- function(celltype,specie) {
  
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(ggthemes)
  library(reshape2)
  library(pheatmap)
  if (!dir.exists(paste0("../output/03.celltype/statistics"))) dir.create(paste0("../output/03.celltype/statistics"))
  
  cellmate <- read.csv("../output/03.celltype/cellmate.csv",header = TRUE)
  cellmate$Tissue <- gsub("-","_",cellmate$Tissue)

 #  计算Ratio O/E
  cell_dist <- calTissueDist(dat.tb = cellmate,colname.cluster = "manual_L3",
                colname.patient = "Tissue",colname.tissue = "Tissue")
  cell_dist <- as.data.frame(cell_dist)
  cell_dist_d <- dcast(cell_dist,Var1~Var2)
  rownames(cell_dist_d) <- cell_dist_d$Var1
  cell_dist_d <- cell_dist_d[,-1]
  bk<-c(seq(min(cell_dist_d),1,length.out = 50),seq(1.01,max(cell_dist_d),length.out = 50))
  
  pdf(paste0("../output/03.celltype/statistics/ratio_0_pheatmap.pdf"),width = 10,height = 10)
  pheatmap(cell_dist_d,cellwidth = 12,cellheight = 10,cluster_cols = FALSE,
           breaks = bk,
           color = colorRampPalette(c("blue", "white", "red"))(length(bk)),
           main = "Ratio_O/E" )
 
 dev.off()
  
  #### all tissues and cells
  p1 <- ggplot(cell_dist,aes(x=Var2,y=Var1,size=Freq,color=Freq))+
            geom_point()+
            scale_color_gradient(low = "blue",high = "red" )+
            theme_minimal()
  p1
  ggsave(paste0("../output/03.celltype/statistics/ratio_1_tissue_cell_heatmap.pdf"),
         plot = p1,width = 12,height = 9)

  #### histogram plot rate
  
  p1 <- ggplot(cell_dist,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat = "identity")+
    facet_wrap(~ Var2,ncol=3)+theme_few()+
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,face = "bold"),
          text= element_text(size = 14,face = "bold"))
  p1
  ggsave(paste0("../output/03.celltype/statistics/ratio_2_bar_single.pdf"),
         plot = p1,width = 16,height = 12)
}
