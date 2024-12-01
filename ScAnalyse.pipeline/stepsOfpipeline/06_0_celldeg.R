#cell <- "CD8"

celldeg <- function(cell,compares,specie) {
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(Seurat)
  library(ggpubr)
  library(tidyverse)
  theme_set(theme_cowplot())
  
  if (!dir.exists(paste0("../output/06.DEG/"))) dir.create(paste0("../output/06.DEG/"))
  if (!dir.exists(paste0("../output/06.DEG/",cell))) dir.create(paste0("../output/06.DEG/",cell))
  celldata=load(paste0("../output/01.QC/7.",projectname,".celltype.check.Rdata") )
  
  seuratdata$Tissue <- gsub("-","_",seuratdata$Tissue)
  Idents(seuratdata) <- seuratdata$Tissue

  #### Identify differential expressed genes across conditions ####
  avg.seuratdata <- AggregateExpression(object = seuratdata, 
                      group.by = c("manual_L1",'Tissue')
                      )$RNA
  # transpose
  avg.seuratdata.t <- t(avg.seuratdata)
  # convert to data.frame
  avg.seuratdata.t <- as.data.frame(avg.seuratdata.t)
  # get values where to split
  splitRows <- gsub('_.*', '', rownames(avg.seuratdata.t))
  # split data.frame
  avg.seuratdata.split <- split.data.frame(avg.seuratdata.t,
                                f = factor(splitRows))
  # fix colnames and transpose
  avg.seuratdata.split.modified <- lapply(avg.seuratdata.split, function(x){
                    rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
                    t(x)
              })
  
  # 1. Get counts matrix
  counts_cell <- avg.seuratdata.split.modified[[1]] %>% as.data.frame
  counts_cell$gene <- rownames(counts_cell)
  write.csv(counts_cell,paste0("../output/06.DEG/",cell,"/avg_tissue.csv"),row.names = FALSE)
  
  #i <- compares[1]
  for (i in compares) {
    print(i)
     comname <- unlist(strsplit(i," VS "))
     out_name <- paste0(comname,collapse = "_")
     comname1 <- gsub("-","_",comname)
     #保存差异表格
     deg <- FindMarkers(seuratdata, ident.1 = comname[1], ident.2 = comname[2],
                            verbose = FALSE,min.pct = 0.25,min.diff.pct = -Inf,
                            min.cells.feature =  3,min.cells.group = 3)
     write.csv(deg,paste0("../output/06.DEG/",cell,"/",out_name,".1.deg.csv"),quote = TRUE)
     #head(deg)
     deg$deg <- ifelse(deg$avg_log2FC>0,"UP", "DOWN") 
     #deg$deg[abs(deg$avg_log2FC)<1] <- "NA" 
     deg$deg[deg$p_val_adj>0.05] <- "NA" 
     deg$SYMBOL <- rownames(deg)
     deg <- deg[order(deg$avg_log2FC,decreasing = TRUE),]
     deg1 <- deg[deg$deg!= "NA",] 
     if(nrow(deg1)>10) {
     top10 <- deg1$SYMBOL[c(1:10,c(nrow(deg1)-9):nrow(deg1))]
     } else (top10 <- deg1$SYMBOL)
     #### 火山图
     p <- ggplot(deg,aes(x=avg_log2FC,y=-log10(p_val_adj),color=deg))+geom_point()+
       ggtitle(out_name)+
       xlab("log2 FC")+ylab("-log10 adjusted p-value")+
       geom_hline(yintercept = 1,linetype =3)+
       geom_vline(xintercept = c(-1,1),linetype=3)+
       scale_color_manual(values =c("deepskyblue4", "grey","red"),
                          name=format("qadj<0.05") )+
       theme(plot.title = element_text(hjust = 0.5))+
       geom_text_repel(data = subset(deg, SYMBOL %in% top10),
                       aes(x=avg_log2FC,y=-log10(p_val_adj),label=SYMBOL),max.overlaps = 50,color="black",
                       box.padding = unit(0.5, "lines"),
                       point.padding = unit(0.5, "lines"), segment.color = "black", show.legend = FALSE)+
        theme_minimal()
     p
     ggsave(paste0("../output/06.DEG/",cell,"/",out_name,".2.deg.Volcano.pdf"),plot=p,width = 8,height =8)
     
     deg.fc1 <- deg[abs(deg$avg_log2FC)>0.59,]
     deg.fc1 <- deg.fc1[deg.fc1$avg_log2FC!="-Inf",]
     deg.fc1 <- deg.fc1[deg.fc1$avg_log2FC!="Inf",]
     deg.fc1 <- deg.fc1[deg.fc1$p_val_adj<0.05,]
     write.csv(deg.fc1,paste0("../output/06.DEG/",cell,"/",out_name,".2.deg.fc1.5.csv"),quote = TRUE)
     DEG_up <- deg.fc1[deg.fc1$avg_log2FC>0,]
     DEG_down <- deg.fc1[deg.fc1$avg_log2FC<0,]
     write.csv(DEG_up,paste0("../output/06.DEG/",cell,"/",out_name,".deg.up.csv"),
               row.names = FALSE,quote = FALSE) 
     write.csv(DEG_down,paste0("../output/06.DEG/",cell,"/",out_name,".deg.down.csv"),
               row.names = FALSE,quote = FALSE) 
     #
     counts_cell.sub <- counts_cell[rownames(counts_cell) %in% rownames(deg.fc1),   ]
     colnames(counts_cell.sub) <- gsub("-","_",colnames(counts_cell.sub))
     p1 <- ggplot(counts_cell.sub, create_aes(list(x=comname1[1],y=comname1[2]))) + 
       geom_point() + ggtitle(out_name)+theme_minimal()
     p1 <- LabelPoints(plot = p1, points = top10, repel = TRUE,col="red")
     p1
     ggsave(paste0("../output/06.DEG/",cell,"/",out_name,".3.genelist.DEG.plot.pdf"),
            plot=p1,width = 8,height =8)
     library(RColorBrewer)
     pdf(paste0("../output/06.DEG/",cell,"/",out_name,".4.top10.DEG.featureplot.pdf"),
         width = 10,height = 5)
     for (j in 1:length(top10)) {
       #j <- 2
       Feature.Plot <- FeaturePlot(seuratdata, features = top10[j],label = TRUE, 
                                   split.by = "Tissue",
                                   #ncol = 3,
                                   #blend=TRUE,
                                   #cols = c("deepskyblue1","white", "red"),
                                   #cols = brewer.pal(3,"Accent"),  
                                   max.cutoff = 3)
       plot(Feature.Plot)
      }
     dev.off()
     # if(length(unique(seuratdata$Tissue))==3) dotcolor=c("deepskyblue1", "bisque2","pink")
     # if(length(unique(seuratdata$Tissue))==2) dotcolor=c("deepskyblue1", "pink")
     # #deg vlnplot
     pdf(paste0("../output/06.DEG/",cell,"/",out_name,".5.top10.DEG.vlnplot.pdf"),
         width = 6.5,height = 4)
     for (z in 1:length(top10)) {
       plots <- VlnPlot(seuratdata, features =top10[z], split.by = "Tissue",
                        #cols = dotcolor,
                        group.by = "manual_L2",ncol = 3,
                        pt.size = 0, combine = FALSE)
       plot(plots[[1]]  )
     }
      dev.off() 
  }
  }

