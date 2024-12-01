#2.通路富集注释####
find_pathway_Profile <- function(){
  if (!dir.exists("../output/04.subcell/07.enrichment")) {
    dir.create("../output/04.subcell/07.enrichment")}
  
#A:Env_setting==================================================================
  
  {#差异分析阈值参数
  DEG_csv='../output/04.subcell/06.DEG/Bcell/Bcell_marker.csv'
  log2FC_cutoff = log2(1)
  pvalue_cutoff = 0.05
  padj_cutoff = 0.05
  }
  
  {#富集分析通路参数
  A <- str_split(a,'\n',simplify = T) %>% as.vector()#a是我们自己挑的通路ID
  B <- str_split(a,'\n',simplify = T) %>% as.vector()
  C <- str_split(a,'\n',simplify = T) %>% as.vector()
  D <- str_split(a,'\n',simplify = T) %>% as.vector()
  E <- str_split(a,'\n',simplify = T) %>% as.vector()
  GSVA_pathway_list=list(
    '00_Bn_IL4R+'=A,
    '01_AtMB_ITGAX+'=B,
    '07_B1_SPN+'=C,
    '08_Plasma_JCHAIN+'=D,
    '09_BINF_IFIT3+'=E
  )
  }

#B:所需包以及函数载入===========================================================
  options(stringsAsFactors = F) 
  library(clusterProfiler)
  #BiocManager::install("clusterProfiler")
  library(enrichplot)
  library(tidyverse)
  library(org.Hs.eg.db)
  #BiocManager::install("org.Hs.eg.db")
  #library(org.Mm.eg.db)
  library(DOSE)
  library(pathview)
  library(msigdbr)
  library(Seurat)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(paletteer)
  library(tidyr)
  library(cowplot)
  library(ggplotify)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(clusterProfiler)
  library(GseaVis)
  library(cowplot)
  library(enrichplot)
  library(ggplot2)
  library(ggplotify)
  library(ggradar)
  library(tidyverse)
  library(patchwork)
  {#library(GSVA),别想了，用docker吧
  #这个GSVA，我可以在终端进行，具体步骤如下（因为我有这些依赖库，只是不在依赖路径里面）：
  #conda activate Rstudio-server
  #export PATH=$PATH:/home/jifanghan/miniconda3/envs/Rstudio-server/lib/
  #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jifanghan/miniconda3/envs/Rstudio-server/lib/
  #R
  #Sys.getenv()
  #.libPaths('~/R/4.4.1/library/')
  #library(GSVA)
  
  #也可以在Rstudio中实现,但不起作用
  #PATH <- Sys.getenv('PATH')
  #LD_LIBRARY_PATH <- Sys.getenv('LD_LIBRARY_PATH')
  #Sys.setenv(PATH=paste0(PATH,':/home/jifanghan/miniconda3/envs/Rstudio-server/lib'))
  #Sys.setenv(LD_LIBRARY_PATH=paste0(LD_LIBRARY_PATH,':/home/jifanghan/miniconda3/envs/Rstudio-server/lib'))
  }
  # library(tidyverse)
  library(org.Hs.eg.db)
  
  #对针对基因集进行富集分析
  trsymbol <- function(genelist,OrgDb="org.Hs.eg.db",organism  = "hsa",ont="ALL"){
    gene_entrez <- as.character(na.omit(bitr(genelist, #数据集
                                             fromType="SYMBOL", #输入格式
                                             toType="ENTREZID", # 转为ENTERZID格式
                                             OrgDb=OrgDb)[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
    kegg_enrich_results <- enrichKEGG(gene = gene_entrez,
                                      organism  = organism, #物种人类 hsa 小鼠mmu
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.2)
    kegg_enrich_results <- DOSE::setReadable(kegg_enrich_results, 
                                             OrgDb=OrgDb, 
                                             keyType='ENTREZID')#ENTREZID to gene Symbol
    go_enrich_results <- enrichGO(gene = gene_entrez,
                                  OrgDb = OrgDb,
                                  ont   = ont  ,     #One of "BP", "MF"  "CC"  "ALL" 
                                  pvalueCutoff  = 0.05,
                                  qvalueCutoff  = 0.2,
                                  readable      = TRUE)
    save(kegg_enrich_results,go_enrich_results,file = "../output/04.subcell/07.enrichment/enrichment_results_temp.Rdata")
  }
  #对针对人工选择的通路获得其基因集
  get_genesets <- function(pathway_list,
                           species = "Homo sapiens",
                           category = "C5",subcategory = "GO:BP"){
    GO_df_ALL <- msigdbr(species = species,
                        category = category ,
                        subcategory = subcategory)  
    GO_df <- dplyr::select(GO_df_ALL,gene_symbol,gs_exact_source,gs_name)
    GO_match <- lapply(pathway_list, function(x){
      GO_list <- x
      GO_results <-lapply(GO_list, function(i){
        GO_SUB <- GO_df %>% filter(gs_exact_source==i)
        GO_SUB <- GO_SUB[!duplicated(GO_SUB$gene_symbol),]
        go_list <- split(GO_SUB$gene_symbol, GO_SUB$gs_name)
      })
    })
    flattened_list <- Reduce(c, unlist(GO_match, recursive = FALSE))#扁平列表
  }


#C:正式pipline==================================================================
  #00 首先先进行Go或者Kegg分析
  Deg_df <- read_csv(DEG_csv) %>% 
    column_to_rownames('...1')
  need_DEG <- Deg_df[,c(2,5:7)]
  need_DEG <- need_DEG[with(need_DEG,
                            avg_log2FC>log2FC_cutoff & 
                              p_val_adj < pvalue_cutoff),]
  path_list <- aggregate(need_DEG$gene,by=list(need_DEG$cluster),function(x){
    genelist <- x
    if (length(genelist) >= 10) {
      trsymbol(genelist,ont = "BP")
      enrich_res <- load("../output/04.subcell/07.enrichment/enrichment_results_temp.Rdata")
      go_results<- go_enrich_results@result
      kegg_results <- kegg_enrich_results@result
      return(list(go_results,kegg_results))
#      write.csv(go_results,file = paste0("../output/04.subcell/07.enrichment/path",Name,"_go_up.csv"))
#      write.csv(kegg_results,file=paste0("../output/04.subcell/07.enrichment/path",Name,"_kegg_up.csv"))
    }
})
  path_list2 <- lapply(path_list$x,function(x){print(x)})
  names(path_list2) <- path_list$Group.1
  
  for (i in names(path_list2)) {
    go_results <-  path_list2[[i]][[1]]
    kegg_results <- path_list2[[i]][[2]]
    write.csv(go_results,file = paste0("../output/04.subcell/07.enrichment/",i,"_go_up.csv"))
    write.csv(kegg_results,file=paste0("../output/04.subcell/07.enrichment/",i,"_kegg_up.csv"))
  }
  #接下来挑选符合你生物学背景的通路吧（还是觉得很麻烦，自己挑5个左右吧）
  
  ###01 各个细胞亚群的Go条形图
  pathway_list <- list.files('../output/04.subcell/07.enrichment',full.names = T,pattern = "go_up.csv$") 
  lapply(pathway_list, function(x){
    df <- read_csv(x) 
    df$qscore <- -log(df$p.adjust,base = 10)
    df <- df[,c("Description","p.adjust","Count","qscore")] %>% arrange(desc(df$qscore))
    df$Description <- factor(df$Description,levels=rev(df$Description))
    
    ggplot(df[1:20,]) +
      geom_col(mapping = aes(x=Description,y=qscore,fill= Count),width = 0.8,)+
      geom_text(aes(x=Description,y=qscore),
                label = df[1:20,]$Description,hjust = 1.04,colour = "white",size = 2,fontface = 'bold')+
      theme_bw()+
      scale_fill_gradient(name="Counts",
                          low='#466983FF',
                           high='#CE3D32FF')+
      theme(axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            #legend.position = "none",
      )+
      scale_y_continuous(expand=c(0,0.01))+
      coord_flip()

    ggsave(paste0("../output/04.subcell/07.enrichment/",str_split(x,"/")[[1]][5],"_go_barplot.pdf"),
           width = 6.5,
           height = 6.6)
  })
  
###02 各个细胞亚群前5-10通路的GSEA图
  if(!dir.exists("../output/04.subcell/07.enrichment/GSEA_res/")){
    dir.create("../output/04.subcell/07.enrichment/GSEA_res/")
  }
  #001准备好差异基因
  genelist <- Deg_df$avg_log2FC
  names(genelist) <- Deg_df$gene
  genelist <- sort(genelist,decreasing = T)
  genelist <- genelist[genelist != 0]
  #002准备好需要富集的通路
  lapply(pathway_list, function(x){
    df <- read_csv(x)
    name <- str_split(x,"/")[[1]][5]
    df$qscore <- -log(df$p.adjust,base = 10)
    df <- df[,c("Description","qscore","ID")] %>% arrange(desc(df$qscore))
    df$ID <- factor(df$ID,levels=df$ID)
    #取qscore前20的通路ID
    pathway_ID <- df$ID[1:20]
    #获得前20的通路具体基因
    lapply(pathway_ID, function(i){
      GO_SUB <- GO_df %>% filter(gs_exact_source==i)
      pathname <- unique(GO_SUB$gs_name)
      if(nrow(GO_SUB)!=0){
        GO_SUB <- GO_SUB[!duplicated(GO_SUB$gene_symbol),] 
        GO_SUB <- GO_SUB[,c(3,1)]
        GO_SUB <- rename(GO_SUB,'gs_name'='term','gene_symbol'='gene')
        set.seed(666)
        egmt <- GSEA(genelist,
                     TERM2GENE = GO_SUB,
                     verbose = F,
                     nPerm=10000,
                     minGSSize = 10,#自定义的geneset某个基因通路很容易没有那么多基因，
                     maxGSSize = 1000,#如果设置太大size，egmt就不会有这个通路，这个时候
                     pvalueCutoff = 1)
        gsea_results <- egmt@result#如果实在想绘制这个通路的gsea图片，就使用gseaplot2吧
        gsea_results2 <- gsea_results[order(gsea_results$enrichmentScore,decreasing = F),]
        #画图
        if(nrow(gsea_results2)!=0){
        terms <- egmt@result$ID#如果terms里面的通路，与egmt的result$ID不一致，就容易报错
        gsea_plots<-gseaNb(object = egmt,
               geneSetID = terms,
               subPlot = 2,
               addPval = T,
               pvalX = 0.75,
               pvalY = 0.75,
               pCol = 'black',
               pvalSize = 4,
               pHjust = 0,)
        ggsave(paste0("../output/04.subcell/07.enrichment/GSEA_res/",name,"_",pathname,"_gsearesult.pdf"),
               plot = gsea_plots,
               width = 8.5,
               height = 4.7)}
      }
    })
})

###03 GSVA显示转录谱
if(F){
  genesets <- get_genesets(GSVA_pathway_list)#得到通路序列
  load('../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata')
  #bulk看看宏观差异
  Idents(seuratdata)<-"manual_L2"
  dat <- AverageExpression(seuratdata, assays = "SCT",  slot = "data")[[1]]
  #选取非零基因
  dat <- dat[rowSums(dat)>0,] 
  #选取我想比较的几种细胞类型
  dat <- dat[,c(1,2,9,10,8)]
  dat <- as.matrix(dat)
  save(genesets,file = '../output/04.subcell/07.enrichment/gsva_genesets.Rdata')
  save(dat,file = '../output/04.subcell/07.enrichment/gsva_dats.Rdata')
  #docker中运行：
  if(F){
    load("gsva_genesets.Rdata")
    load("gsva_dats.Rdata")
    gsva_res <- gsva(expr=dat,
                     gset.idx.list=genesets,
                     method="gsva",
                     mx.diff=FALSE,
                     kcdf="Gaussian",#Gaussian,Poisson
                     verbose=T,
                     parallel.sz = parallel::detectCores())
    gsva.df <- data.frame(Genesets=rownames(gsva_res), gsva_res, check.names = F)
    write.csv(gsva.df,"bulk_gsva.csv")
  }
  bulk_gsva <- read_csv('../output/04.subcell/07.enrichment/bulk_gsva.csv')
  bulk_gsva <- bulk_gsva[,-1] %>% filter(!duplicated(bulk_gsva$Genesets)) %>% column_to_rownames(.,var="Genesets")
  colnames(bulk_gsva) <- gsub("g[0-9]+","",colnames(bulk_gsva))
  colnames(bulk_gsva) <- gsub("^-","",colnames(bulk_gsva))
  df <- rownames_to_column(bulk_gsva)
  write_csv(df,'../output/04.subcell/07.enrichment/bulk_gsva_heatmap.csv',)
  bulk_gsva <- read_csv('../output/04.subcell/07.enrichment/bulk_gsva_heatmap.csv') %>% 
    column_to_rownames(var='rowname')
  rownames(bulk_gsva) <- rownames(bulk_gsva) %>%gsub("^GOBP_","",.)%>%str_to_title()
  #热图
  ##设置颜色：
  col_fun <- colorRamp2(c(-1,0,1),c("navy","#e8f6e5","#d43627"))
  pdf("../output/04.subcell/07.enrichment/Bulk_gsva_heatmap.pdf",
      width =3.5,
      height = 6.7)
  Heatmap(bulk_gsva,
          cluster_columns = F,
          cluster_rows = F,
          show_row_names = T,
          row_names_gp = gpar(fontface="bold",fontsize=4),
          col = col_fun,
          show_heatmap_legend = F,
          border = F,
          gap = unit(1, "mm"),
          rect_gp = gpar(col="white",lwd=1),
          column_names_gp = gpar(fontface="bold",fontsize=5.5),
          column_names_centered = F,
          column_names_rot = 90,
          left_annotation = c(),
          #split = matrix$gene_10.cluster,
          row_title_rot = 90,
          column_title_rot = 0,
          row_title_gp = gpar(fontsize=3,fontface="bold"))
  dev.off()
}#bulk看看宏观差异
if(F){
  subsets <- subset(seuratdata,subset=manual_L2 %in% c("00_Bn_IL4R+","01_AtMB_ITGAX+",
                                                       "07_B1_SPN+","08_Plasma_JCHAIN+",
                                                       "09_BINF_IFIT3+"))
  subsets@meta.data$manual_L2 <- droplevels(subsets@meta.data$manual_L2,
                                            exclude=setdiff(levels(subsets@meta.data$manual_L2),
                                                            unique(subsets@meta.data$manual_L2)))
  scdata <- LayerData(subsets,assay = 'SCT',layer = 'data') %>% as.matrix()
  save(scdata,file = '../output/04.subcell/07.enrichment/sc_GSVA_matrix.Rdata')
  #docker中运行：
  if(F){
    load("gsva_genesets.Rdata")
    load("sc_GSVA_matrix.Rdata")
    gsva_res <- gsva(expr=scdata,
                     gset.idx.list=genesets,
                     method="gsva",
                     mx.diff=FALSE,
                     kcdf="Gaussian",
                     verbose=T,
                     parallel.sz = parallel::detectCores())
    gsva.df <- data.frame(Genesets=rownames(gsva_res), gsva_res, check.names = F)
    write.csv(gsva.df,"bulk_gsva.csv")
  }
  sc_gsva <- read_csv('../output/04.subcell/07.enrichment/sc_gsva.csv')
  scdata <- load('../output/04.subcell/07.enrichment/sc_GSVA_matrix.Rdata')
  sc_gsva <- sc_gsva[,-1] %>% 
    filter(!duplicated(sc_gsva$Genesets)) %>% 
    column_to_rownames(.,var="Genesets")
  write_csv(sc_gsva,'../output/04.subcell/07.enrichment/sc_gsva_dedup.csv')
  #获得gsva评价的meta
  a <- AddMetaData(subsets,t(sc_gsva))
  gsva_df <- a@meta.data[,22:ncol(a@meta.data)] %>% arrange(manual_L2)
  #设置注释colbar
  barnumber <- gsva_df$manual_L2
  #设置注释colbar颜色
  col_anno <- as.character(paletteer_d("ggsci::default_igv",n=5))#n是注释有n种
  col_anno<- rep(col_anno,as.vector(table(barnumber)),each=T)
  names(col_anno) <- barnumber
  col_anno <- HeatmapAnnotation(Module=gsva_df$manual_L2,
                                  col=list(Module=col_anno),
                                  border = F,
                                  show_annotation_name = F,
                                  annotation_legend_param = list(title="",
                                                                 title_position = "topcenter",
                                                                 title_gp = gpar(fontsize=8,fontface="bold")),
                                  show_legend = T
  )
  
}#singlecell看看单细胞差异
if(F){colnames(bulk_gsva)
  bulk_gsva %>% filter(.,`BINF-IFIT3+`>0) %>% rownames() -> a
  #GO_df_ALL <- msigdbr(species = species,
  #                     category = category ,
  #                     subcategory = subcategory)  
  #GO_df <- dplyr::select(GO_df_ALL,gene_symbol,gs_exact_source,gs_name)
  GO_df%>% filter(GO_df$gs_name %in% a) -> b
  b$gs_exact_source %>% unique() -> x
  a <- x}#对通路精挑细选的代码，a再带回14行
}