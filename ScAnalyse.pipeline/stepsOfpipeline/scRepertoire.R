#===============================================================================
#             scRepertoire:免疫组库分析
#===============================================================================
#input_path='../data/Blood_IR/'
#RE='\\d+_B$'
runIR <- function(input_path,threshold=0.85,RE){
  library(scRepertoire)
  library(immunarch)
  library(stringr)
  library(Seurat)
  library(tibble)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
  library(ggraph)
  library(ggpubr)
  library(ggsci)
  library(ggplot2)
  library(circlize)
  library(scales)
  
  if (!dir.exists("../output/04.subcell/05.scRepertoire")) 
    dir.create("../output/04.subcell/05.scRepertoire")
  
#===============================================================================
#             1、scRepertoire:单细胞免疫组库分析标准流程
#===============================================================================
  ####step1:load data-----------------------------------------------------------
  folders=list.files(input_path)#文件夹目录
  folders
  BCR_list = lapply(folders,function(x){
    if(!startsWith(x,"H"))
      {
      matrix_folder<-paste0(input_path,x,'/','Matrix')
      matrix<-list.files(matrix_folder,pattern = ".csv$",full.names = T)
      IR_matrix <- read.csv(matrix)
      }else{
        matrix_folder<-paste0(input_path,x)
        matrix<-list.files(matrix_folder,pattern = "contig_annotations.csv$",full.names = T)
        IR_matrix <- read.csv(matrix) %>% apply(.,2,function(x){gsub("None","",x)})
        vdjc<-IR_matrix[,c(7:10)]
        vdjc<-apply(vdjc,2,function(x){gsub(RE,"",x)})
        IR_matrix[,c(7:10)] <- vdjc
        IR_matrix <- IR_matrix %>% data.frame()
        }
  })
  #数据清洗有点麻烦捏
  ####step2:数据处理------------------------------------------------------------ 
  # 由于CellRanger的输出是两条链（Tcell：TCRA、TCRB。Bcell：IGH，IGK/IGL）的定量，
  # 下一步是通过cell barcode创建TCR/BCR基因(由VDJC基因组成)和CDR3序列的单个列表对象。
  # 使用TCR:combineTCR()/BCR:combineBCR函数。同时根据样品和ID信息对cell barcode重新贴标签，以防止重复。
  # 这里我们的数据是BCR,所以用combineBCR()函数，TCR数据用combineTCR()函数，这两个函数参数不一样
  #如果barcode多，样本多的话，这一步执行还是需要花费一段时间的！
  data_bcr <- combineBCR(BCR_list,
                         ID = folders,#分组
                         samples=gsub('\\d+_B$','',folders),
                         threshold=threshold)
  #data_bcr用于初步画图
  save(data_bcr,file = '../output/04.subcell/05.scRepertoire/data_bcr.RData')
  #虽然在combineBCR函数运行的时候，我们已经添加了样本的一些分组信息，但是如果我们还需要添加更多的内容，
  #例如年龄，性别，其他分组等等，可以使用addVariable函数。
  data_bcr <- addVariable(data_bcr, variable.name = "condition", 
                        variables =c(rep('DLE',4),rep('HC',3),rep('SLE',4)))
  
  save(data_bcr,file = '../output/04.subcell/05.scRepertoire/data_bcr.RData')
  
  #一些其他操作，假设多个样本需要提取几个特定样本，可以使用subsetContig，不过没什么用，既然是list直接取list即可
  # subset <- subsetContig(data_bcr, name = "ID", variables = c("WT1_bcr", "K3_bcr"))

  ####step3:分析可视化----------------------------------------------------------
  #scRepertoire的可视化函数是基于ggplot2，所以图片的修饰可以直接跟theme主题。
  #可视化函数中，有一个共同参数cloneCall，选择不同的参数，可视化clone变化
  #也即是根据不同的指标/关注指标统计
  #cloneCall有四个选项，分别的意思如下：
  #---“gene”：使用含有TCR/Ig的VDJC基因
  #---“nt”：使用CDR3区域的核苷酸序列
  #---“aa”：使用CDR3区域的氨基酸序列
  #---“strict”：使用包含TCR/Ig + CDR3区域核苷酸序列的VDJC基因。
  
  #注意：作者提到，The gene approach will be the most sensitive, 
  # while the use of nt or aa moderately so, 
  # and the most specific for clonotypes being strict. 
  #所以所有的可视化函数中默认的cloneCall参数是strict，但是实际中，我们需要以其他为clontype的时候根据实际情况和自己的目的
  
   #3.1 量化Clonotypes，使用clonalQuant函数,可以返回total Clonotypes或者特异性的Clonotypes####
  load('../output/04.subcell/05.scRepertoire/data_bcr.RData')
  p1 <- clonalQuant(data_bcr, #前期combineBCR或者combineTCR得到的数据
                    cloneCall="strict", 
                    scale = T,#T表示unique clone占比，F表示unique clone的number
                    chain = "both")+#both表示组合链可视化
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(plot = p1,filename = "../output/04.subcell/05.scRepertoire/total_chain_barlot.pdf")
  
  #只想可视化特异性的chain，那么在chain这里选择参数即可。
  #TCR可以选择“TRA”, “TRB”, “TRD”, “TRG”， BCR则是“IGH” or “IGL”
  p2 <- clonalQuant(data_bcr, #前期combineBCR或者combineTCR得到的数据
                    cloneCall="strict", 
                    scale = T,#T表示unique clone占比，F表示unique clone的number
                    chain = "IGH")+
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(plot = p2,filename = "../output/04.subcell/05.scRepertoire/H_chain_barlot.pdf")
  
  p3 <- clonalQuant(data_bcr, #前期combineBCR或者combineTCR得到的数据
                    cloneCall="strict", 
                    scale = T,#T表示unique clone占比，F表示unique clone的number
                    chain = "IGL")+
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(plot = p3,filename = "../output/04.subcell/05.scRepertoire/L_chain_barlot.pdf")
  #我们想要返回数据，而不是图表，则增加参数exportTable=T
  both <- clonalQuant(data_bcr, 
              cloneCall="strict", 
              scale = T,
              chain = "both",
              exportTable=T)
  IGH <- clonalQuant(data_bcr, 
                    cloneCall="strict", 
                    scale = T,
                    chain = "IGH",
                    exportTable=T)
  IGL <- clonalQuant(data_bcr, 
                     cloneCall="strict", 
                     scale = T,
                     chain = "IGL",
                     exportTable=T)
  bycondition_total <- clonalQuant(data_bcr,
                             cloneCall="strict",
                             scale = T,
                             chain = "both",
                             group.by = 'sample',exportTable=T)
  bycondition_H <- clonalQuant(data_bcr,
                                   cloneCall="strict",
                                   scale = T,
                                   chain = "IGH",
                                   group.by = 'sample',exportTable=T)
  df_list <- list(both,IGH,IGL,bycondition_total,bycondition_H )
  names(df_list)<-c('both','IGH','IGL','bycondition_total','bycondition_H ')
  sink(file='../output/04.subcell/05.scRepertoire/unique_clones.txt',append = T)
  df_list;sink()
  if(F){
  IGH$condition<-gsub('_\\w+_*$','',IGH$values) %>% factor(.,level=c('HC','DLE','SLE'))
  
  ggplot(IGH, aes(fill=condition, y=scaled, x=condition))+
    geom_bar(position=position_dodge(),
             stat="summary",
             width=0.7,
             colour = "black",
             size=1)+
    geom_jitter(width = 0.1)+
    stat_summary(fun.data = 'mean_se', 
                 geom = "errorbar", 
                 colour = "black",
                 width = 0.2,
                 position=position_dodge(0.7))+
    theme_classic()+
    stat_compare_means(mapping = aes(condition), 
                       comparisons = list(c("DLE", "HC"),
                                          c("DLE","SLE"),
                                          c("HC","SLE")), 
                       #label = 'p.signif',
                       method = "t.test",
                       hide.ns=T,step.increase=0.2)
  
  }#简单统计看一下
  
   #3.2 Clonotype丰度，使用clonalAbundance函数，通过丰度来检查克隆型的相对分布####
  p1 <- clonalAbundance(data_bcr, cloneCall = "strict", scale = F,
                        #exportTable = T,
                        chain = "both",
                        group.by = 'sample')
  ggsave(plot = p1,filename = "../output/04.subcell/05.scRepertoire/Abundance_lineplot_bygroup.pdf")
  #scale:Converts the graphs into density plots in order to show relative distributions
  p2 <- clonalAbundance(data_bcr, cloneCall = "strict", scale = F,
                        #exportTable = T,
                        )
  ggsave(plot = p2,filename = "../output/04.subcell/05.scRepertoire/Abundance_lineplot.pdf")
  
  
  
   #3.3 比较Clonotype，使用clonalCompare()函数来查看样本之间的克隆类型和动态变化#####
  clonalCompare(data_bcr, 
                top.clones = 5, #The top number clonotype sequences per group
                cloneCall="gene+nt", 
                graph = "alluvial")+ 
    scale_colour_brewer(palette = "Greens")+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'none')
  #可以看任意样本之间的比较
  clonalCompare(data_bcr, 
                          cloneCall="gene+nt", 
                          top.clones = 5,
                          graph = "alluvial",
                          samples = c("HC_HC1_B","SLE_SLE1_B")) + scale_colour_brewer(palette = "Greens")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   #3.4 Clonal Space Homeostasis####
  #通过检查Clonal Space，我们可以有效地查看特定比例的克隆所占据的样本全部的相对比例。
  #函数里面cloneSize比例的切点可以自定义，这里使用的时默认参数
  clonalHomeostasis(data_bcr, cloneCall = "strict", 
                          cloneSize  = c(Rare = 1e-04, 
                                         Small = 0.001, 
                                         Medium = 0.01, 
                                         Large = 0.1, 
                                         Hyperexpanded = 1))
  clonalHomeostasis(data_bcr, cloneCall = "strict", 
                    cloneSize  = c(Rare = 1e-04, 
                                   Small = 0.001, 
                                   Medium = 0.01, 
                                   Large = 0.1, 
                                   Hyperexpanded = 1),
                    group.by = 'sample')
  
   #3.5 克隆比例####
  #克隆比例使用clonalProportion()函数，按总数对克隆进行排序，然后看特定克隆型占总数比例
  #函数参数clonalSplit表示按拷贝数或出现频率排列的克隆型，即1:10为每个样本中排名前10位的克隆型。此参数可自行调整
  clonalProportion(data_bcr, cloneCall = "gene",
                   exportTable =T,
                   clonalSplit = c(10, 100, 1000, 10000, 30000, 1e+05))
  
   #3.6 样本相似性（其实就类似于相关性），clonalOverlap函数，这里不再展示####
  # clonalOverlap(data_bcr, cloneCall = "aa", method = "overlap")
  
   #3.7 样本聚类 clonalSizeDistribution()，按克隆大小分布对样本进行聚类，和样本相似性差不多的意思####
  # clonalSizeDistribution(data_bcr, cloneCall = "strict", method="ward.D2")
  
  
  
  
  
  
  
  
  
  
  
  
  
   #3.7 Diversity Analysis 多样性分析####
   #利用Shannon、Inverse Simpson、Chao和基于丰度的覆盖估计(ACE)指数对样本类型进行基于克隆型的多样性度量。
   #这几个指数什么意思可自行百度学习。总之这些指数高代表克隆多样性高。
   #从这个结果上我们可以看出，HC组的克隆多样性高于disease组。
  p1<- clonalDiversity(data_bcr, 
                         cloneCall = "aa", 
                         n.boots = 1000, 
                         x.axis='sample', #这里的sample就是我们前面data_bcr的分组
                         group='ID')#ID就是每个样本
  ggsave(plot = p1,filename = "../output/04.subcell/05.scRepertoire/Diversity_jiter.pdf")
  
  ####step4:scRepertoire:免疫组库结合单细胞转录组分析####
  load('../output/04.subcell/02.cluster/1.Bcell.cluster.Rdata')
  load('../output/04.subcell/05.scRepertoire/data_bcr.RData')
  #1：修改barcode
  #合并当然需要注意的问题，combineBCR/combineTCR数据中的barcode和单细胞Seurat的barcode要一致
  data_bcr_rename<- lapply(data_bcr, function(x){
    if (unique(x[,14]) %in% c('DLE','SLE')) {
      barcode<-x[,1]
      barcode <- gsub('^[DSLEHC]+_','',barcode)
      barcode <- gsub('-1$','',barcode)
      x[,'barcode'] <- barcode
    }else{
      barcode<-x[,1]
      barcode <- gsub('^[DSLEHC]+_','',barcode)
      x[,'barcode'] <- barcode
    }
    return(x)
  })
#head(Cells(seuratdata))
#[1] "DLE1_B_AAACCTGAGTTTCCTT" "DLE1_B_AAACCTGCACCCTATC" "DLE1_B_AAACGGGTCTGGAGCC"
#[4] "DLE1_B_AAAGATGCAAACCTAC" "DLE1_B_AAAGATGTCTTAGCCC" "DLE1_B_AAAGCAATCCTCAACC"  
  
#head(data_bcr[[1]][,1])
#[1] "DLE_DLE1_B_AGCGGTCCATGGAATA-1" "DLE_DLE1_B_ACGCCAGCACAGACTT-1"
#[3] "DLE_DLE1_B_TTCTACAGTCTCACCT-1" "DLE_DLE1_B_CACATTTAGCGGCTTC-1"
  
  #2：合并并添加一些clone信息
  scBCR_RNA <- combineExpression(data_bcr_rename, #CombineTCR()或者 CombineBCR()产生的对象，或者CombineBCR()(CombineTCR()) list
                                 seuratdata,#barcode与BCR或者TCR对应的单细胞Seurat或者SingleCellExperiment对象
                                 cloneCall="aa", 
                                 cloneSize =c('1'=1, '2'=3, '3'=3, '4'=4, '5'=5,
                                              '6'=6, '7'=7,'8'=8,'9'=9,'>9'=10),#克隆型的频率，分割区间，可自行设置
                                 proportion = FALSE)#如果是T，则cloneTypes分割数值使用百分比
 
  meta <-scBCR_RNA@meta.data
   #合并之后你会发现，metadata中有很多的NA，这是因为并不是所有的细胞中都测到了BCR、TCR
  #况且，在质控的时候，还筛选掉一些细胞，所以才会出现这种情况，不必担忧！
  #添加cell是否有BCR信息
  scBCR_RNA@meta.data$BCR.known <- "yes"
  scBCR_RNA@meta.data$BCR.known[is.na(scBCR_RNA@meta.data$cloneSize)] <-  "no"
  
  #添加cell是否有双链
  scBCR_RNA@meta.data$both_chains <- "yes"
  scBCR_RNA@meta.data$both_chains[grep("NA",scBCR_RNA@meta.data$CTgene)] <- "no"
  scBCR_RNA@meta.data$both_chains[is.na(scBCR_RNA@meta.data$CTaa)]<- "no"
    
  #添加Isotype信息
  scBCR_RNA@meta.data$TYPEs <- scBCR_RNA@meta.data$CTgene
  scBCR_RNA@meta.data$TYPEs[is.na(scBCR_RNA@meta.data$TYPEs)] <- "NA_NA"
  
  heavychains <- unlist(strsplit(scBCR_RNA@meta.data$TYPEs, "[_]"))[seq(1, length(unlist(strsplit(scBCR_RNA@meta.data$TYPEs, "[_]"))), 2)]
  heavychains <- as.data.frame(heavychains)
  heavychains$Isotype <- "unknown"
  heavychains$Isotype <- str_split(heavychains$heavychains,'\\.',simplify = T)[,4]
  heavychains$Isotype[which(heavychains$heavychains=='NA')] <- "unknown"
  scBCR_RNA@meta.data$Isotype <- heavychains$Isotype
  remove(heavychains)
  save(scBCR_RNA,file="../output/04.subcell/02.cluster/scBCR.RData")
  
  DimPlot(scBCR_RNA, group.by = "cloneSize",label = F)
  DimPlot(scBCR_RNA, group.by = "Isotype",label=F)
  
  sink('../output/04.subcell/05.scRepertoire/manual_l2_Isotype.txt',append = F)
  print('这个可以看出每个细胞亚群的克隆特点');
  table(meta$Isotype,meta$manual_L2);
  print('这个可以看出每个样本的克隆特点');
  table(meta$Isotype,meta$orig.ident);
  print('每个细胞亚群每个样本的克隆点');
  table(meta$Isotype,meta$orig.ident,meta$manual_L2);sink()
  
# scBCR_RNA@meta.data$group <- mapvalues(scBCR_RNA@meta.data$orig.ident,
#  from=c("HC1","HC2","HC3","ICB1","ICB2","IEI"),
#  to=c(rep("HC",3), rep("Disease",3))) 这个函数有点意思

  #3：可视化
  
  #3.1 clonalOverlay---可视化克隆扩增细胞在UMAp的位置
  p1 <- clonalOverlay(scBCR_RNA, reduction = "umap", 
                       cutpoint = 10, 
                       bins = 100, 
                       facet = "Isotype") 
  p1
  ggsave('../output/04.subcell/05.scRepertoire/clone_umap_by_isotype.pdf',
         p1,width =8,height = 6
        )
  p2 <- clonalOverlay(scBCR_RNA, reduction = "umap", 
                      cutpoint = 10, 
                      bins = 100, 
                      facet = "condition") 
  p2
  ggsave('../output/04.subcell/05.scRepertoire/clone_umap_by_condition.pdf',
         p2,width = 8,height =7
  )
  #3.2 clonalNetwork---与clonalOverlay()类似，使用clonalNetwork()观察单细胞集群之间共享的克隆型的网络相互作用
  #该函数显示来自开始节点的克隆的相对比例，结束节点由箭头表示。
  
  #所有细胞的互作
  p2 <- clonalNetwork(scBCR_RNA, reduction = "umap", 
                group.by= "ident",
                filter.clones = 1000,
                filter.identity = NULL,
                filter.graph=T,
                cloneCall = "aa",
              #  exportClones = T,
                #exportTable = T
              )
  p2
  ggsave('../output/04.subcell/05.scRepertoire/clone_reaction.pdf',
         p2,width = 5,height = 6)
         
  df1 <- clonalNetwork(scBCR_RNA, reduction = "umap", 
                group.by= "ident",
                filter.clones = 1000,
                filter.identity = NULL,
                filter.graph=T,
                cloneCall = "aa",
                exportClones = T
              )
  df2 <- clonalNetwork(scBCR_RNA, reduction = "umap", 
                       group.by= "ident",
                       filter.clones = 1000,
                       filter.identity = NULL,
                       filter.graph=T,
                       cloneCall = "aa",
                       exportTable = T
                       )
  sink('../output/04.subcell/05.scRepertoire/clonalNetwork.txt',append = F)
  print('这个可以看出克隆性的拷贝多少');
  print(df1,n=100);
  print('这个可以看出克隆的分布特点');
  print(df2,n=100);sink()
  
  #指定细胞类型
  p3 <- clonalNetwork(scBCR_RNA, 
                      reduction = "umap",
                      filter.clones = 100,
                      group.by  = "ident",
                      filter.identity = "07_B1_SPN+",
                      cloneCall = "aa")
  
  
  p3
  ggsave('../output/04.subcell/05.scRepertoire/clonalNetwork_B1.pdf',
         p3,width = 5,height = 6)
  
  #3.3 查看不同细胞亚群中clone type的分布
  p4 <- clonalOccupy(scBCR_RNA, x.axis = "ident",label = F,proportion=T)+
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  p4
  ggsave('../output/04.subcell/05.scRepertoire/clonesiZe_stacked_barplot.pdf',
         p4,width = 5,height = 6)
  #将样本与celltype结合可视化
  scBCR_RNA$condition_celltype <- paste0(scBCR_RNA$condition, "_", scBCR_RNA$manual_L2)
  p5 <- clonalOccupy(scBCR_RNA, x.axis = "condition_celltype",label = F,proportion=T)+
    scale_fill_npg()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  #3.4 alluvialClonotypes,桑葚图展示celltype与Clontype的关系
  p6 <- alluvialClones(scBCR_RNA, cloneCall = "gene", 
                            y.axes = c("manual_L2", "cloneSize"), 
                            color = "manual_L2")
  p6
  ggsave('../output/04.subcell/05.scRepertoire/alluvial_plot.pdf',
         p6,width = 5,height = 6)
  #3.5 getCirclize
  #getCirclize使用弦图来展示celltype之间的互联。它的作用和上面的alluvialClonotypes是差不多的
  circles <- getCirclize(scBCR_RNA, group.by = "ident")
  grid.cols <- scales::hue_pal()(length(unique(scBCR_RNA@active.ident)))#颜色设置
  names(grid.cols) <- levels(scBCR_RNA@active.ident)
  p7 <- circlize::chordDiagram(circles,
                                self.link = 1, 
                                grid.col = grid.cols) 
  
  #3.6 StartracDiversity
  #这里是联合使用startrac包的功能，
  p8 <- StartracDiversity(scBCR_RNA, 
                           type = "manual_L2", 
                           group.by = "condition",
                          #exportTable = T
                          )+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave('../output/04.subcell/05.scRepertoire/StartracDiversity.pdf',
         p8,width = 5,height = 6)
  #输出的图上参数意义
  # expa - Clonal Expansion，克隆扩增：指免疫细胞（如T细胞和B细胞）在识别到特定抗原后，通过有丝分裂进行繁殖，从而增加特异性免疫细胞数量的过程。
  # migr - Cross-tissue Migration：跨组织迁移
  # tran - State Transition：状态转化
  
  #可以画不同细胞亚群的免疫组库谱
  df_mix_celltype_isotype <- table(scBCR_RNA$Isotype,scBCR_RNA@meta.data$manual_L2)
  write.csv(df_mix_celltype_isotype,'../output/04.subcell/05.scRepertoire/mix_celltype_isotype.csv')
  }










