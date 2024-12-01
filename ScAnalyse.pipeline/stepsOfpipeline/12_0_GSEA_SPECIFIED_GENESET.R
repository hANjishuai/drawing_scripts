#这个脚本用于利用自制的基因集进行GSEA
#2024-10-25 Fanghan_Jay Version::0.0.1

#开始之前，想说一说自己的思路：、
#1、为什么要用自己的基因集？
#因为普通的富集分析，导致同一个通路的基因有上调的，和下调的，
#那么怎么评价这个通路到底是抑制，还是激活呢？
#2、该怎么解决上述的问题？
#首先，我人为的把差异基因分为上调的和下调的（根据差异分析的FC值）。
#其次，我对这两部分进行普通的富集分析，本人偏爱GO分析。
#然后，对其中的通路按照自己的生物学知识，注释为几类功能通路集
#毫无疑问，可以看见自己注释的几类功能通路集同时出现在上调区，和下调区
#为解决这个问题，我把相同的功能通路集的基因挑出来，作为gmt文件，进行GSEA分析。


#A、载入包
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
#B、设置环境变量
if(!dir.exists("../output/01HC-DLE-SLE_blood-result/04.subcell/04.GSEA_SPECIFED")){
  dir.create("../output/01HC-DLE-SLE_blood-result/04.subcell/04.GSEA_SPECIFED")
}

{#1
outname="CD5+SPN-_B_cells"
uppath="../output/01HC-DLE-SLE_blood-result/07.enrichment/C08_Bc_CD5+SPN+/common_up_df_among_condition_anno.csv"
downpath="../output/01HC-DLE-SLE_blood-result/07.enrichment/C08_Bc_CD5+SPN+/common_down_df_among_condition_anno.csv"
output='../output/01HC-DLE-SLE_blood-result/04.subcell/04.GSEA_SPECIFED'
}

{#2
  deg_path='../output/01HC-DLE-SLE_blood-result/06.DEG/CD5+_splited2SPN+or-/'#差异集因路径
  pattern='^B1'#保存差异集因的文件集
}

{#3
egmt_paths <- "../output/01HC-DLE-SLE_blood-result/04.subcell/04.GSEA_SPECIFED/"
input_name <- "CD5+SPN-_B_cells"
col_anno <- c("#F0E685FF","#CE3D32FF","#466983FF","#749B58FF","#5050FFFF")
gene_label_list <- c("HLA-C","B2M","CTSD","HLA-DRB5","CD74",#apc
                     "PRF1","KLRD1","IL2RB","JAK1","TGFB1",#bc
                     "LAT","FYB1","CCL5","ITGAL","ZYX",#che
                     "NKG7","KLRK1","ISG15","SYK","BTK",#cyto
                     "ZAP70","BCL11B","SPN","TBX21","STAT4")
}

#C、自定义函数
get_gmt <- function(outname,#保存时的basename
                    uppath,#加载的up通路表格路径
                    downpath,#加载的down通路表格路径
                    output){#输出的路径
  #1、载入之前注释好的富集通路表格
  #因为我这个是不同疾病状态下共同通路的变化，所以我只注释一种疾病状态下的通路注释集
  #进行数据清洗（根据自己的数据进行清洗）：
  common_pathway_UP <- read.csv(uppath,row.names = "X")
  common_pathway_UP <- common_pathway_UP %>% filter(common_pathway_UP$condtion=="DLE")# 
  common_pathway_DOWN <- read.csv(downpath,row.names = "X") %>% select(-1) 
  common_pathway_DOWN <- common_pathway_DOWN %>%  filter(common_pathway_DOWN$condtion=="DLE")#
  #合并两个文件
  common_pathway <- rbind(common_pathway_UP,common_pathway_DOWN)
  #因子化简易注释
  common_pathway$Annotation <- factor(common_pathway$Annotation,
                                      levels = c("A","B","C","D","E","O"),
                                      labels = c("Tc","Bc",
                                                 "Cytokine","APC",
                                                 "Chemotaxis","Others"))#
  common_pathway <- common_pathway %>% arrange(.,common_pathway$Annotation) %>%
    select(5,4)#
  #
  dat1 <- aggregate(common_pathway$geneID,
                    by = list(common_pathway$Annotation),
                    FUN = function(x){paste(x,collapse = "/")}) %>% #得到的是一个dataframe
    rename(term=Group.1,geneset=x)
  
  #将字符串分解为每个基因
  dat2 <- aggregate(dat1$geneset,
                    list(dat1$term),
                    FUN=function(x){
                      genelist <- str_split(x,"/",simplify = T) 
                    }) 
  
  list3 <- lapply(dat2$x,FUN = function(x){
    termset <- x %>% unlist() %>% matrix(ncol = 1) %>% unique()
  })
  
  #注释每一个基因集
  names(list3) <- dat2$Group.1
  #这里有个小技巧，给列表命名后，unlist得到的rowname也有命名哦
  df_geneset <- list3 %>% unlist %>% data.frame() %>% 
    rownames_to_column(var = "term") %>% rename(.,gene=`.`)
  df_geneset$term <- gsub("[0-9]$","",df_geneset$term) 
  df_geneset$term <- gsub("[0-9]$","",df_geneset$term) 
  
  #终于得到了gmt文件格式了，保存吧
  write.csv(df_geneset,file = paste0(output,"/",outname,"_gmt.csv"))
  
}

get_gsea <- function(deg_path,pattern){
  listpath <- list.files(deg_path,pattern = pattern,full.names = T)
  list_deg <- lapply(listpath,function(x){
    deg_file<- read.csv(x,row.names = "X")
    condition <- deg_file$condition %>% unique()#用于后面保存
    genelist <- deg_file$avg_log2FC
    names(genelist) <- deg_file$gene
    genelist <- sort(genelist,decreasing = T)
    genelist <- genelist[genelist != 0]
    set.seed(666)
    egmt <- GSEA(genelist,
                 TERM2GENE = geneset,
                 verbose = F,
                 nPerm=10000,
                 minGSSize = 10,#自定义的geneset某个基因通路很容易没有那么多基因，
                 maxGSSize = 1000,#如果设置太大size，egmt就不会有这个通路，这个时候
                 pvalueCutoff = 1)#gseaNB就得把term改成result$ID,
    gsea_results <- egmt@result#如果实在想绘制这个通路的gsea图片，就使用gseaplot2吧
    gsea_results2 <- gsea_results[order(gsea_results$enrichmentScore,decreasing = F),]
    write.csv(gsea_results2,file = paste0(output,"/",condition,"_gsearesult.csv"))
    save(egmt,file = paste0(output,"/",condition,"_egmt.RData"))
    })
}

get_scancode_plots1<- function(geneset){
  condition <- list.files(path=output,pattern = "_egmt.*$") %>% gsub("_.*$","",.)
  egmt_list <- list.files(path = output,pattern = "_egmt.*$",full.names = T)
  for (i in 1:length(egmt_list)) {
    print(i)
    name<-condition[i]
    load(egmt_list[[i]])
    terms <- egmt@result$ID#如果terms里面的通路，与egmt的result$ID不一致，就容易报错
    gsea_plots<-lapply(terms,function(x){
      gseaNb(object = egmt,
             geneSetID = x,
             subPlot = 2,
             addPval = T,
             pvalX = 0.75,
             pvalY = 0.75,
             pCol = 'black',
             pvalSize = 4,
             pHjust = 0,)
      #
    })
    save(gsea_plots,file = paste0(output,"/",name,"_gseaplot.RData"))
  }
}


#D、正式的pipeline
#1、自制GSEA数据集
get_gmt(outname,uppath,downpath,output)#使用该函数的时候看看数据格式，能要修改一些
geneset <- read.csv(paste0(output,"/",outname,"_gmt.csv"),row.names = "X")

#2、GSEA实现
get_gsea(deg_path = deg_path,pattern = pattern)

#3、可视化
##gesa图####
get_scancode_plots1(geneset = geneset)

#手动移到特定文件夹后
plot_list <- list.files(paste0(output,"/",outname),pattern = "_gseaplot.*$",full.names = T)
name<-list.files(paste0(output,"/",outname),pattern = "_gseaplot.*$") %>%
  gsub("_.*$","",.)

#提取RData后保存为pdf
for (i in 1:length(plot_list)) {
  Sys.sleep(1)
  load(plot_list[i])
  Name<-name[i]
  lapply(gsea_plots, function(x){
    plot<- x
    plotname <- plot[[1]]$labels$title
    polt_gg <- as.ggplot(plot = plot)
    ggsave(filename = paste0(output,"/",outname,"/",Name,"_",plotname,"_gseaplot.pdf"),
           plot = polt_gg,
           width = 7.5,
           height = 4.5)
  })
}

##雷达图####
#（1）构建数据：
egmt_list <- list.files(paste0(egmt_paths,input_name),
                        pattern = "_egmt.*$",full.names = T)

radar_df <- NULL
for (i in 1:length(egmt_list)) {
  egmtname <- egmt_list[i]
  condition_name <- basename(egmtname) %>% gsub("_.*","",.)
  load(egmtname)
  res_df <- egmt@result %>% filter(Description!="Others") %>% 
    select("Description","setSize","qvalue","p.adjust") %>% mutate(condition=condition_name,
                                                        celltype=input_name)
  rownames(res_df) <- NULL 
  radar_df <- rbind(radar_df,res_df)
}
write.csv(radar_df,file = paste0(output,"/",outname,"/","radar_df.csv"))
radar_df_sum <- aggregate(radar_df$setSize,by = list(radar_df$Description),sum)%>%
  mutate(celltype=input_name) %>% rename(.,Description="Group.1",score='x')
write.csv(radar_df_sum,file = paste0(output,"/",outname,"/","radar_df_sum.csv"))
assign(paste0(outname,"_radar_df"),radar_df_sum)

vs_radar_df <- rbind(`CD5+SPN-_B_cells_radar_df`,`CD5+SPN+_B_cells_radar_df`)
vs_radar_df$score <- scale(vs_radar_df$score,center = F,scale = T)
data <- merge(vs_radar_df[which(vs_radar_df$celltype==unique(vs_radar_df$celltype)[1]),],
              vs_radar_df[which(vs_radar_df$celltype==unique(vs_radar_df$celltype)[2]),],
              by="Description")
colnames(data)[2] <- data$celltype.x %>% unique()
colnames(data)[4] <- data$celltype.y %>% unique()
data <- data %>% select(c(1,2,4)) 
rownames(data) <- data$Description 
data <- data[,-1] %>% t() %>% data.frame()%>%rownames_to_column("celltype")#ggradar会把第一列作为axis_label
#作图
p2 <- ggradar(
  data[1,],
  axis.labels = rep(NA, 5),
  grid.min = 0, grid.mid = 1, grid.max = 2,
  # 雷达图线的粗细和颜色：
  group.line.width = 1,
  group.point.size = 2,
  group.colours = "#91CDC8",
  # 背景边框线颜色：
  background.circle.colour = "white",
  gridline.mid.colour = "#7491b7",
  gridline.max.colour = "black",
  legend.position = "none",
  # 不加坐标轴标签：
  label.gridline.min = F,
  label.gridline.mid = F,
  label.gridline.max = F
)+
  theme(plot.background = element_blank(),
        panel.background = element_blank())
p2

############## 圆环注释+文本注释
# 注释数据：
tmp <- data.frame(x = rep(1, 5),#几个通路：5
                  y = rep(1, 5),#不变：1
                  group = colnames(data)[-1])#要注释的几个通路名，记得取消第一列

# 绘图：
#设置注释rowbar颜色
#col_anno <- as.character(paletteer_d("ggsci::default_igv",n=5))

p1 <- ggplot()+
  # 圆环：
  geom_bar(data = tmp, aes(x, y, fill = group), stat = "identity", position = "dodge")+
  # 文本注释：
  geom_text(aes(x = rep(1,25), y = rep(2, 25),
                label = gene_label_list, group = 1:25,fontface = "bold"),
            color = "black", size = 2.5,
            position = position_dodge(width = 0.9))+
  # geom_text(aes(x, y, label = gsub("[.]", " ", df$group), group = group),
  #           color = "white",
  #           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = col_anno )+
  ylim(-5.5,2)+
  # 0.63是计算得来，5个色块，第一个色块的正中心要对准0的位置，
  # 所以2pi/10=0.628即为第一个色块左边界的位置
  coord_polar(start = -0.63)+
  theme_void()+
  theme(legend.position = "none")
p1

p1 + inset_element(p2, left = 0, bottom = 0, right = 0.99, top = 0.99)

ggsave(paste0(output,"/radar_plot.pdf"), height = 5, width = 5)


