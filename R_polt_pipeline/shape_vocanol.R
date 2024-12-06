#Funcation:make a shape volcano plot
#Owner: Fanghan Ji
#History:
#Version:0.0.1 2024/12/06 make a tyoical shape volcano plot

# loading needed packages
library(ggvolcano)
library(latex2exp)
library(ggrepel)

# loading the data
deg_data <- read.csv('~/Documents/text/01_colorful_vocano_plot.csv',row.names = 1)

# 任意指定一组分组变量 - 无实际意义
deg_data$group <- ifelse(sample(1:nrow(deg_data),replace = T)>nrow(deg_data)/2,
                         "Enriched in young",
                         "Enriched in elderly")

#-log10(padj)小于2的设置为灰色
deg_data$color <- ifelse(-log10(deg_data$padj)<2,"#bfc0c1",
                         ifelse(deg_data$log2FoldChange>0,'#e88182','#6489b2'))
head(deg_data)

# 绘图
ggplot(deg_data)+
  geom_point(aes(log2FoldChange,-log10(padj),
                 shape = group,
                 size=-log10(padj)),
             color=deg_data$color,
             alpha=0.7
             )+
  geom_vline(xintercept = 0, linetype="longdash")+
  geom_hline(yintercept = 2, linetype="longdash")+
  scale_size_continuous(range=c(0.2,3))+
  theme_classic()+
  theme(legend.position = 'none')+
  labs(x=latex2exp("$Log_2 \\textit{FC}$"),
       y=latex2exp("$-Log_{10} \\textit{FDR} $"))
ggsave("plot.pdf",height = 5,width = 6)

#添加标签
deg_data$label <- rep('',nrow(deg_data))
deg_data$label[order(deg_data$padj)[1:20]] <- rownames(deg_data)[order(deg_data$padj)[1:20]]

#加标签
ggplot(deg_data,aes(log2FoldChange,-log10(padj)))+
  geom_point(aes(log2FoldChange,-log10(padj),
                 shape = group,
                 size=-log10(padj)),
             color=deg_data$color,
             alpha=0.7
  )+
  geom_vline(xintercept = 0, linetype="longdash")+
  geom_hline(yintercept = 2, linetype="longdash")+
  geom_text_repel(aes(label=label),size=2,color=deg_data$color,max.overlaps = 100,
                 # key_glyph=draw_key_point
                  )+
  scale_size_continuous(range=c(0.2,3))+
  theme_classic()+
  theme(legend.position = 'none')+
  labs(x=latex2exp("$Log_2 \\textit{FC}$"),
       y=latex2exp("$-Log_{10} \\textit{FDR} $"))
ggsave("plot.pdf",height = 5,width = 6)





















