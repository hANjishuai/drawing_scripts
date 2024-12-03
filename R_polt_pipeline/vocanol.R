#Funcation:make a colorful volcano plot
#Owner: Fanghan Ji
#History:
#Version:0.0.1 2024/12/03 make a tyoical volcano plot

#loading packages
library(ggplot2)

# loading data
data <- read.csv('../../text/data02.csv',header = T,row.names = 1)

# adding meta columns
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))

ggplot(data,aes(log2FoldChange,-log10(padj)))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#999999")+
  geom_vline(xintercept = c(-1.2,1.2),linetype="dashed",color="#999999")+
  geom_point(aes(size = -log10(padj),color = -log10(padj)))+
  scale_colour_gradientn(
    #name="q.value",
    values=seq(0,1,0.2),
    colours=c("#39489f","#39bbec","#f9ed36",
                                                   "#f38466","#b81f25"))+
  scale_size_continuous(range=c(1,3))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.8),
        legend.justification = c(0,1))+
  guides(colour=guide_colourbar(title="-Log10_q-value"),
         size="none")+
  geom_text(aes(label = label,colour = -log10(padj)),
            size=2,
            vjust=1,
            hjust=1,
            nudge_x = -0.5,
            check_overlap = T)+
  labs(x="Log2FC",
      y="-Log10(FDR q-value)")
  
  ggsave("vocanol.pdf",width = 10,height = 9)
  
  
  
  
  