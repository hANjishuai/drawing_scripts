#Funcation:make a typical barplot+scatterplot
#Owner: Fanghan Ji
#History:
#Version:0.0.1 2024/12/02 make a tyoical barplot & scatterplot

# create a vitural dataset
Bacteroides_A <- runif(23,min = -17, max = -12)
Bacteroides_B <- runif(23,min = -15.5, max = -10)
BMI_A <- runif(23, min = 3.45, max = 4.2)
BMI_B <- runif(23, min = 3.2, max = 4.0)
data <- cbind(Bacteroides_A,BMI_A,Bacteroides_B,BMI_B)

# load the packages you need
library(ggplot2)
library(tidyr)
library(dplyr)

# simply process on the data
new_data <- as.data.frame(rbind(data[1:23,3:4],data[1:23,1:2]))
colnames(new_data) <- c("Bacteroides","BMI")
new_data$group <- rep(c("group_0M","group_3M"),c(23,23))
new_data$group2 <- rep(paste0('sample',1:23),2)

# draw a scatter plot
ggplot(new_data) +
  geom_point(aes(Bacteroides,BMI,colour = group))+
  scale_color_manual(values=c(group_0M='#ff00ff',group_3M='#8ac53e'))+
  theme_classic()+
  scale_x_continuous(breaks = c(-17:11))+
  scale_y_continuous(breaks = seq(3.2,4.0,0.2))+
  labs(x="Bacteriods (11021)",y="BMI(lg)")+
  theme(legend.position = 'none')
ggsave(filename = "scatter_plot.pdf",width = 5,height = 5)

# get a boxplot leftside of the scatter plot
ggplot(new_data,aes(group,BMI))+
  stat_boxplot(geom = "errorbar",width=0.15,aes(color=group))+
  geom_boxplot(aes(color=group),fill='white')+
  geom_line(aes(group = group2),color='black',linetype="dashed",size=0.2,alpha=0.8)+
  geom_jitter(aes(color=group,fill=group),width = 0.05,shape=21)+
  scale_color_manual(values=c(group_0M='#ff00ff',group_3M='#8ac53e'))+
  scale_fill_manual(values=c(group_0M='#ff00ff',group_3M='#8ac53e'))+
  theme_classic()+
  scale_x_discrete(labels=c("0M","3M"))+
  scale_y_continuous(breaks = seq(3.2,4.0,0.2))+
  labs(x="",y="")+
  theme(legend.position = 'none')
ggsave(filename = "left_boxplot.pdf",width = 2,height = 5)

# get the boxplot which below the scatter plot
ggplot(new_data,aes(group,Bacteroides))+
  stat_boxplot(geom ='errorbar',width=0.05,aes(color=group))+
  geom_boxplot(aes(color=group),fill="white")+
  geom_line(aes(group = group2),color='black',linetype="dashed",size=0.2,alpha=0.8)+
  geom_jitter(aes(color=group,fill=group),shape=21,width = 0.05)+
  scale_color_manual(values=c(group_0M='#ff00ff',group_3M='#8ac53e'))+
  scale_fill_manual(values=c(group_0M='#ff00ff',group_3M='#8ac53e'))+
  theme_classic()+
  scale_x_discrete(labels=c("0M","3M"))+
  scale_y_continuous(breaks = c((-17):(-11)))+
  labs(x="",y="")+
  theme(legend.position = "none")+
  coord_flip()
ggsave(filename = "down_boxplot.pdf",width = 5,height = 2)
  
  
  
  

  




