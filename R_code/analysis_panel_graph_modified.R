library(ggplot2)
library(readxl)
library(stringr)

library(grid)
library(gridExtra)
library(reshape2)
library("igraph")

## FOR REGION NUM
panel1.data <- read.delim("~/Desktop/FINAL_FYT/region_num/plot_1.modified.data.txt")
panel1.data <- panel1.data[,1:4]
panel2.data <- read.delim("~/Desktop/FINAL_FYT/region_num/plot_2.modified.data.txt")
panel3.data <- read.delim("~/Desktop/FINAL_FYT/region_num/plot_3.modified.data.txt")
panel4.data <- read.delim("~/Desktop/FINAL_FYT/region_num/plot_4.modified.data.txt")
data_interaction <- read.delim("~/Desktop/FINAL_FYT/region_num/traj.matrix.txt")


# #### FOR CCF
# panel1.data <- read.delim("~/Desktop/FINAL_FYT/CCF/plot_1.modified.data.txt")
# panel1.data <- panel1.data[,1:4]
# panel2.data <- read.delim("~/Desktop/FINAL_FYT/CCF/plot_2.modified.data.txt")
# panel3.data <- read.delim("~/Desktop/FINAL_FYT/CCF/plot_3.modified.data.txt")
# panel4.data <- read.delim("~/Desktop/FINAL_FYT/CCF/plot_4.modified.data.txt")
# data_interaction <- read.delim("~/Desktop/FINAL_FYT/CCF/traj.matrix.txt")


# #### First panel pre work
gene_interaction <- colnames(data_interaction)
gene.panel1 <- rep('N',length(gene_interaction))
for(i in 2:length(unique(panel1.data[,2]))){
  gene.panel1 <- cbind(gene.panel1,rep('N',length(gene_interaction)))
}
colnames(gene.panel1) <- c(unique(panel1.data[,2]))
rownames(gene.panel1) <- c(gene_interaction)

fetB <- unique(panel1.data[,2])

for(i in 1:length(gene_interaction)){
  for(j in 1:length(fetB)){
    normal_red <- panel1.data[which(panel1.data$FeatureA == gene_interaction[i] & panel1.data$FeatureB == fetB[j] & panel1.data$dotColor == 'A_red'),]
    normal_grey <- panel1.data[which(panel1.data$FeatureA == gene_interaction[i] & panel1.data$FeatureB == fetB[j] & panel1.data$dotColor == 'B_grey'),]
    if(nrow(normal_red)	> 0 ){
      gene.panel1[i,j] <- 'A_red'
    }
    if(nrow(normal_grey)	> 0 ){
      gene.panel1[i,j] <- 'B_grey'
    }
    normal_red <- 0
    normal_grey <- 0
  }
}

# #### second panel pre work
raw_data_gene <- read.delim("~/Desktop/secondsem/dnds.result.txt")
sel_gene <- raw_data_gene$geneID
gene.panel2 <- rep('N',length(sel_gene))
for(i in 2:length(unique(panel2.data[,2]))){
  gene.panel2 <- cbind(gene.panel2,rep('N',length(sel_gene)))
}
colnames(gene.panel2) <- c(unique(panel2.data[,2]))
rownames(gene.panel2) <- c(sel_gene)

fetB2 <- unique(panel2.data[,2])

for(i in 1:length(sel_gene)){
  for(j in 1:length(fetB2)){
    normal_red <- panel2.data[which(panel2.data$FeatureA == sel_gene[i] & panel2.data$FeatureB == fetB2[j] & panel2.data$dotColor == 'A_red'),]
    normal_grey <- panel2.data[which(panel2.data$FeatureA == sel_gene[i] & panel2.data$FeatureB == fetB2[j] & panel2.data$dotColor == 'B_grey'),]
    if(nrow(normal_red)	> 0 ){
      gene.panel2[i,j] <- 'A_red'
    }
    if(nrow(normal_grey)	> 0 ){
      gene.panel2[i,j] <- 'B_grey'
    }
    normal_red <- 0
    normal_grey <- 0
  }
}

#### third panel pre work
gene_interaction <- colnames(data_interaction)
gene.panel3 <- rep('N',length(gene_interaction))
for(i in 2:length(unique(panel3.data[,2]))){
  gene.panel3 <- cbind(gene.panel3,rep('N',length(gene_interaction)))
}
colnames(gene.panel3) <- c(unique(panel3.data[,2]))
rownames(gene.panel3) <- c(gene_interaction)

fetB3 <- unique(panel3.data[,2])


for(i in 1:length(gene_interaction)){
  for(j in 1:length(fetB3)){
    normal_red <- panel3.data[which(panel3.data$listA == gene_interaction[i] & panel3.data$listB == fetB3[j] & panel3.data$dotColor == 'A_red'),]
    normal_grey <- panel3.data[which(panel3.data$listA == gene_interaction[i] & panel3.data$listB == fetB3[j] & panel3.data$dotColor == 'B_grey'),]
    if(nrow(normal_red)	> 0 ){
      gene.panel3[i,j] <- 'A_red'
    }
    if(nrow(normal_grey)	> 0 ){
      gene.panel3[i,j] <- 'B_grey'
    }
    normal_red <- 0
    normal_grey <- 0
  }
}


#### Fourth panel pre work
gene.panel4 <- rep('N',length(sel_gene))
for(i in 2:length(unique(panel4.data[,2]))){
  gene.panel4 <- cbind(gene.panel4,rep('N',length(sel_gene)))
}
colnames(gene.panel4) <- c(unique(panel4.data[,2]))
rownames(gene.panel4) <- c(sel_gene)

fetB4 <- unique(panel4.data[,2])


for(i in 1:length(sel_gene)){
  for(j in 1:length(fetB4)){
    normal_red <- panel4.data[which(panel4.data$listA == sel_gene[i] & panel4.data$listB == fetB4[j] & panel4.data$dotColor == 'A_red'),]
    normal_grey <- panel4.data[which(panel4.data$listA == sel_gene[i] & panel4.data$listB == fetB4[j] & panel4.data$dotColor == 'B_grey'),]
    if(nrow(normal_red)	> 0 ){
      gene.panel4[i,j] <- 'A_red'
    }
    if(nrow(normal_grey)	> 0 ){
      gene.panel4[i,j] <- 'B_grey'
    }
    normal_red <- 0
    normal_grey <- 0
  }
}


#### First panel drawing
plot_1.data <- melt(t(gene.panel1))
temp_arr <- panel1.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel1.data[,3])
plot_1.data <- cbind(plot_1.data, dotSize = 0)
for(i in 1:nrow(plot_1.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_1.data[i,1] == temp_arr[j,2] & plot_1.data[i,2] == temp_arr[j,1]){
      plot_1.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_3.data <- melt(t(gene.panel3))
temp_arr <- panel3.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel3.data[,3])
plot_3.data <- cbind(plot_3.data, dotSize = 0)
for(i in 1:nrow(plot_3.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_3.data[i,1] == temp_arr[j,2] & plot_3.data[i,2] == temp_arr[j,1]){
      plot_3.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_1.data <- rbind(plot_1.data, plot_3.data)
x.scale.111 <- as.character(unique(plot_1.data$Var1))


F1.plot<-ggplot()+theme_classic()

F1.plot<-F1.plot+geom_point(data = plot_1.data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
F1.plot<-F1.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                       legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                       axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1.plot<-F1.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                     labels=c(A_red='Significant',B_grey='Not significant'),
                                     guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
F1.plot<-F1.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)

F1_1.plot<-ggplot()+theme_classic()

F1_1.plot<-F1_1.plot+geom_point(data = plot_1.data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
F1_1.plot<-F1_1.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                       legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                       axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1_1.plot<-F1_1.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                     labels=c(A_red='Significant',B_grey='Not significant'),
                                     guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
F1_1.plot<-F1_1.plot+scale_x_discrete(breaks=x.scale.111,position = "top")+xlab(NULL)+ylab(NULL)
#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-1_merge_landscape.pdf", plot=F1_1.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8A_analysis_landscape_CCF.pdf", plot=final_figure,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

############################
#### transform
x.scale.11 <<- as.character(gene_interaction)
plot_1_1.data <- melt(gene.panel1)
temp_arr <- panel1.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel1.data[,3])
plot_1_1.data <- cbind(plot_1_1.data, dotSize = 0)
for(i in 1:nrow(plot_1_1.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_1_1.data[i,1] == temp_arr[j,1] & plot_1_1.data[i,2] == temp_arr[j,2]){
      plot_1_1.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_3_1.data <- melt(gene.panel3)
temp_arr <- panel3.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel3.data[,3])
plot_3_1.data <- cbind(plot_3_1.data, dotSize = 0)
for(i in 1:nrow(plot_3_1.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_3_1.data[i,1] == temp_arr[j,1] & plot_3_1.data[i,2] == temp_arr[j,2]){
      plot_3_1.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_1_1.data <- rbind(plot_1_1.data, plot_3_1.data)

F1_2.plot<-ggplot()+theme_classic()

F1_2.plot<-F1_2.plot+geom_point(data = plot_1_1.data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
F1_2.plot<-F1_2.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                           legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                           axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1_2.plot<-F1_2.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                         labels=c(A_red='Significant',B_grey='Not significant'),
                                         guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
F1_2.plot<-F1_2.plot+scale_x_discrete(breaks=x.scale.11,position = "top")+xlab(NULL)+ylab(NULL)
#### FOR REGION SUM
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-1_transformed_landscape.pdf", plot=F1_2.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-1_transformed_landscape_CCF.pdf", plot=F1_2.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)



#### Second panel drawing
plot_2.data <- melt(t(gene.panel2))
temp_arr <- panel2.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel2.data[,3])
plot_2.data <- cbind(plot_2.data, dotSize = 0)
for(i in 1:nrow(plot_2.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_2.data[i,1] == temp_arr[j,2] & plot_2.data[i,2] == temp_arr[j,1]){
      plot_2.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_4.data <- melt(t(gene.panel4))
temp_arr <- panel4.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel4.data[,3])
plot_4.data <- cbind(plot_4.data, dotSize = 0)
for(i in 1:nrow(plot_4.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_4.data[i,1] == temp_arr[j,2] & plot_4.data[i,2] == temp_arr[j,1]){
      plot_4.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_2.data <- rbind(plot_2.data, plot_4.data)
x.scale.121 <- as.character(unique(plot_2.data$Var1))


F2.plot<-ggplot()+theme_classic()

F2.plot<-F2.plot+geom_point(data = plot_2.data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
F2.plot<-F2.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                       legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                       axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F2.plot<-F2.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                     labels=c(A_red='Significant',B_grey='Not significant'),
                                     guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
F2.plot<-F2.plot+scale_x_discrete(breaks=x.scale.121,position = "top")+xlab(NULL)+ylab(NULL)
#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-2_merge_landscape.pdf", plot=F2.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8A_analysis_landscape_CCF.pdf", plot=final_figure,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

############################
#### transform
x.scale.12 <<- as.character(sel_gene)
plot_2_1.data <- melt(gene.panel2)
temp_arr <- panel2.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel2.data[,3])
plot_2_1.data <- cbind(plot_2_1.data, dotSize = 0)
for(i in 1:nrow(plot_2_1.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_2_1.data[i,1] == temp_arr[j,1] & plot_2_1.data[i,2] == temp_arr[j,2]){
      plot_2_1.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_4_1.data <- melt(gene.panel4)
temp_arr <- panel4.data[,1:2]
temp_arr <- cbind(temp_arr, V3=panel4.data[,3])
plot_4_1.data <- cbind(plot_4_1.data, dotSize = 0)
for(i in 1:nrow(plot_4_1.data)){
  for(j in 1:nrow(temp_arr)){
    if(plot_4_1.data[i,1] == temp_arr[j,1] & plot_4_1.data[i,2] == temp_arr[j,2]){
      plot_4_1.data[i,4] <- temp_arr[j,3]
    }
  }
}
temp_arr <- 0
plot_2_1.data <- rbind(plot_2_1.data, plot_4_1.data)


F2_1.plot<-ggplot()+theme_classic()

F2_1.plot<-F2_1.plot+geom_point(data = plot_2_1.data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
F2_1.plot<-F2_1.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                       legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                       axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F2_1.plot<-F2_1.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                     labels=c(A_red='Significant',B_grey='Not significant'),
                                     guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
F2_1.plot<-F2_1.plot+scale_x_discrete(breaks=x.scale.12,position = "top")+xlab(NULL)+ylab(NULL)
#### FOR REGION SUM
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-2_transformed_landscape.pdf", plot=F2_1.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-2_transformed_landscape_CCF.pdf", plot=F2_1.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

### merge all graphs
final_figure<-rbind(ggplotGrob(F2.plot),ggplotGrob(F1.plot),size="last")
panels <- final_figure$layout$t[grep("panel", final_figure$layout$name)]

#### FOR REGION SUM
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8_no.._analysis_landscape.pdf", plot=final_figure,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8_no.._analysis_landscape_CCF.pdf", plot=final_figure,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

### FINAL GRAPH-1
final_plot_data <- rbind(plot_2.data, plot_1.data)
x.scale.2 <<- as.character(unique(final_plot_data$Var1))

final.plot<-ggplot()+theme_classic()

final.plot<-final.plot+geom_point(data = final_plot_data,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
final.plot<-final.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                       text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                       legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                       axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
final.plot<-final.plot+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                     labels=c(A_red='Significant',B_grey='Not significant'),
                                     guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
final.plot<-final.plot+scale_x_discrete(breaks=x.scale.2,position = "top")+xlab(NULL)+ylab(NULL)

#### FOR REGION SUM
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8A_analysis_landscape.pdf", plot=final.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8A_analysis_landscape_CCF.pdf", plot=final.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)


### FINAL GRAPH-2
final_plot_data_1 <- rbind(plot_2_1.data, plot_1_1.data)
x.scale.21 <<- as.character(unique(final_plot_data_1$Var1))

final.plot_1<-ggplot()+theme_classic()

final.plot_1<-final.plot_1+geom_point(data = final_plot_data_1,aes(x=Var1,y=Var2,fill=value,size=dotSize,color=as.factor(value)))
final.plot_1<-final.plot_1+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                             text=element_text(size=10,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.1,'cm'),legend.key.height=unit(0.1,'cm'),legend.position='right',
                             legend.direction="horizontal",legend.text=element_text(size=5,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                             axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
final.plot_1<-final.plot_1+scale_colour_manual(name=NULL,values = c(B_grey='grey',A_red='red'),
                                           labels=c(A_red='Significant',B_grey='Not significant'),
                                           guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
final.plot_1<-final.plot_1+scale_x_discrete(breaks=x.scale.21,position = "top")+xlab(NULL)+ylab(NULL)

#### FOR REGION SUM
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8B_analysis_landscape.pdf", plot=final.plot_1,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8B_analysis_landscape_CCF.pdf", plot=final.plot_1,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)


x.scale <<- NULL
x.scale.2 <<- NULL
