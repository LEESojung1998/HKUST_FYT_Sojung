library(ggplot2)
library(readxl)

library(grid)
library(gridExtra)
library(reshape2)
library("igraph")
library(stringr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}


#### Import necessary datasets
raw_data_clinical <- read_excel("~/Desktop/FYT/EGAS00001002247_clinical_data.xlsx")
#raw_data_gene <- read.delim("~/Desktop/firstsem/finding/myowndriver.11281.txt")
raw_data_gene <- read.delim("~/Desktop/secondsem/dnds.result.txt")
#raw_data_interaction <- read.delim("~/Desktop/secondsem/result/FINAL_TEDGedge.modified_CPN.txt")
raw_data_interaction <- read.delim("~/Desktop/FINAL_FYT/CCF/TEDGedge.modified.CCF.txt")
#raw_data_interaction <- read.delim("~/Desktop/secondsem/final+result/FINAL_TEDGedge.modified.txt")
raw_data_all <- read.delim("~/Desktop/secondsem/landscape.input.txt")
raw_data_all_info <- raw_data_all[which(raw_data_all$exonic.func != "NA"  & raw_data_all$exonic.func != "synonymousSNV" & raw_data_all$exonic.func != "unknown" & raw_data_all$exonic.func != "synonymous"),]
raw_subclonal_clonal <- read_excel("~/Desktop/secondsem/subclonal_clonal.xlsx")


sel_gene <- raw_data_gene$geneID
case <- raw_data_clinical$TRACERxID


#### Mutation Table based on mutation types
mut.stopgainSNV <- rep(0,length(case))
mut.stoplossSNV <- rep(0,length(case))
mut.nonsynon <- rep(0,length(case))
mut.nonsynonSNV <- rep(0, length(case))
mut.insert <- rep(0, length(case))
mut.nonframe.substit <- rep(0, length(case))
mut.nonframe.insert <- rep(0, length(case))
mut.sub <- rep(0, length(case))

for (i in 1:length(case)){
  mut.stopgainSNV[i] <- length(which(raw_data_all_info$SampleID == case[i] & (raw_data_all_info$exonic.func == 'stopgainSNV' | raw_data_all_info$exonic.func == 'immediate-stopgain')))
  mut.stoplossSNV[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'stopgainSNV'))
  mut.nonsynon[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'nonsynonymousSNV'))  
  mut.nonsynonSNV[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'nonsynonymous'))  
  mut.insert[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'frameshiftinsertion'))  
  mut.nonframe.substit[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'nonframeshiftsubstitution'))  
  mut.nonframe.insert[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'nonframeshiftinsertion'))  
  mut.sub[i] <- length(which(raw_data_all_info$SampleID == case[i] & raw_data_all_info$exonic.func == 'frameshiftsubstitution'))
}

mutation_num_table <- data.frame(case, mut.stoplossSNV, mut.stopgainSNV, mut.nonsynon, mut.nonsynonSNV, mut.insert, mut.nonframe.substit, mut.nonframe.insert, mut.sub)
colnames(mutation_num_table)<-c('Patients','stoploss', 'stopgain','nonsynonymous_DNV','nonsynonymous_SNV', 'frameshift_insertion', 'frameshift_substitution', 'nonframeshift_substitution', 'nonframeshift_insertion')

#### Gene Matrix
gene.Matrix <- rep('N',length(case))
for(i in 2:length(sel_gene)){
  gene.Matrix <- cbind(gene.Matrix,rep('N',length(case)))
}

data.Sel <- raw_data_all_info[which(raw_data_all_info$Hugo_Symbol %in% sel_gene),]

colnames(gene.Matrix) <- sel_gene
rownames(gene.Matrix) <- case 

for(i in 1:length(case)){
  for(j in 1:length(sel_gene)){
    temp.nonsynonSNV <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'nonsynonymousSNV'),]
    temp.stoplossSNV <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'stoplossSNV'),]
    temp.stopgainSNV <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & (data.Sel$exonic.func == 'stopgainSNV' | data.Sel$exonic.func == 'immediate-stopgain')),]
    temp.nonsynon <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'nonsynonymous'),]
    temp.insert <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'frameshiftinsertion'),]
    temp.nonframe.substit <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'nonframeshiftsubstitution'),]
    temp.nonframe.insert <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'nonframeshiftinsertion'),]
    temp.sub <- data.Sel[which(data.Sel$SampleID == case[i] & data.Sel$Hugo_Symbol == sel_gene[j] & data.Sel$exonic.func == 'frameshiftsubstitution'),]
    if(nrow(temp.nonsynonSNV)	> 0 ){
      gene.Matrix[i,j] <- 'n.s.s'
    }
    if(nrow(temp.stopgainSNV)	> 0 ){
      gene.Matrix[i,j] <- 'sg'
    }
    if(nrow(temp.nonsynon)	> 0 ){
      gene.Matrix[i,j] <- 'n.s'
    }
    if(nrow(temp.insert)	> 0 ){
      gene.Matrix[i,j] <- 'f.insert'
    }
    if(nrow(temp.nonframe.substit)	> 0 ){
      gene.Matrix[i,j] <- 'nf.sub'
    }
    if(nrow(temp.nonframe.insert)	> 0 ){
      gene.Matrix[i,j] <- 'nf.insert'
    }
    if(nrow(temp.sub)	> 0 ){
      gene.Matrix[i,j] <- 'f.sub'
    }
    if(nrow(temp.stoplossSNV)	> 0 ){
      gene.Matrix[i,j] <- 'sl'
    }
    temp.nonsynonSNV <- 0
    temp.stopgainSNV <- 0
    temp.nonsynon <- 0
    temp.insert <- 0
    temp.nonframe.substit <- 0
    temp.nonframe.insert <- 0
    temp.sub <- 0
    temp.stoplossSNV <- 0
  }
}

#### Gene Matrix_clonal/subclonal
clonal.Matrix <- t(raw_subclonal_clonal)
clonal.Matrix <- clonal.Matrix[1,which(clonal.Matrix)]
colnames(clonal.Matrix) <- clonal.Matrix[1,]
clonal.Matrix <- clonal.Matrix[-1,]
clonal.Matrix.2 <- rep('N',length(case))
for(i in 2:length(sel_gene)){
  clonal.Matrix.2 <- cbind(clonal.Matrix.2,rep('N',length(case)))
}
colnames(clonal.Matrix.2) <- sel_gene

colnames(clonal.Matrix.2) <- colnames(clonal.Matrix)
for(i in 1:ncol(clonal.Matrix.2)){
  for(j in 1:ncol(clonal.Matrix)){
    if(colnames(clonal.Matrix.2)[i] == colnames(clonal.Matrix)[j]){
      clonal.Matrix.2[,i] <- clonal.Matrix[,j]
    }
  }
}
for(i in 1:ncol(clonal.Matrix.2)){
  for(j in 1:nrow(clonal.Matrix.2)){
    if(isTRUE(str_detect(clonal.Matrix.2[j,i], "sub"))){
      clonal.Matrix.2[j,i] <- "sub"
    }
    else if(isTRUE(str_detect(clonal.Matrix.2[j,i], "clonal"))){
      clonal.Matrix.2[j,i] <- "clonal"
    }
    else{
      clonal.Matrix.2[j,i] <- "N"
    }
  }
}
### dealing with exception
clonal.Matrix.2[13,2] <- "clonal"

rownames(clonal.Matrix.2) <- case

#### Clinical Matrix
clinic.Matrix <- raw_data_clinical[,2]
for(i in 3:ncol(raw_data_clinical)){
  clinic.Matrix <- cbind(clinic.Matrix, raw_data_clinical[,i])
}
temo <- colnames(raw_data_clinical)
colnames(clinic.Matrix) <- temo[2:17]
rownames(clinic.Matrix) <- case 


#### Trajectory Matrix
traj.Matrix <- rep('N',length(case))
for(i in 2:nrow(raw_data_interaction)){
  traj.Matrix <- cbind(traj.Matrix,rep('N',length(case)))
}
temo2 <- matrix(nrow=nrow(raw_data_interaction), ncol=1)

for(x in 1:nrow(raw_data_interaction)){
  temo2[x] <- paste(raw_data_interaction[x,"geneA"], raw_data_interaction[x,"geneB"], sep = "->")
}
rownames(traj.Matrix) <- case 
colnames(traj.Matrix) <- c(temo2)


for(x in 1:nrow(raw_data_interaction)){
  temo1 <- unlist(strsplit(raw_data_interaction[x,"label"], ";"))
  for(y in 1:length(temo1)){
    for(z in rownames(traj.Matrix)){
      if(z == temo1[y]){
        traj.Matrix[z,x] <- 'Y'
      }
    }
  }
  temo1 <- NULL
}

write.table(gene.Matrix,"~/Desktop/secondsem/gene.matrix.txt",row.names = F, sep = '\t')
write.table(traj.Matrix,"~/Desktop/secondsem/traj.matrix.txt",row.names = F, sep = '\t')
write.table(traj.Matrix,"~/Desktop/FINAL_FYT/CCF/traj.matrix.txt",row.names = F, sep = '\t')


#### Plot the mutation landscape
plot_1a.data <- data.frame(rep(mutation_num_table$Patients,8), c( c(mutation_num_table$stopgain), c(mutation_num_table$nonsynonymous_DNV),c(mutation_num_table$nonsynonymous_SNV),c(mutation_num_table$frameshift_insertion),c(mutation_num_table$frameshift_substitution),c(mutation_num_table$nonframeshift_substitution),c(mutation_num_table$stopgain),c(mutation_num_table$stoploss)), 
                           c( rep('a_stopgain',length(mutation_num_table[,1])), rep('b_nonDNV',length(mutation_num_table[,1])), rep('c_nonSNV',length(mutation_num_table[,1])), rep('d_frameinsert',length(mutation_num_table[,1])), rep('e_framesub',length(mutation_num_table[,1])), rep('f_nonframesub',length(mutation_num_table[,1])), rep('g_nonframeinsert',length(mutation_num_table[,1])), rep('h_stoploss',length(mutation_num_table[,1])) ) )


colnames(plot_1a.data)<-c('Patients','Mutations','Groups')
plot_1a.data$Patients<-as.character(plot_1a.data$Patients)
plot_1a.data$Mutations<-as.numeric(plot_1a.data$Mutations)
plot_1a.data$Groups<-as.factor(plot_1a.data$Groups)

x.scale <<- as.character(mutation_num_table$Patients)


orderID <<- c(1:nrow(plot_1a.data))


F1a.plot<-ggplot()+theme_classic()
F1a.plot<-F1a.plot+geom_bar(data=plot_1a.data,aes(x=reorder(Patients,orderID),y=Mutations,fill=Groups),width=0.7,color="black",stat='identity',position=position_stack())
F1a.plot<-F1a.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.4,2,0.4,2),'lines'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.8,'cm'),legend.position=c(0.92,0.8),legend.text=element_text(size=16,face='bold.italic'),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='bold',color='black'),
                         axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.9,vjust=2,face='bold',color='black'))
F1a.plot<-F1a.plot+ggtitle(NULL)+xlab(NULL)+ylab('Number of somatic mutations')+scale_fill_manual(name=NULL,values=c(gg_color_hue(8)[7], gg_color_hue(8)[8],gg_color_hue(8)[6],gg_color_hue(8)[5],gg_color_hue(8)[4],gg_color_hue(8)[3],gg_color_hue(8)[2],gg_color_hue(8)[1]),labels=c('nonsynonymous SNV','stopgain','nonsynonymous', 'frameshift substitution', 'frameshift insertion', 'nonframeshift substitution', 'nonframeshift insertion', 'stoploss'))
F1a.plot<-F1a.plot+scale_y_continuous(expand=c(0,0),limits=c(0,15000),breaks=seq(0,1500,100))+scale_x_discrete(breaks=x.scale)

yaxis <- mutation_num_table[,2] + mutation_num_table[,3] + mutation_num_table[,4] + mutation_num_table[,5] +mutation_num_table[,6] + mutation_num_table[,7] +  mutation_num_table[,8] +  mutation_num_table[,9]
maxy <- max(yaxis)

if(maxy > 1000){
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
    theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
  split2 <- F1a.plot + coord_cartesian(ylim = c(1000, 1200)) + ggtitle(NULL)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
}else if(maxy > 800){
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
    theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
  split2 <- F1a.plot + coord_cartesian(ylim = c(800, 1000)) + ggtitle(NULL)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
}else if(maxy > 600){
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 390)) +
    theme(legend.position='none',plot.margin=unit(c(0.2,0.6,0.1,2),'lines')) +ylab(NULL)
  split2 <- F1a.plot + coord_cartesian(ylim = c(600, 800)) + ggtitle(NULL)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none',plot.margin=unit(c(1,0.6,0.2,2),'lines'))+ylab(NULL)
}else if(maxy > 400){
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 600)) +
    theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
}else if(maxy > 200){
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 400)) +
    theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
}else{
  split1 <- F1a.plot + coord_cartesian(ylim = c(0, 200)) +
    theme(legend.position='none',plot.margin=unit(c(1,0.6,0.1,2),'lines')) +ylab(NULL)
}


#### First panel - genomic data
plot_1b.data <- melt(gene.Matrix)

F1b.plot<-ggplot()+theme_classic()
F1b.plot<-F1b.plot+geom_tile(data = plot_1b.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
F1b.plot<-F1b.plot+scale_fill_manual(name=NULL,values=c(N='white',n.s.s = gg_color_hue(8)[8],sg = gg_color_hue(8)[7],n.s = gg_color_hue(8)[6], f.sub = gg_color_hue(8)[5], f.insert = gg_color_hue(8)[4], nf.sub = gg_color_hue(8)[3], nf.insert = gg_color_hue(8)[2], sl = gg_color_hue(8)[1]),labels=c(N = 'None', n.s.s = 'Nonsynonymous SNV',sg = 'Stopgain',n.s = 'Nonsynonymous DNV', f.sub = 'Frameshift Substitution', f.insert = 'Frameshift Insertion', nf.sub = 'Nonframeshift Substitution', nf.insert = 'Nonframeshift Insertion', sl = 'Stoploss'))
F1b.plot<-F1b.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='bottom',
                         legend.direction="horizontal",legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1b.plot<-F1b.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)

gene.scale = 1 + (ncol(gene.Matrix) - 12)/12

if(maxy > 600){
  figure_1<-rbind(ggplotGrob(split2),ggplotGrob(split1),ggplotGrob(F1b.plot),size="last")
  panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
  figure_1$heights[panels][1] <- unit(0.513,'null')
  figure_1$heights[panels][3] <- unit(gene.scale,'null')
}else{
  figure_1<-rbind(ggplotGrob(split1),ggplotGrob(F1b.plot),size="last")
  panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
  figure_1$heights[panels][3] <- unit(gene.scale,'null')
}

#### second panel - trajectory data
plot_1c.data <- melt(traj.Matrix)

F1c.plot<-ggplot()+theme_classic()
F1c.plot<-F1c.plot+geom_tile(data = plot_1c.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
F1c.plot<-F1c.plot+scale_fill_manual(name=NULL,values=c(N='white',Y='#8b0000'),labels=c(N = 'None',  Y = 'Observed'))
F1c.plot<-F1c.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='bottom',
                         legend.direction="horizontal",legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                         axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1c.plot<-F1c.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)

traj.scale = 1 + (ncol(traj.Matrix) - 12)

if(maxy > 600){
  figure_1<-rbind(ggplotGrob(split2),ggplotGrob(split1),ggplotGrob(F1b.plot),ggplotGrob(F1c.plot),size="last")
  panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
  figure_1$heights[panels][1] <- unit(0.513,'null')
  figure_1$heights[panels][3] <- unit(gene.scale,'null')
  figure_1$heights[panels][10] <- unit(traj.scale,'null')
}else{
  figure_1<-rbind(ggplotGrob(split1),ggplotGrob(F1b.plot),ggplotGrob(F1c.plot),size="last")
  panels <- figure_1$layout$t[grep("panel", figure_1$layout$name)]
  figure_1$heights[panels][3] <- unit(gene.scale,'null')
  figure_1$heights[panels][10] <- unit(traj.scale,'null')
}


grid.draw(figure_1)
ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure7-2_mutation_landscape_simplfied_CCF.pdf", plot=figure_1,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)


#### remove barplot
F1d.plot<-ggplot()+theme_classic()
F1d.plot<-F1d.plot+geom_tile(data = plot_1b.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
F1d.plot<-F1d.plot+scale_fill_manual(name=NULL,values=c(N='white',n.s.s = gg_color_hue(8)[8],sg = gg_color_hue(8)[7],n.s = gg_color_hue(8)[6], f.sub = gg_color_hue(8)[5], f.insert = gg_color_hue(8)[4], nf.sub = gg_color_hue(8)[3], nf.insert = gg_color_hue(8)[2], sl = gg_color_hue(8)[1]),labels=c(N = 'None', n.s.s = 'Nonsynonymous SNV',sg = 'Stopgain',n.s = 'Nonsynonymous DNV', f.sub = 'Frameshift Substitution', f.insert = 'Frameshift Insertion', nf.sub = 'Nonframeshift Substitution', nf.insert = 'Nonframeshift Insertion', sl = 'Stoploss'))
F1d.plot<-F1d.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='bottom',
                         legend.direction="horizontal",legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1d.plot<-F1d.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)


figure_2<-rbind(ggplotGrob(F1d.plot),ggplotGrob(F1c.plot),size="last")
panels <- figure_2$layout$t[grep("panel", figure_2$layout$name)]

grid.draw(figure_2)
ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure7-1_mutation_landscape_remove_bar_CCF.pdf", plot=figure_2,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

#### subclonal barplot

plot_1e.data <- melt(clonal.Matrix.2)
F1e.plot<-ggplot()+theme_classic()
F1e.plot<-F1e.plot+geom_tile(data = plot_1e.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
F1e.plot<-F1e.plot+scale_fill_manual(name=NULL,values=c(N='white',sub=gg_color_hue(2)[2],clonal=gg_color_hue(2)[1]),labels=c(N = 'None', sub='subclonal', clonal='clonal'))
F1e.plot<-F1e.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                         text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='bottom',
                         legend.direction="horizontal",legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
F1e.plot<-F1e.plot+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)


figure_3<-rbind(ggplotGrob(F1e.plot),ggplotGrob(F1c.plot),size="last")
panels <- figure_3$layout$t[grep("panel", figure_3$layout$name)]

grid.draw(figure_3)
ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure7_mutation_landscape_clonal_CCF.pdf", plot=figure_3,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)




x.scale <<- NULL
orderID <<- NULL 

#### Checking whether number of label are right
data.Sel2 <- matrix(ncol = 1, nrow = 1)
for(y in 1:length(sel_gene)){
  data.Sel2 <- append(data.Sel2, data.Sel[which(data.Sel$Hugo_Symbol == sel_gene[y]),3])
}



