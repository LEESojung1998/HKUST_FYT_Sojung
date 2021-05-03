library(ggplot2)
library(readxl)

#setwd("/Users/sandylee/Desktop/secondsem/final+result")

#### FOR REGION_NUM
setwd("/Users/sandylee/Desktop/FINAL_FYT/region_num")

#### FOR CCF
#setwd("/Users/sandylee/Desktop/FINAL_FYT/CCF")


#### Import necessary datasets
raw_data_clincal <- read_excel("~/Desktop/FYT/EGAS00001002247_clinical_data.xlsx")
raw_data_gene <- read.delim("~/Desktop/secondsem/dnds.result.txt")

#### FOR REGION_NUM
raw_data_interaction <- read.delim("~/Desktop/FINAL_FYT/region_num/TEDGedge.modified.txt")

#### FOR CCF
#raw_data_interaction <- read.delim("~/Desktop/FINAL_FYT/CCF/TEDGedge.modified.CCF.txt")

#### FOR REGION_NUM
data_interaction <- read.delim("~/Desktop/FINAL_FYT/region_num/traj.matrix.txt")

#### FOR CCF
#data_interaction <- read.delim("~/Desktop/FINAL_FYT/CCF/traj.matrix.txt")

data_gene <- read.delim("~/Desktop/FINAL_FYT/gene.matrix.txt")
data_gene[data_gene != "N"] <- "Y"

sel_gene <- raw_data_gene$geneID
case <- raw_data_clincal$TRACERxID
interaction_case <- colnames(data_interaction)


##### panel 1:clinical categorical vs direction
#step1: define categorical features and generate dataframe
categorical_features <- c("Stage","Gender","Ethnicity","Histology","Resection_margins","Vascular_invasion","Pleural_invasion","Adjuvant_therapy","ECOG",
                          "Smoking_status","Genome_doubled","Recurrence_or_death")
cate.table <- raw_data_clincal[,which(colnames(raw_data_clincal) == categorical_features[1])]
for(i in 2:length(categorical_features)){
  tmp <- raw_data_clincal[,which(colnames(raw_data_clincal) == categorical_features[i])]
  cate.table <- cbind(cate.table,tmp)
} 
colnames(cate.table) <- categorical_features

#step2: conduct chisq.test for every-pair features and summarize the pvalue results into a table
plot_1.data <- data.frame(rep(NA,length(categorical_features)*length(interaction_case)),rep(NA,length(categorical_features)*length(interaction_case)),rep(NA,length(categorical_features)*length(interaction_case)),rep(NA,length(categorical_features)*length(interaction_case)),rep(NA,length(categorical_features)*length(interaction_case)))
n<-1
#countfile <- 1
nfeaA <- length(interaction_case)
nfeaB <- length(categorical_features)
for(i in 1:nfeaA){
  for(j in 1:nfeaB){
    plot_1.data[n,2] <- colnames(cate.table)[j]
    plot_1.data[n,1] <- colnames(data_interaction)[i]
    Table.test <- table(cate.table[,j], data_interaction[,i])
    pValue <- chisq.test(Table.test)$p.value
    
    if(pValue < 0.05){
      plot_1.data[n,3] <- -log10(pValue) #for dot size, but if pValue > 0.05 , filling with an instinct number
      plot_1.data[n,4] <- 'A_red'
    }
    else{
      plot_1.data[n,3] <- 0.5 #for dot size, but if pValue > 0.05 , filling with an instinct number
      plot_1.data[n,4] <- 'B_grey'
    }
    plot_1.data[n,5] <- "myTag"
    n <- n+1
  }
}

##### panel 2: gene vs clinical categrical
#step2: conduct chisq.test for every-pair features and summarize the pvalue results into a table
plot_2.data <- data.frame(rep(NA,length(categorical_features)*length(sel_gene)),rep(NA,length(categorical_features)*length(sel_gene)),rep(NA,length(categorical_features)*length(sel_gene)),rep(NA,length(categorical_features)*length(sel_gene)),rep(NA,length(categorical_features)*length(sel_gene)))
n<-1
nfeaA <- length(sel_gene)
nfeaB <- length(categorical_features)
for(i in 1:nfeaA){
  for(j in 1:nfeaB){
    plot_2.data[n,2] <- colnames(cate.table)[j]
    plot_2.data[n,1] <- colnames(data_gene)[i]
    Table.test <- table(cate.table[,j], data_gene[,i])
    pValue <- chisq.test(Table.test)$p.value
    
    if(pValue < 0.05){
      plot_2.data[n,3] <- -log10(pValue) #for dot size, but if pValue > 0.05 , filling with an instinct number
      plot_2.data[n,4] <- 'A_red'
      if(j == 4){
        outfile <- paste0("~/Desktop/FINAL_FYT/", colnames(data_gene)[i])
        write.table(Table.test,file = outfile, sep = '\t')
      }
    }
    else{
      plot_2.data[n,3] <- 0.5 #for dot size, but if pValue > 0.05 , filling with an instinct number
      plot_2.data[n,4] <- 'B_grey'
    }
    plot_2.data[n,5] <- "myTag"
    n <- n+1
  }
}

data_interaction[data_interaction != "N"] <- 1
data_interaction[data_interaction == "N"] <- 0

##### panel 3: gene-interaction vs continuous value
#step1: define continuous features and generate dataframe
continuous_features <- c("Age","Tumor_size","Pack_years","Time_to_recurrence_or_death (months)")
cont.table <- raw_data_clincal[,which(colnames(raw_data_clincal) == continuous_features[1])]
for(i in 2:length(continuous_features)){
  tmp <- raw_data_clincal[,which(colnames(raw_data_clincal) == continuous_features[i])]
  cont.table <- cbind(cont.table,tmp)
} 
cont.table[which(cont.table$Pack_years == "NA"),] <- 0
colnames(cont.table) <- continuous_features

data_interaction[data_interaction == "1"] <- "T"
data_interaction[data_interaction == "0"] <- "F"
#step2: conduct t.test for every-pair features and summarize the pvalue results into a table
plot_3.data <- data.frame(rep(NA,length(continuous_features)*length(interaction_case)),rep(NA,length(continuous_features)*length(interaction_case)),rep(NA,length(continuous_features)*length(interaction_case)),rep(NA,length(continuous_features)*length(interaction_case)))
n<-1
nfeaA <- length(interaction_case)
nfeaB <- length(continuous_features)
for(i in 1:nfeaA){
  for(j in 1:nfeaB){
    plot_3.data[n,2] <- colnames(cont.table)[j]
    plot_3.data[n,1] <- colnames(data_interaction)[i]
    
    Table.test <- cbind(data_interaction[,i],as.numeric(cont.table[,j]))
    colnames(Table.test) <- c('categor','numeric')
    if(length(unique(Table.test[,1]))==1){
      plot_3.data[n,4] <- 'B_grey'
    }
    else{
      pValue <- t.test(as.numeric(Table.test[,2])~Table.test[,1],data=Table.test, var.equal = TRUE)$p.value
      if(pValue < 0.05){
        plot_3.data[n,3] <- -log10(pValue) #for dot size, but if pValue > 0.05 , filling with an instinct number
        plot_3.data[n,4] <- 'A_red'
      }
      else{
        plot_3.data[n,3] <- 0.5 #for dot size, but if pValue > 0.05 , filling with an instinct number
        plot_3.data[n,4] <- 'B_grey'
      }
    }
    n <- n+1
  }
}

##### panel 4: gene vs continuous value
#step1: define continuous features and generate dataframe

#step2: conduct t.test for every-pair features and summarize the pvalue results into a table
plot_4.data <- data.frame(rep(NA,length(continuous_features)*length(sel_gene)),rep(NA,length(continuous_features)*length(sel_gene)),rep(NA,length(continuous_features)*length(sel_gene)),rep(NA,length(continuous_features)*length(sel_gene)))
n<-1
nfeaA <- length(sel_gene)
nfeaB <- length(continuous_features)
for(i in 1:nfeaA){
  for(j in 1:nfeaB){
    plot_4.data[n,2] <- colnames(cont.table)[j]
    plot_4.data[n,1] <- colnames(data_gene)[i]
    
    Table.test <- cbind(data_gene[,i],as.numeric(cont.table[,j]))
    colnames(Table.test) <- c('categor','numeric')
    if(length(unique(Table.test[,1]))==1){
      plot_4.data[n,4] <- 'B_grey'
    }
    else{
      pValue <- t.test(as.numeric(Table.test[,2])~Table.test[,1],data=Table.test, var.equal = TRUE)$p.value
      if(pValue < 0.05){
        plot_4.data[n,3] <- -log10(pValue) #for dot size, but if pValue > 0.05 , filling with an instinct number
        plot_4.data[n,4] <- 'A_red'
      }
      else{
        plot_4.data[n,3] <- 0.5 #for dot size, but if pValue > 0.05 , filling with an instinct number
        plot_4.data[n,4] <- 'B_grey'
      }
    }
    n <- n+1
  }
}


#step3: draw the figure
#### panel 1
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}
plot_1.data<-plot_1.data[which(!(is.na(plot_1.data[,1]))),]
colnames(plot_1.data) <- c('FeatureA','FeatureB','dotSize','dotColor','myTag')
orderID<-c(1:nrow(plot_1.data))
coMu.plot<-ggplot()+theme_classic()
coMu.plot<-coMu.plot+geom_point(data = plot_1.data,aes(x=reorder(FeatureB,orderID),y=reorder(FeatureA,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot<-coMu.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.direction="horizontal",
                           legend.text=element_text(size=16,vjust=0.5,hjust=0,face='plain'),axis.text.y=element_text(size=16,vjust=0.5,hjust=1,face='bold',color='black'),
                           axis.text.x=element_text(size=16,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))

coMu.plot<-coMu.plot+scale_colour_manual(name=NULL,values=c(A_red='red',B_grey='grey'),
                                         labels=c(A_red='Significant',B_grey='Not significant'),
                                         guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
coMu.plot<-coMu.plot+scale_size(name=NULL,  breaks = c(1, 2, 10, 20, 30, 40),labels = expression(10^-1, 10^-2,10^-10, 10^-20, 10^-30,10^-40) ,guide = guide_legend(nrow=6))

coMu.plot<-coMu.plot+scale_x_discrete(position = "top")

plotAll<-rbind(ggplotGrob(coMu.plot),size="first")
plotAll

#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-1_stat_test.pdf", plot=plotAll,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)

#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-1_stat_test_CCF.pdf", plot=plotAll,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)

#### panel 2
plot_2.data<-plot_2.data[which(!(is.na(plot_2.data[,1]))),]
colnames(plot_2.data) <- c('FeatureA','FeatureB','dotSize','dotColor','myTag')
orderID<-c(1:nrow(plot_2.data))
coMu.plot2<-ggplot()+theme_classic()
coMu.plot2<-coMu.plot2+geom_point(data = plot_2.data,aes(x=reorder(FeatureB,orderID),y=reorder(FeatureA,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMu.plot2<-coMu.plot2+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.direction="horizontal",
                           legend.text=element_text(size=16,vjust=0.5,hjust=0,face='plain'),axis.text.y=element_text(size=16,vjust=0.5,hjust=1,face='bold',color='black'),
                           axis.text.x=element_text(size=16,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3)))

coMu.plot2<-coMu.plot2+scale_colour_manual(name=NULL,values=c(A_red='red',B_grey='grey'),
                                         labels=c(A_red='Significant',B_grey='Not significant'),
                                         guide = guide_legend(override.aes=list(size=4),nrow=4),na.translate = F) 
coMu.plot2<-coMu.plot2+scale_size(name=NULL,  breaks = c(1, 2, 10, 20, 30, 40),labels = expression(10^-1, 10^-2,10^-10, 10^-20, 10^-30,10^-40) ,guide = guide_legend(nrow=6))

coMu.plot2<-coMu.plot2+scale_x_discrete(position = "top")

figure_2<-rbind(ggplotGrob(coMu.plot2),size="first")


#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-2_stat_test.pdf", plot=figure_2,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)
#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-1_stat_test_CCF.pdf", plot=plotAll,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)

#### panel 3
colnames(plot_3.data) <- c('listA','listB','dotSize','dotColor')

orderID <<- c(1:nrow(plot_3.data))

coMut.plot3<-ggplot()+theme_classic()
coMut.plot3<-coMut.plot3+geom_point(data = plot_3.data,aes(x=reorder(listB,orderID),y=reorder(listA,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMut.plot3<-coMut.plot3+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.7,0.7,0.7,0.7),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                               text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='bottom',legend.direction="horizontal",
                               legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                               axis.text.x=element_text(size=14,angle=90,vjust=0.5,hjust=1,face='bold.italic',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),
                               axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
coMut.plot3<-coMut.plot3+scale_colour_manual(name=NULL,values=c(A_red='red',B_grey='grey'),
                                             labels=c(B_grey='Not significant',A_red='Significant'),
                                             guide = guide_legend(override.aes=list(size=4),nrow=3),na.translate = F)+guides(size=FALSE)

figure_3<-rbind(ggplotGrob(coMut.plot3),size="first")

#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-3_stat_test.pdf", plot=figure_3,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)
#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-1_stat_test_CCF.pdf", plot=plotAll,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)


#### panel 4
colnames(plot_4.data) <- c('listA','listB','dotSize','dotColor')

orderID <<- c(1:nrow(plot_4.data))

coMut.plot4<-ggplot()+theme_classic()
coMut.plot4<-coMut.plot4+geom_point(data = plot_4.data,aes(x=reorder(listB,orderID),y=reorder(listA,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMut.plot4<-coMut.plot4+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.7,0.7,0.7,0.7),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                               text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='bottom',legend.direction="horizontal",
                               legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                               axis.text.x=element_text(size=14,angle=90,vjust=0.5,hjust=1,face='bold.italic',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),
                               axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))

coMut.plot4<-coMut.plot4+scale_colour_manual(name=NULL,values=c(A_red='red',B_grey='grey'),
                                             labels=c(B_grey='Not significant',A_red='Significant'),
                                             guide = guide_legend(override.aes=list(size=4),nrow=3),na.translate = F)+guides(size=FALSE)

figure_4<-rbind(ggplotGrob(coMut.plot4),size="first")

#### FOR REGION SUM
#ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure8-4_stat_test.pdf", plot=figure_4,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)
#### FOR CCF
#ggsave(file="~/Desktop/FINAL_FYT/CCF/Figure8-1_stat_test_CCF.pdf", plot=plotAll,bg = 'white', width = 34, height = 25, units = 'cm', dpi = 600)


## FOR REGION_SUM
write.table(plot_1.data,"~/Desktop/FINAL_FYT/region_num/plot_1.modified.data.txt",row.names = F, sep = '\t')
write.table(plot_2.data,"~/Desktop/FINAL_FYT/region_num/plot_2.modified.data.txt",row.names = F, sep = '\t')
write.table(plot_3.data,"~/Desktop/FINAL_FYT/region_num/plot_3.modified.data.txt",row.names = F, sep = '\t')
write.table(plot_4.data,"~/Desktop/FINAL_FYT/region_num/plot_4.modified.data.txt",row.names = F, sep = '\t')

# #### FOR CCF
# write.table(plot_1.data,"~/Desktop/FINAL_FYT/CCF/plot_1.modified.data.txt",row.names = F, sep = '\t')
# write.table(plot_2.data,"~/Desktop/FINAL_FYT/CCF/plot_2.modified.data.txt",row.names = F, sep = '\t')
# write.table(plot_3.data,"~/Desktop/FINAL_FYT/CCF/plot_3.modified.data.txt",row.names = F, sep = '\t')
# write.table(plot_4.data,"~/Desktop/FINAL_FYT/CCF/plot_4.modified.data.txt",row.names = F, sep = '\t')

x.scale <<- NULL
orderID <<- NULL


##### binomial testing



