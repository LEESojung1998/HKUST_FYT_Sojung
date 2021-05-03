setwd("/Users/sandylee/Desktop/FINAL_FYT/copy_num")

library(ggplot2)
library(readxl)

library(grid)
library(gridExtra)
library(reshape2)
library("igraph")
library(stringr)

raw.data <- read.delim("~/Desktop/FINAL_FYT/copy_num/one_gene_per_one.txt", header=TRUE)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}

arm.data <- read.delim("~/Desktop/FINAL_FYT/copy_num/sort.cytoband.txt", header=FALSE)
patients_region <- read.delim("~/Desktop/FINAL_FYT/patients_region.txt")
for(u in 1:nrow(arm.data)){
  if(grepl(".", arm.data[u,4]) == TRUE){
    temo4 <- unlist(strsplit(arm.data[u,4], "[.]"))
    arm.data[u,4] <- temo4[1]
  }
}
arm.data <- cbind(arm.data, chr_arm = rep("N"))
for(i in 1:nrow(arm.data)){
  if(grepl("p", arm.data[i,4]) == TRUE && grepl("q", arm.data[i,4]) == TRUE){
    arm.data[i,6] <- paste0(arm.data[i,1], "pq")
  }else if(grepl("p", arm.data[i,4]) == TRUE){
    arm.data[i,6] <- paste0(arm.data[i,1], "p")
  }else{
    arm.data[i,6] <- paste0(arm.data[i,1], "q")
  }
}
arm.val <- unique(arm.data$chr_arm)
nec.data <- raw.data[,c(1,5,10,19,20)]

nec.data <- cbind(nec.data, determine = rep("N"))
#nec.data <- cbind(nec.data, region = rep("N"))


for(i in 1:nrow(nec.data)){
  #### patient id extraction
  #temo <- unlist(strsplit(nec.data[i,2],"-"))
  #nec.data[i,2] <- temo[1]
  #nec.data[i,7] <- temo[2]
  temo1 <- unlist(strsplit(nec.data[i,5],","))
  tem.final <- rep("N", length(temo1))
  for(j in 1:length(temo1)){
    if(grepl(".", temo1[j]) == TRUE){
      temo2 <- unlist(strsplit(temo1[j], "[.]"))
      tem.final[j] <- temo2[1]
    }else{
      tem.final[j] <- temo1[j]
    }
  }
  tem.final <- na.omit(tem.final)
  tem.final <- unique(tem.final)
  finalword <- ""
  for(k in 1:length(tem.final)){
    finalword <- paste0(finalword, tem.final[k], ";")
  }
  nec.data[i,5] <- finalword
  finalword <- ""
  #### determine cpn
  if(nec.data[i,3] == 2){
    nec.data[i,6] = "normal"
  }else if(nec.data[i,3] == 1){
    nec.data[i,6] = "loss"
  }
  else if(nec.data[i,3] == 0){
    nec.data[i,6] = "deletion"
  }
  else if(nec.data[i,3] == 3){
    nec.data[i,6] = "gain"
  }else{
    nec.data[i,6] = "amp"
  }
}
nec.data <- cbind(nec.data, chr_arm = rep("N"))
for(i in 1:nrow(nec.data)){
  if(grepl("p", nec.data[i,5]) == TRUE && grepl("q", nec.data[i,5]) == TRUE){
    nec.data[i,7] <- paste0(paste0(nec.data[i,1], "p;"),paste0(nec.data[i,1], "q"))
    
  }else if(grepl("p", nec.data[i,5]) == TRUE){
    nec.data[i,7] <- paste0(nec.data[i,1], "p")
  }else{
    nec.data[i,7] <- paste0(nec.data[i,1], "q")
  }
}
### only top three?
# topone <- nec.data[nec.data$V19=="LRP1B",]
# topone <- rbind(topone, c(nec.data[nec.data$V19=="EPHA5",]))
# topone <- rbind(topone, c(nec.data[nec.data$V19=="FAT2",]))

#### based on gene
candidates.new <- read.delim("~/Desktop/FINAL_FYT/copy_num/candidates.new.txt", header=FALSE)
gene.list <- candidates.new$V1

nec.data.sum <- cbind(nec.data, raw.data[,c(2,3)])

case <- sort(unique(nec.data.sum$V5))
arm.case <- sort(unique(arm.data$chr_arm))
copynum.Matrix.sum <- rep('N',length(case))
for(i in 2:length(arm.case)){
  copynum.Matrix.sum <- cbind(copynum.Matrix.sum,rep('N',length(case)))
}

rownames(copynum.Matrix.sum) <- c(case)
colnames(copynum.Matrix.sum) <- c(arm.case)


for(i in 1:nrow(nec.data.sum)){
  if(grepl(";", nec.data.sum[i,7])==TRUE){
    temo10 <- unlist(strsplit(nec.data.sum[i,7],";"))
    for(j in 1:length(temo10)){
      if(copynum.Matrix.sum[nec.data.sum[i,2],temo10[j]] != "N" && copynum.Matrix.sum[nec.data.sum[i,2],temo10[j]] != nec.data.sum[i,6]){
        copynum.Matrix.sum[nec.data.sum[i,2], temo10[j]] <- paste(copynum.Matrix.sum[nec.data.sum[i,2], temo10[j]], nec.data.sum[i,6],sep=";")
      }else{
        copynum.Matrix.sum[nec.data.sum[i,2],temo10[j]] <- nec.data.sum[i,6]
      }
    }
  }else{
    if(copynum.Matrix.sum[nec.data.sum[i,2],nec.data.sum[i,7]] != "N" && copynum.Matrix.sum[nec.data.sum[i,2],nec.data.sum[i,7]] != nec.data.sum[i,6] ){
      copynum.Matrix.sum[nec.data.sum[i,2],nec.data.sum[i,7]] <- paste(copynum.Matrix.sum[nec.data.sum[i,2],nec.data.sum[i,7]], nec.data.sum[i,6], sep=";")
    }else{
      copynum.Matrix.sum[nec.data.sum[i,2],nec.data.sum[i,7]] <- nec.data.sum[i,6]
    }
  }
}
copynum.Matrix.sum[copynum.Matrix.sum=="N"] <- NA
x.scale <<- as.character(case)
copynum.Matrix.sum <- copynum.Matrix.sum[,colSums(is.na(copynum.Matrix.sum))<nrow(copynum.Matrix.sum)]
plot.data.sum <- melt(copynum.Matrix.sum)

for(vale.check.num in 1:nrow(plot.data.sum)){
  if(grepl(";", plot.data.sum[vale.check.num,"value"]) == TRUE && grepl("amp", plot.data.sum[vale.check.num,"value"]) == TRUE && grepl("gain", plot.data.sum[vale.check.num,"value"]) == TRUE){
    plot.data.sum[vale.check.num,"value"] <- "amp"
  }
  else if(grepl(";", plot.data.sum[vale.check.num,"value"]) == TRUE && grepl("deletion", plot.data.sum[vale.check.num,"value"]) == TRUE && grepl("loss", plot.data.sum[vale.check.num,"value"]) == TRUE){
    plot.data.sum[vale.check.num,"value"] <- "deletion"
  }
  else if(grepl(";", plot.data.sum[vale.check.num,"value"]) == TRUE){
    temp.seg1 <- unlist(strsplit(plot.data.sum[vale.check.num,"value"], ";"))
    temp.seg1[temp.seg1 == "normal"] <- NA
    temp.seg1 <- na.omit(temp.seg1)
    if(length(temp.seg1) ==1){
      plot.data.sum[vale.check.num,"value"] <- temp.seg1[1]
    }
    else{
      if(length(names(which.max(table(temp.seg1)))) != 1){
        print("warning")
      }else{
        plot.data.sum[vale.check.num,"value"] <- names(which.max(table(temp.seg1)))
      }
    }
  }
}
plot.data.sum[is.na(plot.data.sum)] = "normal"
plot.name.final<-ggplot()+theme_classic()
plot.name.final<-plot.name.final+geom_tile(data = plot.data.sum,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
plot.name.final<-plot.name.final+scale_fill_manual(name=NULL,values=c(amp = gg_color_hue(4)[4],normal = 'white', loss = gg_color_hue(4)[3], deletion = gg_color_hue(4)[2], gain = gg_color_hue(4)[1]),labels=c(gain = "gain", loss = "loss", deletion = "deletion", amp = "amplification", normal = "normal"))
plot.name.final<-plot.name.final+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='right',
                           legend.direction="horizontal",legend.text=element_text(size=10,face='bold.italic'),axis.text.y=element_text(size=10,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                           axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
plot.name.final<-plot.name.final+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)

ggsave(file="~/Desktop/FINAL_FYT/copy_num/copynum.matrix/final_sum_graph.pdf", plot=plot.name.final,bg = 'white', width = 120, height = 40, units = 'cm', dpi = 600)


whole_table <- rep(0,length(gene.list))
for(i in 2:5){
  whole_table <- cbind(whole_table, rep(0, length(gene.list)))
}
########### divide by gene
count <- 0
for(gene in gene.list){
  count <- count+1
  nec.data <- cbind(nec.data, raw.data[,c(2,3)])
  nec.data.mod <- nec.data[nec.data$V19 == gene,]
  
  deltcount<- length(which(nec.data.mod$determine=="deletion"))
  gaincount<- length(which(nec.data.mod$determine=="gain"))
  ampcount<- length(which(nec.data.mod$determine=="amp"))
  losscount<- length(which(nec.data.mod$determine=="loss"))
  normcount<- length(which(nec.data.mod$determine=="normal"))
  
  if(deltcount == 0){
    whole_table[count, 4] <- 0
  }else{
    whole_table[count, 4] <- as.numeric((deltcount/length(nec.data.mod$determine))*100)
  }
  if(gaincount == 0){
    whole_table[count, 2] <- 0
  }else{
    whole_table[count, 2] <- as.numeric((gaincount/length(nec.data.mod$determine))*100)
  }
  if(ampcount == 0){
    whole_table[count, 1] <- 0
  }else{
    whole_table[count, 1] <- as.numeric((ampcount/length(nec.data.mod$determine))*100)
  }
  if(losscount == 0){
    whole_table[count, 3] <- 0
  }else{
    whole_table[count, 3] <- as.numeric((losscount/length(nec.data.mod$determine))*100)
  }
  if(normcount == 0){
    whole_table[count, 5] <- 0
  }else{
    whole_table[count, 5] <- as.numeric((normcount/length(nec.data.mod$determine))*100)
  }
  
  case <- sort(unique(nec.data.mod$V5))
  arm.case <- sort(unique(arm.data$chr_arm))
  copynum.Matrix <- rep('N',length(case))
  for(i in 2:length(arm.case)){
    copynum.Matrix <- cbind(copynum.Matrix,rep('N',length(case)))
  }
  
  rownames(copynum.Matrix) <- c(case)
  colnames(copynum.Matrix) <- c(arm.case)
  
  
  for(i in 1:nrow(nec.data.mod)){
    if(grepl(";", nec.data.mod[i,7])==TRUE){
      temo10 <- unlist(strsplit(nec.data.mod[i,7],";"))
      for(j in 1:length(temo10)){
        if(copynum.Matrix[nec.data.mod[i,2],temo10[j]] != "N" && copynum.Matrix[nec.data.mod[i,2],temo10[j]] != nec.data.mod[i,6]){
          copynum.Matrix[nec.data.mod[i,2], temo10[j]] <- paste(copynum.Matrix[nec.data.mod[i,2], temo10[j]], nec.data.mod[i,6],sep=";")
        }else{
          copynum.Matrix[nec.data.mod[i,2],temo10[j]] <- nec.data.mod[i,6]
        }
      }
    }else{
      if(copynum.Matrix[nec.data.mod[i,2],nec.data.mod[i,7]] != "N" && copynum.Matrix[nec.data.mod[i,2],nec.data.mod[i,7]] != nec.data.mod[i,6] ){
        copynum.Matrix[nec.data.mod[i,2],nec.data.mod[i,7]] <- paste(copynum.Matrix[nec.data.mod[i,2],nec.data.mod[i,7]], nec.data.mod[i,6], sep=";")
      }else{
        copynum.Matrix[nec.data.mod[i,2],nec.data.mod[i,7]] <- nec.data.mod[i,6]
      }
    }
  }
  copynum.Matrix[copynum.Matrix=="N"] <- NA
  plot.name <- paste0(gene,".plot")
  x.scale <<- as.character(case)
  if(length(unique(nec.data.mod$chr_arm)) == 1){
    copynum.Matrix <- as.data.frame(copynum.Matrix[,unique(nec.data.mod$chr_arm)])
    plot.data <- as.data.frame(rep("N", length(case)))
    plot.data <- cbind(plot.data,rep(unique(nec.data.mod$chr_arm), length(case)))
    plot.data <- cbind(plot.data,rep("N", length(case)))
    colnames(plot.data) <- c("Var1", "Var2", "value")
    plot.data[,1] <- c(case)
    plot.data[,3] <- c(copynum.Matrix)
    
  }else{
    copynum.Matrix <- copynum.Matrix[,colSums(is.na(copynum.Matrix))<nrow(copynum.Matrix)]
    plot.data <- melt(copynum.Matrix)
  }
  
  for(vale.check.num in 1:nrow(plot.data)){
    if(grepl(";", plot.data[vale.check.num,"value"]) == TRUE && grepl("amp", plot.data[vale.check.num,"value"]) == TRUE && grepl("gain", plot.data[vale.check.num,"value"]) == TRUE){
      plot.data[vale.check.num,"value"] <- "amp"
    }
    else if(grepl(";", plot.data[vale.check.num,"value"]) == TRUE && grepl("deletion", plot.data[vale.check.num,"value"]) == TRUE && grepl("loss", plot.data[vale.check.num,"value"]) == TRUE){
      plot.data[vale.check.num,"value"] <- "deletion"
    }
    else if(grepl(";", plot.data[vale.check.num,"value"]) == TRUE){
      temp.seg <- unlist(strsplit(plot.data[vale.check.num,"value"], ";"))
      temp.seg[temp.seg == "normal"] <- NA
      temp.seg <- na.omit(temp.seg)
      if(length(temp.seg) ==1){
        plot.data[vale.check.num,"value"] <- temp.seg[1]
      }
      else{
        if(length(names(which.max(table(temp.seg)))) != 1){
          print("warning")
        }else{
          plot.data[vale.check.num,"value"] <- names(which.max(table(temp.seg)))
        }
      }
    }
  }
  
  if(count == 23){
    plot.name<-ggplot()+theme_classic()
    plot.name<-plot.name+geom_tile(data = plot.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
    plot.name<-plot.name+scale_fill_manual(name=NULL,values=c(amp = gg_color_hue(4)[4],normal = 'white', loss = gg_color_hue(4)[3], deletion = gg_color_hue(4)[2], gain = gg_color_hue(4)[1]),labels=c(gain = "gain", loss = "loss", deletion = "deletion", amp = "amplification", normal = "normal"))
    plot.name<-plot.name+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                               text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='right',
                               legend.direction="horizontal",legend.text=element_text(size=10,face='bold.italic'),axis.text.y=element_text(size=20,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                               axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
    plot.name<-plot.name+scale_x_discrete(breaks=x.scale,position = "bottom")+xlab(NULL)+ylab(NULL)+ggtitle(gene)
  
    outfile2 <- paste0("~/Desktop/FINAL_FYT/copy_num/copynum.matrix/",paste0(gene,"_copynum.Matrix.pdf"))
    ggsave(file=outfile2, plot=plot.name,bg = 'white', width = 70, height = 10, units = 'cm', dpi = 600)
    
    plot1 <- plot.name
   }else{
    plot.name<-ggplot()+theme_classic()
    plot.name<-plot.name+geom_tile(data = plot.data,aes(x=Var1,y=Var2,fill=value),color='black',width=0.7,height=0.7,stat='identity')
    plot.name<-plot.name+scale_fill_manual(name=NULL,values=c(amp = gg_color_hue(4)[4],normal = 'white', loss = gg_color_hue(4)[3], deletion = gg_color_hue(4)[2], gain = gg_color_hue(4)[1]),labels=c(gain = "gain", loss = "loss", deletion = "deletion", amp = "amplification", normal = "normal"))
    plot.name<-plot.name+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.15,1,1,2),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                               text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.68,'cm'),legend.position='right',
                               legend.direction="horizontal",legend.text=element_blank(),axis.text.y=element_text(size=20,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                               axis.text.x=element_blank(),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))
    plot.name<-plot.name+scale_x_discrete(breaks=x.scale,position = "top")+xlab(NULL)+ylab(NULL)+ggtitle(gene)
    outfile2 <- paste0("~/Desktop/FINAL_FYT/copy_num/copynum.matrix/",paste0(gene,"_copynum.Matrix.pdf"))
    ggsave(file=outfile2, plot=plot.name,bg = 'white', width = 70, height = 5, units = 'cm', dpi = 600)
    if(count == 2){
      plot2 <- plot.name
    }else if(count == 3){
      plot3 <- plot.name
    }else if(count == 4){
      plot4 <- plot.name
    }else if(count == 5){
      plot5 <- plot.name
    }else if(count == 6){
      plot6 <- plot.name
    }else if(count == 7){
      plot7 <- plot.name
    }else if(count == 8){
      plot8 <- plot.name
    }else if(count == 9){
      plot9 <- plot.name
    }else if(count == 10){
      plot10 <- plot.name
    }else if(count == 11){
      plot11 <- plot.name
    }else if(count == 12){
      plot12 <- plot.name
    }else if(count == 13){
      plot13 <- plot.name
    }else if(count == 14){
      plot14 <- plot.name
    }else if(count == 15){
      plot15 <- plot.name
    }else if(count == 16){
      plot16 <- plot.name
    }else if(count == 17){
      plot17 <- plot.name
    }else if(count == 18){
      plot18 <- plot.name
    }else if(count == 19){
      plot19 <- plot.name
    }else if(count == 20){
      plot20 <- plot.name
    }else if(count == 21){
      plot21 <- plot.name
    }else if(count == 22){
      plot22 <- plot.name
    }else{
      plot23 <- plot.name
    }
  }
  # outfile <- paste0("~/Desktop/FINAL_FYT/copy_num/copynum.matrix/",paste0(gene,"_copynum.Matrix.txt"))
  # write.table(copynum.Matrix,file = outfile, sep = '\t')
}
rownames(whole_table) <- c(gene.list)
write.table(whole_table,file = "~/Desktop/FINAL_FYT/copy_num/wholedata.txt", sep = '\t')

figure_2 <- rbind(ggplotGrob(plot23),ggplotGrob(plot2),ggplotGrob(plot3),ggplotGrob(plot4),ggplotGrob(plot5),ggplotGrob(plot6),ggplotGrob(plot7),ggplotGrob(plot8),ggplotGrob(plot9),ggplotGrob(plot10),ggplotGrob(plot11),ggplotGrob(plot12),ggplotGrob(plot13),ggplotGrob(plot14),ggplotGrob(plot15),ggplotGrob(plot16),ggplotGrob(plot17),ggplotGrob(plot18),ggplotGrob(plot19),ggplotGrob(plot20),ggplotGrob(plot21),ggplotGrob(plot22),ggplotGrob(plot1),size="last")

ggsave(file="~/Desktop/FINAL_FYT/copy_num/copynum.matrix/final_graph.pdf", plot=figure_2,bg = 'white', width = 120, height = 80, units = 'cm', dpi = 600)

# plot.name <- paste0(gene.list[1],".plot")
# figure_2 <- plot.name
# 
# for(j in 2:length(gene.list)){
#   plot.name <- paste0(gene.list[j],".plot")
#   figure_2<-rbind(ggplotGrob(figure_2),ggplotGrob(F1c.plot),size="last")
#   panels <- figure_2$layout$t[grep("panel", figure_2$layout$name)]
# }
# 
# 
# 



