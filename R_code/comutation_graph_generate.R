library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library("igraph")

#### input variable
mutGene <- read.delim("~/Desktop/secondsem/gene.matrix.txt")

mutGene[mutGene != "N"] <- 1
mutGene[mutGene == "N"] <- 0

plot_2.data <- data.frame(rep(0,ncol(mutGene)^2),rep(0,ncol(mutGene)^2),rep(0,ncol(mutGene)^2),rep(0,ncol(mutGene)^2))

n<-1
for(i in 1:ncol(mutGene)){
  for(j in 1:ncol(mutGene)){
    plot_2.data[n,1] <- colnames(mutGene)[i]
    plot_2.data[n,2] <- colnames(mutGene)[j]
    Mut.FEtest <- cbind(c(0,0),c(0,0))
    
    if(i > j){ #Primary Tumor
      Mut.FEtest[1,1] <- length(which(mutGene[,i] == 1 & mutGene[,j] == 1))
      Mut.FEtest[1,2] <- length(which(mutGene[,i] == 1 & mutGene[,j] == 0))
      Mut.FEtest[2,1] <- length(which(mutGene[,i] == 0 & mutGene[,j] == 1))
      Mut.FEtest[2,2] <- length(which(mutGene[,i] == 0 & mutGene[,j] == 0))
      coMutation <- ((Mut.FEtest[1,1]+1)*(Mut.FEtest[2,2]+1)) / ((Mut.FEtest[1,2]+1)*(Mut.FEtest[2,1]+1))
      
      pValue <- fisher.test(Mut.FEtest,alternative ="two.sided")$p.value
      
      if(pValue < 0.1){
        plot_2.data[n,3] <- -log10(pValue) + 1 #for dot size, but if pValue > 0.1 , filling with an instinct number
        if(coMutation > 1){
          plot_2.data[n,4] <- 'D_red'
        }
        else{
          plot_2.data[n,4] <- NA
        }
      }
      else{
        plot_2.data[n,3] <- 1 #for dot size, but if pValue > 0.1 , filling with an instinct number
        plot_2.data[n,4] <- 'C_grey'
      }
    }
    else{ 
      plot_2.data[n,3] <- NA 
      plot_2.data[n,4] <- NA 
    }
    n <- n+1
  }
}

#plot_2.data <- na.omit(plot_2.data)


#### Figuring Square
colnames(plot_2.data) <- c('listA','listB','dotSize','dotColor')

orderID <<- c(1:nrow(plot_2.data))

coMut.plot<-ggplot()+theme_classic()
coMut.plot<-coMut.plot+geom_point(data = plot_2.data,aes(x=reorder(listA,orderID),y=reorder(listB,orderID),size=dotSize,color=as.factor(dotColor)))+ylab(NULL)+xlab(NULL)+ggtitle(NULL)
coMut.plot<-coMut.plot+theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin=unit(c(0.7,0.7,0.7,0.7),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                             text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='bottom',legend.direction="horizontal",
                             legend.text=element_text(size=16,face='bold.italic'),axis.text.y=element_text(size=14,vjust=0.5,hjust=1,face='bold.italic',color='black'),
                             axis.text.x=element_text(size=14,angle=90,vjust=0.5,hjust=1,face='bold.italic',color='black'),axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),
                             axis.title.y=element_text(size=24,hjust=0.5,vjust=2,face='bold',color='black'))

coMut.plot<-coMut.plot+scale_colour_manual(name=NULL,values=c(E_black='black',C_grey='grey',D_red='red'),
                                           labels=c(E_black='Co-occurrence in recurrent',C_grey='Not significant',D_red='Co-occurrence Mutation'),
                                           guide = guide_legend(override.aes=list(size=4),nrow=3),na.translate = F)+guides(size=FALSE)

figure_2<-rbind(ggplotGrob(coMut.plot),size="first")
grid.draw(figure_2)
ggsave(file="~/Desktop/secondsem/final+result/Figure6_CoMutation.pdf", plot=figure_2,bg = 'white', width = 20, height = 22, units = 'cm', dpi = 600)
ggsave(file="~/Desktop/FINAL_FYT/Figure2_CoMutation.pdf", plot=figure_2,bg = 'white', width = 20, height = 22, units = 'cm', dpi = 600)


x.scale <<- NULL
orderID <<- NULL
