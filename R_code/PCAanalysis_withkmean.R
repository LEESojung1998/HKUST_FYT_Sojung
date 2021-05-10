#setwd('~/')
setwd("~/Users/sandylee/Desktop/secondsem")


#install.packages("factoextra")
#install_github("vqv/ggbiplot")

library(ggplot2)
library(dplyr)
library(devtools)
require(ggbiplot)
library(factoextra)
library(glmpca)

#### get required datasets
#myowndriver <- read.delim("~/Desktop/firstsem/finding/myowndriver.11281.txt")
myowndriver <- read.delim("~/Desktop/secondsem/dnds.result.txt")


#### initialize a dataset countfreq of 100 patient
wholedata <- matrix(nrow = nrow(myowndriver), ncol = 100) 
colnames(wholedata) <- c(paste0("p", 1:100))
rownames(wholedata) <- c(myowndriver$geneID)

#### loop through all the dataset
for(x in 1:100){
  infile <- paste0(paste0(paste0(paste0("p0", x),"_muttable"),".txt"))
  #infilepath <- paste0("~/Desktop/secondsem/result/muttable/",infile)
  infilepath <- paste0("~/Desktop/FINAL_FYT/CCF/mut_table/",infile)
  
  muttable <- read.csv(infilepath, sep="")
  
  #### fill up the final dataset
  countrow <- 0
  for(i in rownames(wholedata)){
    countrow <- countrow+1
    if(length(which(rownames(muttable)==i)) != 0){
      #### normalized the dataset
      totalregion <- ncol(muttable)-1
      wholedata[countrow, x] <- as.numeric(muttable[which(rownames(muttable)==i), ncol(muttable)]/totalregion)
    }
    ### fill with 0 for NA
    else{
      wholedata[countrow, x] <- 0
    }
  }
}

#### Add one total

#### draw a PCA data - preprocess the dataset
wholedata[is.na(wholedata)] <- 0
mdwholedata <- data.frame(na.omit(wholedata))
mdwholedata <- mdwholedata[,apply(mdwholedata, 2, var, na.rm=TRUE) != 0]

### PCA for patients

pca <- prcomp(mdwholedata,nrow(mdwholedata)) # principal components analysis using correlation matrix
pc <- prcomp(mdwholedata)
comp <- data.frame(pc$x[,1:4])

k <- kmeans(comp, 4, nstart=25, iter.max=1000)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 120)[1:n]
}
comp <- cbind(comp, value = as.character(c(k$cluster)))

plot1<-ggplot()+theme_classic()
plot1<-plot1+geom_point(data = comp,aes(x=comp[,1],y=comp[,2],fill = value, color = value,size = 5))+ylab("PCA2")+xlab("PCA1")+ggtitle(NULL)
plot1<-plot1+scale_fill_manual(name=NULL,values=c("1"=gg_color_hue(4)[4],"2"=gg_color_hue(4)[3],"3"=gg_color_hue(4)[2],"4"=gg_color_hue(4)[1]),labels=c("1"="cluster 1","2"="cluster 2","3"="cluster 3","4"="cluster 4"))
plot1<-plot1+theme(plot.margin=unit(c(1,1,1,1),'lines'),plot.title=element_text(size=24,vjust=0.5,hjust=0.5,face='bold.italic'),
                           text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.3,'cm'),legend.position='right',legend.direction="horizontal",
                           legend.text=element_text(size=16,vjust=0.5,hjust=0,face='plain'),axis.text.y=element_text(size=16,vjust=0.5,hjust=1,face='bold',color='black'),
                           axis.text.x=element_text(size=16,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.title.x=element_text(size=20,vjust=0,hjust=0.5,face='plain',color='black'),
                           axis.title.y=element_text(size=20,hjust=0.5,vjust=2,face='plain',color='black'),
                           strip.text = element_text(size=18,face='bold',vjust=0.5,hjust=0.5),strip.background = element_rect(colour="black", fill=gg_color_hue(3))) + geom_text(aes(comp[,1], comp[,2],label=rownames(mdwholedata),size=5),nudge_y = -0.05, nudge_x = -0.05)
ggsave(file="~/Desktop/FINAL_FYT/sample.pdf", plot=plot1,bg = 'white', width = 60, height = 50, units = 'cm', dpi = 600)



