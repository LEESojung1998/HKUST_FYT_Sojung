#setwd('~/')
#setwd("~/Users/sandylee/Desktop/secondsem")

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
  #infilepath <- paste0("~/Desktop/secondsem/result/muttable_CPN/",infile)
  infilepath <- paste0("~/Desktop/FINAL_FYT/CCF/mut_table/",infile)
  muttable <- read.csv(infilepath, sep="")
  
  #### fill up the final dataset
  countrow <- 0
  for(i in rownames(wholedata)){
    countrow <- countrow+1
    if(length(which(rownames(muttable)==i)) != 0){
      #### already normalized so no need to normalized again the dataset
      totalregion <- ncol(muttable)-1
      wholedata[countrow, x] <- as.numeric(muttable[which(rownames(muttable)==i), ncol(muttable)])
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

#### wrap into table
write.table(mdwholedata, file = "wholedata_ccf.txt")

### PCA for genes
#pca <- glmpca(mdwholedata,20)
#factors <- pca$factors
#plot1 <- plot(factors[,1],factors[,2],pch=19)
#plot1 <- ggbiplot(factors, ellipse=TRUE, xlab = "Dimension1", ylab = "Dimension2")+theme_classic()+ scale_color_discrete(name = '')
#plot1_1 <- plot1 + ggtitle("PCA of all genes")
#print(plot1_1)

#ggsave(file="Figure1_PCA_graphforpatients_th10.pdf", plot=plot1, bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)

### PCA for patients
pca1 <- glmpca(t(mdwholedata),nrow(t(mdwholedata)))
factors1 <- pca1$factors
factors1new <- cbind(factors1, colorname = 0)

count <- 0
for(gene in myowndriver$geneID){
  count <- count + 1
  if(gene == "KRAS"){
    factors1new[count, "colorname"] <- 1
  }
  else if(gene == "FAM135B"){
    factors1new[count, "colorname"] <- 2
  }
  else if(gene == "EGFR"){
    factors1new[count, "colorname"] <- 3
  }
  else if(gene == "TP53"){
    factors1new[count, "colorname"] <- 4
  }
  else if(gene == "PIK3CA"){
    factors1new[count, "colorname"] <- 5
  }
  else if(gene == "KEAP1"){
    factors1new[count, "colorname"] <- 6
  }
  else{
    factors1new[count, "colorname"] <- 7
  }
}

plot1 <- plot(factors1new[,1], factors1new[,2], bty="n",main="PCA of patients data", col=c("#00AFBB", "#E7B800", "#FC4E07", "#912742",'#66FFFF','#9999FF', '#00FF00',"#169c47")[factors1new$colorname],xlab = "PCA1", ylab = "PCA2", pch = 20)
plot1_1 <- plot1 + text(factors1new[,1], factors1new[,2], labels=rownames(mdwholedata), cex=0.3, pos=4,cex.lab=1.5, cex.axis=0.3)
print(plot1_1)


### clustering analysis
#res.km <- eclust(mdwholedata, "hclust", nboot = 2)
#plot(res.km)
#groups <- cutree(res.km, k=4) 
#rect.hclust(res.km, k=4, border="red")
#fit <- kmeans(mdwholedata, 4)
#library(fpc)
#plotcluster(mdwholedata, fit$cluster)



#plot(factors1[,1],factors1[,2],pch=20,xlab="PCA1", ylab="PCA2")
#plot2_1 <- plot2 + ggtitle("PCA of all patients")
#plot2 <- ggbiplot(factors1, ellipse=TRUE, xlab = "Dimension1", ylab = "Dimension2")+theme_bw()+ scale_color_discrete(name = '')
#plot2_1 <- plot2 + ggtitle("PCA of all patients")
#print(plot2_1)

#ggsave(file="Figure2_PCA_graphforgenes_th10.pdf", plot=plot2, bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)
