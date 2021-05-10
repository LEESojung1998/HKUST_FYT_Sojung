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

#### wrap into table
write.table(mdwholedata, file = "wholedata_th05.txt")
write.table(mdwholedata, file = "~/Desktop/FINAL_FYT/region_num/wholedata_th05.txt")

