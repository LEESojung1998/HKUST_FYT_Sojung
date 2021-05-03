#setwd("/Users/sandylee/Desktop/secondsem/result/")
setwd("/Users/sandylee/Desktop/FINAL_FYT/CCF")

library(dplyr)
library(stringr)
library(grid)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(randomcoloR)
library("igraph")
library(data.table)
#myowndriver <- read.delim("~/Desktop/firstsem/finding/myowndriver.11281.txt")
myowndriver <- read.delim("~/Desktop/secondsem/dnds.result.txt")

#### Initialization of the tables
#node.table <- data.frame(0, 0, 0, 0)
#colnames(node.table) <- c("Gene", "P_CDF", "FC", "Occurrence")
node.table <- data.frame(0, 0)
colnames(node.table) <- c("Gene", "Occurrence")
edge.table <- data.frame(0, 0, 0, 0)
colnames(edge.table) <- c("geneA", "geneB", "weight", "label")

#### Node Matrix
for(x in myowndriver$geneID){
  node.table <- rbind(node.table, c(x, 0, 0, 0))
}
node.table <- node.table[-1,]


#### Edge Matrix
for(x in myowndriver$geneID){
  for(y in myowndriver$geneID){
    if(x != y){
      edge.table <- rbind(edge.table, c(x, y, 0, 0))
    }
  }
}
edge.table <- edge.table[-1,]

#### Initialized Edge Node
countfile <- 0

for(patient in 1:100){
  countfile <- countfile + 1
  infile <- paste0(paste0(paste0(paste0("p0", patient),"_muttable"),".txt"))
  infilepath <- paste0("~/Desktop/FINAL_FYT/CCF/mut_table/",infile)
  muttable <- read.csv(infilepath, sep="")
  muttable <- cbind(muttable, total = c(0))
  if(is.na(muttable["countfreq"]) == FALSE){
    ### Remove thing below 10
    for(j in 1:nrow(muttable)){
      sumtemp <- 0
      for(i in 3:ncol(muttable)-2){
        if(muttable[j,i] < 0){
          muttable[j,i] <- 0
        }
        else{
          print(sumtemp)
          sumtemp <- sumtemp + muttable[j,i]
        }
      }
      muttable[j,"total"] <- sumtemp
    }
    
    ### Remove row with all 0 + sort the file
    muttable <- muttable[apply(muttable[,-1], 1, function(x) !all(x==0)),]
    muttable.sort.temp <- muttable[order(muttable$countfreq, decreasing = TRUE),]
    muttable.sort.temp <- muttable.sort.temp[order(muttable.sort.temp$total, decreasing = TRUE),]
    muttable.sort <- muttable.sort.temp[order(muttable.sort.temp$countfreq, decreasing = TRUE),]
    
    temp <- row.names(muttable.sort)
    
    #### node.table occurrence
    for(x in 1:nrow(muttable.sort)){
      temp.occur <- as.numeric(node.table[which(node.table$Gene == temp[x]),"Occurrence"])
      node.table[which(node.table$Gene == temp[x]),"Occurrence"] <- temp.occur + 1
    }
    
    #### Determine time line of the mutation with difference threshold 0
    if(length(temp) > 1){
      for(x in 2:nrow(muttable.sort)-1){
        for(y in (x+1):nrow(muttable.sort)){
          if(((muttable.sort[x,"total"] - muttable.sort[y, "total"]) > 0.1) && muttable.sort[x,1] >= muttable.sort[y,1]){
            temp.value <- as.numeric(edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"weight"])
            edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"weight"] <- temp.value + 1
            if(edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] == 0){
              if(patient == 100){
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste0("CRUK0", patient)
              }
              else if(patient < 10){
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste0("CRUK000", patient)
              }
              else{
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste0("CRUK00", patient)
              }
            }
            else{
              if(patient == 100){
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste(edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"], paste0("CRUK0", patient), sep=";")
              }
              else if(patient < 10){
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste(edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"], paste0("CRUK000", patient), sep=";")
              }
              else{
                edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"] <- paste(edge.table[which(edge.table$geneA == temp[x] & edge.table$geneB == temp[y]),"label"], paste0("CRUK00", patient), sep=";")
              }
            }
          }
        }
      }
    }
  }
}


#### remove 0 in the label part
edge.table <- edge.table[-c(which(edge.table["weight"] == 0)),]


### wrap into table or pdf
write.table(edge.table,"TEDGedge.original.CCF.txt",row.names = F, sep = '\t')
write.table(node.table,"TEDGnode.CCF.txt",row.names = F, sep = '\t')

#edge.table.sorted <- edge.table[order(edge.table$weight, decreasing = TRUE),]

#### Draw TEDG Graph
node <- data.frame(node.table)
edge <- data.frame(edge.table)

write.table(edge[which(as.numeric(edge$weight) > 2),],"TEDGedge.modified.CCF.txt",row.names = F, sep = '\t')
#write.table(edge.table.sorted[which(as.numeric(edge.table.sorted$weight)),],"TEDGedge.sorted.txt",row.names = F, sep = '\t')


net <- graph_from_data_frame(d=edge[which(as.numeric(edge$weight) > 1),], vertices=node, directed=T)

color.gradient <- function(x, colors=c("darkgreen","yellow","red"), colsteps=10000) {
  colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] 
}

colrs <- randomColor(count = length(V(net)$name))
V(net)$color <- colrs[1:length(V(net)$name)]
#V(net)$color <- color.gradient(V(net)$Name)

V(net)$size <- as.numeric(V(net)$Occurrence)
E(net)$width <- as.numeric(E(net)$weight)/2

E(net)$label <- NA

E(net)$arrow.size <- 0.4
E(net)$edge.color <- "gray90"

graph_attr(net, "layout") <- layout_with_lgl
#pdf(file="Figure2_TEDG_CPN.pdf", bg = 'white', width = 8, height = 6)
plot(net)
dev.off()
plot(net)





