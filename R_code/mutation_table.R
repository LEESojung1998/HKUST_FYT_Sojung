#setwd('~/')
#install.packages("stringr")
#install.packages("rlist")
#install.packages ("purrr")

setwd("/Users/sandylee/Desktop/secondsem/result/muttable")
#setwd("/Users/sandylee/Desktop/FINAL_FYT/region_num/mut_table")

#### required libraries
library(stringr)
library(ggplot2)
library ("purrr")
library(MASS)

##### download the dataset
#myowndriver <- read.delim("~/Desktop/firstsem/finding/myowndriver.11281.txt")
myowndriver <- read.delim("~/Desktop/secondsem/dnds.result.txt")
patientfileraw <- read.delim("~/Desktop/secondsem/region.filtered.txt")
patientname <- sort(unique(patientfileraw$SampleID))
region.Matrix <- matrix(ncol = 2, nrow = 1)


countfile <- 0

for(patientele in patientname){
  countfile <- countfile + 1
  patientfile <- patientfileraw[which(patientfileraw$SampleID == patientele),2:3]
  ## count what genes exist among driver genes
  existgene <- intersect(patientfile$Hugo_Symbol, myowndriver$geneID)
  
  ## first row of the mutation table
  a<-patientfile[1, "RegionSum"]
  regionnum <- str_count(a, ";")+1
  region.temp <- c(patientele,as.numeric(regionnum))
  region.Matrix <- rbind(region.Matrix, c(region.temp))
  
  ### Lymph node metastatsis
  #tempcheck <- unlist(strsplit(unlist(strsplit(a, ";")),":"))
  #if((which(tempcheck == "LN1") != 0)){
    #regionnum <- regionnum - 1
  #}
  #if((which(tempcheck == "LN1") != 0)){
    #regionnum <- regionnum - 1
  #}
  
  #### when there is no commongene in the patients
  if(length(existgene) == 0 || patientele == "CRUK0099" || patientele == "CRUK0013"){
    muttablefinal <- matrix(nrow = length(myowndriver$geneID), ncol = 1) 
    rownames(muttablefinal) <- c(myowndriver$geneID)
    colnames(muttablefinal) <- c("countfreq")
    print(patientele)
  }
  else{
    #### make a final mutation table for each patients and set the column and row names
    muttable<- matrix(nrow = length(existgene), ncol = regionnum+10) 
    rownames(muttable) <- c(existgene)
    muttablefinal<- matrix(nrow = length(existgene), ncol = regionnum) 
    colnames(muttablefinal) <- c(paste0("R", 1:regionnum))
    rownames(muttablefinal) <- c(existgene)
    muttablefinal <- cbind(muttablefinal, countfreq = c("patient"))
    
    #### create a temporary mutation table
    count <- 0
    for(gene in existgene){
      count <- count+1
      b <- grep(gene, patientfile$Hugo_Symbol)
      counttemp <- 0
      for(i in b){
        counttemp <- counttemp+1
        temp <- patientfile[i, "RegionSum"]
        muttable[count, counttemp] <- temp
      }
    }
    #### fill out the region in the final table
    for(ele in 1:length(existgene)){
      temo <- unlist(strsplit(muttable[ele,], ";"))
      na.omit(temo)
      for(element in 1:regionnum){
        for(element2 in 1:length(grep(ele,patientfile$Hugo_Symbol))){
          na.omit(temo[element+regionnum*element2])
          temo[element] <- paste(temo[element],temo[element + regionnum*element2], sep=":")
        }
      }
      
      #### calculate the tumor region percentage
      temofinal <- character(regionnum)
      for(m in 1:regionnum){
        tempstr <- unlist(strsplit(temo[m], ":"))
        na.omit(tempstr)
        tempstr <- str_subset(tempstr, "/")
        tempstrlist <- unlist(strsplit(tempstr, "/"))
        index1 <- 0
        index2 <- 0
        for(k in 1:length(tempstrlist)){
          if(k %% 2 != 0){
            index1 <- index1 + as.numeric(tempstrlist[k])
          }
          else{
            index2 <- index2 + as.numeric(tempstrlist[k])
          }
        }
        #### calculate the final percentage 
        muttablefinal[ele,m] <- fractions(index1/index2) * 100
      }
    }
    
    #### calculate the frequency based on the percentage
    for(m in 1:length(existgene)){
      finalcount <- 0
      for(n in 1:regionnum){
        if(as.integer(muttablefinal[m,n]) >= 5){
          finalcount <- finalcount + 1
        }
      }
      muttablefinal[m,ncol(muttablefinal)] <- finalcount
    }
  }
  
  ### wrap into table or pdf
  outfile <- paste0(paste0(paste0(paste0("p0", countfile),"_muttable"),".txt"))
  write.table(muttablefinal, file = outfile)
}

#### number of region for each patients
region.Matrix <- na.omit(region.Matrix)
F1.plot<-ggplot()+theme_classic()
F1.plot<-F1.plot+geom_bar(data=as.data.frame(region.Matrix),aes(x=V1,y=V2),width=0.7,color="black",stat='identity',position=position_stack())
F1.plot<-F1.plot+theme(panel.background=element_rect(fill='transparent',color='#7E7E7E'),plot.margin=unit(c(0.4,2,0.4,2),'lines'),
                         plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=24,vjust=1.4,hjust=0.5,face='bold'),
                         legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.8,'cm'),legend.position=c(0.92,0.8),legend.text=element_text(size=16,face='bold.italic'),
                         axis.text.x=element_text(size=12,angle=90,vjust=0.5,hjust=0,face='bold',color='black'),axis.text.y=element_text(size=12,vjust=0.5,hjust=1,face='bold',color='black'),
                         axis.title.x=element_text(size=24,vjust=0,hjust=0.5,face='bold',color='black'),axis.title.y=element_text(size=24,hjust=0.9,vjust=2,face='bold',color='black'))
F1.plot<-F1.plot+ggtitle(NULL)+xlab(NULL)+ylab('Number of Region')
write.table(region.Matrix,file = "patients_region.txt",row.names = F, sep = '\t')

ggsave(file="~/Desktop/secondsem/final+result/Figure6_number_of_region.pdf", plot=F1.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)
ggsave(file="~/Desktop/FINAL_FYT/Figure3_number_of_region.pdf", plot=F1.plot,bg = 'white', width = 40, height = 33, units = 'cm', dpi = 600)


