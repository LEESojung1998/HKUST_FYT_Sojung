#setwd('~/')
#install.packages("stringr")
#install.packages("rlist")
#install.packages ("purrr")

#setwd("/Users/sandylee/Desktop/secondsem/result/muttable_CPN")
setwd("/Users/sandylee/Desktop/FINAL_FYT/CCF/mut_table")

#### required libraries
library(stringr)
library(ggplot2)
library ("purrr")
library(MASS)

##### download the dataset
#myowndriver <- read.delim("~/Desktop/firstsem/finding/myowndriver.11281.txt")
myowndriver <- read.delim("~/Desktop/secondsem/dnds.result.txt")
#### cpn value dataset
patientfileraw_cpn <- read.delim("~/Desktop/secondsem/analysis.filtered.ws.txt")
patientfileraw_cpn <- patientfileraw_cpn[,c(2,8:10,18)]
patientfileraw_cpn <- patientfileraw_cpn[which(patientfileraw_cpn$exonic.func != "synonymousSNV" & patientfileraw_cpn$exonic.func != "NA" & patientfileraw_cpn$exonic.func != "unknown" & patientfileraw_cpn$exonic.func != "synonymous"),]
unique(patientfileraw_cpn$exonic.func)
patientfileraw <- patientfileraw_cpn[,c(1:2,5)]

patientname <- sort(unique(patientfileraw$SampleID))
region.Matrix <- matrix(ncol = 2, nrow = 1)


countfile <- 0

for(patientele in patientname){
  countfile <- countfile + 1
  patientfile <- patientfileraw[which(patientfileraw$SampleID == patientele),2:3]
  patientfile <- na.omit(patientfile)
  ## count what genes exist among driver genes
  existgene <- intersect(patientfile$Hugo_Symbol, myowndriver$geneID)
  
  #### just for checking
  print(patientele)
  
  ## first row of the mutation table
  a<-patientfile[1, "PhyloCCF"]
  regionnum <- str_count(a, ";")+1
  region.temp <- c(patientele,as.numeric(regionnum))
  region.Matrix <- rbind(region.Matrix, c(region.temp))
  #### when there is no commongene in the patients
  if(length(existgene) == 0 || patientele == "CRUK0099" || patientele == "CRUK0013"){
    muttablefinal <- matrix(nrow = length(myowndriver$geneID), ncol = 1) 
    rownames(muttablefinal) <- c(myowndriver$geneID)
    colnames(muttablefinal) <- c("countfreq")
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
        temp <- patientfile[i, "PhyloCCF"]
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
        countnum <- 0
        totalval <- 0
        for(k in 1:length(tempstr)){
          if(k %% 2 == 0 & tempstr[k] != "NA"){
            countnum <- countnum+1
            totalval <- totalval+as.numeric(tempstr[k])
          }
        }
        muttablefinal[ele,m] <- as.numeric(as.numeric(totalval)/as.numeric(countnum))
      }
    }
    
    #### calculate the frequency by taking average of all percentage in one region
    for(m in 1:length(existgene)){
      finalcount <- 0
      for(n in 1:regionnum){
        if(muttablefinal[m,n] >= 0.01){
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

