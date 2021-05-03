#setwd("~/Users/sandylee/Desktop/secondsem/result/threshold")

library(ggplot2)

threshold.table <- matrix(nrow = 2, ncol = 10)
threshold.table[1,] <- c(1:10)
temo <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

for(x in 1:10){
  if(x != 10){
    temo[x]<- paste0("th0", x)
  }
  else{
    temo[x]<- paste0("th", x)
  }
}
threshold.table[1,] <- c(temo)

for(x in 1:10){
  infile <- paste0(paste0(paste0("wholedata_th0", x),".txt"))
  infilepath <- paste0("~/Desktop/secondsem/result/threshold/",infile)
  threshold <- read.csv(infilepath, sep="")
  threshold.table[2,x] <- ncol(threshold)
}

threshold.table <- t(threshold.table)
colnames(threshold.table) <- c("thresholdnum","percentage")

#### scatter plot
p <- ggplot(as.data.frame(threshold.table), aes(x = thresholdnum, y = percentage, group=1)) + 
   geom_line(color = "#2874A6") + geom_point(color = "#5DADE2", size = 2.5) + theme_light() + ggtitle("Change in Patients' coverage over threshold") + xlab("threshold") + ylab("Percentage(%)") + theme(plot.title = element_text(hjust=0.5, face="bold"))

ggsave(file="~/Desktop/secondsem/final+result/Figure3_threshold coverage.pdf", plot=p, bg = 'white', width = 17, height = 10, units = 'cm', dpi = 600)
ggsave(file="~/Desktop/FINAL_FYT/region_num/Figure5_threshold coverage.pdf", plot=p, bg = 'white', width = 17, height = 10, units = 'cm', dpi = 600)

