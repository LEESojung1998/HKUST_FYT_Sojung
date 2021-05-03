#library(devtools); install_github("im3sanger/dndscv")
library("dndscv")

dnds.raw <- read.delim("~/Desktop/secondsem/dnds.raw.txt")
dnds.raw[dnds.raw$var == "#NAME?","var"] <- '-'
dnds.raw[dnds.raw$var=="-C","var"] <- '-'
colnames(dnds.raw) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations2<-unique(dnds.raw)
mutations2 <- as.data.frame(mutations2)

#data("dataset_simbreast", package="dndscv")
#mutations<-unique(mutations)
#head(mutations)

dndsout = dndscv(mutations2, max_muts_per_gene_per_sample = 1000)
sel_cv = dndsout$sel_cv 
print(head(sel_cv), digits = 3)
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

signif_genes <- signif_genes[signif_genes$gene_name != "CDKN2A.p16INK4a",]
signif_genes[signif_genes$gene_name == "CDKN2A.p14arf", 'gene_name'] <- 'CDKN2A'
colnames(signif_genes) <- c("geneID", "qvalue")
print(dndsout$nbreg$theta)

write.table(signif_genes,"~/Desktop/secondsem/final+result/dnds.result.txt",row.names = F, sep = '\t')
write.table(signif_genes,"~/Desktop/FINAL_FYT/dnds.result.txt",row.names = F, sep = '\t')


