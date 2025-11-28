rm(list = ls());options(stringsAsFactors=FALSE)
wd = "/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/Gene_Counts"
setwd(wd);system(paste("cd",wd))
readLines("./4365.counts.txt")[1:6]
##获取SampleID
filenames=read.table(file="/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/stringtie/mergelist.txt",head=F)[,1]
#filenames=filenames[1:129]
#sample_names = c("1", "2", "4760","4764")
sample_names=as.character(filenames)
for(sample in sample_names) {
  count_table_each <- read.delim(sprintf("./%s.counts.txt", sample), header=FALSE)
  names(count_table_each) <- c("geneID", sample)
  if(!exists("count_table")){ count_table <- count_table_each } else {
count_table <- merge(count_table, count_table_each, by="geneID", all=TRUE)
  }
}
#write.table(count_table,"Count_table_merge.txt",row.names=F,col.names=T,quote=F)
#2.将基因长度的信息加入到这个数据框中
gene_length <- read.delim("./4365.geneLength.txt", header=FALSE)
names(gene_length) <- c("geneID", "geneLength")
count_table <- merge(gene_length, count_table, by="geneID", all.y=TRUE)
#3.只保留达到一定阈值的基因作为有表达的基因，比如：9个样本中总的reads数>=10
#f<-function(x) sum(x==0)
#count_table_expressed <- count_table[apply(count_table[ ,-c(1,2)],1,f) < 26, ]

f<-function(x) length(which(x==0))/length(x)

count_table_expressed <- count_table[apply(count_table[ ,-c(1,2)],1,f) < 0.5, ]

#count_table_expressed=count_table
#count_table_expressed <- subset(count_table,count_table[,3]>5&count_table[,33]>5)
nrow(count_table) - nrow(count_table_expressed)
#二、DESeq2 差异分析
library(DESeq2)
#1.构建一个DESeqData对象，准备相应文件
countData <- as.matrix(count_table_expressed[ , -c(1,2)])#删去第一列和第二列
row.names(countData) <- count_table_expressed[ , 1]
condition <- factor(substr(sample_names,1,length(filenames)))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)

mcols(dds)$basepairs <- count_table_expressed$geneLength  #添加基因长度信息
#进行差异分析

dds <- DESeq(dds)
#进行差异分析

#dds <- DESeq(dds) #用 ?DESeq 查看各个参数的设置，这里使用的全部是默认参数

#提取fpkm值，就是RPKM
rpkm_table_normalized<- fpkm(dds, robust = TRUE)  #robust = TRUE选项提取经过sizeFactors归一化的
#如果选FALSE，则只是单纯的rpkm归一化，推荐选TRUE
#将这个表达矩阵转化为数据框
rpkm_table_normalized <- data.frame(rpkm_table_normalized)
#colnames(rpkm_table_normalized) <- sub("X","",colnames(rpkm_table_normalized))
head(rpkm_table_normalized)

write.table(rpkm_table_normalized,"/home/ZYFLAB/projects/BackFat_IMF_RNAseq/01_Gene_expression/ALL_Tissues/ALL_Genes/rpkm_norm_ALL_555.txt",quote=F,row.name=T,col.name=T,sep="\t")
