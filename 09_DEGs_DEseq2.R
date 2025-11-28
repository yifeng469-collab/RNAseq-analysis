rm(list = ls());options(stringsAsFactors=FALSE)
wd = "/Users/billzhang/Desktop/"
setwd(wd)
count_table=read.table(file="GenesExp_geneCounts.txt",head=T)
f<-function(x) length(which(x==0))/length(x)
count_table_expressed <- count_table[apply(count_table[ ,-c(1,2)],1,f) < 0.5, ]
count_table_expressed=na.omit(count_table_expressed)
count_table_expressed=count_table
#count_table_expressed=count_table
#count_table_expressed <- subset(count_table,count_table[,3]>5&count_table[,33]>5)
nrow(count_table) - nrow(count_table_expressed)
#二、DESeq2 差异分析
library(DESeq2)
#1.构建一个DESeqData对象，准备相应文件
countData <- as.matrix(count_table_expressed[ , -c(1,2)])#删去第一列和第二列
row.names(countData) <- count_table_expressed[ , 1]
filenames=colnames(count_table_expressed)[-c(1,2)]
countData=countData[,-4]
#condition <- factor(substr(sample_names,1,length(filenames)))
condition=factor(c("HFD_B4","HFD_B4","HFD_B4","HFD","HFD","HFD"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
mcols(dds)$basepairs <- count_table_expressed$GeneLength  #添加基因长度信息
#进行差异分析
dds <- DESeq(dds)
#提取fpkm值，就是RPKM
rpkm_table_normalized<- fpkm(dds, robust = TRUE)  #robust = TRUE选项提取经过sizeFactors归一化的
#如果选FALSE，则只是单纯的rpkm归一化，推荐选TRUE
#将这个表达矩阵转化为数据框
rpkm_table_normalized <- data.frame(rpkm_table_normalized)
#write.table(rpkm_table_normalized,file="B4_rpkm_table_normalized.txt",col.names=T,row.names=T,quote=F)
############################################################# differential gene expression 
res1 <- results(dds, contrast=c("condition", "HFD_B4","HFD"),pAdjustMethod = "fdr")
results_frame <- as.data.frame(res1) 
names(results_frame) <- paste(names(results_frame), "HFD_B4_vs_HFD", sep="_")
DEG=results_frame
library(qvalue)
DEG[,6]=qvalue(as.numeric(DEG[,5]))$qvalues
#logFC_cutoff=with(DEG,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff=0
k1=(DEG$padj_HFD_B4_vs_HFD<0.05)&(DEG$log2FoldChange_HFD_B4_vs_HFD < logFC_cutoff)
k2=(DEG$padj_HFD_B4_vs_HFD<0.05)&(DEG$log2FoldChange_HFD_B4_vs_HFD > logFC_cutoff)
DEG$change=ifelse(k1,"DR",ifelse(k2,"UR","NS"))
table(DEG$change)
DESeq2_DEG=DEG
write.table(DESeq2_DEG,"HFD_B4_vs_HFD_DESeq2_DEG_2022-1-14.txt",col.names=T,row.names=T,quote=F)
DESeq2_DEG$change <- factor(DESeq2_DEG$change,levels=c("UR","DR","NS")) 
library(ggplot2)
png(paste("HFD_B4_vs_HFD_DiffGeneExpFC",".png",sep=""),height=1600,width=2600,res=400)
g <-ggplot() + geom_point(data=DESeq2_DEG,aes(x=log2FoldChange_HFD_B4_vs_HFD,y=-log10(padj_HFD_B4_vs_HFD),colour=change))+ 
#geom_vline(aes(xintercept=0)+#按照up_down分组来显示颜色
#  scale_y_continuous(limits=c(0,31),breaks=c(5,27)) +#设置纵坐标极限,以及刻度的位置
 # scale_x_continuous(limits=c(-2,2),breaks=-4:4*0.5,
 #labels=c("-2.00","-01.5", "-1", -0.5, 0, 0.5, 1, 1.5, "2.00000")) +  
  #这里用来设置很坐标刻度的标记，可以很自由地根据自己的需要进行改变
  scale_color_manual(breaks=c("UR","DR","NS"), ##选择图例中要显示的组以及显示顺序(从上到下或者从左到右)
 values=c("#0000FF","#FF0000","#bdbdbd"),  #设置各组的颜色定义，这个颜色顺序设置与因子顺序一致，如果不设置为因子，则按组名字母顺序
 labels=c("Up Regulated","Down Regulated","NoSig"),
 name="") + #设置图例的标题，这里冒号之间为空，就不显示标题
  theme(axis.text.x= element_text(size=10, family="myFont", color="black", face= "bold", vjust=0.5, hjust=0.5))+
  theme(axis.text.y= element_text(size=10, family="myFont", color="black",face= "bold", vjust=0.5, hjust=0.5))+#纵轴标题
  xlab("log2 Fold Change") +  #设置横轴标题
  ylab("-log10 Pvalue") + 
  theme_bw() + 
  guides(colour=FALSE)+
  xlim(-10,10)+
  theme(axis.title.y=element_text(face="bold",size=10, angle=90, vjust=0.5, hjust=0.5))+
  theme(axis.title.x=element_text(face="bold",size=10, vjust=0.5, hjust=0.5)) 
g
dev.off()


##c("#FF0000","#0000FF","#bdbdbd")





