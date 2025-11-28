######################################################

rm(list=ls());options(stringAsFactors=FALSE)
wd="/Users/billzhang/Desktop"
setwd(wd)
#install.packages("GUniFrac")
#library(GUniFrac) #用于计算Unifrac距离
library(ape) # 用于pcoa分析
library(ggplot2) #用于画图####GO KEGG
dat=read.csv(file="DEGs.csv",head=T)
UR=subset(dat,dat[,3]>0)[,1]

UR_Genes=toupper(UR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
GOGenes=toupper(UR_Genes)

gene.df <- bitr(GOGenes, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENSEMBL", "ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db)

ENTREZIDID=gene.df[,3]
ENTREZIDID=unique(ENTREZIDID)
# GO 分析：
ego <- enrichGO(gene          = ENTREZIDID,
                # universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",   ###GO:BP, GO:CC and GO:MF 三种富集类型
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
write.csv(summary(ego),file=paste("DEGs_GOBP_Pathway_","Genes",".csv",sep=""),row.names =F,col.names=T,quote=F)
 
#ekk <- enrichKEGG(gene= ENTREZIDID,organism  = 'hsa', qvalueCutoff = 0.05) #KEGG富集分析
 
#write.csv(summary(ekk),"KEG_enrichment.csv",row.names =F) 
png(paste("DEGs_GOBP_Pathway_","Genes",".png",sep=""),height=2200,width=2600,res=500)

a=dotplot(ego,showCategory=10,title="Enrichment GO Top10") #plot

print(a)

dev.off()

##KEGG
ekk=enrichKEGG(ENTREZIDID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
pAdjustMethod = "BH",qvalueCutoff = 0.05, use_internal_data = FALSE)
png(paste("DEGs_KEGG_Pathways_","Genes",".png",sep=""),height=2200,width=2600,res=400)
a=dotplot(ekk,showCategory=10,title="Enrichment KEGG Top10") #plot
print(a)
dev.off()
z = setReadable(ekk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(summary(z),"DEGs_KEGG_Pathways.csv",row.names =F) 


