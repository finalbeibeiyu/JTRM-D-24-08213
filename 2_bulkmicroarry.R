###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("COMBAT")
BiocManager::install("ComBat")
BiocManager::install("genefilter",force = TRUE)
###去除批次效应
library(sva)
library(limma)
setwd("C:\\Users\\micro\\Desktop\\ischemic\\2_batch")
dta=read.table("merge.txt",sep="\t",header=T)
dta=as.matrix(dta)
rownames(dta)=dta[,1]
geneexp=dta[,2:ncol(dta)]
dta2=matrix(as.numeric(as.matrix(geneexp)),nrow=nrow(geneexp),dimnames=list(rownames(geneexp),colnames(geneexp)))
group1=c(rep(1,10),rep(2,12))##分组一定要修改 分三组 每组多少
group2=c(rep("t1",4),rep("t2",6),rep("t1",4),rep("t2",8))###对照组与实验组修改
md = model.matrix(~as.factor(group2))
cmb=ComBat(dta2,group1, md, par.prior=T)
cmb=rbind(id=colnames(cmb),cmb)
write.table(cmb,"input.txt",sep="\t",quote=F,col.names=F)



library(limma)
setwd("C:\\Users\\micro\\Desktop\\ischemic\\1_merge")
data0=read.table("GSE5296.txt",header = T,sep = "\t")
data1=read.table("GSE47681.txt",header = T,sep = "\t")
data4=merge(data0,data1,by="symbol")
#第三个数据集
data3=read.table("GSE16561.txt",header = T,sep = "\t")
data4=merge(data2,data3,by="symbol")
#相同基因取平均值
data4=as.matrix(data4)
rownames(data4)=data4[,1]
geneE=data4[,2:ncol(data4)]
newgene=matrix(as.numeric(as.matrix(geneE)),nrow=nrow(geneE),
               dimnames=list(rownames(geneE),colnames(geneE)))
newgene=avereps(newgene)
newgene=rbind(id=colnames(newgene),newgene)
write.table(newgene,"merge.txt",sep = "\t",quote = F,col.names = F)




#library
library(limma)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
setwd("D:\\analysis\\ischemic\\3_diff")
mydata<-read.table("input.txt",header = T,sep="\t")
mydata<-as.matrix(mydata)
rownames(mydata)=mydata[,1]
a=grep("///",mydata[,1])#一个探针对应多个基因
mydata=mydata[-a,]
geneExp=mydata[,2:ncol(mydata)]
matrixNa=list(rownames(geneExp),colnames(geneExp))
geneExp=matrix(as.numeric(as.matrix(geneExp)),nrow=nrow(geneExp),dimnames=matrixNa)
#WGCA分析以及免疫浸润分析
wgcna=rbind(id=colnames(geneExp),geneExp)
write.table(wgcna,"wgcna.txt",sep = "\t",quote = F,col.names = F)

group2=c(rep("t1",4),rep("t2",6),rep("t1",4),rep("t2",8))
design <- model.matrix(~ 0+factor(group2))
colnames(design) <- c("t1", "t2")
fit <- lmFit(geneExp, design)
con_matrix<-contrast.matrix <- makeContrasts(t2-t1, levels=design)
myfit <- contrasts.fit(fit, con_matrix)
myfit <- eBayes(myfit)
allgene<-topTable(myfit, coef=1, adjust="fdr",number = "all")
allgene=allgene[,-7]
allgene=cbind(id=rownames(allgene),allgene)
write.table(allgene,"allgene.txt",sep = "\t",quote = F,row.names = F)

diffgene<-allgene[abs(allgene$logFC)>=1 & allgene$P.Val<0.05,]
write.table(diffgene,"diffgene.txt",sep = "\t",quote = F,row.names = F)
upgene<-allgene[allgene$logFC>=1 & allgene$P.Value<0.05,]
write.table(upgene,"upgene.txt",sep = "\t",quote = F,row.names = F)
downgene<-allgene[allgene$logFC<= (-1) & allgene$P.Value<0.05,]
write.table(downgene,"downgene.txt",sep = "\t",quote = F,row.names = F)
#小提琴图
pdf("1_差异基因数量.pdf",6,6)
library(vioplot)
x1 <- c(441,1)
x2 <- c(956,1)
vioplot(x1,x2,
        names=c("Down genes","Up genes"),
        col=c("#7FABCC","#F9918C"))
title("Healthy_VS_SCI",ylab = "number of genes",xlab = "DEG")






###########################################多彩火山图1####################################
##3.1 设定阈值（FC和adjP)
threshold_logFC <- log2(2)
threshold_P.Value <- 0.05

##3.2 明确标签：标记自己挑选的差异基因（下面只是个例子）
# 选择logFC最高的前10个基因和最低的前10个基因
top_genes <- diffgene$id[order(diffgene$logFC, decreasing = TRUE)[1:10]]   # 最高的前10个基因
bottom_genes <- diffgene$id[order(diffgene$logFC, decreasing = FALSE)[1:10]] # 最低的前10个基因

# 合并这两个基因列表
gene_show <- unique(c(top_genes, bottom_genes))

# 根据gene_show创建label列，只有这些选定的基因会被标记
allgene$label <- ifelse(allgene$id %in% gene_show, allgene$id, '')
pdf("1_火山图.pdf",8,6)

# 重新构建volcano plot
volcano <- ggplot(data = allgene) +
  #点
  geom_point(aes(x = logFC,y = -log10(P.Value),color = logFC,size = -log10(P.Value))) +
  #线
  geom_vline(xintercept = c(-threshold_logFC,threshold_logFC), linetype = "dashed", color = "black", size = 0.2) +
  geom_hline(yintercept = threshold_P.Value, linetype = "dashed", color = "black", size = 0.2) +
  #标签
  geom_text_repel(aes(x = logFC,y = -log10(P.Value),label = label),size = 6,
                  box.padding = unit(0.2,"lines"),point.padding = unit(0.3,"lines"),
                  segment.color = "black",show.legend = FALSE)+
  #颜色
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  #scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),values = seq(0, 1, 0.2)) +
  #scale_color_continuous(low = '#08519c', high = '#de2d26')+
  #scale_color_continuous(values = c('#08519c','grey','#de2d26'))+
  #坐标轴
  labs(x = "log2(FoldChange)",y = "-log10( P value)") +
  #主题
  theme_bw()+
  theme(axis.text = element_text(size = 8,colour = 'black'),axis.title = element_text(size = 10), # 坐标轴标签和标题
        plot.title = element_text(hjust = 0.5,size = 6,face = "bold"), # 标题
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.key.size = unit(0.1, "cm"),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

# 显示图形
volcano

dev.off()


data0=read.table("input.txt",header = T,sep = "\t")
data1=read.table("2.txt",header = T,sep = "\t")

table(data0$transcript_id)
data2=merge(data0,data1,sort = F,by="id")
rownames(data2)=data2[,1]
write.table(data2,"差异基因.txt",sep = "\t",quote = F,row.names = F)







###Volcano
library(ggplot2)
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))

ggplot(allgene,aes(logFC, -log10(P.Value)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(P.Value), color= -log10(P.Value)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank())




head(allgene, 10)
Q <- ifelse(
  allgene$logFC<=(-0.4)&allgene$P.Value<0.05,'#0073C2CC',
  ifelse(allgene$logFC>=(0.4)&allgene$P.Value<0.05,'#EFC000CC',
         '#868686CC'))
#Q[is.na(Q)]<-'#868686CC'
names(Q)[Q=='#0073C2CC']<-'Up'
names(Q)[Q=='#868686CC']<-'Nodiff'
names(Q)[Q=='#EFC000CC']<-'Down'

#downvals<-c('CXCR2','MNDA')
#upvals<-c('PRG4','DNER')
a <- EnhancedVolcano(allgene,
                     x="logFC",
                     y="P.Value",
                     lab=allgene$id,
                     pCutoff=10e-1/20,#y轴阈值线(水平)
                     FCcutoff=0.4,#x轴阈值线（垂直）
                     pointSize=3,#散点大小
                     labSize=3.5,#标签大小
                     xlim=c(-2, 2),#限制X轴范围
                     ylim=c(0,10),#限制Y轴范围
                     colCustom=Q,#用group覆盖默认配色方案
                     colAlpha=0.6,#调整透明度
                     #selectLab=c(downvals,upvals),#使用selectLab参数选定所关注的标签
                     cutoffLineType='longdash',#阈值线类型，可选“blank”、“solid”、“dashed”、“dotted”、“dotdash”、“longdash”和“twodash”
                     cutoffLineCol='black',#阈值线颜色
                     cutoffLineWidth=0.38,#阈值线粗细
                     title="Volcano Plot Exp",#主标题
                     subtitle="Differential expression",#副标题
                     caption=bquote(~Log[2]~"fold change cutoff,0.4;p-value cutoff,0.05"),#注释说明
                     legendPosition='right',#图例位置
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE,
                     legendLabSize=12,#图例文字大小
                     legendIconSize=6)+coord_flip()#图例符号大小
#+1
a
#heatmap
library(ggsci)
library("scales")

#show_col(pal_jco("default", alpha = 0.8)(10))
heatmapexp<-geneExp[row.names(diffgene),]
Geo=c(rep("GSE5296",10),rep("GSE47681",12))
Type=c(rep("Normal",4),rep("SCI",6),rep("Normal",4),rep("SCI",8))
names(Type)=colnames(heatmapexp)
annotation=cbind(Geo,Type)
annotation=as.data.frame(annotation)
coul <- colorRampPalette(brewer.pal(9, "RdBu"))(50)
#coul <- scale_fill_jco()

#连续型
display.brewer.all(type = "seq")
#离散型
display.brewer.all(type = "div")
#极端型
display.brewer.all(type = "qual")

pdf("pheatmap.pdf",7,6)
#par(mar=c(5,4,10,15)+1)
pheatmap(heatmapexp, 
         annotation=annotation, 
         angle_col = "45",
         scale="row",
         #color =coul,
         color = colorRampPalette(c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"))(50),
         #gaps_row = c(3,10),#分行
         gaps_col = c(10),
         cluster_cols =F,fontsize_row=4,fontsize_col=6,show_rownames = F)
dev.off()

############################GSEA#########################
library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)
# load diff

diff <- allgene %>%
  arrange(desc(logFC))

genelist <- diff$logFC
names(genelist) <- diff$id

# check
head(genelist,3)
#    Aplnr    Foxf1     Bmp5
# 13.45176 13.35322 12.02845

# load tpm
expr <- read.table("input.txt",header = T,sep = "\t")
# check
head(expr,3)

#富集分析
ego <- gseGO(geneList     = genelist,
             OrgDb        = org.Mm.eg.db,
             ont          = "ALL",
             keyType      = 'SYMBOL',
             minGSSize    = 5,
             maxGSSize    = 500,
             pvalueCutoff = 1,
             verbose      = FALSE)
b <- ego@result
write.table(b,"BP.txt",sep = "\t",quote = F,row.names = F)
save(ego,file = 'GSEA.RData')
load("GSEA.RData")
input <- read.table("select.txt",header = T,sep="\t")
dotplotGsea(data = input,
            order.by = 'GeneRatio')
pdf("top10基因比例图.pdf",width = 10, height = 8)
dotplotGsea(data = ego,topn = 10,
            #order.by = 'NES',
            add.seg = T)
dev.off()

pdf("top10NES比例图.pdf",width = 10, height = 8)
dotplotGsea(data = ego,topn = 10,
            order.by = 'NES',
            add.seg = T)
dev.off()
#########多个GSEA

ego
terms <- c('GO:0006909',
           'GO:0050727',
           'GO:0042116',
           'GO:1905517',"GO:0097191","GO:0007409","GO:0006979","GO:0007269")
colorRampPalette(C("#CD9DC5", "white", "#1C84BC"))
pdf("top6单独NES比例图.pdf",width = 16, height = 6)

lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         subPlot=2,
         addPval = T,
         pvalX = 0.5, pvalY = 0.7,
         pCol = 'black',
         curveCol=c("#1C84BC", "grey","#CD9DC5"),
         rankCol=c("#1C84BC", "grey","#CD9DC5"),
         #rankSeq=600,
         pHjust = 0)
}) -> gseaList


# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 4,align = 'hv')
dev.off()

pdf("top6单独NES比例图新样式.pdf",width = 8, height = 8)

######多样式
lapply(terms, function(x){
  gseaNb(object = ego,
         geneSetID = x,
         newGsea = T,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         markTopgene=T,
         topGeneN=5,
         subPlot = 2)
}) -> gseaList1

# combine
cowplot::plot_grid(plotlist = gseaList1,ncol = 2,align = 'hv')
dev.off()
######一张图多个通路 最多三个三个
pdf("一张图多个通路.pdf",width = 6, height = 4)

gseaNb(object = ego,
       geneSetID = c('GO:0097190',
                     'GO:0006955','GO:0008219'),
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8))
dev.off()


