###############################普通流程 不去除双细胞 不进行细胞周期基因去除
getwd()
setwd("/mnt/data/home/tycloud/SCImac/0_data/")
options(stringsAsFactors = F)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)

load("sce1.RData")
table(sce$Group)
#读取数据
#1、批量读取单细胞的数据
dir_name=c('Sham','3dpi','7dpi')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0('/mnt/data/home/tycloud/SCImac/0_data/',dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 200)
}
#修改名称
names(datalist)=dir_name
#2、细胞质控####
# 批量计算线粒体和rRNA占比
##############过滤红细胞基因

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern = "^mt")# 计算线粒体占比 ^MT- human
  datalist[[i]] <- sce
  rm(sce)}


#质控前的
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')


#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$percent.mt)

hb_genes <- rownames(datalist$Sham)[grep("^Hb[^(p)]", rownames(datalist$Sham),ignore.case = T)]
hb_genes <- c("Hbb-bt","Hbb-bs","Hbb-y","Hbs1l","Hba-a1","Hbq1b","Hba-a2","Hbegf")
sce[["HB_percent"]] <- PercentageFeatureSet(sce, features=hb_genes) 

pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt","HB_percent"),
                         pt.size = 0.01, 
                         ncol = 4,
                         cols=c('#AB3282', '#23452F', '#BD956A', '#CCC9E6', '#585658', '#58A4C3'))
pearplot_befor1
ggsave(filename = 'QC_before1.pdf',plot = pearplot_befor1,he=7,wi=15)
rm(sce)


#过滤

datalist1 <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 300
            & nFeature_RNA < 10000  
            &percent.mt < 15)
})

#合并数据
sce <- merge(datalist1[[1]],y=datalist1[2:length(datalist1)])

sce[["HB_percent"]] <- PercentageFeatureSet(sce, features=hb_genes) 
sce <- subset(sce,HB_percent<0.1)



clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)


#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$Group <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","Group")
counts[,1] <- as.numeric(counts[,1])
library(ggsci)
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = Group, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + 
  scale_fill_jama() +
  theme(text=element_text(size=25),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 1),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = 'fit_cell_count.pdf',plot = fit_cell_count,width = 8,height = 5)
#质控后的小提琴图
pearplot_after1 <- VlnPlot(sce,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt","HB_percent"), 
                           pt.size = 0.01,
                           ncol = 4,
                           cols=c('#AB3282', '#23452F', '#BD956A', '#CCC9E6', '#585658', '#58A4C3'))
pearplot_after1
ggsave(filename = 'QC_after1.pdf',plot = pearplot_after1,he=6,wi=12)



#保存datalist文件

#3、数据预处理####
#合并数据

sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
### 可视化前20个高变基因
top20 <- head(VariableFeatures(sce), 20)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE, size=3.0)

feat_20 <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
feat_20
ggsave(filename = 'feat_20.pdf',plot = feat_20,he=10,wi=15)
#ScaleData
scale.genes <-  rownames(sce)
sce <- ScaleData(sce, features = scale.genes)
#样本的分组
meta1<-data.frame(matrix(nrow=length(sce@meta.data$orig.ident), ncol=2)) 
colnames(meta1)=c('orig.ident','Group')
meta1$orig.ident=sce@meta.data$orig.ident
table(sce$orig.ident)

### Group Tumor 为原发性肿瘤；Normal：正常
meta1[grep("Sham",meta1$orig.ident),]$Group="Sham"
meta1[grep("3dpi",meta1$orig.ident),]$Group="SCI"
meta1[grep("7dpi",meta1$orig.ident),]$Group="SCI"

sce <- AddMetaData(sce, meta1$Group,col.name = "Group")

table(sce@meta.data$Group)


sce$Group <- factor(sce$Group,levels = c("Sham",
                                         "SCI"))
sce$orig.ident <- factor(sce$orig.ident,levels = c("Sham",
                                         "3dpi",
                                         "7dpi"))

sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca") 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
ggsave(filename = 'sc_pca.pdf',plot = sc_pca,he=10,wi=15)
#可视化前2个PC的top20个基因
VizDimLoadings(sce, dims = 1:10, nfeatures = 20, reduction = "pca")
before_pca <- DimPlot(sce, reduction = "pca", group.by = "orig.ident",cols=my36colors)
ggsave('before_pca.pdf',before_pca,he=4,wi=5)
ElbowPlot(sce, ndims=50)


# 不一定要聚类，但modelHomotyp、

library(harmony)
sce <- RunHarmony(sce, group.by.vars = "orig.ident")
ElbowPlot(sce,reduction = 'harmony')
dim.usage<- 15
sce <- FindNeighbors(sce, dims = 1:dim.usage, reduction = "harmony") 
table(sce$Group)


library(clustree)

sce <- FindClusters(sce,resolution = c(0.2,0.3,0.4,0.5,0.6))

q <- clustree(sce@meta.data, prefix = "RNA_snn_res.")+ ggsci::scale_color_frontiers()
ggsave('树状图.pdf',q,he=10,wi=15)


sce <- RunUMAP(sce, dims = 1:15,reduction = 'harmony')
sce <- RunTSNE(sce, 
               dims=1:15, 
               reduction="harmony",
               perplexity=30,
               max_iter=1000)
# 去批次成功
after_pca <- DimPlot(sce, reduction = "pca", group.by = "orig.ident",cols=my36colors)
ggsave('after_pca.pdf',after_pca,he=4,wi=5)
table(sce$Group)
#颜色
Idents(sce) <- sce$RNA_snn_res.0.2


load("")
############### 继续做柱状图###################

#可视化
library(SCP)
library(BiocParallel)
library(stats)
library(ggsci)
color_cluster <- pal_locuszoom()(10)

table(sce$cell_type)
f <- CellDimPlot(sce, 
                 group.by = "RNA_snn_res.0.2", 
                 reduction = "umap", 
                 label =T,
                 #hex = TRUE, 
                 hex.bins = 30, 
                 palcolor=my36colors)

ggsave('RNA0.2UMAP图.pdf',f,he=6,wi=8)


#marker基因的筛选################
#需要修改#############
#寻找差异基因时的差异倍数
Logfc = 1
#差异基因时最小的表达比例
Minpct = 0.15
DefaultAssay(sce) <- "RNA"
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, 
                              min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)

write.table(sce.markers,'sce_marker_gene_0.2.txt',quote = F,row.names = F,sep='\t')

DimPlot(sce,
        group.by ="RNA_snn_res.0.2", 
        reduction="umap",
        #reduction="tsne",
        label = "T", 
        pt.size = 0.2,
        label.size = 5)


##############################注释###################  
library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scRNAtoolVis)

####################anno##################
hpca.se=get(load("ref_Mouse_all.RData"))
#hpca.se <- HumanPrimaryCellAtlasData()
table(sce$Group)

#获取基因的表达谱的count数据
testdata <- GetAssayData(sce, slot="data")
#获取聚类的亚群
clusters <- sce@meta.data$RNA_snn_res.0.2
pred.sce <- SingleR(test =  testdata, 
                    ref = hpca.se, 
                    labels = hpca.se$label.fine,
                    method = "clusters",
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

plotScoreHeatmap(pred.sce)

celltype = data.frame(ClusterID=rownames(pred.sce), celltype=pred.sce$labels
                      , stringsAsFactors = F)

celltype
write.table(celltype,'celltype.txt',quote = F,sep = '\t',row.names = F)
#手工注释
celltype=read.table("celltype.txt", header=T, sep="\t", check.names=F)
celltype = data.frame(celltype, stringsAsFactors = F)

sce$cell_type=sce$RNA_snn_res.0.2
sce$seurat_clusters <- sce$cell_type
for (i in 1:nrow(celltype)){
  sce$cell_type=gsub(paste0('^',celltype[i,1],'$'),celltype[i,2],sce$cell_type)
}
unique(sce$cell_type)
table(sce1$cell_type)

length(table(sce$cell_type))

table(sce1$Group)
Idents(sce) <- sce$cell_type
###去除基质细胞################

pdf("1_Umap图注释总.pdf", width=5, height=4)
CellDimPlot(sce, 
            group.by = "cell_type", 
            label =T,
            #hex = TRUE, 
            hex.bins = 30,
            palcolor=color_cluster)
dev.off()


library(parallel)
mc.cores <- detectCores()
options(mc.cores = 24)
register(SnowParam(workers = 10, progressbar = TRUE))
library(SCP)
library(BiocParallel)
library(stats)

sce <-  RunDEtest(srt = sce, group_by = "cell_type", fc.threshold = 1.2, min.pct = 0.15,only.pos = FALSE)
M <- sce@tools$DEtest_cell_type$AllMarkers_wilcox
write.table(M,'细胞sce_marker_gene.txt',quote = F,row.names = F,sep='\t')

CellStatPlot(sce,   
             stat.by = "cell_type", 
             group.by = "orig.ident",
             individual = F, 
             label = T,
             plot_type = "trend",
             palcolor=color_cluster)
ggsave("修改1_桑基图.pdf", height = 5, width = 4)


load('sce20240113.RData')

table(sce$cell_type)

myMarkers <- c("Ms4a7",	"Fabp5",	 "Gpnmb",	"Lyz2",	"Lgals3", #巨噬细胞
               "Cx3cr1",	"Tmem119",	"Siglech",	"Csf1r",	"P2ry12", #小胶质细胞
               "S100a9",	"Il1b",	"Ly6g",	"Ltf",#	中性粒
               "Nkg7",	"Gzma"	,"Ccl5","Klrd1","Cd3d", #MK-cell 
               "Ly6d",	"Cd79a",	"Cd79b",	"Iglc2" )	

library(ggsci)
CellDimPlot(sce,
            group.by = "cell_type", 
            split.by = "orig.ident", 
            reduction = "umap", 
            add_density = T,
            label = T,
            show_stat = FALSE,
            cells.highlight = TRUE, 
            theme_use = "theme_blank", 
            palcolor=color_cluster) 
CellDimPlot(sce, group.by = "cell_type", split.by = "orig.ident", 
reduction = "UMAP", palcolor=color_cluster)

ggsave("修改2_分组umap图.pdf", height = 4, width = 10)

table(sce$cell_type)
sce$cell_type <- factor(sce$cell_type,levels = c("Microglia","Macrophages","Neutrophils",
                                                 "NK cells","B cells"))

myMarkers2 <- c("Ms4a7",#巨噬细胞
                "Siglech",#小胶
                "Ltf",	#	中性粒 
                "Nkg7"	, #NK细胞 
                "Ly6d")	#T细胞
pdf("3_点热图.pdf",width = 10, height = 5)
jjDotPlot(object = sce,
          gene  =  myMarkers,
          anno = T,
          plot.margin = c(3,1,1,1),
          point.geom = T,
          id = 'celltype',
          tile.geom = T,
          x.text.angle=60,
          dot.col=c('White','#4182A4'),
          cluster.order=c("Macrophages","Microglia","Neutrophils",
                          "NK cells","B cells"))
dev.off()

pdf("3_umapmarker.pdf",width = 12, height = 3)
FeatureCornerAxes(object = sce,
                  reduction = 'umap',
                  show.legend = T,
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = myMarkers2)
dev.off()


####################Macrophages######################
load("注释Macrophages.RData")
Macrophages <- subset(sce,cell_type=="Macrophages")

Macrophages<- FindVariableFeatures(Macrophages, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(Macrophages)
Macrophages <- ScaleData(Macrophages, features = scale.genes)
Macrophages<- RunPCA(Macrophages, features = VariableFeatures(Macrophages))
pdf("1_去批次前PCA.pdf",width = 4, height = 3)
DimPlot(Macrophages, reduction = "pca", group.by = "orig.ident")
dev.off()
ElbowPlot(Macrophages)
## 重新harmony
library(harmony)
set.seed(1000)
Macrophages <- RunHarmony(Macrophages, group.by.vars = "orig.ident")
pdf("1_去批次后PCA.pdf",width = 4, height = 3)
DimPlot(Macrophages, reduction = "harmony", group.by = "orig.ident")
dev.off()
ElbowPlot(Macrophages,reduction = 'harmony')

library(clustree)
Macrophages <- FindNeighbors(Macrophages, reduction = 'harmony',dims = 1:15)
Macrophages <- FindClusters(Macrophages,resolution = c(0.1,0.2,0.3))
q <- clustree(Macrophages@meta.data, prefix = "RNA_snn_res.")+ ggsci::scale_color_frontiers()
ggsave('树状图2.pdf',q,he=10,wi=15)

Macrophages <- RunUMAP(Macrophages, reduction = 'harmony',dims = 1:15)

Idents(Macrophages)="RNA_snn_res.0.1"

# 瞄一眼，小亚群,0,2,4!!
DimPlot(Macrophages, label=TRUE,group.by = 'RNA_snn_res.0.1')

pdf("2_Umap图注释总Macrophages.pdf", width=5, height=4)
CellDimPlot(Macrophages, 
            group.by = "RNA_snn_res.0.1", 
            reduction = "umap", 
            label =T,
            #hex = TRUE, 
            hex.bins = 30)
dev.off()

library(parallel)
mc.cores <- detectCores()
options(mc.cores = 24)
register(SnowParam(workers = 10, progressbar = TRUE))
library(SCP)
library(BiocParallel)
library(stats)


#marker基因的筛选################
#需要修改#############
#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.15
DefaultAssay(Macrophages) <- "RNA"
Macrophages.markers <- FindAllMarkers(object = Macrophages,logfc.threshold = Logfc, 
                                 min.pct = Minpct,only.pos = T)
Macrophages.markers["pct.diff"]=Macrophages.markers$pct.1-Macrophages.markers$pct.2
length(unique(Macrophages.markers$gene))
head(Macrophages.markers)
write.table(Macrophages.markers,'巨噬细胞分群细胞marker.txt',quote = F,row.names = F,sep='\t')



hpca.se=get(load("ref_Mouse_all.RData"))
#hpca.se <- HumanPrimaryCellAtlasData()
#获取基因的表达谱的count数据
testdata <- GetAssayData(Macrophages, slot="data")
#获取聚类的亚群

clusters <- Macrophages@meta.data$RNA_snn_res.0.1
pred.Macrophages <- SingleR(test =  testdata, 
                       ref = hpca.se, 
                       labels = hpca.se$label.fine,
                       #method = "clusters",
                       clusters = clusters, 
                       assay.type.test = "logcounts", 
                       assay.type.ref = "logcounts")

plotScoreHeatmap(pred.Macrophages,
                 color = (c("#44BDED ", "white", "#EFD595")(50)),
)

celltype = data.frame(ClusterID=rownames(pred.Macrophages), celltype=pred.Macrophages$labels
                      , stringsAsFactors = F)

celltype
write.table(celltype,'巨噬细胞celltype.txt',quote = F,sep = '\t',row.names = F)
#手工注释
celltype=read.table("巨噬细胞celltype.txt", header=T, sep="\t", check.names=F)
celltype = data.frame(celltype, stringsAsFactors = F)
########修改
Macrophages$cell_type=Macrophages$RNA_snn_res.0.1
for (i in 1:nrow(celltype)){
  Macrophages$cell_type=gsub(paste0('^',celltype[i,1],'$'),celltype[i,2],Macrophages$cell_type)
}
unique(Macrophages$cell_type)
table(Macrophages1$cell_type)

length(table(Macrophages$cell_type))
Idents(Macrophages) <- Macrophages$cell_type
Macrophages1 <- Macrophages1[, Macrophages1$cell_type != "Dendritic"]

Macrophages <- Macrophages1

Macrophages$cell_type <- factor(Macrophages$cell_type, levels=c("T cells", "Macrophages","SCs",
                                                "FBs","Macrophagess","OLCs","SGCs"))
color_cluster <- c("#2E2A2B","#CE4E9C","#8B57A2")

pdf("修改2_巨噬细胞Umap图注释.pdf", width=8, height=4)
CellDimPlot(Macrophages, 
            group.by = "cell_type", 
            reduction = "umap", 
            label =F,
            #hex = TRUE, 
            hex.bins = 30,
            hex = F, 
            #hex.bins = 40, 
            palcolor=color_cluster)
dev.off()
library(parallel)
mc.cores <- detectCores()
options(mc.cores = 24)
register(SnowParam(workers = 10, progressbar = TRUE))
library(SCP)
library(BiocParallel)
library(stats)

Macrophages <-  RunDEtest(srt = Macrophages, group_by = "cell_type", fc.threshold = 1.2, min.pct = 0.15,only.pos = FALSE)

Macrophages <- RunEnrichment(
  srt = Macrophages, group_by = "cell_type", db = "GO_BP", species = "Mus_musculus",
  DE_threshold = "avg_log2FC > log2(1.2) & p_val < 0.05")
M <- Macrophages@tools$Enrichment_cell_type_wilcox$enrichment
M <- Macrophages@tools$DEtest_cell_type$AllMarkers_wilcox
write.table(M,'巨噬细胞marker基因.txt',quote = F,sep = '\t',row.names = F)


pdf("3_功能比较1.pdf", width=12, height=10.5)

EnrichmentPlot(
  srt = Macrophages, group_by = "cell_type", 
  palcolor = c("#4471CC","white","#DF6A67"),
  #group_use ="0",
  plot_type = "comparison",
  topTerm = 6)
dev.off()

pdf("3_火山图.pdf", width=10, height=6)

VolcanoPlot(srt = Macrophages, group_by = "cell_type")
dev.off()

table(Macrophages$cell_type)

pdf("3_Borderenrichmap.pdf", width=12, height=10)
EnrichmentPlot(
  srt = Macrophages, group_by = "cell_type", group_use = "Border-associated Macrophages",
  plot_type = "enrichmap")
dev.off()


pdf("3_Chemotaxisenrichmap.pdf", width=12, height=10)
EnrichmentPlot(
  srt = Macrophages, group_by = "cell_type", group_use = "Chemotaxis-inducing Macrophages",
  plot_type = "enrichmap")
dev.off()

pdf("3_Inflammatory Macrophagesenrichmap.pdf", width=12, height=8)
EnrichmentPlot(
  srt = Macrophages, group_by = "cell_type", group_use = "Inflammatory Macrophages",
  plot_type = "enrichmap")
dev.off()

EnrichmentPlot(
  srt = Macrophages, group_by = "cell_type", group_use = "Inflammatory Macrophages",
  words_excluded ="class",
  plot_type = "network")

Macrophages <- RunGSEA(
  srt = Macrophages, group_by = "cell_type", db = "GO_BP", species = "Mus_musculus",
  DE_threshold = " p_val < 0.05")

Top5 <- c("Mrc1","Maf","Ccl12","Trf","Selenop",#border
          "Selplg","Btg2","Ly86","Cd81","Jun",#Che
          "Fabp5","Arg1","Spp1","Lgals3","Il1b"
)
pdf("3_Macrophages表达热图.pdf",width = 8, height = 7)
AverageHeatmap(object = Macrophages,
               markerGene = Top5,
               group.by ="cell_type",
               clusterAnnoName = T,
               annoCol = F)
dev.off()
pdf("3_Macrophages分组点图.pdf",width = 10, height =8)
jjDotPlot(object = Macrophages,
          gene = Top5,
          id = 'cell_type',
          scale = T,
          #anno = T,
          ytree=F,
          plot.margin = c(3,1,1,1),
          point.shape = 21,
          dot.col=c('White','#DF6A67'),
          x.text.angle=60,
          #split.by="Sample",
          flip = F)
dev.off()

save(Macrophages,file = '注释Macrophages.RData')
##################巨噬细胞##################################

library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(GSEABase)

library(readxl)
SenMayo <- openxlsx::read.xlsx("./phagocytosis_regulators.xlsx")
unique(SenMayo$ID)
a <- SenMayo[SenMayo$ID %in% rownames(Macrophages),]
# AUCell_buildRankings
cells_rankings <- AUCell_buildRankings(Macrophages@assays$RNA@counts,featureType = "genes",plotStats=TRUE,splitByBlocks=TRUE)  
## 制作基因集
geneSets <- a
extraGeneSets <- c(GeneSet(sample(geneSets),setName="PCD_gene"))

geneSets1 <- GeneSetCollection(extraGeneSets)
names(geneSets1) 
cells_AUC <- AUCell_calcAUC(geneSets1, cells_rankings)
cells_AUC
pdf("小胶质细胞Macrophages后cells_assignment.pdf",width = 8, height = 6)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=T, assignCells=TRUE) 
dev.off()
thresholds <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE)
getThresholdSelected(cells_assignment)
cells_assignment$ABPC_gene$aucThr$selected <- cells_assignment$ABPC_gene$aucThr$thresholds
cells_assignment$ABPC_gene$aucThr$thresholds

##set gene set of interest here for plotMacrophagesng
geneSet <- "PCD_gene"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Macrophages$AUC <- aucs
df<- data.frame(Macrophages@meta.data, Macrophages@reductions[["umap"]]@cell.embeddings)
####################提取特定细胞#######################
meta1<-data.frame(matrix(nrow=length(Macrophages@meta.data$AUC), ncol=2)) 
colnames(meta1)=c('AUC1','AUC2')
meta1$AUC1=Macrophages@meta.data$AUC
unique(meta1$AUC1)
meta1$AUC2 <- ifelse(meta1$AUC1>0.03455305,"SenMayo","No_SenMayo")
Macrophages <- AddMetaData(Macrophages, meta1$AUC2,col.name = "AUC2")
table(Macrophages$AUC2)

colnames(df)

class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

FeatureDimPlot(
  Macrophages, features = paste0("AUC"), 
  reduction = "umap",
  palette="BluGrn")
show_palettes()
ggsave("3_AUC评分UMAP图.pdf", height = 4, width = 5)
comparisons = list(c("Monocyte","B cells"),
                   c( "Monocyte","Dendritic cells"),
                   c( "Monocyte","Erythroid cells"),
                   c( "Monocyte","Granulocytic cells"),
                   c( "Monocyte","Plasma cells"),
                   c( "Monocyte","T cells")), 

FeatureStatPlot(Macrophages, 
                stat.by = "AUC", 
                group.by = "cell_type", 
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = T,
                bg.by = "cell_type")

ggsave("3_AUC评分细胞比较.pdf", height = 5, width = 8)


FeatureStatPlot(Macrophages, 
                stat.by = "AUC", 
                group.by = "Group", 
                #split.by = "cell_type",
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = TRUE,
                #comparisons = list(c("ABPC","FBMSC"),
                                 #  c( "ABPC","Old")), 
                bg.by = "Group")
ggsave("3_AUC评分分组比较.pdf", height = 4, width = 5)


Macrophages$G2M.Score
FeatureStatPlot(Macrophages, 
                stat.by = "G2M.Score", 
                group.by = "Group", 
                #split.by = "cell_type",
                add_trend = TRUE,
                palcolor=color_cluster, 
                add_box = TRUE,
                comparisons = list(c("ABPC","FBMSC"),
                                   c( "ABPC","Old")), 
                bg.by = "Group")
ggsave("3_G2M评分分组比较.pdf", height = 2.3, width = 4)

table(Macrophages$Phase)
CellStatPlot(Macrophages, 
             #split.by = "cell_type",
             stat.by = c("AUC2"), 
             plot_type = "trend",
             group.by = "Group",
             label = T,
             palcolor=c("#EB7467","#79AF97"))
ggsave("3_AUC评分柱状图比较.pdf", height = 4, width = 4)

table(Macrophages$Group)
M <- Macrophages@tools$DEtest_cell_type$AllMarkers_wilcox
write.table(M,'celltypemarker.txt',quote = F,sep = '\t',row.names = F)


save(sce,file = 'sce1.RData')
save(Macrophages,file = '注释Macrophages.RData')
load('注释Macrophages.RData')
#################细胞轨迹#############################
table(Macrophages$cell_type)
FeatureDimPlot(Macrophages, features = paste0("Lineage", 1), reduction = "UMAP", theme_use = "theme_blank")
Macrophages <- RunPAGA(
  srt = Macrophages, group_by = "cell_type",
  linear_reduction = "PCA", nonlinear_reduction = "UMAP")
Macrophages <- RunPAGA(
  srt = Macrophages, group_by = "cell_type", linear_reduction = "PCA", nonlinear_reduction = "UMAP",
  embedded_with_PAGA = TRUE, infer_pseudotime = TRUE, root_group = "Inflammatory Macrophages"
)
PAGAPlot(srt = Macrophages, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE,palcolor=color_cluster)
ggsave("4_细胞轨迹PAGA图.pdf", height = 4, width = 4)

CellDimPlot(Macrophages, group.by = "cell_type", reduction = "UMAP", lineages = paste0("Lineage", 1), lineages_span = 0.1,palcolor=color_cluster)
ggsave("4_细胞轨迹图.pdf", height = 5, width = 8)
Macrophages <- RunDynamicFeatures(srt = Macrophages, lineages = c("Lineage1"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = Macrophages, lineages = c("Lineage1"),
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", cell_annotation = "cell_type",
   separate_annotation_palette = c("cosmic", "jama"),
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "blue",
  height = 5, width = 2
)
print(ht$plot)
ggsave("4_细胞轨迹热图.pdf", height = 10, width = 12)
load("注释Macrophages.RData")
DynamicPlot(
  srt = Macrophages, lineages = c("Lineage1"), group.by = "cell_type",
  features = c("Fcrls", "Olfml3", "Ly86", "Csf1r", "Selplg", "Cfh", "Cx3cr1", "Hexb"),
  compare_lineages = TRUE, compare_features = FALSE,  point_palette = "ucscgb"
)
ggsave("4_细胞轨迹热图.pdf", height = 10, width = 10)

FeatureDimPlot(
  srt = pancreas_sub, features = c("Fcrls", "Olfml3", "Ly86", "Csf1r", "Selplg", "Cfh", "Cx3cr1", "Hexb"),
  reduction = "UMAP", theme_use = "theme_scp")

theme_use？
  