# 推荐本地安装
#devtools::install_local('hdWGNCA.zip')
#BiocManager::install('harmony',update=F,ask=F)
library(hdWGCNA)
#加载单细胞分析包
library(Seurat)
#加载作图包
library(tidyverse)
library(cowplot)
library(patchwork)
#加载共表达网络分析包
#BiocManager::install('WGCNA',update=F,ask=F)
library(WGCNA)
gc()

#设置随机种子
set.seed(12345)


## 读取上节课的数据
load("注释Macrophages.RData")


#过滤出至少在5%的细胞中表达的基因
Macrophages <- SetupForWGCNA(
  Macrophages,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
)

table(Macrophages$cell_type)
table(Macrophages$orig.ident)
#构建metacells!!这一步非常重要，WGCNA对数据的稀疏性非常敏感，与普通转录组的WGCNA分析相比
# 单细胞的稀疏矩阵的解决方法是WGCNA流程的问题核心
table(Macrophages$dissociationMethod)
# construct metacells  in each group
Macrophages<- MetacellsByGroups(
  seurat_obj = Macrophages,k=20,
  max_shared = 10,
  # group.by一般关注的是组织类型和细胞类型!
  group.by = c("orig.ident",'cell_type'), # 也可以选择别的groupby
  ident.group = 'orig.ident' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
Macrophages <- NormalizeMetacells(Macrophages)
metacell_obj <- GetMetacellObject(Macrophages)


#转置表达矩阵
# 安全起见，另起一个对
seurat_obj  <- SetDatExpr(
  Macrophages
)


#选择softpower
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # 这边不用转置了，前面已转置
)


# plot the results:
plot_list <- PlotSoftPowers(seurat_obj,selected_power = 8)
# assemble with patchwork

pdf("1softpower.pdf", width=6, height=4)
wrap_plots(plot_list,ncol=2)
dev.off()


#查看POWER table
power_table <- GetPowerTable(seurat_obj)


#构建共表达网络
softpower=8

seurat_obj <- ConstructNetwork(seurat_obj,soft_power = softpower,
                               setDatExpr = F,overwrite_tom = T)

dev.off()


pdf("2WGCNA.pdf", width=6, height=4)

PlotDendrogram(seurat_obj, main='Macrophages hdWGCNA Dendrogram')
dev.off()


#(可选)获取TOM矩阵，可用来进行其他高级分析
TOM <- GetTOM(seurat_obj)


#计算模块协调特征
#记得scale一下 or else harmony throws an error:
seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = GetWGCNAGenes(seurat_obj),
  
)
# 计算ME，根据组织类型分组
# harmony必须biocManager安装，不可以用github安装！！！
library(harmony)

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use = "linear", #  choices are "linear", "poisson", or "negbinom"
  assay = NULL, # 默认:DefaultAssay(seurat_obj)
  pc_dim = 1)
  #group.by.vars="Group1" # 根据样本去批次化 harmonize

seurat_obj <- ModuleConnectivity(seurat_obj)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# plot genes ranked by kME for each module
#可视化每个模块中，按照kME打分的基因
PlotKMEs(seurat_obj, ncol=5)

# 获取hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 5000)
write.table(hub_df,'hub基因.txt',quote = F,sep = '\t',row.names = F)

head(hub_df)
pdf("3单独模块.pdf", width=5, height=6)
PlotKMEs(seurat_obj, ncol=3)
dev.off()

#记得保存上面hdWGNCA关键分析过程！
saveRDS(seurat_obj, file='hdWGCNA_object.rds')




####------------一些可视化-----------------------
## 模块间的相关性
library(igraph)
library(qgraph)
# 载入保存的

seurat_obj=readRDS('hdWGCNA_object.rds')

# 画模块间相关性图
dev.off()
pdf("3相关性.pdf", width=6, height=6)
addcol <- colorRampPalette(c("#F9918C", "white", "#7FABCC"))

ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2,col = addcol(100),
                  type ="lower",method = "ellipse")

ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2,col = addcol(100),
                        type ="upper",method = "number",add = T)

dev.off()

## 1，2，6 microglia到底与哪些模块相关
# 由于识别到了每个模块的hub基因，可以去计算hub基因打分
# compute gene scoring for the top 25 hub genes by kME for each module
# (方法一)with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# (方法二)with UCell method #推荐这种方法
# 由于Ucell刚刚更新，所以4.1.x的同学请用本地安装,依赖包自行安装
#devtools::install_local("UCell-1.3.zip")
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# featureplot
# 瞄一眼
Idents(Macrophages)="cell_type"
DimPlot(Macrophages, label=TRUE,split.by = 'Group1') 

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE ,# order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf("4模块UMAP.pdf", width=6, height=4)
wrap_plots(plot_list)
dev.off()

### dotplot
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
pdf("5模块相关性.pdf", width=8, height=8)
p <- DotPlot(seurat_obj, features=mods)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='#F9918C', mid='grey95', low='#7FABCC')

# plot output
p

dev.off()

library(igraph)
ModuleNetworkPlot(seurat_obj)

# hubgene network
pdf("5network.pdf", width=6, height=6)

HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
pdf("5network2.pdf", width=8, height=8)
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE)
dev.off()


# convert sex to factor
seurat_obj$msex <- as.factor(seurat_obj$msex)

# convert age_death to numeric
seurat_obj$age_death <- as.numeric(seurat_obj$age_death)

# 需要与模块进行关联的特征，该特征从Seurat_obj@meta.data抽取
cur_traits <- c('braaksc', 'pmi', 'msex', 'age_death', 'doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'total_counts_mt')

# 将基因模块与特征进行关联
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,  # 需要与模块进行关联的特征
  group.by = 'cell_type'  # 该参数表示以cell_type进行分组，以cell_type为单位计算每个cell_type的模块-特征相关性矩阵
)

# 查看相关性矩阵
mt_cor <- GetModuleTraitCorrelation(seurat_obj)

names(mt_cor)
# 'cor''pval''fdr'

names(mt_cor$cor)
# 'all_cells''INH''EX''OPC''ODC''ASC''MG'

head(mt_cor$cor$all_cells[, 1:5])
# 模块-特征相关性热图
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',  # 可选pval、fdr作为显著性分数
  label_symbol = 'stars',  # 以*号作为显著性标记，numeric则显示数值
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)