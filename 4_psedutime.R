## 读取上节课的数据
load("注释Macrophages.RData")
load('注释Macrophages.RData')
#################细胞轨迹#############################


############################monocle2####################
devtools::load_all("monocle")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(monocle)
data <- as(as.matrix(Microglia_harmony@assays$RNA@counts), 'sparseMatrix')
?new
pd <- new('AnnotatedDataFrame', data = Microglia_harmony@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
########过滤低质量细胞
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM)
orderCells(HSMM, root_state=5)
table(Microglia_harmony$RNA_snn_res.0.1)

load(file = "HSMM.RData")
#State轨迹分布图

library(ggsci)
plot1 <- plot_cell_trajectory(HSMM, color_by = "State")+ scale_color_locuszoom(10)
plot1
pdf("6拟时序分析细胞排序.pdf", width=6, height=4)
plot_cell_trajectory(HSMM, color_by = "cell_type")+scale_color_manual(values = c("#2E2A2B", "#CE4E9C", "#8B57A2"))

dev.off()
##Pseudotime轨迹图
pdf("6拟时序分析时间排序.pdf", width=6, height=4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",root_state = 1)+viridis::scale_color_viridis(option="D")

dev.off()

HSMM=orderCells(HSMM,root_state = 3) 

table(Microglia_harmony$cell_type)
pdf("6拟时序分析分组.pdf", width=6, height=4)
plot_cell_trajectory(HSMM, color_by = "Group")+  scale_color_locuszoom()
dev.off()

table(HSMM@phenoData)
pdf("6拟时序分析年龄分组.pdf", width=6, height=4)
plot_cell_trajectory(HSMM, color_by = "Age")+ scale_color_simpsons()
dev.off()

#######添加monocle结果到cds对象中
Microglia_harmony@meta.data$Pseudotime <- HSMM@phenoData@data$Pseudotime#添加monocle结果到cds对象中
Microglia_harmony@meta.data$cytotrace<- as.numeric(Microglia_harmony@meta.data$cytotrace)
library(dplyr)
df<- data.frame(Microglia_harmony@meta.data, Microglia_harmony@reductions[["umap"]]@cell.embeddings)
class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p1 <- ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = Pseudotime)) +
  scale_color_viridis(option="G") +
  ggrepel::geom_text_repel(aes(label = cell_type),
                           data = class_avg,
                           size = 6,
                           segment.color = NA) + 
  #scale_colour_manual(values = my9color)+
  theme(legend.position = "none") + theme_classic()
p1 <- p1 & theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                 legend.position = "none")
pdf("4_monocle.pdf",width = 4, height = 4)

p1
dev.off()
p1 <- plot_cell_trajectory(HSMM, x = 1, y = 2, color_by = "cell_type") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_manual(values = c("#2E2A2B", "#CE4E9C", "#8B57A2"))
pdf("6拟时序分析复杂分组.pdf", width=6, height=4)

plot_complex_cell_trajectory(HSMM, x = 1, y = 2,
                                   color_by = "cell_type") +
  theme(legend.title = element_blank()) + scale_color_manual(values = c("#2E2A2B", "#CE4E9C", "#8B57A2"))
dev.off()
p1 / p2

root_state <- get_pseudotime(HSMM, type = "pseudotime")

print(root_state)
save(HSMM,file = "HSMM.Rdata")
