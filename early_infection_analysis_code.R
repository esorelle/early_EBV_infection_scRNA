library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Signac)
library(cellranger)
library(monocle3)
library(magrittr)
library(data.table)
library(patchwork)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(GGally)
library(gplots)
library(plotly)
library(cowplot)
library(scales)


#### PREP ####

setwd(paste0(getwd()))

# QC filter control panel
min_cell = 3
mito_pct = 20
min_feat = 200
max_count = 65000

# get gene name lists for cell cycle scoring and regression (s and g2m phases)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


#### READ AND MERGE TIMECOURSE DATA ####

### TX1241 ###

## Day 0
# load and process data
TX1241_0.data <- Read10X(data.dir = './TX1241_0/outs/filtered_feature_bc_matrix/')
TX1241_0 <- CreateSeuratObject(counts = TX1241_0.data$'Gene Expression', project = "TX1241_0", assay = 'RNA')
TX1241_0[["percent.mt"]] <- PercentageFeatureSet(TX1241_0, pattern = "^MT-")
rm(TX1241_0.data)

# QC check
VlnPlot(
  object = TX1241_0,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1241_0 <- subset(
  x = TX1241_0, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1241_0,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1241_0)


## Day 2
# load and process data
TX1241_2.data <- Read10X(data.dir = './TX1241_2/outs/filtered_feature_bc_matrix/')
TX1241_2 <- CreateSeuratObject(counts = TX1241_2.data$'Gene Expression', project = "TX1241_2", assay = 'RNA')
TX1241_2[["percent.mt"]] <- PercentageFeatureSet(TX1241_2, pattern = "^MT-")
rm(TX1241_2.data)

# QC check
VlnPlot(
  object = TX1241_2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1241_2 <- subset(
  x = TX1241_2, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1241_2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1241_2)


## Day 5
# load and process data
TX1241_5.data <- Read10X(data.dir = './TX1241_5/outs/filtered_feature_bc_matrix/')
TX1241_5 <- CreateSeuratObject(counts = TX1241_5.data$'Gene Expression', project = "TX1241_5", assay = 'RNA')
TX1241_5[["percent.mt"]] <- PercentageFeatureSet(TX1241_5, pattern = "^MT-")
rm(TX1241_5.data)

# QC check
VlnPlot(
  object = TX1241_5,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1241_5 <- subset(
  x = TX1241_5, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1241_5,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1241_5)


## Day 8
# load and process data
TX1241_8.data <- Read10X(data.dir = './TX1241_8/outs/filtered_feature_bc_matrix/')
TX1241_8 <- CreateSeuratObject(counts = TX1241_8.data$'Gene Expression', project = "TX1241_8", assay = 'RNA')
TX1241_8[["percent.mt"]] <- PercentageFeatureSet(TX1241_8, pattern = "^MT-")
rm(TX1241_8.data)

# QC check
VlnPlot(
  object = TX1241_8,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1241_8 <- subset(
  x = TX1241_8, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1241_8,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1241_8)


# merge timecourse data for donor
tx1241.big <- merge(TX1241_0, y = c(TX1241_2, TX1241_5, TX1241_8), add.cell.ids = c("d0", "d2", "d5", "d8"), project = "tx1241")
## optional: remove individual timepoint Seurat objects to free up memory
#rm(TX1241_0)
#rm(TX1241_2)
#rm(TX1241_5)
#rm(TX1241_8)
tx1241.big <- NormalizeData(tx1241.big, normalization.method = "LogNormalize", scale.factor = 10000)
tx1241.big <- CellCycleScoring(tx1241.big, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tx1241.big <- ScaleData(tx1241.big, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(tx1241.big))
#tx1241.big <- ScaleData(tx1241.big, features = rownames(tx1241.big))         ### use instead for scaling without cell cycle marker regression
tx1241.big <- FindVariableFeatures(tx1241.big, selection.method = "vst", nfeatures = 2000)
tx1241.big <- RunPCA(tx1241.big, features = VariableFeatures(object = tx1241.big))
tx1241.big <- FindNeighbors(tx1241.big, dims = 1:25)
tx1241.big <- RunUMAP(tx1241.big, dims=1:25)
tx1241.big <- FindClusters(tx1241.big, resolution = 0.28) #0.3/0.28 for 1241, 0.2/0.22 for 1242





### TX1242 ###

## Day 0
# load and process data
TX1242_0.data <- Read10X(data.dir = './TX1242_0/outs/filtered_feature_bc_matrix/')
TX1242_0 <- CreateSeuratObject(counts = TX1242_0.data$'Gene Expression', project = "TX1242_0", assay = 'RNA')
TX1242_0[["percent.mt"]] <- PercentageFeatureSet(TX1242_0, pattern = "^MT-")
rm(TX1242_0.data)

# QC check
VlnPlot(
  object = TX1242_0,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1242_0 <- subset(
  x = TX1242_0, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1242_0,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1242_0)


## Day 2
# load and process data
TX1242_2.data <- Read10X(data.dir = './TX1242_2/outs/filtered_feature_bc_matrix/')
TX1242_2 <- CreateSeuratObject(counts = TX1242_2.data$'Gene Expression', project = "TX1242_2", assay = 'RNA')
TX1242_2[["percent.mt"]] <- PercentageFeatureSet(TX1242_2, pattern = "^MT-")
rm(TX1242_2.data)

# QC check
VlnPlot(
  object = TX1242_2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1242_2 <- subset(
  x = TX1242_2, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1242_2,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1242_2)


## Day 5
# load and process data
TX1242_5.data <- Read10X(data.dir = './TX1242_5/outs/filtered_feature_bc_matrix/')
TX1242_5 <- CreateSeuratObject(counts = TX1242_5.data$'Gene Expression', project = "TX1242_5", assay = 'RNA')
TX1242_5[["percent.mt"]] <- PercentageFeatureSet(TX1242_5, pattern = "^MT-")
rm(TX1242_5.data)

# QC check
VlnPlot(
  object = TX1242_5,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1242_5 <- subset(
  x = TX1242_5, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1242_5,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1242_5)


## Day 8
# load and process data
TX1242_8.data <- Read10X(data.dir = './TX1242_8/outs/filtered_feature_bc_matrix/')
TX1242_8 <- CreateSeuratObject(counts = TX1242_8.data$'Gene Expression', project = "TX1242_8", assay = 'RNA')
TX1242_8[["percent.mt"]] <- PercentageFeatureSet(TX1242_8, pattern = "^MT-")
rm(TX1242_8.data)

# QC check
VlnPlot(
  object = TX1242_8,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# QC filtering
TX1242_8 <- subset(
  x = TX1242_8, subset = nFeature_RNA > min_feat & nCount_RNA < max_count & percent.mt < mito_pct)

# post-filter QC
VlnPlot(
  object = TX1242_8,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

print(TX1242_8)


# merge timecourse data for donor
tx1242.big <- merge(TX1242_0, y = c(TX1242_2, TX1242_5, TX1242_8), add.cell.ids = c("d0", "d2", "d5", "d8"), project = "tx1242")
## optional: remove individual timepoint Seurat objects to free up memory
#rm(TX1242_0)
#rm(TX1242_2)
#rm(TX1242_5)
#rm(TX1242_8)
tx1242.big <- NormalizeData(tx1242.big, normalization.method = "LogNormalize", scale.factor = 10000)
tx1242.big <- CellCycleScoring(tx1242.big, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tx1242.big <- ScaleData(tx1242.big, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(tx1242.big))
#tx1242.big <- ScaleData(tx1242.big, features = rownames(tx1241.big))         ### use instead for scaling without cell cycle marker regression
tx1242.big <- FindVariableFeatures(tx1242.big, selection.method = "vst", nfeatures = 2000)
tx1242.big <- RunPCA(tx1242.big, features = VariableFeatures(object = tx1242.big))
tx1242.big <- FindNeighbors(tx1242.big, dims = 1:25)
tx1242.big <- RunUMAP(tx1242.big, dims=1:25)
tx1242.big <- FindClusters(tx1242.big, resolution = 0.22)



#### DATA VISUALIZATION (MODIFY FOR DONOR OF INTEREST) ####

# color maps settings for expression data and clusters (tune accordingly)
levels(Idents(tx1241.big))
custom_fill <- scale_fill_viridis(option='mako', discrete = FALSE)
custom_fill_2 <- scale_fill_viridis(option='viridis', discrete = FALSE)
custom_fill_3 <- scale_fill_viridis(option='plasma', discrete = FALSE)
custom_pal <- colorRampPalette(c("gray20", "royalblue2", "aquamarine3", "greenyellow", 'lemonchiffon', 'seashell1'), alpha = TRUE)(7)
custom_cmap <- scale_color_gradientn(colors = custom_pal, limits=NULL)
heatmap_custom_pal <- colorRampPalette(c("gray15", "midnightblue", "darkmagenta", "mediumvioletred", 'firebrick1', 'coral1', 'goldenrod2', 'gold'), alpha = TRUE)(100)
heatmap_custom_pal_2 <- colorRampPalette(c("gray15", "darkgreen", "forestgreen", "lightgreen", "olivedrab3", 'gold', 'goldenrod1'), alpha = TRUE)(100)

mycol <- rgb(0.1, 0.1 , 0.1, 0.05)
test <- colorRampPalette(brewer.pal(11, "Spectral"))(length(levels(Idents(tx1241.big))))
sample_id_cols <- c(test[10], test[8], test[4], test[2])
phase_cols <- c(test[10], test[5], test[3])

# pick which one to use for ordered cluster coloring here:
cluster_cols <- test
my_cols2 <- cluster_cols[order(as.integer(names(cluster_cols)))]
cluster_id_cols <- my_cols2
scales::show_col(my_cols2)

# optional: manually assign ordered anotations based on clustering (levels = identified Seurat clusters; labels = new custom annotation)
tx1241_combo@meta.data$OrderedClusters <- factor(tx1241_combo@meta.data$seurat_clusters, 
                                                 levels = c('8', '3', '4', '7', '6', '0', '1', '2', '5'), 
                                                 labels = c('Naive', 'Mem', 'Antiviral', 'Arrest', 'Hyperprolif', 'AP-eMBC', 'Naive Act Intermediate', 'NF-kB Act', 'Diff MBC/pre-PB'))


### PLOTTING ###

# basic sample identity and global characterization
DimPlot(tx1241.big, cols = cluster_id_cols)
DimPlot(tx1241.big, group.by = 'orig.ident', cols = sample_id_cols)
DimPlot(tx1241.big, group.by = 'Phase', cols = phase_cols)

s1 <- DimPlot(tx1241.big, group.by = 'orig.ident', cols = c(test[10], mycol, mycol, mycol)) & NoLegend() & ggtitle('tx1241_0')
s2 <- DimPlot(tx1241.big, group.by = 'orig.ident', cols = c(mycol, test[8], mycol, mycol)) & NoLegend() & ggtitle('tx1241_2')
s3 <- DimPlot(tx1241.big, group.by = 'orig.ident', cols = c(mycol, mycol, test[4], mycol)) & NoLegend() & ggtitle('tx1241_5')
s4 <- DimPlot(tx1241.big, group.by = 'orig.ident', cols = c(mycol, mycol, mycol, test[2])) & NoLegend() & ggtitle('tx1241_8')

VlnPlot(tx1242.big, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'S.Score', 'G2M.Score'),
        pt.size = 0, ncol = 5, 
        cols=cluster_id_cols) & geom_boxplot(width=0.1, fill="gray", color='black')

VlnPlot(tx1242.big, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'S.Score', 'G2M.Score'),
        pt.size = 0, ncol = 5, group.by='orig.ident',
        cols=sample_id_cols) & geom_boxplot(width=0.1, fill="gray", color='black')

VlnPlot(tx1242.big, c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'S.Score', 'G2M.Score'),
        pt.size = 0, ncol = 5, group.by='Phase',
        cols=phase_cols) & geom_boxplot(width=0.1, fill="gray", color='black')


# violin plots for gene(s) of interest (note: can plot by orig.ident, phase, seurat_clusters, or OrderedClusters)
vln_feature = 'EBI3'

VlnPlot(tx1241.big, vln_feature, group.by='seurat_clusters', 
        pt.size = 0, cols = sample_id_cols) & geom_boxplot(width=0.1, fill="gray", color='black')


# UMAPs for gene(s) of interest
umap_feature = 'EBI3'

FeaturePlot(tx1241.big, umap_feature, pt.size = 1, cols = cluster_id_cols) & custom_cmap & theme_minimal()


# find one vs all cluster markers
big.markers.clust.1 <- FindAllMarkers(tx1241.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
big.markers.clust.2 <- FindAllMarkers(tx1242.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
top40_by_clust.1 <- big.markers.clust.1 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)
top40_by_clust.2 <- big.markers.clust.2 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)
top15_by_clust.1 <- big.markers.clust.1 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top15_by_clust.2 <- big.markers.clust.2 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
# write.csv(top25_by_clust.1, './_merged_top_cluster_markers_tx1241.csv')
# write.csv(top25_by_clust.2, './_merged_top_cluster_markers_tx1242.csv')
print(top40_by_clust.1)
print(top40_by_clust.2)

# find specific cluster markers vs other clusters or groups
cluster_n_vs_m.markers <- FindMarkers(tx1241.big, ident.1 = c(3), ident.2 = c(0, 1, 2, 5, 6), min.pct = 0.25)

# single cell heatmap by cluster and sample ID for top differential genes
DoHeatmap(subset(tx1241.big, downsample = 100), features = c(top15_by_clust.1$gene), size = 3) & custom_fill


# modify for highlighting specific clusters with specific colors
gr = "gray20"
color_vector = c(gr, gr, "tomato1", "goldenrod1", "royalblue1",gr, gr, gr, gr, gr, gr)
DimPlot(tx1241.big) + scale_color_manual(values = color_vector)
# write.csv(cluster_n_vs_m.markers, './c3_vs_c01256_markers.csv')



# feature correlation plotting
n = 200
feat_corr_set = head(VariableFeatures(tx1241.big), n)

av.exp <- AverageExpression(tx1241.big, features = feat_corr_set)$RNA
clust.corr = cor(av.exp)
# correlation by cluster ([0:9] to exclude T cell and CD14 monocytes when 11 phenotypes are identified)
heatmap.2(clust.corr[0:9, 0:9], col = viridis(50), tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "Pearson's R", key.ylab = "", 
          key.ytickfun = "")

bulk.var.exp = as.matrix(rowMeans(av.exp))
bulk.eucl.dist <- as.matrix(dist(bulk.var.exp, method = "manhattan"))
dist.log <- log(bulk.eucl.dist)
dist.log <- replace(dist.log, is.infinite(dist.log), 0)
# euclidean distance among top n variable genes *with diagonal -Inf replaced*
heatmap.2(dist.log, col = magma(50), tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "log distance", key.ylab = "", 
          key.ytickfun = "")


# Pearson's R for top features
n = 200
feat_set = head(VariableFeatures(tx1241.big), n)
# empty matrix
gene_corr_mat <- matrix(0, nrow = n, ncol = n)
cells = NULL
cells <- cells %||% colnames(tx1241.big)
group.by <- NULL
tx1241.big[['ident']] <- Idents(object = tx1241.big)
group.by <- group.by %||% 'ident'
top.feat.data <-  FetchData(
  object = tx1241.big,
  vars = feat_set,
  cells = cells,
  slot = 'data'
)

heatmap.2(cor(top.feat.data, method='pearson'), col = inferno(50), tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "Pearson's R", key.ylab = "", 
          key.ytickfun = "",cexRow=0.6, cexCol=0.6)








##### PSEUDOTIME ANALYSIS #####

# monocle3 (modify for donor)
cds <- as.cell_data_set(tx1241.big) 
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = TRUE)
cds <- learn_graph(cds, use_partition = TRUE)
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  label_principal_points = TRUE,
)

# ordering and pseudotime plot
root = c('Y_3') # manually define root node for graph
cds <- order_cells(cds, root_pr_nodes = root)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
)

cds <- order_cells(cds, root_cells = root)

# assign calculated pseudotime value to original Seurat object (merged infection timecourse data)
tx1241.big[['pseudotime']] <- cds@principal_graph_aux$UMAP$pseudotime
FeaturePlot(tx1241.big, 'pseudotime', pt.size = 1) & custom_cmap_2 & theme_void()

VlnPlot(tx1241.big, 'pseudotime', 
        group.by = 'orig.ident', 
        cols = sample_id_cols,
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black')
VlnPlot(tx1241.big, 'pseudotime', 
        group.by = 'OrderedClusters', 
        cols = my_cols1,
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black')


VlnPlot(tx1241.big, 'pseudotime', pt.size = 0, 
              cols = sample_id_cols) & geom_boxplot(width=0.1, fill="gray", color='black') & NoLegend()
# optional: remove cds object after assigning pseudotime to original Seuart object to free up memory
# rm(cds)


# scatterplot genes of interest in pseudotime with cluster/sample annotation in color
target_genes <- "CCR7"
FeatureScatter(tx1242.big, 'pseudotime', target_genes, shuffle = TRUE, seed = 42)
FeatureScatter(tx1242.big, 'pseudotime', target_genes, shuffle = TRUE, seed = 42, group.by = 'orig.ident')

# get genes from FindMarkers output for plotting in pseudotime (pseudotime heatmap)
top_markers_csv <- read.csv('./_merged_top_cluster_markers_tx1241.csv')
count.data <- GetAssayData(tx1241.big[['RNA']], slot = 'counts')
genes_for_pseudotime <- top_markers_csv[1:200, ]$gene
genes_for_pseudotime <- unique(genes_for_pseudotime)
genes_for_pseudotime <- c(genes_for_pseudotime)
pseudotime.count.data <- count.data[rownames(count.data) %in% genes_for_pseudotime, ]
test_set <- as.matrix(pseudotime.count.data)

mini_set <- test_set[,sample(ncol(test_set), size = 50000), drop = FALSE]
mini_set <- LogNormalize(mini_set)
mini_set <- as.matrix(mini_set)
mini_set <- mini_set[, order(as.numeric(colnames(mini_set)), decreasing = TRUE)]
mini_set <- t(mini_set)
n <- 20
mini_avg <- aggregate(mini_set, list(rep(1:(nrow(mini_set) %/% n + 1), each = n, len = nrow(mini_set))), mean)[-1];

mini_avg <- as.matrix(mini_avg)
heatmap.2(mini_avg, col = heat_pal, tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "expression", key.ylab = "", 
          key.ytickfun = "", dendrogram = 'col' , Rowv = FALSE)



# clusterwise differential gene expression
csv_path <- "./plots_and_figure_output/tx1241_big/_clusterwise_dge/"

clust_dge_mat <- matrix(0, 9, 9)

ind <- 1
for(i in 0:7){
  for(j in ind:8){
    print(paste(i, j))
    test <- read.csv(paste0(csv_path, 'c', i, '_vs_c', j, '_markers.csv'))
    num_dge <- dim(test)[1]
    clust_dge_mat[i+1, j+1] <- num_dge
    clust_dge_mat[j+1, i+1] <- num_dge
  }
  ind <- ind + 1
}

log_clust_dge_mat <- log(clust_dge_mat)
for(d in 1:9){
  log_clust_dge_mat[d, d] <- 0
}

heatmap.2(clust_dge_mat, col = inferno(50), tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "# differential genes", key.ylab = "", 
          key.ytickfun = "", cexRow=1.1, cexCol=1.1, 
          labRow = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8'),
          labCol = c('c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8'))

# set clusters for desired comparison for dge
c1 <- 4
c2 <- 6
c3 <- 7

# load marker csv and find state-transition unique genes
a_genes <- read.csv(paste0(csv_path, 'c', c1, '_vs_c', c2, '_markers.csv'))
b_genes <- read.csv(paste0(csv_path, 'c', c1, '_vs_c', c3, '_markers.csv'))

#cluster_n_vs_m.markers <- FindMarkers(tx1241.big, ident.1 = c(3, 8), ident.2 = c(4), min.pct = 0.25)
#write.csv(cluster_n_vs_m.markers, './plots_and_figure_output/tx1241_big/_clusterwise_dge/c38_vs_c4_markers.csv')
a_genes <- read.csv(paste0(csv_path, 'c2_vs_c5_markers.csv'))
a_genes <- a_genes[order(a_genes[,3],decreasing=TRUE),]
mtx <- head(a_genes, 7)
mtx <- rbind(mtx, tail(a_genes, 7))
MTX <- cbind(as.matrix(mtx[,3]), -as.matrix(log10(mtx[,6])))
marker_vector <- c(mtx[['X']])
new.levels <- c("", "", "", "", "", "", "", "2", "5", "", "")
names(new.levels) <- levels(tx1241.big)
tx1241.big <- RenameIdents(tx1241.big, new.levels)
# downsampled single cells
DoHeatmap(subset(tx1241.big, downsample = 200, ident = c('2', '5')), 
          features = rev(marker_vector), 
          size = 5,
          disp.min = -3,
          disp.max = 3,
          group.colors = c("royalblue1", "tomato1"), group.bar.height = 0.1
          ) & scale_fill_distiller(type='div', palette='RdYlBu') & NoLegend() + theme(text = element_text(size = 20))

# reset active indentity for next test
Idents(tx1241.big) <- 'OrderedClusters'
# tx1241.big <- FindClusters(tx1241.big, resolution = 0.3)

# cluster averaged expression
# feat_corr_set <- rev(marker_vector)
av.exp.mat <- AverageExpression(tx1241.big, features = feat_corr_set, return.seurat = FALSE)
avg.exp.mat <- as.matrix(av.exp.mat$RNA)
log.avg.mat <- log10(avg.exp.mat)
heatmap.2(log.avg.mat[,1:2], tracecol=NA, col = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(50)),
          keysize = 1, key.title= "", key.xlab = "log10(1+Avg Expresion)", key.ylab = "",
          key.ytickfun = "", cexRow=.7, cexCol=1, srtCol = 45, srtRow = 0, Colv = FALSE, Rowv = FALSE,
          dendrogram = c("both"))



DoHeatmap(subset(tx1241.big, downsample = 100), 
          features = c(top15_by_clust.1$gene), 
          size = 5,
          disp.min = -2,
          disp.max = 2, group.by = 'OrderedClusters',
          group.colors = donor_2_colors, group.bar.height = 0.02
          ) & scale_fill_distiller(type='div', palette='RdYlBu') # + theme(text = element_text(size = 8))



# empty gene list
ab_genes <- c()
ab_genes <- paste0(a_genes[['X']], b_genes[['X']])

length(ab_genes)
# ab_genes

c_not_ab_genes <- c()
c_genes <- read.csv(paste0(csv_path, 'c', c2, '_vs_c', c3, '_markers.csv'))
unique_map_2 <- !c_genes[['X']] %in% ab_genes
for(row in 1:length(unique_map_2)){
  if(unique_map_2[row]){
    c_not_ab_genes <- append(c_not_ab_genes, c_genes[['X']][row])
  }
}

length(c_not_ab_genes)
# c_not_ab_genes


color_vector = c(gr, gr, gr, gr, gr ,gr, gr, gr, gr, gr, gr)
color_vector[c1 + 1] <- "tomato1"
color_vector[c2 + 1] <- "royalblue1"
color_vector[c3 + 1] <- "goldenrod1"

DimPlot(tx1241.big) + scale_color_manual(values = color_vector)
head(c_not_ab_genes, 100)

FeaturePlot(tx1241.big, c_not_ab_genes[1]) & custom_cmap




# GET CLUSTER MEMBERSHIP (# OF CELLS) BY SAMPLE IN MERGED DATA
# As described in: https://www.biostars.org/p/399234/
md <- tx1241.big@meta.data %>% as.data.table
## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md[, .N, by = c("orig.ident", "OrderedClusters")]
## with additional casting after the counting

membership_table <- md[, .N, by = c("orig.ident", "OrderedClusters")] %>% dcast(., orig.ident ~ OrderedClusters, value.var = "N")
membership_table[is.na(membership_table)] <- 0

x <- colSums(membership_table[1:4, 2:10])
y <- rbind(x, x)
z <- rbind(y, y)
pct <- as.matrix(membership_table[1:4, 2:10]) / z
pct <- round(pct, 3)

membership_table2 <- md[, .N, by = c("Phase", "ident")] %>% dcast(., Phase ~ ident, value.var = "N")
test2 <- membership_table2

membership_table3 <- md[, .N, by = c("percent.mt", "ident")] %>% dcast(., percent.mt ~ ident, value.var = "N")
test3 <- membership_table3


heatmap.2(as.matrix(test[1:4, 2:12]), col = heatmap_custom_pal, tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "# of cells", key.ylab = "", 
          key.ytickfun = "", cexRow=1.1, cexCol=1.1, 
          xlab = 'Suerat Cluster ID',
          labRow = c('tx1242_0', 'tx1242_2', 'tx1242_5', 'tx1242_8'),
          labCol = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'),
          cellnote = pct, notecol = 'dodgerblue', notecex = 1.4)


heatmap.2(as.matrix(test2[1:3, 2:12]), col = heatmap_custom_pal_2, tracecol=NA, 
          keysize = 1, key.title= "", key.xlab = "# of cells", key.ylab = "", 
          key.ytickfun = "", cexRow=1.1, cexCol=1.1, 
          xlab = 'Suerat Cluster ID',
          labRow = c("G1", "G2M", "S"),
          labCol = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'),
          cellnote = as.matrix(test2[1:3, 2:12]), notecol = 'dodgerblue', notecex = 1.4)


# tx1241.big
sorted_pct <- cbind(pct[ , 1], pct[ , 2], pct[ , 3], pct[ , 4], pct[ , 5], pct[ , 6], pct[ , 7], pct[ , 8], pct[ , 9])
# tx1242.big
sorted_pct <- cbind(pct[ , 4], pct[ , 8], pct[ , 3], pct[ , 10], pct[ , 7], pct[ , 11], pct[ , 5], pct[ , 9], pct[ , 6], pct[ , 2], pct[ , 1])


x <- 1:9
y_0 <- sorted_pct[1, ]
y_2 <- sorted_pct[2, ]
y_5 <- sorted_pct[3, ]
y_8 <- sorted_pct[4, ]


plot(x, y_0, col = sample_id_cols[1], pch=19, xaxt="none", ylab="Relative Enrichment", xlab="Seurat Cluster ID") & 
  axis(1, at = seq(1, 9, by = 1), las=2,labels=c('c8', 'c3', 'c4', 'c7', 'c6', 'c0', 'c1', 'c2', 'c5'))
points(x, y_2, col = sample_id_cols[2], pch=19)
points(x, y_5, col = sample_id_cols[3], pch=19)
points(x, y_8, col = sample_id_cols[4], pch=19)
lo <- loess(y_0~x)
lines(predict(lo), col=sample_id_cols[1], lwd=4)
lo <- loess(y_2~x)
lines(predict(lo), col=sample_id_cols[2], lwd=4)
lo <- loess(y_5~x)
lines(predict(lo), col=sample_id_cols[3], lwd=4)
lo <- loess(y_8~x)
lines(predict(lo), col=sample_id_cols[4], lwd=4)
legend('topright', legend = c('tx1242_0', 'tx1242_2', 'tx1242_5', 'tx1242_8'),
       col = sample_id_cols, pch=19)



# NEW CMAP QC VIOLIN PLOTS
VlnPlot(tx1242.big, 'percent.mt', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'percent.mt', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 


VlnPlot(tx1242.big, 'nFeature_RNA', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'nFeature_RNA', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 


VlnPlot(tx1242.big, 'nCount_RNA', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'nCount_RNA', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 


VlnPlot(tx1242.big, 'G2M.Score', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'G2M.Score', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 


VlnPlot(tx1242.big, 'S.Score', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'S.Score', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 


VlnPlot(tx1242.big, 'pseudotime', 
        cols = cluster_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 
VlnPlot(tx1242.big, 'pseudotime', 
        group.by = 'orig.ident', 
        cols = sample_id_cols, 
        pt.size = 0) & geom_boxplot(width=0.1, fill="gray", color='black') 



# Interactive 3d UMAP plotting
# Reference: https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0/blob/master/3D%20UMAP%20Plotting%20v1.3.R

tx1241.big <- RunUMAP(tx1241.big, dims = 1:30, n.components = 3L)
umap_1 <- tx1241.big[["umap"]]@cell.embeddings[,1]
umap_2 <- tx1241.big[["umap"]]@cell.embeddings[,2]
umap_3 <- tx1241.big[["umap"]]@cell.embeddings[,3]
Embeddings(object = tx1241.big, reduction = "umap")


gene_to_3d_plot = "EBI3"

plot.data <- FetchData(object = tx1241.big, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters", gene_to_3d_plot))
plot.data$label <- paste(rownames(plot.data))
# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~EBI3,
               colors = "inferno",
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 2, width=1.5), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube 







# merging with other scRNA datasets (requires other data -- example here is with GM12878 and Day 8 from this resource)
# load and merge with other datasets as shown previously for timecourse merging
# (assuming gm12878_and_day8 Seurat object has been created via merging -- customize levels and labels based on clustering resolution and phenotype by gene expression):
gm12878_and_day8@meta.data$OrderedClusters <- factor(gm12878_and_day8@meta.data$seurat_clusters, 
                                               levels = c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26'),
                                               labels = c('d8_2','d8_0','d8_6','d8_1','d8_0','lcl_5',
                                                          'lcl_int','d8_5','lcl_5','lcl_5','d8_2','d8_1','d8_4',
                                                          'lcl_int','lcl_6','d8_1','d8_5','lcl_5','d8_7','d8_7',
                                                          'lcl_int','lcl_2','d8_0','lcl_int','d8_9','lcl_int','d8_10'))



# labeled scatterplot of top variable features
options(ggrepel.max.overlaps = 100)
top50 <- head(VariableFeatures(tx1241.big), 50)
plot1 <- VariableFeaturePlot(tx1241.big)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE) + theme(legend.position = "none") 
plot2



# dotplot generation
d <- DotPlot(tx1241.big, features = c('EBI3', 'NFKB2', 'ICAM1', 'TNFAIP3', 'CD27', 'CD38', 'PRDM1', 'XBP1'), 
        col.min = -3, col.max = 3, dot.scale = 7, idents = c('0','1','2','3','4','5','6','7','8'),
        group.by="OrderedClusters") & custom_cmap
d + theme(axis.text.x = element_text(angle = 90))
d

# heatmap generation for average expression
h <- DoHeatmap(avg.exp.mat, 
          features = feat_corr_set,
          size = 5,
          disp.min = -2,
          disp.max = 4,
          draw.lines = FALSE,
          group.colors = lcl_map, group.bar.height = 0.02
) & scale_fill_gradientn(colors = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(11))) # + theme(text = element_text(size = 8))
h



# biological replicate QC -- correlation plotting (requires individually loaded timepoints for each donor)
TX1242_0.data <- Read10X(data.dir = './TX1242_0/outs/filtered_feature_bc_matrix/')
TX1242_0 <- CreateSeuratObject(counts = TX1242_0.data$'Gene Expression', project = "TX1242_0", assay = 'RNA')
TX1242_0[["percent.mt"]] <- PercentageFeatureSet(TX1242_0, pattern = "^MT-")
rm(TX1242_0.data)
TX1242_2.data <- Read10X(data.dir = './TX1242_2/outs/filtered_feature_bc_matrix/')
TX1242_2 <- CreateSeuratObject(counts = TX1242_2.data$'Gene Expression', project = "TX1242_2", assay = 'RNA')
TX1242_2[["percent.mt"]] <- PercentageFeatureSet(TX1242_2, pattern = "^MT-")
rm(TX1242_2.data)
TX1242_5.data <- Read10X(data.dir = './TX1242_5/outs/filtered_feature_bc_matrix/')
TX1242_5 <- CreateSeuratObject(counts = TX1242_5.data$'Gene Expression', project = "TX1242_5", assay = 'RNA')
TX1242_5[["percent.mt"]] <- PercentageFeatureSet(TX1242_5, pattern = "^MT-")
rm(TX1242_5.data)
TX1242_8.data <- Read10X(data.dir = './TX1242_8/outs/filtered_feature_bc_matrix/')
TX1242_8 <- CreateSeuratObject(counts = TX1242_8.data$'Gene Expression', project = "TX1242_8", assay = 'RNA')
TX1242_8[["percent.mt"]] <- PercentageFeatureSet(TX1242_8, pattern = "^MT-")
rm(TX1242_8.data)


rep1 <- AverageExpression(TX1242_5, assays = 'RNA', return.seurat = TRUE, group.by = 'orig.ident')
rep2 <- AverageExpression(TX1242_8, assays = 'RNA', return.seurat = TRUE, group.by = 'orig.ident')
rep1_rna <- rep1@assays$RNA@data
rep2_rna <- rep2@assays$RNA@data
rep1_rna <- rep1_rna[!is.na(rep1_rna) & !is.infinite(rep1_rna), ]
rep2_rna <- rep2_rna[!is.na(rep2_rna) & !is.infinite(rep2_rna), ]
rep1_rna <- as.data.frame(rep1_rna)
rep2_rna <- as.data.frame(rep2_rna)
shared_features <- intersect(rownames(rep1_rna), rownames(rep2_rna))
rep1_rna <- rep1_rna[shared_features, ]
rep2_rna <- rep2_rna[shared_features, ]
rep1_rna <- log2(rep1_rna + 1)
rep2_rna <- log2(rep2_rna + 1)
rep_corr <- round(cor(rep1_rna, rep2_rna), digits = 3)
model <- lm(rep1_rna ~ rep2_rna)
r_sq <- round(summary(model)$adj.r.squared, digits = 3)
# d0_color <- rgb(.1,.5,1,0.02)
# d2_color <- rgb(.2,.8,.5,0.02)
# d5_color <- rgb(.8,.4,.1,0.02)
# d8_color <- rgb(1,.2,0,0.02)
plot(rep1_rna, rep2_rna, log = 'xy', 
     pch = 19, col = rgb(0,0,0,.02),
     main = paste0('Pearson R = ', rep_corr, '\n', 'Adj. r^2 = ', r_sq))
abline(model, col = 'firebrick3', lwd = 2)




# example usage to generate zero-preserving imputation with ALRA (sets default assay to ALRA)
TX1241_8 <- RunALRA(TX1241_8)

# compare differences with and without imputation
DefaultAssay(TX1241_8) <- 'RNA'
FeaturePlot(TX1241_8, 'TBX21', reduction = 'umap', pt.size = 2) & custom_cmap & theme_void()

DefaultAssay(TX1241_8) <- 'alra'
FeaturePlot(TX1241_8, 'TBX21', reduction = 'umap', pt.size = 2) & custom_cmap & theme_void()



