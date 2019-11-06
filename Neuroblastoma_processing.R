library(Seurat)
sample_paths_neurobl1 <- c("Neuroblastoma1/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/",
       "Neuroblastoma1/4602STDY7685341/outs/filtered_gene_bc_matrices/GRCh38/",
       "Neuroblastoma1/4602STDY7685342/outs/filtered_gene_bc_matrices/GRCh38/")
names(sample_paths_neurobl1) <- c("STDY7685340", "STDY7685341", "STDY7685342")

sample_paths_neurobl2 <- c("Neuroblastoma2/4602STDY7733084/outs/filtered_gene_bc_matrices/GRCh38/",
                           "Neuroblastoma2/4602STDY7733085/outs/filtered_gene_bc_matrices/GRCh38/",
                           "Neuroblastoma2/4602STDY7733086/outs/filtered_gene_bc_matrices/GRCh38/")
names(sample_paths_neurobl2) <-c("STDY7733084","STDY7733085", "STDY7733086")

neuroblastoma1_mtx <- Read10X(data.dir = sample_paths_neurobl1)
neuroblastoma2_mtx <- Read10X(data.dir = sample_paths_neurobl2)
neuroblastoma1 <- CreateSeuratObject(counts=neuroblastoma1_mtx,project = "Neuroblastoma", min.cells = 3, min.features = 200)
neuroblastoma2 <- CreateSeuratObject(counts=neuroblastoma2_mtx,project = "Neuroblastoma", min.cells = 3, min.features = 200)

###run through the pipeline w/ NB1

neuroblastoma1[["percent.mt"]] <- PercentageFeatureSet(neuroblastoma1, pattern = "^MT-")
levels(neuroblastoma1@active.ident) <- rep("nb1", 3)
VlnPlot(neuroblastoma1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(neuroblastoma1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(neuroblastoma1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#rm cells with >10% mito reads and <200 features
neuroblastoma1 <- subset(neuroblastoma1, subset = nFeature_RNA > 200 & percent.mt < 10)
#log normalize 
neuroblastoma1 <- NormalizeData(neuroblastoma1, normalization.method = "LogNormalize", scale.factor = 10000)
#find 2k most variable genes
neuroblastoma1 <- FindVariableFeatures(neuroblastoma1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(neuroblastoma1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(neuroblastoma1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#scale the data
all.genes <- rownames(neuroblastoma1)
neuroblastoma1 <- ScaleData(neuroblastoma1, features = all.genes)
#run PCA
neuroblastoma1 <- RunPCA(neuroblastoma1, features = VariableFeatures(object = neuroblastoma1), npcs = 100)
print(neuroblastoma1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(neuroblastoma1, dims = 40:50, reduction = "pca")
DimPlot(neuroblastoma1, reduction = "pca", group.by = "orig.ident")
DimHeatmap(neuroblastoma1, dims = 1:15, cells = 500, balanced = TRUE)
#jackstraw
neuroblastoma1 <- JackStraw(neuroblastoma1, num.replicate = 200, dims = 100)
neuroblastoma1 <- ScoreJackStraw(neuroblastoma1, dims = 1:100)
JackStrawPlot(neuroblastoma1, dims = 1:100) +theme(legend.position = "none")
JackStrawPlot(neuroblastoma1, dims = 1:100)
#elbow
ElbowPlot(neuroblastoma1, ndims = 100)

neuroblastoma1_10 <- FindNeighbors(neuroblastoma1, dims = 1:10)
neuroblastoma1_20 <- FindNeighbors(neuroblastoma1, dims = 1:20)
neuroblastoma1_50 <- FindNeighbors(neuroblastoma1, dims = 1:50)
neuroblastoma1_100 <- FindNeighbors(neuroblastoma1, dims = 1:100)

neuroblastoma1_10 <- FindClusters(neuroblastoma1_10, resolution = 0.5)
neuroblastoma1_10 <- FindClusters(neuroblastoma1_10, resolution = 0.8)
neuroblastoma1_10 <- FindClusters(neuroblastoma1_10, resolution = 1)

neuroblastoma1_20 <- FindClusters(neuroblastoma1_20, resolution = 0.5)
neuroblastoma1_20 <- FindClusters(neuroblastoma1_20, resolution = 0.8)
neuroblastoma1_20 <- FindClusters(neuroblastoma1_20, resolution = 1)

neuroblastoma1_50 <- FindClusters(neuroblastoma1_50, resolution = 0.5)
neuroblastoma1_50 <- FindClusters(neuroblastoma1_50, resolution = 0.8)
neuroblastoma1_50 <- FindClusters(neuroblastoma1_50, resolution = 1)

neuroblastoma1_100 <- FindClusters(neuroblastoma1_100, resolution = 0.5)
neuroblastoma1_100 <- FindClusters(neuroblastoma1_100, resolution = 0.8)
neuroblastoma1_100 <- FindClusters(neuroblastoma1_100, resolution = 1)



neuroblastoma1_10 <- RunUMAP(neuroblastoma1_10, dims = 1:10)
neuroblastoma1_20 <- RunUMAP(neuroblastoma1_20, dims = 1:20)
neuroblastoma1_50 <- RunUMAP(neuroblastoma1_50, dims = 1:50)
neuroblastoma1_100 <- RunUMAP(neuroblastoma1_100, dims = 1:100)

umap_10_05 <-DimPlot(neuroblastoma1_10, reduction = "umap", group.by = "RNA_snn_res.0.5") + ggtitle("UMAP 10 PCs, res 0.5")
umap_10_08 <-DimPlot(neuroblastoma1_10, reduction = "umap", group.by = "RNA_snn_res.0.8") + ggtitle("UMAP 10 PCs, res 0.8")
umap_10_1 <-DimPlot(neuroblastoma1_10, reduction = "umap",group.by = "RNA_snn_res.1") + ggtitle("UMAP 10 PCs, res 1")

umap_20_05 <-DimPlot(neuroblastoma1_20, reduction = "umap", group.by = "RNA_snn_res.0.5")+ ggtitle("UMAP 20 PCs, res 0.5")
umap_20_08 <-DimPlot(neuroblastoma1_20, reduction = "umap", group.by = "RNA_snn_res.0.8")+ ggtitle("UMAP 20 PCs, res 0.8")
umap_20_1 <-DimPlot(neuroblastoma1_20, reduction = "umap", group.by = "RNA_snn_res.1")+ ggtitle("UMAP 20 PCs, res 1")

umap_50_05 <-DimPlot(neuroblastoma1_50, reduction = "umap", group.by = "RNA_snn_res.0.5")+ ggtitle("UMAP 50 PCs, res 0.5")
umap_50_08 <-DimPlot(neuroblastoma1_50, reduction = "umap", group.by = "RNA_snn_res.0.8")+ ggtitle("UMAP 50 PCs, res 0.8")
umap_50_1 <-DimPlot(neuroblastoma1_50, reduction = "umap", group.by = "RNA_snn_res.1")+ ggtitle("UMAP 50 PCs, res 1")

umap_100_05 <-DimPlot(neuroblastoma1_100, reduction = "umap", group.by = "RNA_snn_res.0.5")+ ggtitle("UMAP 100 PCs, res 0.5")
umap_100_08 <-DimPlot(neuroblastoma1_100, reduction = "umap", group.by = "RNA_snn_res.0.8")+ ggtitle("UMAP 100 PCs, res 0.8")
umap_100_1 <-DimPlot(neuroblastoma1_100, reduction = "umap", group.by = "RNA_snn_res.1")+ ggtitle("UMAP 100 PCs, res 1")
CombinePlots(plots=list(umap_10_05, umap_10_08, umap_10_1,
                        umap_20_05, umap_20_08, umap_20_1,
                        umap_50_05, umap_50_08, umap_50_1,
                        umap_100_05, umap_100_08, umap_100_1), ncol = 3)

neuroblastoma1_50$seurat_clusters <- neuroblastoma1_50$RNA_snn_res.0.5
neuroblastoma1_50@active.ident <- neuroblastoma1_50$RNA_snn_res.0.5
#will be using 50PCs, res 0.5 clusters for all downstream stuff

saveRDS(neuroblastoma1_50, file="")


neuroblastoma1_50 <-readRDS("Neuroblastoma/Neuroblastoma1/neuroblastoma_1_50pcs.rds")
#find all markers 
neuroblastoma1_50@active.ident
neuroblastoma1_50$seurat_clusters <- neuroblastoma1_50$RNA_snn_res.0.5
neuroblastoma1_50@active.ident <- neuroblastoma1_50$RNA_snn_res.0.5
all.markers <- FindAllMarkers(neuroblastoma1_50)
a <-all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) 
features_0_5 <- a$gene[1:12]
features_6_10 <- a$gene[13:22]
features_11_12 <- a$gene[23:26]
VlnPlot(neuroblastoma1_50, features = c("MS4A1", "IGKC"))
VlnPlot(neuroblastoma1_50, features = c("GNLY", "GZMB"))
FeaturePlot(neuroblastoma1_50, features = features_0_5)
FeaturePlot(neuroblastoma1_50, features = features_6_10)
FeaturePlot(neuroblastoma1_50, features = features_11_12)
FeaturePlot(neuroblastoma1_50, features = "PTPRC")

genexp <-sort(Matrix::rowMeans(neuroblastoma1_50@assays$RNA@counts), decreasing = T)
names(genexp)[1:10]
FeaturePlot(neuroblastoma1_50, features = names(genexp)[1:12])     
