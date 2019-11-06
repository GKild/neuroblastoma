###run through the pipeline w/ NB2

neuroblastoma2[["percent.mt"]] <- PercentageFeatureSet(neuroblastoma2, pattern = "^MT-")
levels(neuroblastoma2@active.ident) <- rep("nb2", 3)
VlnPlot(neuroblastoma2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(neuroblastoma2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(neuroblastoma2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#rm cells with >10% mito reads and <200 features
neuroblastoma2 <- subset(neuroblastoma2, subset = nFeature_RNA > 200 & percent.mt < 10)
#log normalize 
neuroblastoma2 <- NormalizeData(neuroblastoma2, normalization.method = "LogNormalize", scale.factor = 10000)
#find 2k most variable genes
neuroblastoma2 <- FindVariableFeatures(neuroblastoma2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(neuroblastoma2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(neuroblastoma2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#scale the data
all.genes <- rownames(neuroblastoma2)
neuroblastoma2 <- ScaleData(neuroblastoma2, features = all.genes)
#run PCA
neuroblastoma2 <- RunPCA(neuroblastoma2, features = VariableFeatures(object = neuroblastoma2))
print(neuroblastoma2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(neuroblastoma2, dims = 1:10, reduction = "pca")
DimPlot(neuroblastoma2, reduction = "pca", group.by = "orig.ident")
DimHeatmap(neuroblastoma2, dims = 1:15, cells = 500, balanced = TRUE)
#jackstraw
neuroblastoma2 <- JackStraw(neuroblastoma2, num.replicate = 100)
neuroblastoma2 <- ScoreJackStraw(neuroblastoma2, dims = 1:20)
JackStrawPlot(neuroblastoma2, dims = 1:20)

#elbow
ElbowPlot(neuroblastoma2, ndims = 20)
neuroblastoma2 <- FindNeighbors(neuroblastoma2, dims = 1:10)
neuroblastoma2 <- FindClusters(neuroblastoma2, resolution = 0.5)
neuroblastoma2 <- RunUMAP(neuroblastoma2, dims = 1:10)
DimPlot(neuroblastoma2, reduction = "umap")
nb2_markers <-FindAllMarkers(neuroblastoma2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
b <-nb2_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
nb2_top <- b$gene
FeaturePlot(neuroblastoma2, features = nb2_top)
