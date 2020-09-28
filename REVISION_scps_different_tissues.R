library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
scp_mtx=readMM('all_scp_matrix.mtx')

obs=read.csv("all_scps/obs.csv")
cells=read.csv("all_scps/var.csv")
rownames(cells)=cells$index
rownames(scp_mtx)=as.character(obs$index)
colnames(scp_mtx)=as.character(cells$index)

scp_seurat=CreateSeuratObject(scp_mtx, meta.data = cells)
process_10x =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


scp_seurat=process_10x(scp_seurat)

scp_seurat@meta.data$new_labs="x"
scp_seurat@meta.data$new_labs[grep("bone",scp_seurat@meta.data$organ)]="bone SCPs"
scp_seurat@meta.data$new_labs[grep("skin",scp_seurat@meta.data$organ)]="skin SCPs"
scp_seurat@meta.data$new_labs[grep("gut",scp_seurat@meta.data$organ)]="gut SCPs"
scp_seurat@meta.data$new_labs[grep("adrenal",scp_seurat@meta.data$organ)]="adrenal SCPs"
FeaturePlot(scp_seurat, c("SOX2", "SOX10", "ERBB3", "MPZ", "PLP1", "ASCL1", "DLL3"))
scp_seurat=subset(scp_seurat, idents=c(0,6,17,14), invert=T)

scp_pred=predictSimilarity(train_with_podos, scp_seurat@assays$RNA@counts, scp_seurat@meta.data$new_labs, minGeneMatch = 0.8)

similarityHeatmap(scp_pred, column_order=lr_col_ord)
 

FeaturePlot(scp_seurat, c("SOX2", "SOX10", "ERBB3", "MPZ", "PLP1", "ASCL1", "DLL3"))
