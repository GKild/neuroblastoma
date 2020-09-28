#process rms3

paths=c("/lustre/scratch119/realdata/mdt1/team274/ek12/RMS/RMS3rdRun/cellranger302_count_32349_CG_SB_NB8608056_GRCh38-1_2_0/filtered_feature_bc_matrix/",
        "/lustre/scratch119/realdata/mdt1/team274/ek12/RMS/RMS3rdRun/cellranger302_count_32349_CG_SB_NB8608057_GRCh38-1_2_0/filtered_feature_bc_matrix/",
        "/lustre/scratch119/realdata/mdt1/team274/ek12/RMS/RMS3rdRun/cellranger302_count_32349_CG_SB_NB8608058_GRCh38-1_2_0/filtered_feature_bc_matrix/")
names(paths)=c("NB8608056", "NB8608057", "NB8608058")

rms3_mtx=Read10X(data.dir = paths)
dim(rms3_mtx)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#MT cut-off
mtCut = 0.2
#nGene Cut
minNumGenes = 300
#nUMI cut
minNumUMI = 1000
#How many doublets in a cluster before you throw the entire cluster
maxScrubFrac = 0.3
#UMAP parameters to try in grid
minDists = c(0,.01,.1,.3,.5)
NNs = c(5,10,20,50,100)
#Final UMAP parameters
finalMinDist = 0.5
finalNN = 50
#Clustering resolution
clusterRes=0.8
#npcs
nPCs=50

process_10X=function(mtx, file_name){
  srat = CreateSeuratObject(mtx)
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
  pdf(paste0(file_name, "_vln_plots.pdf"), height = 10, width = 10)
  gg=VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot(gg)
  dev.off()
  srat <- subset(srat, subset = nFeature_RNA > minNumGenes & percent.mt < mtCut*100 & nCount_RNA>minNumUMI)
  keep_features <- rownames(srat@assays$RNA@data)[-grep("^MT-",rownames(srat@assays$RNA@data))]
  srat <- subset(srat,features= keep_features)
  srat <- NormalizeData(srat)
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  nPCs <- 75
  srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat <- ScaleData(srat, features = rownames(srat))
  srat <- RunPCA(srat, npcs = nPCs)
  srat <- FindNeighbors(srat, dims=1:nPCs)
  srat <- FindClusters(srat, resolution = clusterRes)
  srat <- RunUMAP(srat, dims=1:nPCs, min.dist = finalMinDist, n.neighbors = finalNN)
  pdf(paste0(file_name, "_umap.pdf"), height = 14, width = 14)
  gg=DimPlot(srat)
  plot(gg)
  dev.off()
  return(srat)
}


srat_rms3=process_10X(rms3_mtx, "rms3")
VlnPlot(srat_rms3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

FeaturePlot(srat_rms3, c('HBB','HBA1','PECAM1','PTPRC','PDGFRB'))
FeaturePlot(srat_rms3, c('MYOD1', "MYOG", "PAX7", "FOXO1", "MYF5", "DES"))
DimPlot(srat_rms3, group.by = "Phase")   
dim(srat_rms3)


DimPlot(srat_rms3)



srat_rms3@meta.data$idents_for_plot="x"



srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(0))]="tumour cluster 1"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(1))]="tumour cluster 2"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(2))]="tumour cluster 3"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(3))]="tumour cluster 4"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(4))]="tumour cluster 5"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(5))]="tumour cluster 6"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(6))]="tumour cluster 7"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(7))]="tumour cluster 8"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(10))]="tumour cluster 9"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(13))]="tumour cluster 10"

srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(8,9))]="immune"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(12,15))]="vascular endothelium"
srat_rms3@meta.data$idents_for_plot[which(srat_rms3@meta.data$seurat_clusters%in%c(11,14))]="mesenchyme"

DimPlot(srat_rms3, group.by = "idents_for_plot")

srat_rms3@meta.data$idents_for_plot=factor(srat_rms3@meta.data$idents_for_plot,
                                              levels=c("tumour cluster 1", "tumour cluster 2",
                                                       "tumour cluster 3", "tumour cluster 4",
                                                       "tumour cluster 5", "tumour cluster 6",
                                                       "tumour cluster 7", "tumour cluster 8",
                                                       "tumour cluster 9", "tumour cluster 10",
                                                       "mesenchyme","immune","vascular endothelium"))


DimPlot(srat_rms3, group.by = "idents_for_plot",
        cols = c(colfunc(10),"#117733", "#DDCC77", "#CC6677"),
        pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank())+ylab("UMAP 2") +xlab("UMAP 1")


