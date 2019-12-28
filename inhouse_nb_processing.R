#in-house NB samples

library(Matrix)
library(MASS)
library(dplyr)
library(Seurat)
library(ggplot2)
library(glmnet)
source("gene_sets.R")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
#MT cut-off
#DEFINE PARAMS
exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
#MT cut-off
mtCut = 20
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
clusterRes=1.0
#Populations of genes to exclude from the analysis

#Number of PCs to use for quick analysis during QC
defNumPCs=75

all_paths = c("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843576/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843577/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843578/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733084/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733085/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733086/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004894/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004902/outs/filtered_gene_bc_matrices/GRCh38/",
               "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004910/outs/filtered_gene_bc_matrices/GRCh38/")
names(all_paths) = c("STDY7843576", "STDY7843577", "STDY7843578","STDY7733084","STDY7733085", "STDY7733086", "STDY8004894",
                     "STDY8004902", "STDY8004910")

neuroblastoma_mtx = Read10X(data.dir = all_paths)

process_10x_nb =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > minNumGenes & nCount_RNA > minNumUMI & percent.mt < mtCut)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
  keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = defNumPCs)
  srat = FindNeighbors(srat, dims=1:defNumPCs)
  srat = FindClusters(srat, resolution = clusterRes)
  srat = RunUMAP(srat, dims=1:defNumPCs, min.dist = finalMinDist, n.neighbors = finalNN)
  return(srat)
}

nb_srat=process_10x_nb(neuroblastoma_mtx)




inhouse_nb = process_10x_nb_keepgenes(neuroblastoma_mtx)
dim(inhouse_nb)
nb_srat=FindClusters(nb_srat, resolution = 1)
DimPlot(nb_srat, group.by = R)
nb_srat@meta.data
inhouse_nb@meta.data$seurat_clusters

nb_srat.markers <- FindAllMarkers(nb_srat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x=nb_srat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
