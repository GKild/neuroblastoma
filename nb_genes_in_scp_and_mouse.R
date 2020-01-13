#Neuroblastoma genes in SCP and mouse

#read baby adr ref
adr=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
#add mature adrenal neuronal dataset
mature_count_table =read.table("matureAdrenalNeuronalCounts.tsv", header = T)
mature_metadata = read.table("matureAdrenalNeuronalMetadata.tsv", header = T)
#make it into seurat object
mature_matrix = as.matrix(mature_count_table)
colnames(mature_matrix) = rownames(mature_metadata)
mature_srat=CreateSeuratObject(mature_matrix, meta.data = mature_metadata)

#mouse data, three timepoints

e10_5_trunk=read.table("mouse_neural_crest/GSE129114_E10.5_trunk_Wnt1_counts.txt", header = T)
e10_5_trunk=as.matrix(e10_5_trunk)

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


process_10x_mouse =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > minNumGenes & nCount_RNA > minNumUMI & percent.mt < mtCut)
  #srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  #srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  #exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
  #keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  #srat=subset(srat, features=keep_features)
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

m1=process_10x_mouse(e10_5_trunk)
m1

DimPlot(m1, label = T, label.size=6)
FeaturePlot(m1, features=c("Sox2", "Sox10", "Phox2a", "Phox2b", "Isl1", "Mycn"))
m1.markers <- FindAllMarkers(m1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x=m1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
