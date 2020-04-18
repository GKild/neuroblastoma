library(Seurat)
source("Logistic_regression.R")
source("gene_sets.R")
adr_all <- readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)

mtCut = 20
#nGene Cut
minNumGenes = 200
#nUMI cut
minNumUMI = 500
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
#Number of PCs to use for quick analysis during QC
defNumPCs=75
source("quickmarkers.R")
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb <- data
genes_10x <- read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
#main matrix 
mtx <- Jan_nb$data$counts
#column metadata 
meta_data <- Jan_nb$col.annot
#calculate fraction of spikeins
meta_data$spikeIn_fraction<- colSums(mtx[grep("^ERCC-",rownames(mtx)),])/colSums(mtx)
#remove spikeins
mtx <- mtx[-grep("^ERCC-",rownames(mtx)),]
#calculate fraction of reads lost by using 10X ensembl IDs, and fraction of reads lost from  MT genes that are NOT in 10X. 
mtx_10x <-mtx[which(rownames(mtx)%in%genes_10x$V1),]
not_10x<-mtx[which(!rownames(mtx)%in%genes_10x$V1),]
meta_data$reads_lost <-(colSums(mtx)-colSums(mtx_10x))/colSums(mtx)
rownames(not_10x)<-Jan_nb$row.annot$gene[match(rownames(not_10x),Jan_nb$row.annot$ensemblId)]
meta_data$mt_reads_lost <-colSums(not_10x[grep("^MT-",rownames(not_10x)),])/colSums(mtx)

patient_ids <-unique(Jan_nb$col.annot$PatientId)
cells_each_sample <-sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})

cells_list <- list("NB039_Tumour_organoids_in_20perc_HP_3D"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO_TGFb_inh"=cells_each_sample$NB039$`Tumour organoids in OM_WNT_RSPO_TGFb inh`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO"=cells_each_sample$NB039$`Tumour organoids in OM _WNT_RSPO`, 
                   "NB060_Fresh_biopsy"=as.character(cells_each_sample$NB060), 
                   "NB067_Fresh_biopsy"=as.character(cells_each_sample$NB067),
                   "NB086_Fresh_biopsy"=as.character(cells_each_sample$NB086),
                   "NB098_Fresh_biopsy"=cells_each_sample$NB098$`Fresh biopsy`,
                   "NB098_peripheral_blood_from_fresh_biopsy"=cells_each_sample$NB098$`peripheral blood from fresh biopsy`,
                   "NB098_TILs_from_fresh_biopsy"=cells_each_sample$NB098$`TILs from fresh biopsy`,
                   "NB098_Tumour_organoids"=cells_each_sample$NB098$`Tumour organoids`, 
                   "NB106_Fresh_biopsy"=as.character(cells_each_sample$NB106),
                   "NB107_Fresh_biopsy"=as.character(cells_each_sample$NB107),
                   "NB123_Fresh_biopsy"=as.character(cells_each_sample$NB123),
                   "NB124_Fresh_biopsy"=as.character(cells_each_sample$NB124),
                   "NB125_Fresh_biopsy"=as.character(cells_each_sample$NB125),
                   "NB130_Fresh_biopsy"=as.character(cells_each_sample$NB130),
                   "NB132_Fresh_biopsy"=as.character(cells_each_sample$NB132),
                   "NB138_Fresh_biopsy"=as.character(cells_each_sample$NB138),
                   "NB151_Fresh_biopsy"=as.character(cells_each_sample$NB151),
                   "NB152_Fresh_biopsy"=cells_each_sample$NB152$`Fresh biopsy`,
                   "NB152_Tumour_organoids_P1"=cells_each_sample$NB152$`Tumour organoids P1`)


#add the unique sample to metadata 

meta_data$unique_sample <- "sample"
for (x in 1:length(cells_list)) {
  meta_data$unique_sample[which(meta_data$CellId%in%cells_list[[x]])] <-names(cells_list[x])
}

# going to be working with only 10X from now on, so use mtx_10x,
# convert to gene names, and force uniqueness on each rowname. 

rownames(mtx_10x)<-genes_10x$V2[match(rownames(mtx_10x),genes_10x$V1)]
rownames(mtx_10x) <- make.unique(rownames(mtx_10x))

process_10x =function(mtx, meta_data){
  srat=CreateSeuratObject(mtx, meta.data = meta_data)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > minNumGenes & nCount_RNA > minNumUMI & percent.mt < mtCut)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
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

fresh_biopsies =unique(meta_data$unique_sample)[grep("Fresh",unique(meta_data$unique_sample))]
organoids =unique(meta_data$unique_sample)[grep("organoids",unique(meta_data$unique_sample))]
organoids=organoids[1:4]
tumour_only_mtx =mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%fresh_biopsies)]]
tumour_only_meta =meta_data[which(meta_data$unique_sample%in%fresh_biopsies),]

org_only_mtx=mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%organoids)]]
org_only_meta=meta_data[which(meta_data$unique_sample%in%organoids),]

srat_tumour =process_10x(tumour_only_mtx,tumour_only_meta)
srat_org=process_10x(org_only_mtx, org_only_meta)


all_paths = c("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685341/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685342/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843576/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843577/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843578/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733084/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733085/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733086/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004894/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004902/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004910/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787237/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787238/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787239/outs/filtered_gene_bc_matrices/GRCh38/")
names(all_paths) = c("STDY7685340","STDY7685341","STDY7685342","STDY7843576", "STDY7843577", "STDY7843578","STDY7733084",
                     "STDY7733085", "STDY7733086", "STDY8004894",
                     "STDY8004902", "STDY8004910", "STDY7787237", "STDY7787238", "STDY7787239")
inhouse_srats=Read10X(all_paths)

process_inhouse =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 300 & nCount_RNA > 1000 & percent.mt < mtCut)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
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

srat_inhouse=process_inhouse(inhouse_srats)

DimPlot(srat_inhouse)


new_levels=c("SCPs","bridge","right_clust","left_clust", "11", "18", "23", 
                  "5","22","28","0","1", "2", "10", "8", "6", "9", "21",
                  "27","30","7","4","12", "31", "29", "26", "16", "3", "15", "14", "17")
levels(adr_all)=new_levels

dim(adr_all)
dim(srat_tumour)
dim(srat_inhouse)
dim(srat_org)

adr_common_features=subset(adr_all, features=rownames(srat_tumour))
inhouse_common_features=subset(srat_inhouse, features=rownames(srat_tumour))
adr_common_features@meta.data$new_clust=as.character(adr_common_features@meta.data$seurat_clusters)

adr_common_features@meta.data$new_clust[which(adr_common_features@meta.data$new_clust%in%c("25"))]="SCPs"
adr_common_features@meta.data$new_clust[which(adr_common_features@meta.data$new_clust%in%c("13"))]="left_clust"
adr_common_features@meta.data$new_clust[which(adr_common_features@meta.data$new_clust%in%c("24"))]="bridge"
adr_common_features@meta.data$new_clust[which(adr_common_features@meta.data$new_clust%in%c("19", "20"))]="right_clust"
fit = trainModel(adr_common_features@assays$RNA@counts,adr_common_features@meta.data$new_clust,maxCells=5000)

save(fit, file="adr_all_lr_train.RData")


ps_dutch = predictSimilarity(fit, srat_tumour@assays$RNA@counts,srat_tumour@active.ident)
ps_inhouse=predictSimilarity(fit, inhouse_common_features@assays$RNA@counts, inhouse_common_features@active.ident)
ps_org=predictSimilarity(fit, srat_org@assays$RNA@counts, srat_org@active.ident)
dutch_row_ord=c(12,15,23,9, 11,7,0,26,25,16,20,10,14,2,1,8,24,21,19,17,3,5,6,13,4,18,22)
inhouse_row_ord=c(25,12,8,16,5,11,21,22,4,6,3,18,17,15,13,22,9,2,24,7,0,23,19,14,20,10)
print(similarityHeatmap(ps_dutch, column_order=new_levels, row_order=as.character(dutch_row_ord)))
print(similarityHeatmap(ps_inhouse, column_order=new_levels, row_order=as.character(inhouse_row_ord)))
print(similarityHeatmap(ps_org, column_order=new_levels))


saveRDS(srat_inhouse, "srat_inhouse.rds")
saveRDS(srat_tumour, "srat_dutch_primary.rds")



