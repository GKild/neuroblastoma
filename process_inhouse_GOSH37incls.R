
source("/home/jovyan/Neuroblastoma/gene_sets.R")
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
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787239/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch117/casm/team274/my4/oldScratch/HCA/Neuroblastoma3/cellranger302_count_33883_CG_SB_NB8759214_GRCh38-1_2_0/filtered_feature_bc_matrix/", 
              "/lustre/scratch117/casm/team274/my4/oldScratch/HCA/Neuroblastoma3/cellranger302_count_33883_CG_SB_NB8759215_GRCh38-1_2_0/filtered_feature_bc_matrix/", 
              "/lustre/scratch117/casm/team274/my4/oldScratch/HCA/Neuroblastoma3/cellranger302_count_33883_CG_SB_NB8759216_GRCh38-1_2_0/filtered_feature_bc_matrix/")
names(all_paths) = c("STDY7685340","STDY7685341","STDY7685342",
                     "STDY7843576", "STDY7843577", "STDY7843578",
                     "STDY7733084","STDY7733085", "STDY7733086",
                     "STDY8004894","STDY8004902", "STDY8004910",
                     "STDY7787237", "STDY7787238", "STDY7787239", 
                     "NB8759214", "NB8759215", "NB8759216")


inhouse_srats=Read10X(all_paths)

process_inhouse =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 300 & nCount_RNA > 1000 & percent.mt < 30)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("/home/jovyan/Neuroblastoma/useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
  keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 75)
  srat = FindNeighbors(srat, dims=1:75)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:75, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


srat_inhouse=process_inhouse(inhouse_srats)
lit_genes=c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
            "CHGB","GAP43", "BCL2", "NPY", "PHOX2A", "PHOX2B","MYCN", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
            "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC")
FeaturePlot(srat_inhouse, c("CHGB","GAP43", "BCL2", "NPY", "PHOX2A", "PHOX2B","MYCN"))

srat_inhouse@meta.data$GOSH_ID=as.character(srat_inhouse@meta.data$orig.ident)
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="GOSH014"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="GOSH023"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="GOSH019"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="GOSH025"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="GOSH021"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("NB8759214", "NB8759215", "NB8759216"))]="GOSH037"
DimPlot(srat_inhouse, group.by = "Phase")

load("/home/jovyan/Neuroblastoma/adr_all_podo_train.RData")

inh_adr_ps=predictSimilarity(fit, srat_inhouse@assays$RNA@counts[fit$podocytes$glmnet.fit$beta@Dimnames[[1]],], srat_inhouse@meta.data$seurat_clusters)
length(which(fit$podocytes$glmnet.fit$beta@Dimnames[[1]]%in%rownames(srat_inhouse)))

similarityHeatmap(inh_adr_ps)

colnames(inh_adr_ps)[5]="chromaffin"
colnames(inh_adr_ps)[7]="sympathoblastic"






