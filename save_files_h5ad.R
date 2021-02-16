library(sceasy)
library(reticulate)
use_condaenv('sceasy_recs', required = T)
loompy <- reticulate::import('loompy')

adr_all=readRDS("/lustre/scratch117/casm/team274/gk14/adr_all_annot.rds")


keep_meta=adr_all@meta.data


keep_meta=keep_meta[,c("sample_name","new_clust", "nCount_RNA", "nFeature_RNA", "mtGenes", "hspGenes", "riboGenes")]

keep_meta$gestational_age=keep_meta$sample_name

keep_meta$gestational_age=sapply(strsplit(keep_meta$gestational_age, "_"), "[", 1)
keep_meta=keep_meta[,c("gestational_age","sample_name","new_clust", "nCount_RNA", "nFeature_RNA", "mtGenes", "hspGenes", "riboGenes")]

colnames(keep_meta)=c("GestationalAge","SampleName","Annotation", "nCount_RNA", "nFeature_RNA", "mtGenes", "hspGenes", "riboGenes")


sceasy::convertFormat(adr_all, from="seurat", to="anndata",
                      outFile='/lustre/scratch117/casm/team274/gk14/adr_all_cellxgene.h5ad')
saveRDS(adr_all, "/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds")

srat_inhouse=readRDS('/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_inhouse.rds')

srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25))]="Tumour cluster 1"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(12))]="Tumour cluster 2"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(5))]="Tumour cluster 3"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                         17,15,18,3,6,4,21,22,24))]="Leukocytes"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="Endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="Mesenchyme"
srat_inhouse@meta.data$idents_for_plot=factor(srat_inhouse@meta.data$idents_for_plot, levels=c("Tumour cluster 1", "Tumour cluster 2",
                                                                                               "Tumour cluster 3","Mesenchyme",
                                                                                               "Leukocytes","Endothelium"))

srat_dutch=readRDS('/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_dutch.rds')
DimPlot(srat_dutch, group.by = "idents_for_plot")

srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(13))]="Tumour cluster 1"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(24))]="Tumour cluster 2"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(15))]="Ambiguous 1"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(10))]="Tumour cluster 3"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(17))]="Ambiguous 2"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(6))]="Schwannian stroma"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(7,29))]="Endothelium"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(28,18,23,5,26,0,8,20,9,3,22,25,4,19,21,27))]="Leukocytes"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(1,2,11,12,14,16))]="Mesenchyme"
srat_dutch@meta.data$idents_for_plot=factor(srat_dutch@meta.data$idents_for_plot, levels=c(
  "Tumour cluster 1","Tumour cluster 2","Tumour cluster 3","Ambiguous 1", "Ambiguous 2",
  "Schwannian stroma","Endothelium", "Mesenchyme", "Leukocytes"
))


saveRDS(srat_dutch, '/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_dutch.rds')
saveRDS(srat_inhouse, '/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_inhouse.rds')



srat_inhouse@meta.data$PD_ID=as.character(srat_inhouse@meta.data$orig.ident)

srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="PD42184"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="PD42752-1"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="PD42752-2"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="PD46693"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="PD43255"

keep_meta=srat_inhouse@meta.data
keep_meta=keep_meta[,c("PD_ID","idents_for_plot", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hspGenes", "percent.riboGenes")]

colnames(keep_meta)=c("SampleName","Annotation", "nCount_RNA", "nFeature_RNA", "mtGenes", "hspGenes", "riboGenes")

srat_inhouse@meta.data=keep_meta

sceasy::convertFormat(srat_inhouse, from="seurat", to="anndata",
                      outFile='/lustre/scratch117/casm/team274/gk14/cellxgene/nb_GOSH_cellxgene.h5ad')
saveRDS(srat_inhouse, "/lustre/scratch117/casm/team274/gk14/cellxgene/nb_GOSH.rds")

keep_meta=srat_dutch@meta.data

keep_meta=keep_meta[,c("unique_sample","idents_for_plot", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hspGenes", "percent.riboGenes")]
keep_meta$unique_sample=sapply(strsplit(keep_meta$unique_sample, "_"), "[", 1)

colnames(keep_meta)=c("SampleName","Annotation", "nCount_RNA", "nFeature_RNA", "mtGenes", "hspGenes", "riboGenes")

srat_dutch@meta.data=keep_meta

sceasy::convertFormat(srat_dutch, from="seurat", to="anndata",
                      outFile='/lustre/scratch117/casm/team274/gk14/cellxgene/nb_PMC_cellxgene.h5ad')
saveRDS(srat_dutch, "/lustre/scratch117/casm/team274/gk14/cellxgene/nb_PMC.rds")