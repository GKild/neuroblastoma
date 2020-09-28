
library(Seurat)
source("logisticRegression.R")
source("gene_sets.R")
adr_all <- readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
fetal_pods=readRDS("/lustre/scratch117/casm/team274/gk14/fetal_podocytes_forGerda.RDS")
genes_10x <- read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)
load("common_genes.rdata")

fetal_pds_counts=fetal_pods@assays$RNA@counts
fetal_pds_counts=fetal_pds_counts[which(rownames(fetal_pds_counts)%in%genes_10x$V1),]

rownames(fetal_pds_counts)=as.character(genes_10x$V2)[match(rownames(fetal_pds_counts), genes_10x$V1)]
rownames(fetal_pds_counts)=make.unique(rownames(fetal_pds_counts))



adr_all@meta.data$new_clust=as.character(adr_all@meta.data$seurat_clusters)

adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("13"))]="Chromaffin"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("24"))]="Bridge"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("19", "20"))]="Sympathoblastic"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("3", "15", "14","17"))]="Vascular_endothelium"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("23", "11", "18"))]="Mesenchyme"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("22", "5", "28", "0", "1", "2", "10", "6", "8", "9","21"))]="Cortex"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("16", "26"))]="Leukocytes"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("30", "7", "4", "12"))]="Erythroblasts"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("27", "31", "29"))]="Other"

comm_genes=intersect(rownames(fetal_pds_counts), obj)

pods_comm_genes=fetal_pds_counts[comm_genes, ]
adr_comm_genes=adr_all@assays$RNA@counts[comm_genes,]

merged_cells=cbind(adr_comm_genes, pods_comm_genes)
cell_labs=c(adr_all@meta.data$new_clust, rep("Podocytes", 278))

fit = trainModel(merged_cells, cell_labs)

save(fit, file="adr_all_podo_train.RData")

load("adr_all_podo_train.RData")


ps_dutch = predictSimilarity(fit, srat_tumour@assays$RNA@counts[comm_genes,],srat_tumour@meta.data$idents_for_plot)
ps_inhouse=predictSimilarity(fit, srat_inhouse@assays$RNA@counts[comm_genes,], srat_inhouse@meta.data$idents_for_plot)
ps_org=predictSimilarity(fit, srat_org@assays$RNA@counts[comm_genes,], srat_org@meta.data$idents_for_plot)
ps_self=predictSimilarity(fit, merged_cells, cell_labs)



lr_col_ord=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme","Cortex","Leukocytes", "Vascular_endothelium", "Erythroblasts","Other", "podocytes")
pdf("dutch_lr.pdf", height = 4, width = 7)
print(similarityHeatmap(ps_dutch, column_order=lr_col_ord,row_order=as.character(c("Tumour cluster 1", "Tumour cluster 2",
                                                                              "Schwannian stroma",
                                                                              "Mesenchyme", "Leukocytes", "Vascular endothelium"))))
dev.off()

pdf("inhouse_lr.pdf", height = 4, width = 7)
print(similarityHeatmap(ps_inhouse, column_order=lr_col_ord, row_order=c("Tumour cluster 1", "Tumour cluster 2",
                                                                   "Tumour cluster 3","Mesenchyme",
                                                                   "Leukocytes", "Vascular endothelium")))
dev.off()
pdf("org_lr.pdf", height = 4, width = 7)
print(similarityHeatmap(ps_org, column_order=lr_col_ord, row_order=c("tumour cluster 1", "tumour cluster 2",
                                                                   "tumour cluster 3", "tumour cluster 4",
                                                               "tumour cluster 5", "tumour cluster 6",
                                                               "tumour cluster 7", "tumour cluster 8",
                                                               "tumour cluster 9", "tumour cluster 10",
                                                               "tumour cluster 11", "tumour cluster 12",
                                                               "mesenchyme","other")))
dev.off()

pdf("self_lr.pdf", height = 4, width = 7)
similarityHeatmap(ps_self, column_order=lr_col_ord, row_order=lr_col_ord)
dev.off()

quick_with_podos=quickMarkers(merged_cells, cell_labs, N=Inf)
quick_with_podos_rel=quick_with_podos[which(quick_with_podos$cluster%in%c("SCPs","bridge","left_clust", "right_clust","pods")),]
podo_marks_filt=dplyr::filter(quick_with_podos_rel, tfidf>1 & geneFrequencySecondBest <0.2)

similarityHeatmap(ps_self, row_order=lr_col_ord, column_order=lr_col_ord)
similarityHeatmap(ps_inhouse,column_order=lr_col_ord)

dev.off()