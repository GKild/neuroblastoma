genes_10x=read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/scData/lateAdrenals/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz")


adr_all@meta.data$orig_ident2=as.character(adr_all@meta.data$orig.ident)
adr_all@meta.data$orig_ident2[which(adr_all@meta.data$orig_ident2=="babyAdrenal1")]="WSSS8012017"
adr_all@meta.data$orig_ident2[which(adr_all@meta.data$orig_ident2=="babyAdrenal2")]="WSSS8011223"
adr_all@meta.data$orig_ident2[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)=="Adr8710632")]="WSSS8710632"
adr_all@meta.data$orig_ident2[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)=="Adr8710633")]="WSSS8710633"
adr_all@meta.data$orig_ident2[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)=="Adr8710634")]="WSSS8710634"
adr_all@meta.data$orig_ident2[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)=="Adr8710635")]="WSSS8710635"

adr_all@meta.data$UMAP_1=Embeddings(adr_all, reduction =  'umap')[,1]
adr_all@meta.data$UMAP_2=Embeddings(adr_all, reduction =  'umap')[,2]

adr_all_cm=adr_all@assays$RNA@counts

colnames(adr_all_cm)=paste0(adr_all@meta.data$orig_ident2, ":", adr_all@meta.data$barcode, ":", adr_all@meta.data$new_clust)
rownames(adr_all_cm)=paste0(as.character(genes_10x$V2), ":", as.character(genes_10x$V1))

writeMM(adr_all_cm, "adrenalTableOfCounts.mtx")

rowLabels=data.frame(RowNumber=1:nrow(adr_all),
                     GeneLabel=paste0(as.character(genes_10x$V2), ":", as.character(genes_10x$V1)),
                     Symbol=as.character(genes_10x$V2), 
                     EnsemblID=as.character(genes_10x$V1))
adr_all@meta.data$barcode=gsub("^.*_", "", colnames(adr_all))

colLabels=data.frame(ColNumber=1:ncol(adr_all),
                     DropletID=paste0(adr_all@meta.data$orig_ident2, ":", adr_all@meta.data$barcode, ":", adr_all@meta.data$new_clust),
                     SangerID=adr_all@meta.data$orig_ident2, Barcode=adr_all@meta.data$barcode, Annotation=adr_all@meta.data$new_clust)

write.table(colLabels, "adrenalTableOfCounts_colLabels.tsv", row.names = F, col.names = T, sep = "\t", quote = F)                    
write.table(rowLabels, "adrenalTableOfCounts_rowLabels.tsv", row.names = F, col.names = T, sep = "\t", quote = F)


adr_all@meta.data$MTfrac=adr_all@meta.data$mtGenes
adr_all@meta.data$nUMI=adr_all@meta.data$nCount_RNA
adr_all@meta.data$nGenes=adr_all@meta.data$nFeature_RNA


adr_all@meta.data$sample_name="x"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]="w21_1"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)%in%c("Adr8710634", "Adr8710635"))]="w21_2"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident=="babyAdrenal1")]="w8"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident=="babyAdrenal2")]="w8d6"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5388STDY7717452","5388STDY7717453",
                                                                    "5388STDY7717454","5388STDY7717455"))]="w10d5_1"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5388STDY7717456","5388STDY7717458",
                                                                    "5388STDY7717459"))]="w10d5_2"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5698STDY7839907","5698STDY7839909",
                                                                    "5698STDY7839917"))]="w11"

adr_cell_manifest=data.frame(barcode=adr_all@meta.data$barcode, SangerID=adr_all@meta.data$orig_ident2, Annotation=adr_all@meta.data$new_clust,
                             DropletID=paste0(adr_all@meta.data$orig_ident2, ":", adr_all@meta.data$barcode, ":", adr_all@meta.data$new_clust),
                             nUMI=adr_all@meta.data$nUMI, nGenes=adr_all@meta.data$nGenes, MTfrac=adr_all@meta.data$MTfrac, UMAP_1=adr_all@meta.data$UMAP_1,
                             UMAP_2=adr_all@meta.data$UMAP_2, Gestation_weeks=adr_all@meta.data$sample_name)

write.table(adr_cell_manifest, "adr_cell_manifest.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


srat_inhouse=readRDS("srat_inhouse.rds")
srat_tumour=readRDS("srat_dutch.rds")


srat_tumour@meta.data$new_idents="x"
srat_tumour@meta.data$idents_for_plot="x"

srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(22,8))]="tumour"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(13))]="schwannian stroma"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(0,2,11,14,15,18,25,10,12,20))]="mesenchyme"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(1,17,19,21,16,5,9,4,6,26,23,3,24))]="immune"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(7))]="vascular endothelium"

srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(8))]="Tumour cluster 1"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(22))]="Tumour cluster 2"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(13))]="Schwannian stroma"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(0,2,11,14,15,18,25,10,12,20))]="Mesenchyme"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(1,17,19,21,16,5,9,4,6,26,23,3,24))]="Leukocytes"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(7))]="Endothelium"


srat_inhouse@meta.data$new_idents="x"
srat_inhouse@meta.data$idents_for_plot="x"

srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))]="tumour"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                    17,15,18,3,6,4,21,22,24))]="immune"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="vascular endothelium"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="mesenchyme"

srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25))]="Tumour cluster 1"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(12))]="Tumour cluster 2"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(5))]="Tumour cluster 3"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                         17,15,18,3,6,4,21,22,24))]="Leukocytes"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="Endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="Mesenchyme"

length(which(genes_10x$V2%in%rownames(srat_tumour)))

paste0(rownames(srat_tumour), ":",as.character(genes_10x$V1)[match(rownames(srat_tumour), as.character(genes_10x$V2))])

CELseq2rowLabels=data.frame(RowNumber=1:nrow(srat_tumour),
                     GeneLabel=paste0(rownames(srat_tumour), ":",as.character(genes_10x$V1)[match(rownames(srat_tumour), as.character(genes_10x$V2))]),
                     Symbol=rownames(srat_tumour), 
                     EnsemblID=as.character(genes_10x$V1)[match(rownames(srat_tumour), as.character(genes_10x$V2))])

srat_tumour@meta.data$


CELseq2colLabels=data.frame(ColNumber=1:ncol(srat_tumour),DropletID=paste0(srat_tumour@meta.data$CellId,":",srat_tumour@meta.data$idents_for_plot),
                            CellID=srat_tumour@meta.data$CellId, Annotation=srat_tumour@meta.data$idents_for_plot)

CELseq2_cell_manifest=data.frame(CellID=srat_tumour@meta.data$CellId, Annotation=srat_tumour@meta.data$idents_for_plot,
                                 DropletID=paste0(srat_tumour@meta.data$CellId,":",srat_tumour@meta.data$idents_for_plot),
                                 nUMI=srat_tumour@meta.data$nCount_RNA, nGenes=srat_tumour@meta.data$nFeature_RNA, 
                                 MTfrac=srat_tumour@meta.data$percent.mt/100, UMAP_1=Embeddings(srat_tumour)[,1],
                                 UMAP_1=Embeddings(srat_tumour)[,2], Tumour_sample=srat_tumour@meta.data$unique_sample)

writeMM(srat_tumour@assays$RNA@counts, "CELseq2TableOfCounts.mtx")
write.table(CELseq2_cell_manifest, "CELseq2_cell_manifest.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(CELseq2colLabels, "CELseq2_colLabels.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(CELseq2rowLabels, "CELseq2_rowLabels.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

writeMM(srat_inhouse@assays$RNA@counts, "tum_10x_TableOfCounts.mtx")

srat_inhouse@meta.data$barcode=sapply(strsplit(colnames(srat_inhouse), "_"), "[", 2)

tum_10x_colLabels=data.frame(ColNumber=1:ncol(srat_inhouse),
                     DropletID=paste0(srat_inhouse@meta.data$orig.ident, ":", srat_inhouse@meta.data$barcode, ":", srat_inhouse@meta.data$idents_for_plot),
                     SangerID=srat_inhouse@meta.data$orig.ident, Barcode=srat_inhouse@meta.data$barcode, Annotation=srat_inhouse@meta.data$idents_for_plot)

tum_10x_rowLabels=data.frame(RowNumber=1:nrow(srat_inhouse),
                     GeneLabel=paste0(as.character(genes_10x$V2[-which(genes_10x$V2%in%exclude_genes$V1)]), ":", as.character(genes_10x$V1[-which(genes_10x$V2%in%exclude_genes$V1)])),
                     Symbol=as.character(genes_10x$V2[-which(genes_10x$V2%in%exclude_genes$V1)]), 
                     EnsemblID=as.character(genes_10x$V1[-which(genes_10x$V2%in%exclude_genes$V1)]))

dim(srat_inhouse)
which(sapply(strsplit(rownames(srat_inhouse)[which(!rownames(srat_inhouse)%in%genes_10x$V2)], "\\."), '[', 1)%in%genes_10x$V2)


exclude_genes=read.table("useful_files/excludeGenes.tsv")
dim(exclude_genes)
exclude_genes$V1[which(!exclude_genes$V1%in%genes_10x$V2)]

genes_10x=read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv")
as.character(genes_10x$V2[-which(genes_10x$V2%in%exclude_genes$V1)])

tum_10x_cell_manifest=data.frame(barcode=srat_inhouse@meta.data$barcode, SangerID=srat_inhouse@meta.data$orig.ident, Annotation=srat_inhouse@meta.data$idents_for_plot,
                                 DropletID=paste0(srat_inhouse@meta.data$orig.ident, ":", srat_inhouse@meta.data$barcode, ":", srat_inhouse@meta.data$idents_for_plot),
                                 nUMI=srat_inhouse@meta.data$nCount_RNA, nGenes=srat_inhouse@meta.data$nFeature_RNA, 
                                 MTfrac=srat_inhouse@meta.data$percent.mt/100, UMAP_1=Embeddings(srat_inhouse)[,1],
                                 UMAP_1=Embeddings(srat_inhouse)[,2], Tumour_sample=srat_inhouse@meta.data$GOSH_ID)

write.table(tum_10x_cell_manifest, "tum_10x_cell_manifest.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(tum_10x_colLabels, "tum_10x_colLabels.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(tum_10x_rowLabels, "tum_10x_rowLabels.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


sapply(unique(adr_all@meta.data$sample_name), function(x){
  table(adr_all@meta.data$new_clust[which(adr_all@meta.data$sample_name==x)])
})
3