adr=readRDS('/lustre/scratch117/casm/team274/gk14/fetalProcessing/Results/processedSeurat.RDS')
DimPlot(adr, label=T)

FeaturePlot(adr, "ISL1")

adr@meta.data$annot="x"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(24))]="SCP/Bridge"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(13,23))]="Chromaffin"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(18,20))]="Sympathoblast"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(16,12,22))]="Mesenchyme"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(9,17,28, 0,4,5,11,2,3,7,26,29))]="Cortex"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(6,10,15,21,25))]="Vascular endothelium"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(1,8,19))]="Erythroblasts"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(32,14,27))]="Leukocytes"
adr@meta.data$annot[which(adr@meta.data$seurat_clusters%in%c(30,31,33))]="Other"





DimPlot(adr, group.by = "annot")

saveRDS(adr,'annotSeurat.RDS')

