srat_inhouse=readRDS("srat_inhouse.rds")
srat_org=readRDS("srat_org.rds")
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
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(7))]="Vascular endothelium"




srat_org@meta.data$new_idents="x"
srat_org@meta.data$idents_for_plot="x"

srat_org@meta.data$new_idents[which(srat_org@meta.data$seurat_clusters%in%c(0,1,2,3,4,6,7,8,9,10,11,13))]="tumour"
srat_org@meta.data$new_idents[which(srat_org@meta.data$seurat_clusters%in%c(12))]="mesenchyme"
srat_org@meta.data$new_idents[which(srat_org@meta.data$seurat_clusters%in%c(5))]="other"

srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(0))]="tumour cluster 1"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(1))]="tumour cluster 2"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(2))]="tumour cluster 3"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(3))]="tumour cluster 4"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(4))]="tumour cluster 5"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(6))]="tumour cluster 6"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(7))]="tumour cluster 7"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(8))]="tumour cluster 8"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(9))]="tumour cluster 9"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(10))]="tumour cluster 10"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(11))]="tumour cluster 11"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(13))]="tumour cluster 12"

srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(12))]="Mesenchyme"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$seurat_clusters%in%c(5))]="Other"




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
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="Vascular endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="Mesenchyme"



srat_tumour@meta.data$idents_for_plot=factor(srat_tumour@meta.data$idents_for_plot,
                                             levels=c("Tumour cluster 1", "Tumour cluster 2","Schwannian stroma",
                                                      "Vascular endothelium","Mesenchyme", "Leukocytes"))
srat_org@meta.data$idents_for_plot=factor(srat_org@meta.data$idents_for_plot,
                                          levels=c("tumour cluster 1", "tumour cluster 2",
                                                   "tumour cluster 3", "tumour cluster 4",
                                                   "tumour cluster 5", "tumour cluster 6",
                                                   "tumour cluster 7", "tumour cluster 8",
                                                   "tumour cluster 9", "tumour cluster 10",
                                                   "tumour cluster 11", "tumour cluster 12",
                                                   "mesenchyme", "other"))
srat_inhouse@meta.data$idents_for_plot=factor(srat_inhouse@meta.data$idents_for_plot,
                                              levels=c("Tumour cluster 1", "Tumour cluster 2",
                                                       "Tumour cluster 3","Mesenchyme",
                                                       "Leukocytes","Vascular endothelium"))

colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))

pdf("inhouse_umap.pdf", height = 3, width=4,useDingbats = F)
gg=DimPlot(srat_inhouse, group.by = "idents_for_plot", label=T,
           cols = c(colfunc(3),"#7E481C", "#DDCC77", "#CC6677"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
pdf("org_umap.pdf", height = 7, width=7,useDingbats = F)
gg=DimPlot(srat_org, group.by = "idents_for_plot",
           cols = c(colfunc(12), "#7E481C", "grey"), label=T,
           pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank() ,legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
pdf("dutch_umap.pdf", height = 3, width=4,useDingbats = F)
gg=DimPlot(srat_tumour, group.by = "idents_for_plot",label=T,
           cols = c("#084594", "#4292C6","#7c5295", "#CC6677", "#7E481C","#DDCC77"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), 
                                 legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")


plot(gg)
dev.off()


