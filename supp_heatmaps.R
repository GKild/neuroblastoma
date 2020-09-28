adr_all = readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")
library(Seurat)

adr_all@meta.data$new_clust=as.character(adr_all@meta.data$seurat_clusters)
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("13"))]="Chromaffin"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("24"))]="Bridge"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("19", "20"))]="Sympathoblastic"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("3", "15", "14","17"))]="Endothelium"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("23", "11", "18"))]="Mesenchyme"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("22", "5", "28", "0", "1", "2", "10", "6", "8", "9","21"))]="Cortex"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("16", "26"))]="Leukocytes"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("30", "7", "4", "12"))]="Erythroblasts"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("27", "31", "29"))]="Other"


adr_all_supp = adr_all

adr_all_supp@active.ident

adr_names=names(adr_all_supp@active.ident)

adr_all_supp@active.ident=factor(adr_all_supp@meta.data$new_clust)

names(adr_all_supp@active.ident)=adr_names

DimPlot(adr_all_supp)
adr_all_supp=ScaleData(adr_all_supp, features = rownames(adr_all_supp))

cluster.averages <- AverageExpression(adr_all_supp, return.seurat = T)

for_HM=cluster.averages@assays$RNA@scale.data[c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
                                                "CHGB","GAP43", "BCL2", "NPY", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
                                                "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC"),]
col_ord=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme", "Vascular Endothelium",
          "Erythroblasts" ,"Cortex", "Leukocytes", "Other")
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
ht1=Heatmap(for_HM, cluster_rows = F, col=col_fun, column_order = col_ord, name = "Normalized expression", row_names_gp = gpar(fontsize = 8))

levels(adr_all_supp)=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme", "Vascular Endothelium",
                       "Erythroblasts" ,"Cortex", "Leukocytes", "Other")

top_5_all_clust=quickMarkers(adr_all_supp@assays$RNA@counts, adr_all_supp@active.ident, N=5)

for_hm_inm=cluster.averages@assays$RNA@scale.data[as.character(top_5_all_clust$gene),]

ht2=Heatmap(for_hm_inm, col=col_fun, cluster_rows=F, column_order = col_ord, name = "Normalized expression", row_names_gp = gpar(fontsize = 8))

ht_list = ht1 %v% ht2 

draw(ht_list)


lit_genes=c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
  "CHGB","GAP43", "BCL2", "NPY", "PHOX2A", "PHOX2B","MYCN", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
  "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC")

DimPlot(srat_inhouse, group.by = "idents_for_plot")
DimPlot(srat_tumour, group.by = "idents_for_plot")


srat_inhouse2=srat_inhouse
srat_tumour2=srat_tumour


inhouse_names=names(srat_inhouse@active.ident)
dutch_names=names(srat_tumour@active.ident)


srat_inhouse2@active.ident=factor(srat_inhouse2@meta.data$idents_for_plot)
names(srat_inhouse2@active.ident)=inhouse_names

srat_tumour2@active.ident=factor(srat_tumour2@meta.data$idents_for_plot)
names(srat_tumour2@active.ident)=dutch_names

levels(srat_inhouse2)=c("Tumour cluster 1", "Tumour cluster 2",
                        "Tumour cluster 3","Mesenchyme",
                        "Leukocytes","Vascular endothelium")
levels(srat_tumour2)=c("Tumour cluster 1", "Tumour cluster 2","Schwannian stroma",
                       "Vascular endothelium","Mesenchyme", "Leukocytes")
inhouse_avg=AverageExpression(srat_inhouse2, return.seurat = T)

dutch_avg=AverageExpression(srat_tumour2, return.seurat = T)
inhouse_markers=quickMarkers(srat_inhouse2@assays$RNA@counts, srat_inhouse2@active.ident, N=5)
dutch_markers=quickMarkers(srat_tumour2@assays$RNA@counts, srat_tumour2@active.ident, N=5)

im1=Heatmap(inhouse_avg@assays$RNA@scale.data[lit_genes,], col = col_fun, cluster_rows = F, cluster_columns = F,
        row_names_gp = gpar(fontsize = 8), name = "Normalized expression")
im2=Heatmap(inhouse_avg@assays$RNA@scale.data[as.character(inhouse_markers$gene),], col=col_fun,cluster_rows = F,
            cluster_columns = F,  row_names_gp = gpar(fontsize = 8), name = "Normalized expression")

dm1=Heatmap(dutch_avg@assays$RNA@scale.data[lit_genes,], col = col_fun, cluster_rows = F, cluster_columns = F,
            row_names_gp = gpar(fontsize = 8), name = "Normalized expression")
dm2=Heatmap(dutch_avg@assays$RNA@scale.data[as.character(dutch_markers$gene),], col=col_fun,cluster_rows = F,
            cluster_columns = F,  row_names_gp = gpar(fontsize = 8), name = "Normalized expression")

im_list = im1 %v% im2 
dm_list= dm1%v% dm2
draw(im_list)
draw(dm_list)

inhouse_all=quickMarkers(srat_inhouse2@assays$RNA@counts, srat_inhouse2@active.ident, N=Inf, FDR=0.01)
dutch_all=quickMarkers(srat_tumour2@assays$RNA@counts, srat_tumour2@active.ident, N=Inf, FDR=0.01)

write.table(inhouse_all, "inhouse_markers.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dutch_all, "dutch_markers.txt", sep = "\t", quote = F, row.names = F, col.names = T)


quick_with_podos=quickMarkers(merged_cells, cell_labs, N=Inf)
quick_with_podos_rel=quick_with_podos[which(quick_with_podos$cluster%in%c("SCPs","Bridge","Sympathoblastic", "Chromaffin","Podocytes")),]
podo_marks_filt=dplyr::filter(quick_with_podos_rel, tfidf>1 & geneFrequencySecondBest <0.2)

write.table(podo_marks_filt, "marks_with_podos.txt", sep = "\t", row.names = F, col.names = T, quote = F)
