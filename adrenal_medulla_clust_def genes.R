adr_all_new =adr_all
adr_all_new@meta.data$new_clust= as.character(adr_all_new@meta.data$seurat_clusters)

adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("13"))]="left_clust"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("24"))]="bridge"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("19", "20"))]="right_clust"

tr2mark=quickMarkers(adr_all_new@assays$RNA@counts, adr_all_new@meta.data$new_clust, N=30)
tr2mark_rel=tr2mark[which(tr2mark$cluster%in%c("SCPs", "left_clust", "bridge", "right_clust")),]

tr2mark_rel$gene=make.unique(as.character(tr2mark_rel$gene))

tr2_unique_forpl=tr2mark_rel[which(tr2mark_rel$gene%in%rownames(adr_all_med)),]
tr2_unique_forpl$cluster=factor(tr2_unique_forpl$cluster, levels = c("SCPs", "bridge", "left_clust", "right_clust"))

adr_all_med@meta.data$new_clust= as.character(adr_all_med@meta.data$seurat_clusters)
adr_all_med@meta.data$new_clust[which(adr_all_med@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all_med@meta.data$new_clust[which(adr_all_med@meta.data$new_clust%in%c("13"))]="left_clust"
adr_all_med@meta.data$new_clust[which(adr_all_med@meta.data$new_clust%in%c("24"))]="bridge"
adr_all_med@meta.data$new_clust[which(adr_all_med@meta.data$new_clust%in%c("19", "20"))]="right_clust"

Heatmap(adr_all_med@assays$RNA@scale.data[tr2_unique_forpl$gene,], cluster_rows = F, cluster_columns = F, show_column_names = F, 
        row_split = tr2_unique_forpl$cluster, column_split = factor(adr_all_med@meta.data$new_clust, levels = c("SCPs", "bridge", "left_clust", "right_clust")))

mouse_big=c("SOX10", "HTR3A", "TH", "ERBB3", "DLL3", "CHGB", "FOXD3", "THBD", "FOXQ1", "CARTPT")
which(mouse_big%in%tr2_unique_forpl$gene)

adr_all_med@active.ident

new_idents=adr_all_med@meta.data$new_clust 
names(new_idents)=names(adr_all_med@active.ident)
new_idents=factor(new_idents, levels = c("SCPs", "bridge", "left_clust","right_clust" ))

new_med=adr_all_med

new_med@active.ident=new_idents

new_avg=AverageExpression(new_med, features = rownames(new_med), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)

new_genes=c(tr2_unique_forpl$gene, "HTR3A","CHGB" ,"THBD", "FOXQ1")
new_type=c(as.character(tr2_unique_forpl$cluster), "other_mouse", "other_mouse", "other_mouse", "other_mouse")

new_df=data.frame(genes=as.character(new_genes), clust=new_type)
new_df$clust=factor(new_df$clust, levels = c("SCPs", "bridge", "left_clust","right_clust", "other_mouse"))
col_fun2=colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(n = 5, name = "RdYlBu")))
Heatmap(new_avg@assays$RNA@scale.data[new_df$genes,], col=col_fun2, cluster_rows = F, cluster_columns = F, show_column_names = T, 
        row_split = new_df$clust)
