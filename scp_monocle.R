library(monocle3)
adr=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")

scp_traj=subset(adr, idents=c(10,17,11,7,19))
nocc_cells=rownames(scp_traj@meta.data)[which(scp_traj@meta.data$Phase=="G1")]
scp_traj_nocc=subset(scp_traj, cells=nocc_cells)
DimPlot(scp_traj, label=T)

counts.data <-GetAssayData(object = scp_traj, slot = "counts", assay="RNA")
pheno.data <- as.data.frame(scp_traj@meta.data, row.names = rownames(scp_traj@meta.data))#babyadr_meta[babyadr_scp,], row.names = rownames(babyadr_meta[babyadr_scp,]))
feature.data <- as.data.frame(scp_traj@assays$RNA@counts@Dimnames[[1]], row.names = scp_traj@assays$RNA@counts@Dimnames[[1]])


cds <- new_cell_data_set(counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- cds@colData$seurat_clusters
names(list_cluster) <- scp_traj@assays[["RNA"]]@data@Dimnames[[2]]


cds@clusters@listData[["UMAP"]][["clusters"]] <- as.factor(list_cluster)


cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds@reducedDims@listData[["UMAP"]] <-Embeddings(scp_traj, reduction =  'umap')
cells_embeddings=Embeddings(scp_traj, reduction =  'umap')
cds <- learn_graph(cds,use_partition = F, learn_graph_control = list(minimal_branch_len = 10))
cds <- order_cells(cds, root_cells = "WSSS8011223_TAAGTGCAGACATAAC")
pseudo_srat = pseudotime(cds, reduction_method = 'UMAP')
length(pseudo_srat)
plot_cells(cds, color_cells_by = 'pseudotime', cell_size = 1.5) +ggtitle("Fetal adrenal SCP trajectory")
colData(cds)$pseudot=pseudo_srat
rowData(cds)$gene_short_name=rownames(scp_traj)
tf_table=read.table("useful_files/Homo_sapiens_TF.txt", sep = "\t", header = T)
tfs=as.character(tf_table$Symbol)
cds_subset = cds[rowData(cds)$gene_short_name %in% tfs,]

scp_graph_test = graph_test(cds_subset, neighbor_graph="principal_graph")
pr_deg_ids =subset(scp_graph_test, q_value < 0.01)
pr_deg_ids = pr_deg_ids[order(pr_deg_ids$q_value),]
pr_deg_ids_nocc=pr_deg_ids[which(!rownames(pr_deg_ids)%in%c("HMGB1", "FOXM1", "CENPA", "MXD3", "YBX1", "HMGB2", "HMGB3","MIS18BP1", "E2F1", "MYBL2", "TCF19")),]
m = match(rownames(pr_deg_ids),rownames(cds_subset))[seq(20)]
pt = pseudotime(cds_subset)
pt_group = data.frame(cells = names(pt),
                      group = cut_number(pt,100))
dat = aggregate_gene_expression(cds_subset[m,],cell_group_df = pt_group,scale_agg_value=FALSE)
#Z-transform rows
dat = t(scale(t(dat)))

dat[abs(dat)>3] = sign(dat[abs(dat)>3]) *3
pheatmap(dat,cluster_cols=FALSE)


m = match(rownames(pr_deg_ids_nocc),rownames(cds_subset))[seq(63)]
pt = pseudotime(cds_subset)
pt_group = data.frame(cells = names(pt),
                      group = cut_number(pt,100))
dat = aggregate_gene_expression(cds_subset[m,],cell_group_df = pt_group,scale_agg_value=FALSE)
#Z-transform rows
dat = t(scale(t(dat)))

dat[abs(dat)>3] = sign(dat[abs(dat)>3]) *3
pheatmap(dat,cluster_cols=FALSE)
col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
r_ord=c("CREB5", "ZFP36L1", "TFAP2A", "SOX2", "HMGA2", "ID1", "SNAI1",
  "ASCL1","TBX2", "TEAD2", "ID3", "ID4", "SOX10","TLX2", "GATA3", "HAND2", "PHOX2A", "GATA2", "EPAS1", "NR4A1")

scp_tfs=c("SFPQ","PHOX2B","ID1","ID2","ASCL1","TBX2","DLX2","DLX1","DLX5","SNAI1","TFDP2","TCF12","SOX4","TEAD2",
  "ID3","ID4","SOX10","SOX6","HOXB9","SOX5","TGIF1","HEYL","MAFF","ETS1","DNAJC1","RXRG","HES1","CREB5",
  "MEF2C","ZFP36L1","TFAP2A","SOX2","HMGA2","ZEB2","GLI3","REST","AHR")
rev_scp_tfs=rev(scp_tfs)
all_tfs=c(rev_scp_tfs, "TLX2","LBX2","FOXS1" ,
          "TCF4", "HAND2","GATA3","TBX20","ZNF467",
          "ETS2","RARA","FOSL2","NR4A1","NPAS4","ZNF331" ,"HEY1",
          "FOS","JUNB","PROX1","EPAS1","GATA2","PHOX2A","INSM1",
          "HAND1","ISL1","POU2F2","FOXN4")
lgd=Legend(col_fun = col_fun, title = "Normalized gene expression")
hm=Heatmap(dat, cluster_columns = F, show_column_names = F, show_row_dend = F, col=col_fun, row_order = all_tfs, show_row_names = F, 
           name="Normalized gene expression", heatmap_legend_param = list(direction="horizontal",heatmap_legend_side = "bottom"))

draw(hm)



DimPlot(scp_traj)
pdf("scp_traj.pdf", height = 5, width=7,useDingbats = F)
gg=FeaturePlot(scp_traj, features = "pseudotime", pt.size = 2) + scale_color_viridis(option="plasma", direction = 1) +
  ggtitle("") + theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none") +ylab("UMAP 2") +xlab("UMAP 1")

plot(gg)
dev.off()
pdf("scp_clusters.pdf", height = 5, width=7,useDingbats = F)
gg=DimPlot(scp_traj, pt.size = 2, cols = viridis(5)[c(4,1,3,2,5)]) +
  ggtitle("") + theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
