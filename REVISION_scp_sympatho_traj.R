just_bridge=subset(adr_all,
                cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$new_clust%in%c("Bridge", "Sympathoblastic"))])


counts.data <-GetAssayData(object = just_bridge, slot = "counts", assay="RNA")
pheno.data <- as.data.frame(just_bridge@meta.data, row.names = rownames(just_bridge@meta.data))
feature.data <- as.data.frame(just_bridge@assays$RNA@counts@Dimnames[[1]], row.names = just_bridge@assays$RNA@counts@Dimnames[[1]])


cds <- new_cell_data_set(counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- cds@colData$seurat_clusters
names(list_cluster) <- just_bridge@assays[["RNA"]]@data@Dimnames[[2]]


cds@clusters@listData[["UMAP"]][["clusters"]] <- as.factor(list_cluster)


cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds@reducedDims@listData[["UMAP"]] <-Embeddings(just_bridge, reduction =  'umap')
cells_embeddings=Embeddings(just_bridge, reduction =  'umap')
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@reducedDims$UMAP) <- NULL
cds <- learn_graph(cds, use_partition = F,learn_graph_control = list(minimal_branch_len = 8,
                                                                     geodesic_distance_ratio = 3))
plot_cells(cds, cell_size = 1.5) +ggtitle("Fetal adrenal SCP trajectory")
cds <- order_cells(cds)
pseudo_srat = pseudotime(cds, reduction_method = 'UMAP')
just_bridge@meta.data$pseudotime=pseudo_srat

FeaturePlot(just_bridge, features = "pseudotime", pt.size = 2) + scale_color_viridis(option="plasma", direction = 1) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank()) +ylab("UMAP 2") +xlab("UMAP 1")
rowData(cds)$gene_short_name=rownames(just_bridge)

pr_test_res <- graph_test(cds, neighbor_graph="principal_graph")
pr_deg_ids <- subset(pr_test_res, q_value < 0.001)
pr_deg_ids_ord=pr_deg_ids[order(pr_deg_ids$morans_I, decreasing = T),]


#make the "left" heatmap
m_traj = match(rownames(pr_deg_ids_ord),rownames(cds))[seq(100)]

pt_group = data.frame(cells = names(pseudotime(cds)),
                           group = cut_number(pseudotime(cds),100))
dat = aggregate_gene_expression(cds[m_traj,],cell_group_df = pt_group,scale_agg_value=FALSE)

scaled_dat = t(scale(t(dat)))
scaled_dat[abs(scaled_dat)>3] = sign(scaled_dat[abs(scaled_dat)>3]) *3

pheatmap(scaled_dat,cluster_cols=FALSE)



