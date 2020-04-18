adr_all=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")

pdf("adr_all_clusters.pdf", height = 5, width=7,useDingbats = F)
gg=DimPlot(adr_all, pt.size = 1.2, label = T, label.size = 5) +
  ggtitle("") + theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
adr_all_top10=quickMarkers(adr_all@assays$RNA@counts, adr_all@active.ident, N=10)
adr_all_med=subset(adr_all, idents = c(13,19,20,24,25))

DimPlot(adr_all_med)


counts.data <-GetAssayData(object = adr_all_med, slot = "counts", assay="RNA")
pheno.data <- as.data.frame(adr_all_med@meta.data, row.names = rownames(adr_all_med@meta.data))#babyadr_meta[babyadr_scp,], row.names = rownames(babyadr_meta[babyadr_scp,]))
feature.data <- as.data.frame(adr_all_med@assays$RNA@counts@Dimnames[[1]], row.names = adr_all_med@assays$RNA@counts@Dimnames[[1]])


cds <- new_cell_data_set(counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- cds@colData$seurat_clusters
names(list_cluster) <- adr_all_med@assays[["RNA"]]@data@Dimnames[[2]]


cds@clusters@listData[["UMAP"]][["clusters"]] <- as.factor(list_cluster)


cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds@reducedDims@listData[["UMAP"]] <-Embeddings(adr_all_med, reduction =  'umap')
cells_embeddings=Embeddings(adr_all_med, reduction =  'umap')
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@reducedDims$UMAP) <- NULL
cds <- learn_graph(cds,use_partition = F)
cds <- order_cells(cds)
pseudo_srat = pseudotime(cds, reduction_method = 'UMAP')
length(pseudo_srat)
plot_cells(cds, color_cells_by = 'pseudotime', cell_size = 1.5) +ggtitle("Fetal adrenal SCP trajectory")
FeaturePlot(adr_all_med, features = "pseudotime")
adr_all_med@meta.data$pseudotime=pseudo_srat
pdf("all_medulla_pseudotime.pdf", height = 5, width=7,useDingbats = F)
gg=FeaturePlot(adr_all_med, features = "pseudotime", pt.size = 2) + scale_color_viridis(option="plasma", direction = 1) +
  ggtitle("") + theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank()) +ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()

pdf("all_medulla_clusters.pdf", height = 5, width=7,useDingbats = F)
gg=DimPlot(adr_all_med, pt.size = 2, label = T, label.size = 5, cols = viridis(5)[c(4,1,3,2,5)]) +
  ggtitle("") + theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
#do DE along the trajectory for TFs
rowData(cds)$gene_short_name=rownames(adr_all_med)
tf_table=read.table("useful_files/Homo_sapiens_TF.txt", sep = "\t", header = T)
tfs=as.character(tf_table$Symbol)
cds_subset = cds[rowData(cds)$gene_short_name %in% tfs,]
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph")
pr_deg_ids <- subset(subset_pr_test_res, q_value < 0.05)
pr_deg_ids = pr_deg_ids[order(pr_deg_ids$q_value),]
# subset cds to cells in "left" and "right" trajectory

left_cells=WhichCells(adr_all_med, idents = c(25,24,13))
right_cells=WhichCells(adr_all_med, idents = c(25,24,19,20))

left_cds=cds_subset[,left_cells]
right_cds=cds_subset[,right_cells]

left_pseudot=pseudotime(left_cds)
right_pseudot=pseudotime(right_cds)
#make the "left" heatmap
m_left = match(rownames(pr_deg_ids),rownames(left_cds))[seq(100)]

pt_group_left = data.frame(cells = names(left_pseudot),
                      group = cut_number(left_pseudot,100))
dat_left = aggregate_gene_expression(left_cds[m_left,],cell_group_df = pt_group_left,scale_agg_value=FALSE)
#Z-transform rows
dat_left = t(scale(t(dat_left)))

dat_left[abs(dat_left)>3] = sign(dat_left[abs(dat_left)>3]) *3
pheatmap(dat_left,cluster_cols=FALSE)

#make the "right" heatmap

m_right = match(rownames(pr_deg_ids),rownames(right_cds))[seq(100)]

pt_group_right = data.frame(cells = names(right_pseudot),
                           group = cut_number(right_pseudot,100))
dat_right = aggregate_gene_expression(right_cds[m_right,],cell_group_df = pt_group_right,scale_agg_value=FALSE)
#Z-transform rows
dat_right = t(scale(t(dat_right)))

dat_right[abs(dat_right)>3] = sign(dat_right[abs(dat_right)>3]) *3
pheatmap(dat_right,cluster_cols=FALSE)

left_reverse=dat_left[,order(ncol(dat_left):1)]
pheatmap(right_reverse,cluster_cols=FALSE)


all_mtx=cbind(left_reverse, dat_right)
pheatmap(all_mtx,cluster_cols=FALSE)
adr_all_med_qm=quickMarkers(adr_all_med@assays$RNA@counts, adr_all_med@active.ident, N=20)

dat_left_reverse=dat_left[,order(ncol(dat_left):1)]

joined_unscaled=cbind(dat_left_reverse, dat_right)

joined_scaled = t(scale(t(joined_unscaled)))
joined_scaled[abs(joined_scaled)>3] = sign(joined_scaled[abs(joined_scaled)>3]) *3

pheatmap(joined_scaled,cluster_cols=FALSE)
kmeans_1 =kmeans(joined_scaled, 10)

split_genes=kmeans_1$cluster

split_genes=factor(split_genes, levels=c(6,8,10,1,3,7,2,9,4,5))
col_fun = colorRamp2(c(-3,-1,0,1,3), rev(brewer.pal(n = 5, name = "RdYlBu")))
Heatmap(joined_scaled,col=col_fun, cluster_columns = F, cluster_rows=F, row_order = names(sort(split_genes)), show_column_names =F)
Heatmap(joined_scaled, cluster_columns = F, cluster_rows=F, row_split = split_genes, show_column_names =F)

cds_subset <- choose_cells(cds)

levels(adr_all_med)= c(25, 24,13,19,20)

adr_all_med=ScaleData(adr_all_med, features = rownames(adr_all_med))
adr_all=ScaleData(adr_all, features = rownames(adr_all))
mouse_genes=c("ZIC3", "OLIG3", "BMP6", "LMX1A",
              "DLX5", "PAK3", "HAPLN1",
              "PHOX2B", "ASCL1", "LCP2",
              "PRRX1", "MEOX1", "TWIST1", "ETV4",
              "PLP1", "ZNF488", "NKAIN4", "EGFLAM", "MSTN",
              "NEUROG2", "POU4F1", "NEUROG1", "NEUROD4", "NEUROD1", "EYA2", "ISL1",
              "SIX1")
mouse_gene_type=c("premigratory","premigratory","premigratory","premigratory", 
                  "delamination","delamination","delamination",
                  "autonomic", "autonomic","autonomic",
                  "mesenchyme", "mesenchyme", "mesenchyme", "mesenchyme", 
                  "glia", "glia", "glia", "glia", "glia",
                  "sensory",  "sensory", "sensory", "sensory", "sensory", "sensory", "sensory", "sensory")

mouse_gene_df=data.frame(gene=mouse_genes, type=mouse_gene_type)
mouse_gene_df$type=factor(mouse_gene_df$type, levels=c("premigratory","delamination", "autonomic","mesenchyme", 
                                                       "glia","sensory"))
col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
Heatmap(adr_all_med@assays$RNA@scale.data[as.character(mouse_gene_df$gene),],col=col_fun, row_split = mouse_gene_df$type, cluster_columns = F, cluster_rows = F, show_column_names = F, column_split = adr_all_med@active.ident)
adr_all_med@meta.data$devTime


FeaturePlot(adr_all_med, features = c("SOX10", "HTR3A", "TH", "ERBB3", "DLL3", "CHGB", "FOXD3", "THBD", "FOXQ1", "CARTPT"))

all_mouse_genes=read.csv("mouse_genes.csv")

mouse_pseudot=all_mouse_genes[order(all_mouse_genes$peak.time),]
mouse_pseudot$gene=as.character(mouse_pseudot$gene)
overlap=mouse_pseudot[which(mouse_pseudot$gene%in%rownames(adr_all_med)),]

musGenes <- c("Hmmr", "Tlx3", "Cpeb4")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  
  # Print the first 6 genes found to the scree
  return(genesV2)
}

mouse_proper=convertMouseGeneList(all_mouse_genes$gene)


col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdBu")))
Heatmap(adr_all_med@assays$RNA@scale.data[overlap$gene,],col=col_fun, cluster_columns = F, cluster_rows = F, show_column_names = F, show_row_names=F,column_split = adr_all_med@active.ident)
