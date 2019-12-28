library(Seurat)
library(ggplot2)
exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
#get baby adrenal 2 reference
adr=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
#get pseudotime of scp traj
adr_pseudot= read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/cellranger220_count_30714_WSSS8011223_GRCh38-1_2_0/monocle_pseudotime.csv", sep = "\t")
#unique cell barcode in adr
adr@meta.data$cell_tags=sapply(strsplit(colnames(adr), "_"), "[", 2)
#unique barcode for SCP cells that are in pseudotime list
scp_cells=rownames(adr@meta.data)[which(adr@meta.data$cell_tags%in%sapply(strsplit(rownames(adr_pseudot), "___"), "[", 2))]
#seurat object that only has scp trajectory cells
just_scp_traj=subset(adr, cells = scp_cells)
#add pseudotime to only SCP seurat object
just_scp_traj@meta.data$pseudot=adr_pseudot$x[match(just_scp_traj@meta.data$cell_tags,
                                                    sapply(strsplit(rownames(adr_pseudot), "___"), "[", 2))]
#plot pseudotime on UMAP
FeaturePlot(just_scp_traj, features = "pseudot") + scale_color_viridis(option="viridis", name="Pseudotime")
+ggtitle("SCP trajectory coloured by pseudotime")

#add mature adrenal neuronal dataset
mature_count_table =read.table("matureAdrenalNeuronalCounts.tsv", header = T)
mature_metadata = read.table("matureAdrenalNeuronalMetadata.tsv", header = T)
#make it into seurat object
mature_matrix = as.matrix(mature_count_table)
colnames(mature_matrix) = rownames(mature_metadata)
mature_srat=CreateSeuratObject(mature_matrix, meta.data = mature_metadata)
#split cells into cycling and non-cycling
mature_non_cycling_cells=rownames(mature_srat@meta.data)[which(mature_srat@meta.data$Phase=="G1")]
mature_cycling_cells=rownames(mature_srat@meta.data)[which(mature_srat@meta.data$Phase!="G1")]
#final mature order
final_mature_order=c(mature_non_cycling_cells, mature_cycling_cells)
#all scp trajectory cells ordered by pseudotime
scp_trajectory_cells=rownames(just_scp_traj@meta.data)[order(just_scp_traj@meta.data$pseudot)]
#split SCPs into cycling and non-cycling
scp_non_cycling=rownames(adr@meta.data[scp_trajectory_cells,])[which(adr@meta.data[scp_trajectory_cells,]$Phase=="G1")]
scp_cycling=rownames(adr@meta.data[scp_trajectory_cells,])[which(adr@meta.data[scp_trajectory_cells,]$Phase!="G1")]
#get cortex cells
cortex_cells=rownames(adr@meta.data)[which(adr@meta.data$seurat_clusters%in%c("0","1","2", "4","6","8", "12", "13"))]
#split into cycling and non-cycling
cortex_non_cycling=rownames(adr@meta.data[cortex_cells,])[which(adr@meta.data[cortex_cells,]$Phase=="G1")]
cortex_cycling=rownames(adr@meta.data[cortex_cells,])[which(adr@meta.data[cortex_cells,]$Phase!="G1")]
#final cortex order
final_cortex_order=c(cortex_non_cycling, cortex_cycling)
#get mesenchymal cells
mesen_cells=rownames(adr@meta.data)[which(adr@meta.data$seurat_clusters%in%c("3", "9", "15"))]
#split into cycling and non-cycling
mesen_non_cycling=rownames(adr@meta.data[mesen_cells,])[which(adr@meta.data[mesen_cells,]$Phase=="G1")]
mesen_cycling=rownames(adr@meta.data[mesen_cells,])[which(adr@meta.data[mesen_cells,]$Phase!="G1")]
#final order
final_mesen_order=c(mesen_non_cycling, mesen_cycling)
# a function that takes in seurat objects relevant for neuroblastoma, scales them together and returns a list of complexHeatmaps objects
#objects have to have the same NGenes. ideally low quality cells are also removed at this point.
list_of_cells=list(scp_ncc=scp_non_cycling, scp_cc=scp_cycling, mature_adr=final_mature_order, cortex_cells=final_cortex_order, mesen_cells=final_mesen_order)

get_dist_matrix = function(srat, ctrl_cells, tumour_cells){
  ref_mtx=srat@assays$RNA@scale.data[,ctrl_cells]
  test_mtx=srat@assays$RNA@scale.data[,tumour_cells]
  tumour_clusters=unique(srat@meta.data[tumour_cells, ]$seurat_clusters)
  clust_avg=sapply(tumour_clusters, function(x){
    cells_in_cl=rownames(srat@meta.data[tumour_cells, ])[which(srat@meta.data[tumour_cells, ]$seurat_clusters==x)]
    rowMeans(test_mtx[,cells_in_cl])
  })
  
  thing=apply(ref_mtx, 2, function(x){
    apply(clust_avg, 2,function(y){
      1-cor(x,y)
    })
  })
  return(thing)
}

get_corrdist_plt=function(tum_obj, scp_ref, mature_ref, ref_cell_list){
  #combine into the same object
  combined_in_nb = merge(tum_obj,y=list(mature_ref, scp_ref), project="corrdist")
  keep_features=rownames(combined_in_nb)[which(!rownames(combined_in_nb)%in%exclude_genes)]
  combined_in_nb=subset(combined_in_nb, features = keep_features)
  combined_in_nb=NormalizeData(combined_in_nb)
  combined_in_nb=ScaleData(combined_in_nb, features = rownames(combined_in_nb))
  plot_width=c(6, 4,4,4,4)
  plot_title=c("SCP trajectory, non-cycling", "SCP trajectory, cycling", "Mature Adrenals", "Cortex", "Mesenchyme")
  hm_objects=lapply(c(1:5), function(x){
    col_marker = colorRamp2(c(0, 2), c("white", "black"))
    cycle_col=colorRamp2(c(-1,0, 1), c("black", "red", "yellow"))
    col_fun=colorRamp2(c(0.75, 1, 1.25), c("green", "white", "brown"))
    markers_mtx=combined_in_nb@assays$RNA@scale.data[c("SOX2", "SOX10", "PHOX2A", "PHOX2B", "ISL1"),ref_cell_list[[x]]]
    dist_mtx=get_dist_matrix(combined_in_nb, ref_cell_list[[x]], colnames(tum_obj))
    ha = HeatmapAnnotation(S_score=combined_in_nb@meta.data[ref_cell_list[[x]],]$S.Score,
                           G2M_score=combined_in_nb@meta.data[ref_cell_list[[x]],]$G2M.Score, exp=t(markers_mtx),
                           annotation_name_side = "left" , col=list(exp=col_marker, S_score=cycle_col, G2M_score=cycle_col),
                           border = T, gap = unit(c(2, 10), "mm"))
    hm=ComplexHeatmap::Heatmap(dist_mtx, name="correlation distance", col=col_fun,
                               column_order = colnames(dist_mtx), row_order = rownames(dist_mtx),
                               show_column_names = F, row_names_side = "left",
                               column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                               column_title_gp = gpar(fontsize=10), width = unit(plot_width[x], "cm"), 
                               column_title = plot_title[x], bottom_annotation = ha)
    
    return(hm)
  })
  return(hm_objects)
}


test1=get_corrdist_plt(nb_srat, adr, mature_srat, list_of_cells)
thing=test1[[1]] + test1[[2]] +test1[[3]]+ test1[[4]] +test1[[5]]
pdf("heatmap_all.pdf", width = 16, height=8)
draw(thing, ht_gap = unit(2, "cm"), auto_adjust=F)
dev.off()

common_genes=rownames(srat_tumour)

mature_dutch_subs=subset(mature_srat, features = common_genes)
adr_dutch_subs=subset(adr, features = common_genes)


test2=get_corrdist_plt(srat_tumour,adr_dutch_subs,mature_dutch_subs, list_of_cells)

thing2=test2[[1]] + test2[[2]] +test2[[3]]+ test2[[4]] +test2[[5]]
pdf("heatmap_all_dutch.pdf", width = 16, height=8)
draw(thing2, ht_gap = unit(2, "cm"), auto_adjust=F)
dev.off()

DimPlot(srat_tumour, group.by = "unique_sample")
DimPlot(nb_srat, label=T)
srat_tumour@meta.data$unique_sample
  

