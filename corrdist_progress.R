adr=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
adr_pseudot= read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/cellranger220_count_30714_WSSS8011223_GRCh38-1_2_0/monocle_pseudotime.csv", sep = "\t")

adr@meta.data$cell_tags=sapply(strsplit(colnames(adr), "_"), "[", 2)

scp_cells=rownames(adr@meta.data)[which(adr@meta.data$cell_tags%in%sapply(strsplit(rownames(adr_pseudot), "___"), "[", 2))]
just_scp_traj=subset(adr, cells = scp_cells)
is.unsorted(match(just_scp_traj@meta.data$cell_tags, sapply(strsplit(rownames(adr_pseudot), "___"), "[", 2)))

just_scp_traj@meta.data$pseudot=adr_pseudot$x[match(just_scp_traj@meta.data$cell_tags,
                                                    sapply(strsplit(rownames(adr_pseudot), "___"), "[", 2))]
FeaturePlot(just_scp_traj, features = "pseudot") + scale_color_viridis(option="viridis", name="Pseudotime")
+ggtitle("SCP trajectory coloured by pseudotime")
rownames(just_scp_traj@meta.data)[order(just_scp_traj@meta.data$pseudot)]

mature_count_table =read.table("matureAdrenalNeuronalCounts.tsv", header = T)
mature_metadata = read.table("matureAdrenalNeuronalMetadata.tsv", header = T)
mature_matrix = as.matrix(mature_count_table)
colnames(mature_matrix) = rownames(mature_metadata)
mature_srat=CreateSeuratObject(mature_matrix, meta.data = mature_metadata)
mature_non_cycling_cells=rownames(mature_srat@meta.data)[which(mature_srat@meta.data$Phase=="G1")]
mature_cycling_cells=rownames(mature_srat@meta.data)[which(mature_srat@meta.data$Phase!="G1")]

scp_trajectory_cells=rownames(just_scp_traj@meta.data)[order(just_scp_traj@meta.data$pseudot)]
scp_non_cycling=rownames(adr@meta.data[scp_trajectory_cells,])[which(adr@meta.data[scp_trajectory_cells,]$Phase=="G1")]
scp_cycling=rownames(adr@meta.data[scp_trajectory_cells,])[which(adr@meta.data[scp_trajectory_cells,]$Phase!="G1")]

cortex_cells=rownames(adr@meta.data)[which(adr@meta.data$seurat_clusters%in%c("1","2", "4","6", "12", "13"))]
cortex_non_cycling=rownames(adr@meta.data[cortex_cells,])[which(adr@meta.data[cortex_cells,]$Phase=="G1")]
cortex_cycling=rownames(adr@meta.data[cortex_cells,])[which(adr@meta.data[cortex_cells,]$Phase!="G1")]

mesen_cells=rownames(adr@meta.data)[which(adr@meta.data$seurat_clusters%in%c("3", "9", "15"))]
mesen_non_cycling=rownames(adr@meta.data[mesen_cells,])[which(adr@meta.data[mesen_cells,]$Phase=="G1")]
mesen_cycling=rownames(adr@meta.data[mesen_cells,])[which(adr@meta.data[mesen_cells,]$Phase!="G1")]



scp_nc_df=data.frame(cell=scp_non_cycling, type="SCP NC")
scp_cc_df=data.frame(cell=scp_cycling, type="SCP CC")

mature_nc_df=data.frame(cell=mature_non_cycling_cells, type="MatureAdrNeur NC")
mature_cc_df=data.frame(cell=mature_cycling_cells, type="MatureAdrNeur CC")

cortex_nc_df=data.frame(cell=cortex_non_cycling, type="Cortex NC")
cortex_cc_df=data.frame(cell=cortex_cycling, type="Cortex CC")

mesen_nc_df=data.frame(cell=mesen_non_cycling, type="Mesenchyme NC")
mesen_cc_df=data.frame(cell=mesen_cycling, type="Mesenchyme CC")

annotated_cells=rbind(scp_nc_df, scp_cc_df, mature_nc_df, mature_cc_df, cortex_nc_df, cortex_cc_df, mesen_nc_df, mesen_cc_df)
scp_df = rbind(scp_nc_df, scp_cc_df)
scp_df$cell= as.character(scp_df$cell)
inhouse_nb$seurat_clusters = nb_srat@meta.data$seurat_clusters
combined_in_nb = merge(inhouse_nb,y=list(mature_srat, adr), project="corrdist")

combined_in_nb[["percent.mt"]] = PercentageFeatureSet(combined_in_nb, pattern = "^MT-")

combined_in_nb=NormalizeData(combined_in_nb)
combined_in_nb=ScaleData(combined_in_nb, features = rownames(combined_in_nb))
combined_in_nb@assays$RNA@scale.data

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
annotated_cells$cell=as.character(annotated_cells$cell)
x=get_dist_matrix(combined_in_nb, annotated_cells$cell, colnames(inhouse_nb))
y=get_dist_matrix(combined_in_nb, scp_non_cycling, colnames(inhouse_nb))
z=get_dist_matrix(combined_in_nb, scp_cycling, colnames(inhouse_nb))
col_fun = colorRamp2(c(0.75, 1, 1.25), c("green", "white", "brown"))
just_scp_corrdist=as.matrix(x)
ComplexHeatmap::Heatmap(x, name="correlation distance", col=col_fun, 
                        column_order = colnames(x), row_order = rownames(x),
                        column_split = annotated_cells$type,
                        show_column_names = F, row_names_side = "left",
                        column_title_rot = 90, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                        column_title_gp = gpar(fontsize=8), 
                        column_gap = unit(c(2, 7,2,7,2,7,2), "mm"),
                        )
dim(y)
hm_nc=ComplexHeatmap::Heatmap(y, name="correlation distance", col=col_fun, 
                        column_order = colnames(y), row_order = rownames(y),
                        show_column_names = F, row_names_side = "left",
                        column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                        column_title_gp = gpar(fontsize=10), width = unit(6, "cm"), 
                        column_title = "SCP non-cycling, ordered by pseudotime", bottom_annotation = ha)

hm_c=ComplexHeatmap::Heatmap(z, name="correlation distance", col=col_fun, 
                        column_order = colnames(z), row_order = rownames(z),
                        show_column_names = F, row_names_side = "left",
                        column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                        column_title_gp = gpar(fontsize=10), width = unit(4, "cm"), column_title = "SCP cycling, ordered by pseudotime")


markers_scp_ncc = combined_in_nb@assays$RNA@scale.data[c("SOX2", "SOX10", "PHOX2A", "PHOX2B", "ISL1"),scp_non_cycling]
markers_scp_cc=combined_in_nb@assays$RNA@scale.data[c("SOX2", "SOX10", "PHOX2A", "PHOX2B", "ISL1"),scp_cycling]
col_marker = colorRamp2(c(0, 2), c("white", "black"))
markers_ncc=ComplexHeatmap::Heatmap(markers_scp_ncc, name="scaled expression", col=col_marker, 
                                    column_order = colnames(markers_scp_ncc), row_order = rownames(markers_scp_ncc),
                                    show_column_names = F, row_names_side = "left",
                                    column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                                    column_title_gp = gpar(fontsize=10), width = unit(6, "cm"), column_title = "SCP cycling, ordered by pseudotime", border = T)
markers_cc=ComplexHeatmap::Heatmap(markers_scp_cc, name="scaled expression", col=col_marker, 
                                    column_order = colnames(markers_scp_cc), row_order = rownames(markers_scp_cc),
                                    show_column_names = F, row_names_side = "left",
                                    column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                                    column_title_gp = gpar(fontsize=10), width = unit(4, "cm"), column_title = "SCP cycling, ordered by pseudotime", border = T)

ht_list=(hm_nc %v% markers_ncc)
list=(hm_c%v%markers_cc)
draw(ht_list, ht_gap = unit(1.5, "cm"))
draw(list, ht_gap = unit(1.5, "cm"))
cycle_col=colorRamp2(c(-1,0, 1), c("black", "red", "yellow"))
total_list=ht_list + list
markers_scp_ncc=as.data.frame(markers_scp_ncc)
markers_scp_ncc=t(markers_scp_ncc)
ha = HeatmapAnnotation(S_score=combined_in_nb@meta.data[scp_non_cycling,]$S.Score,
                       G2M_score=combined_in_nb@meta.data[scp_non_cycling,]$G2M.Score, exp=t(markers_scp_ncc),
    annotation_name_side = "left" , col=list(exp=col_marker, S_score=cycle_col, G2M_score=cycle_col),
    border = T, gap = unit(c(2, 10), "mm"))

hb = HeatmapAnnotation(S_score=combined_in_nb@meta.data[scp_cycling,]$S.Score,
                       G2M_score=combined_in_nb@meta.data[scp_cycling,]$G2M.Score,
                       exp=t(markers_scp_cc),
                       annotation_name_side = "left" , col=list(exp=col_marker, S_score=cycle_col, G2M_score=cycle_col),
                       border = T, gap = unit(c(2, 10), "mm"))
hm_nc=ComplexHeatmap::Heatmap(y, name="correlation distance", col=col_fun,
                              column_order = colnames(y), row_order = rownames(y),
                              show_column_names = F, row_names_side = "left",
                              column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                              column_title_gp = gpar(fontsize=10), width = unit(6, "cm"), 
                              column_title = "SCP non-cycling, ordered by pseudotime", bottom_annotation = ha)
hm_c=ComplexHeatmap::Heatmap(z, name="correlation distance", col=col_fun, 
                             column_order = colnames(z), row_order = rownames(z),
                             show_column_names = F, row_names_side = "left",
                             column_title_rot = 0, row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                             column_title_gp = gpar(fontsize=10), width = unit(4, "cm"),
                             column_title = "SCP cycling, ordered by pseudotime",bottom_annotation = hb)

test_list=hm_nc+hm_c

draw(test_list,  ht_gap = unit(2, "cm"), auto_adjust=F)
