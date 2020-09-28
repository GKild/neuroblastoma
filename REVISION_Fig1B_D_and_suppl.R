library(Seurat)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(SoupX)
library(viridis)
source("logisticRegression.R")
#literature markers in adrenal gland
adr_all=readRDS("fAdrenal/processedSeurat.RDS")
fetal_pods=readRDS("fetal_podocytes_forGerda.RDS")
genes_10x <- read.table("features.tsv.gz", sep = "\t", header = F)

FeaturePlot(adr_all, c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
                       "CHGB","GAP43", "BCL2", "NPY", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
                       "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC"))
#annotate louvain clusters
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

adr_all@meta.data$new_clust=factor(adr_all@meta.data$new_clust, levels = c("SCPs","Bridge","Sympathoblastic",
                                                                           "Chromaffin","Endothelium",
                                                                           "Mesenchyme", "Cortex",
                                                                           "Leukocytes", "Erythroblasts", "Other"))
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

#make figure 1B
pdf("adrenal_gland.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(adr_all, group.by = "new_clust", label=T,label.size = 3,cols=c(viridis(10)[c(1,3,5,7)],c(rep("grey38",5), "lightgrey"))) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()

adr_all_med=subset(adr_all, idents=c(25,24,13,19,20))
pdf("just_medulla.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(adr_all_med, group.by = "new_clust", label=T,label.size = 3,cols=c(viridis(10)[c(1,3,5,7)])) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()

####algorithmic and lit markers in the adrenal gland######
#assign annotation as active.ident
adr_all_supp = adr_all
adr_names=names(adr_all_supp@active.ident)
adr_all_supp@active.ident=factor(adr_all_supp@meta.data$new_clust)
names(adr_all_supp@active.ident)=adr_names
adr_all_supp=ScaleData(adr_all_supp, features = rownames(adr_all_supp))
cluster.averages <- AverageExpression(adr_all_supp, return.seurat = T)
lit_marks=c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
            "CHGB","GAP43", "BCL2", "NPY", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
            "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC", "FOXD3", "HTR3A", "THBD", "CARTPT", "FOXQ1")

levels(adr_all_supp)=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme", "Endothelium",
                       "Erythroblasts" ,"Cortex", "Leukocytes", "Other")

#literature marker cluster average heatmap
col_ord=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme", "Endothelium",
          "Erythroblasts" ,"Cortex", "Leukocytes", "Other")
col_fun = colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(n = 5, name = "RdYlBu")))

#algorithmic markers of each cat
top_all_clust=quickMarkers(adr_all_supp@assays$RNA@counts, adr_all_supp@active.ident, N=Inf)
# get top 5 for each cat that aren't in the lit markers
top_clust_no_lit=top_all_clust[-which(top_all_clust$gene%in%lit_marks), ]
top_5=by(top_clust_no_lit, top_clust_no_lit["cluster"], head, n=5)
top5=Reduce(rbind, top_5)

lit_mark_clust=c("SCPs","SCPs","SCPs","SCPs", 
                 "Bridge",  "Bridge", 
                 "Sympathoblastic", "Sympathoblastic","Sympathoblastic","Sympathoblastic","Sympathoblastic","Sympathoblastic","Sympathoblastic",
                 "Chromaffin", "Mesenchyme", "Mesenchyme", "Endothelium","Endothelium","Endothelium","Endothelium","Erythroblasts", "Erythroblasts",
                 "Erythroblasts", "Cortex", "Cortex", "Leukocytes", "SCPs", "Bridge", "Bridge", "Sympathoblastic", "Chromaffin")

lit_df=data.frame(gene=lit_marks, cluster=lit_mark_clust)
algo_df=data.frame(gene=top5$gene, cluster=top5$cluster)
joint_df=rbind(lit_df, algo_df)
joint_df$cluster=factor(joint_df$cluster, levels=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme", "Endothelium",
                                                   "Erythroblasts" ,"Cortex", "Leukocytes", "Other"))


for_hm_inm=cluster.averages@assays$RNA@scale.data[as.character(joint_df$gene),]
#algorithmic marker heatmap - figure 1C
Heatmap(for_hm_inm, col=col_fun, cluster_rows=F, column_order = col_ord,
            row_split = joint_df$cluster, name = "Normalized expression", row_names_gp = gpar(fontsize = 8), row_gap = unit(1, "mm"), row_title_rot = 0)

####scps in different tissues########

s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
scp_mtx=readMM('all_scp_matrix.mtx')

obs=read.csv("all_scps/obs.csv")
cells=read.csv("all_scps/var.csv")
rownames(cells)=cells$index
rownames(scp_mtx)=as.character(obs$index)
colnames(scp_mtx)=as.character(cells$index)

scp_seurat=CreateSeuratObject(scp_mtx, meta.data = cells)
process_10x =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


scp_seurat=process_10x(scp_seurat)

scp_seurat@meta.data$new_labs="x"
scp_seurat@meta.data$new_labs[grep("bone",scp_seurat@meta.data$organ)]="bone SCPs"
scp_seurat@meta.data$new_labs[grep("skin",scp_seurat@meta.data$organ)]="skin SCPs"
scp_seurat@meta.data$new_labs[grep("gut",scp_seurat@meta.data$organ)]="gut SCPs"
scp_seurat@meta.data$new_labs[grep("adrenal",scp_seurat@meta.data$organ)]="adrenal SCPs"

#exclude clusters that aren't absolutely SCPs
FeaturePlot(scp_seurat, c("SOX2", "SOX10", "ERBB3", "MPZ", "PLP1", "ASCL1", "DLL3"))
scp_seurat=subset(scp_seurat, idents=c(0,6,17,14), invert=T)
#logistic regression to compare adrenal SCPs with all the other SCPs. 
scp_pred=predictSimilarity(train_with_podos, scp_seurat@assays$RNA@counts, scp_seurat@meta.data$new_labs, minGeneMatch = 0.8)
similarityHeatmap(scp_pred, column_order=lr_col_ord)

#comparing bilateral adrenal gland with all the other grands (Fig 1D)
just_med=subset(adr_all,
                cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$new_clust%in%c("SCPs", "Bridge",
                                                                                         "Sympathoblastic", "Chromaffin"))])
just_bilat=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name%in%c("w10d5_1", "w10d5_2"))])
bone_scps=subset(scp_seurat,
                 cells=rownames(scp_seurat@meta.data)[which(scp_seurat@meta.data$new_labs=="bone SCPs")])

bilat_med=subset(just_bilat,
                 cells=rownames(just_bilat@meta.data)[which(just_bilat@meta.data$new_clust%in%c("SCPs", "Bridge",
                                                                                                "Sympathoblastic", "Chromaffin", 
                                                                                                "Endothelium", "Mesenchyme"))])
bone_med_genes=intersect(rownames(bone_scps), rownames(just_med))

med_and_bone=cbind(just_med@assays$RNA@counts[bone_med_genes,], bone_scps@assays$RNA@counts[bone_med_genes,])

med_and_bone_labs=c(as.character(just_med@meta.data$new_clust), bone_scps@meta.data$new_labs)

fetal_pds_counts=fetal_pods@assays$RNA@counts
rownames(fetal_pds_counts)=genes_10x$V2


merged_cells=cbind(adr_all@assays$RNA@counts, fetal_pds_counts)
cell_labs=c(adr_all@meta.data$new_clust, rep("Podocytes", 278))

bilat_and_podo=cbind(bilat_med@assays$RNA@counts, fetal_pds_counts)
bilat_podo_labs=c(as.character(bilat_med@meta.data$new_clust), rep("Podocytes", 278))
bilat_train=trainModel(bilat_and_podo, bilat_podo_labs, workers=NULL)

med_ag_bilat=predictSimilarity(bilat_train, med_and_bone, minGeneMatch = 0.8)
sample_ha=c(just_med@meta.data$sample_name, rep("bone SCPs", 299))
logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')

tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#0c3823", "#999933", "#DDCC77", "#ff9a00","#AA4499")
names(tol8qualitative)=c("w8", "w8d6", "w10d5_1", "w10d5_2", "w11", "w21_1", "w21_2", "bone SCPs")

ha = rowAnnotation(Sample = factor(sample_ha, levels = c("w8", "w8d6", "w10d5_1", "w10d5_2", "w11", "w21_1", "w21_2", "bone SCPs")),
                   col = list(Sample = tol8qualitative))
col_fun_2=colorRamp2(c(-3,-2,-1,0,1,2,3), c("#000000", "#404040", "#bfbfbf","#f8f8f8","#f8f8ba", "#f8f87c", "#f8f800"))
similarityHeatmap(med_ag_bilat,
        show_row_names = F, cluster_rows=F, cluster_columns = F, name = "Predicted\nSimilarity\n(Logit)",
        row_split=factor(med_and_bone_labs, levels=c("bone SCPs", "SCPs", "Bridge", "Sympathoblastic", "Chromaffin")), 
        column_order = c("SCPs", "Bridge", "Sympathoblastic", "Chromaffin", "Mesenchyme", "Endothelium", "Podocytes"), left_annotation = ha)

just_med@meta.data$broad_time="x"

just_med@meta.data$broad_time[which(just_med@meta.data$sample_name%in%c("w8", "w8d6","w10d5_1", "w10d5_2", "w11"))]="early"
just_med@meta.data$broad_time[which(just_med@meta.data$sample_name%in%c("w21_1", "w21_2"))]="late"

just_med@meta.data$time_ct=paste0(just_med@meta.data$broad_time, "_", just_med@meta.data$new_clust)
table(just_med@meta.data$time_ct)

excludeGenes=c(riboGenes, mtGenes, riboGenes)
keep_genes=rownames(just_med)[which(!rownames(just_med)%in%excludeGenes)]
just_med=subset(just_med, features=keep_genes)
just_med=ScaleData(just_med, features=rownames(just_med))
med_cell_names=names(just_med@active.ident)
new_med_ident=just_med@meta.data$time_ct
names(new_med_ident)=med_cell_names
just_med@active.ident=as.factor(new_med_ident)
scp_de = FindMarkers(just_med, ident.1 = "early_SCPs", ident.2 = "late_SCPs", min.pct = 0.25)
bridge_de=FindMarkers(just_med, ident.1 = "early_Bridge", ident.2 = "late_Bridge", min.pct = 0.25)
symp_de=FindMarkers(just_med, ident.1 = "early_Sympathoblastic", ident.2 = "late_Sympathoblastic", min.pct = 0.25)
chromaffin_de=FindMarkers(just_med, ident.1 = "early_Chromaffin", ident.2 = "late_Chromaffin", min.pct = 0.25)

scp_de$category="early_vs_late_SCPs"
bridge_de$category="early_vs_late_Bridge"
symp_de$category="early_vs_late_Sympathoblast"
chromaffin_de$category="early_vs_late_Chromaffin"

de_marks_early_late=rbind(scp_de, bridge_de, symp_de, chromaffin_de)
de_marks_early_late=de_marks_early_late[which(de_marks_early_late$p_val_adj<0.01),]
de_marks_early_late$gene=rownames(de_marks_early_late)
de_marks_early_late=de_marks_early_late[,c(7,1:6)]
write.table(de_marks_early_late, "de_marks_early_late.txt", sep = "\t", quote = F, row.names = F, col.names = T)
#####medulla pseudotime#########
counts.data <-GetAssayData(object = adr_all_med, slot = "counts", assay="RNA")
pheno.data <- as.data.frame(adr_all_med@meta.data, row.names = rownames(adr_all_med@meta.data))
feature.data <- as.data.frame(adr_all_med@assays$RNA@counts@Dimnames[[1]], row.names = adr_all_med@assays$RNA@counts@Dimnames[[1]])

#create a monocle3 object
cds <- new_cell_data_set(counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- cds@colData$new_clust
names(list_cluster) <- adr_all_med@assays[["RNA"]]@data@Dimnames[[2]]


cds@clusters@listData[["UMAP"]][["clusters"]] <- as.factor(list_cluster)


cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


# Assign UMAP coordinate

cds@reducedDims@listData[["UMAP"]] <-Embeddings(adr_all_med, reduction =  'umap')
cells_embeddings=Embeddings(adr_all_med, reduction =  'umap')
#learn principal graph
cds <- learn_graph(cds,use_partition = F)
#order cells with SCP node in shiny app
cds <- order_cells(cds)
#if shiny app doesn't work, use these cells
cds <- order_cells(cds, root_cells = c("babyAdrenal2_GTACTCCAGTCAAGGC","5388STDY7717454_CCAATCCAGTACACCT",
                                       "5388STDY7717456_AACACGTTCGTTACGA","5388STDY7717456_ATCTACTAGTAAGTAC",
                                       "5388STDY7717456_GACGTGCCATGCTGGC","5388STDY7717458_AAAGATGCATTACCTT",
                                       "5388STDY7717458_ACGCCGATCAGCACAT","5388STDY7717458_ACGGGCTTCGTCTGAA",
                                       "5388STDY7717458_CCTACCAAGGACAGAA","5388STDY7717458_TTCTCAAGTCTAGCGC",
                                       "5388STDY7717459_AACTCAGAGACCACGA" ,"5388STDY7717459_TGAAAGACAATGGATA",
                                       "cellranger302_count_32644_WSSS_F_Adr8710635_GRCh38-1_2_0_TTGGGCGAGCGTGAAC",
                                       "5698STDY7839907_GAATAAGCAGGAACGT","5698STDY7839907_TTCTCAATCGCTTGTC",
                                       "5698STDY7839909_AACCTGATCTGGCTGG","5698STDY7839909_ACACCAATCTCATAGG",
                                       "5698STDY7839909_ATACCTTCAGGCTATT","5698STDY7839909_CATAAGCTCCCTCATG",
                                       "5698STDY7839909_CTGCCATCAACTCGAT","5698STDY7839909_CTGGACGTCGTTCTGC",
                                       "5698STDY7839909_CTTACCGAGGTTCATC","5698STDY7839909_GAGGCCTAGACTCTAC",
                                       "5698STDY7839909_GTCTCACGTCACTGAT","5698STDY7839909_TAGTGCAAGAGAACCC",
                                       "5698STDY7839909_TGGATGTTCTTTCCAA" ,"5698STDY7839909_TTCCACGGTACGTAGG",
                                       "5698STDY7839917_AGTAGTCTCATTTGGG","5698STDY7839917_CAAGAAACACTTGGAT",
                                       "5698STDY7839917_CTCGTCATCGCAAGCC","5698STDY7839917_CTCTGGTCACTTAACG",
                                       "5698STDY7839917_CTTAGGAAGGTGTTAA" ,"5698STDY7839917_GGAATAAGTACAGTGG",
                                       "5698STDY7839917_TACAGTGCATACTCTT","5698STDY7839917_TTCCCAGAGCACCGCT"))
pseudo_srat = pseudotime(cds, reduction_method = 'UMAP')
adr_all_med@meta.data$pseudotime=pseudo_srat

pdf("all_medulla_pseudotime.pdf", height = 3, width=3,useDingbats = F)
gg=FeaturePlot(adr_all_med, features = "pseudotime", pt.size = 2) + scale_color_viridis(option="plasma", direction = 1) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank()) +ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()

pdf("all_medulla_pseudotime_no_leg.pdf", height = 3, width=3,useDingbats = F)
gg=FeaturePlot(adr_all_med, features = "pseudotime", pt.size = 2) + scale_color_viridis(option="plasma", direction = 1) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none") +ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
#same but just bridge to sympathoblast
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

pr_deg_ids_ord$gene=rownames(pr_deg_ids_ord)
pr_deg_ids_ord=pr_deg_ids_ord[,c(8,1:4,7)]
pr_deg_ids_ord$category="bridge_to_symp_only"

write.table(pr_deg_ids_ord, "bridge_to_symp_traj.txt", row.names = F, col.names = T, quote = F, sep = "\t")
###analyse the medulla separately to show that when you do it on it's own the bridge isn't there 

redo_adr =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 33)
  srat = FindNeighbors(srat, dims=1:33)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:33, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}
just_med_redo=redo_adr(just_med)
just_med_redo@meta.data$new_clust=factor(just_med_redo@meta.data$new_clust, levels=c("SCPs", "Bridge", "Chromaffin","Sympathoblastic"))
DimPlot(just_med_redo, group.by = "new_clust", label = T, cols=c(viridis(10)[c(1,3,5,7)])) +theme(legend.position = "none") 


#####markers in human and mouse and LR train on inhouse and match against mouse#########

m13=read.delim("GSE99933_E13.5_counts.txt.gz")
m13=as.matrix(m13)

m13
m13_srat=CreateSeuratObject(m13)
m13_srat = NormalizeData(m13_srat)
m13_srat =FindVariableFeatures(m13_srat, selection.method = "vst", nfeatures = 2000)
m13_srat = ScaleData(m13_srat, features = rownames(m13_srat))
m13_srat = RunPCA(m13_srat, npcs = 25)
m13_srat = FindNeighbors(m13_srat, dims=1:25)
m13_srat = FindClusters(m13_srat, resolution = 0.6)
m13_srat = RunUMAP(m13_srat, dims=1:25, min.dist = 0.5, n.neighbors = 30)

DimPlot(m13_srat, label=T)
FeaturePlot(m13_srat, c("Sox10", "Mpz", "Dll3","Tlx2", "Th", "Cartpt", "Pnmt","Gap43", "Mapt"))

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(m13) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
dim(m13_srat)
dim(genesV2)
mouse_proper=convertMouseGeneList(rownames(m13_srat))

genesV2$MGI.symbol[grep("Mpz",genesV2$MGI.symbol)]


m13_srat@meta.data$new_clust="x"

m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==0)]="sympathoblast"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==1)]="chromaffin"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==2)]="scp"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==3)]="bridge"
DimPlot(m13_srat, group.by = "new_clust", label=T, pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                                                         axis.ticks = element_blank() ,legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") 

FeaturePlot(m13_srat, c("Sox10", "Erbb3", "Foxd3","Htr3a","Dll3","Thbd", "Th", "Chgb","Foxq1", "Cartpt", "Slc18a3"))


chromaffin_for_m=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="chromaffin"),]$gene)[1:20]
sympathoblastic_form_m=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="sympathoblastic"),]$gene)[1:20]

right_as_mouse=genesV2$MGI.symbol[which(genesV2$HGNC.symbol%in%chromaffin_for_m)]
left_as_mouse=genesV2$MGI.symbol[which(genesV2$HGNC.symbol%in%sympathoblastic_form_m)]


FeaturePlot(m13_srat, right_as_mouse)
FeaturePlot(m13_srat, left_as_mouse)

length(which(genesV2$MGI.symbol%in%rownames(m13_srat)))

t.first <- genesV2[match(unique(as.character(genesV2$HGNC.symbol)), as.character(genesV2$HGNC.symbol)),]
dim(t.first)

t.second=t.first[match(unique(as.character(t.first$MGI.symbol)), as.character(t.first$MGI.symbol)),]

final_overlap=t.second[which(t.second$HGNC.symbol%in%rownames(adr_all)),]

mouse_subs=subset(m13_srat, features=as.character(final_overlap$MGI.symbol))
adr_subs=subset(adr_all, features=as.character(final_overlap$HGNC.symbol))


adr_mtx=adr_subs@assays$RNA@counts


adr_right_order=adr_mtx[as.character(final_overlap$HGNC.symbol),]

rownames(adr_right_order)=as.character(final_overlap$MGI.symbol)


mouse_mtx=m13_srat@assays$RNA@counts
mouse_right_order=mouse_mtx[as.character(final_overlap$MGI.symbol),]



fit_for_mouse = trainModel(adr_right_order, adr_all@meta.data$new_clust)

ps_mouse = predictSimilarity(fit_for_mouse, mouse_right_order,m13_srat@meta.data$new_clust)
similarityHeatmap(ps_mouse, row_order=c("scp", "bridge", "sympathoblast", "chromaffin"),
                  column_order=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme","Cortex","Leukocytes", "Vascular endothelium", "Erythroblasts","Other"))


FeaturePlot(adr_med, c("SOX10", "ERBB3", "FOXD3","HTR3A","DLL3","THBD", "TH", "CHGB","FOXQ1", "CARTPT", "SLC18A3"), ncol=3) +
  theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank(), legend.position = "none") +ylab("UMAP 2") +xlab("UMAP 1")


FeaturePlot(adr_med, c("CARTPT", "SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH","DBH", "PHOX2B", "CHGB", "GAP43", "BCL2", "NPY", "PNMT"), ncol=5)
FeaturePlot(m13_srat, c( "Sox2", "Sox10", "Mpz", "Plp1", "Erbb3", "Dll3", "Tlx2", "Th","Dbh", "Phox2b", "Chgb", "Gap43", "Bcl2", "Npy", "Pnmt"), ncol = 5)

