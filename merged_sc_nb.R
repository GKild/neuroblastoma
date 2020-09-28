source("gene_sets.R")
#DEFINE PARAMS
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
#MT cut-off
mtCut = 20
#nGene Cut
minNumGenes = 200
#nUMI cut
minNumUMI = 500
#How many doublets in a cluster before you throw the entire cluster
maxScrubFrac = 0.3
#UMAP parameters to try in grid
minDists = c(0,.01,.1,.3,.5)
NNs = c(5,10,20,50,100)
#Final UMAP parameters
finalMinDist = 0.5
finalNN = 50
#Clustering resolution
clusterRes=1.0
#Number of PCs to use for quick analysis during QC
defNumPCs=50
source("new_quickmarkers.R")
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb <- data
genes_10x <- read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
#main matrix 
mtx <- Jan_nb$data$counts
#column metadata 
meta_data <- Jan_nb$col.annot
#calculate fraction of spikeins
meta_data$spikeIn_fraction<- colSums(mtx[grep("^ERCC-",rownames(mtx)),])/colSums(mtx)
#remove spikeins
mtx <- mtx[-grep("^ERCC-",rownames(mtx)),]
#calculate fraction of reads lost by using 10X ensembl IDs, and fraction of reads lost from  MT genes that are NOT in 10X. 
mtx_10x <-mtx[which(rownames(mtx)%in%genes_10x$V1),]
not_10x<-mtx[which(!rownames(mtx)%in%genes_10x$V1),]
meta_data$reads_lost <-(colSums(mtx)-colSums(mtx_10x))/colSums(mtx)
rownames(not_10x)<-Jan_nb$row.annot$gene[match(rownames(not_10x),Jan_nb$row.annot$ensemblId)]
meta_data$mt_reads_lost <-colSums(not_10x[grep("^MT-",rownames(not_10x)),])/colSums(mtx)

patient_ids <-unique(Jan_nb$col.annot$PatientId)
cells_each_sample <-sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})

cells_list <- list("NB039_Tumour_organoids_in_20perc_HP_3D"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO_TGFb_inh"=cells_each_sample$NB039$`Tumour organoids in OM_WNT_RSPO_TGFb inh`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO"=cells_each_sample$NB039$`Tumour organoids in OM _WNT_RSPO`, 
                   "NB060_Fresh_biopsy"=as.character(cells_each_sample$NB060), 
                   "NB067_Tumour_organoids_in_OM"=as.character(cells_each_sample$NB067),
                   "NB086_Fresh_biopsy"=as.character(cells_each_sample$NB086),
                   "NB098_Fresh_biopsy"=cells_each_sample$NB098$`Fresh biopsy`,
                   "NB098_peripheral_blood_from_fresh_biopsy"=cells_each_sample$NB098$`peripheral blood from fresh biopsy`,
                   "NB098_TILs_from_fresh_biopsy"=cells_each_sample$NB098$`TILs from fresh biopsy`,
                   "NB098_Tumour_organoids"=cells_each_sample$NB098$`Tumour organoids`, 
                   "NB106_Fresh_biopsy"=as.character(cells_each_sample$NB106),
                   "NB107_Fresh_biopsy"=as.character(cells_each_sample$NB107),
                   "NB123_Fresh_biopsy"=as.character(cells_each_sample$NB123),
                   "NB124_Fresh_biopsy"=as.character(cells_each_sample$NB124),
                   "NB125_Fresh_biopsy"=as.character(cells_each_sample$NB125),
                   "NB130_Fresh_biopsy"=as.character(cells_each_sample$NB130),
                   "NB132_Fresh_biopsy"=as.character(cells_each_sample$NB132),
                   "NB138_Fresh_biopsy"=as.character(cells_each_sample$NB138),
                   "NB151_Fresh_biopsy"=as.character(cells_each_sample$NB151),
                   "NB152_Fresh_biopsy"=cells_each_sample$NB152$`Fresh biopsy`,
                   "NB152_Tumour_organoids_P1"=cells_each_sample$NB152$`Tumour organoids P1`)


#add the unique sample to metadata 

meta_data$unique_sample <- "sample"
for (x in 1:length(cells_list)) {
  meta_data$unique_sample[which(meta_data$CellId%in%cells_list[[x]])] <-names(cells_list[x])
}

# going to be working with only 10X from now on, so use mtx_10x,
# convert to gene names, and force uniqueness on each rowname. 

rownames(mtx_10x)<-genes_10x$V2[match(rownames(mtx_10x),genes_10x$V1)]
rownames(mtx_10x) <- make.unique(rownames(mtx_10x))

process_10x =function(mtx, meta_data){
  srat=CreateSeuratObject(mtx, meta.data = meta_data)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > minNumGenes & nCount_RNA > minNumUMI & percent.mt < mtCut)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
  keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = defNumPCs)
  srat = FindNeighbors(srat, dims=1:defNumPCs)
  srat = FindClusters(srat, resolution = clusterRes)
  srat = RunUMAP(srat, dims=1:defNumPCs, min.dist = finalMinDist, n.neighbors = finalNN)
  return(srat)
}

fresh_biopsies =unique(meta_data$unique_sample)[grep("Fresh",unique(meta_data$unique_sample))]
organoids =unique(meta_data$unique_sample)[grep("organoids",unique(meta_data$unique_sample))]
organoids=organoids[1:5]
tumour_only_mtx =mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%fresh_biopsies)]]
tumour_only_meta =meta_data[which(meta_data$unique_sample%in%fresh_biopsies),]

org_only_mtx=mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%organoids)]]
org_only_meta=meta_data[which(meta_data$unique_sample%in%organoids),]

srat_tumour =process_10x(tumour_only_mtx,tumour_only_meta)
srat_org=process_10x(org_only_mtx, org_only_meta)


srat_tumour@meta.data$new_idents="x"
srat_tumour@meta.data$idents_for_plot="x"

srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(22,8))]="tumour"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(13))]="schwannian stroma"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(0,2,11,14,15,18,25,10,12,20))]="mesenchyme"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(1,17,19,21,16,5,9,4,6,26,23,3,24))]="immune"
srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$seurat_clusters%in%c(7))]="vascular endothelium"

srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(8))]="tumour cluster 1"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(22))]="tumour cluster 2"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(13))]="schwannian stroma"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(0,2,11,14,15,18,25,10,12,20))]="mesenchyme"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(1,17,19,21,16,5,9,4,6,26,23,3,24))]="immune"
srat_tumour@meta.data$idents_for_plot[which(srat_tumour@meta.data$seurat_clusters%in%c(7))]="vascular endothelium"


#make UMAP pretty

bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
srat_tumour@meta.data$idents_for_plot=factor(srat_tumour@meta.data$idents_for_plot, levels=c("tumour cluster 1", "tumour cluster 2",
                                                                                            "schwannian stroma","vascular endothelium","mesenchyme", "immune"))

DimPlot(srat_tumour, group.by = "idents_for_plot",label=T,
        cols = c("#084594", "#4292C6",viridis(4)[1], "#CC6677", "#7E481C","#DDCC77"),
        pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), 
                              legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")

freq_tab=as.data.frame(table(srat_tumour@meta.data$unique_sample[which(srat_tumour@meta.data$new_idents=="tumour")]))
freq_tab$type="tumour"

freq_tab2=as.data.frame(table(srat_tumour@meta.data$unique_sample[which(srat_tumour@meta.data$new_idents=="mesenchyme")]))
freq_tab2$type="mesenchyme"

freq_tab_all=rbind(freq_tab, freq_tab2)
freq_tab_all$Var1=factor(freq_tab_all$Var1, levels = ord$Var1)
freq_tab_all$type=factor(freq_tab_all$type, levels=c("tumour", "mesenchyme"))
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
ggplot(freq_tab_all, aes(x=type, y=Freq, fill=Var1)) +geom_bar(position="fill", stat="identity") + scale_fill_manual(values=tol14rainbow) + theme(panel.background = element_blank())
sort(freq_tab$Freq, decreasing = T)

adr_nocc=subset(adr, cells=rownames(adr@meta.data)[which(adr@meta.data$Phase=="G1")])
adr_quickmarkers=quickMarkers(adr@assays$RNA@counts, as.character(adr@active.ident), N=30)
scp_markers=adr_quickmarkers[which(adr_quickmarkers$cluster%in%c("17","11","7")),]

adr.markers <- FindAllMarkers(adr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
adr.markers_top=adr.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
scp.markers=adr.markers_top[which(adr.markers_top$cluster%in%c("17","11","7")),]

adr_scp_only=subset(adr, idents = c(17,11,7))
scp_only_quickmarkers=quickMarkers(adr_scp_only@assays$RNA@counts, adr_scp_only@active.ident, N=30)


adr_new@meta.data$new_clust= as.character(adr_new@meta.data$seurat_clusters)

adr_new@meta.data$new_clust[which(adr_new@meta.data$new_clust%in%c("17", "10"))]="3_SCPs"
adr_new@meta.data$new_clust[which(adr_new@meta.data$new_clust%in%c("19", "7"))]="1_neuroblasts"
adr_new@meta.data$new_clust[which(adr_new@meta.data$new_clust%in%c("11"))]="2_bridge"

tr2mark=quickMarkers(adr_new@assays$RNA@counts, adr_new@meta.data$new_clust, N=30)
tr2mark_rel=tr2mark[which(tr2mark$cluster%in%c("3_SCPs", "1_neuroblasts", "2_bridge")),]

final_list=tr2mark_rel[which(tr2mark_rel$gene%in%unique_genes),]

just_tum_cells=subset(srat_tumour, idents = c(9,12,15,23))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[final_list$gene,]
final_list$cluster=factor(final_list$cluster, levels = c("3_SCPs", "2_bridge", "1_neuroblasts"))
col_fun=col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
names(tol14rainbow)= freq_tab$Var1
just_tum_cells@meta.data$unique_sample=factor(just_tum_cells@meta.data$unique_sample, levels = as.character(freq_tab$Var1[order(freq_tab$Freq)]))
ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = tol14rainbow))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = final_list$cluster, col=col_fun,bottom_annotation = ha)

#now do the same for organoids and in-house NB. 
srat_org=FindClusters(srat_org, resolution = 0.5)
srat_org@meta.data$new_idents="x"
srat_org@meta.data$idents_for_plot="x"

srat_org@meta.data$new_idents[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(0,1,2,4,5,6,7,8,9,10))]="tumour"
srat_org@meta.data$new_idents[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(3))]="mesenchyme"


srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(0))]="tumour cluster 1"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(1))]="tumour cluster 2"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(2))]="tumour cluster 3"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(4))]="tumour cluster 4"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(5))]="tumour cluster 5"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(6))]="tumour cluster 6"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(7))]="tumour cluster 7"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(8))]="tumour cluster 8"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(9))]="tumour cluster 9"
srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(10))]="tumour cluster 10"

srat_org@meta.data$idents_for_plot[which(srat_org@meta.data$RNA_snn_res.0.5%in%c(3))]="mesenchyme"




DimPlot(srat_org, group.by = "idents_for_plot")

bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")

colfunc <- colorRampPalette(c("#084594", "#C6DBEF"))
colfunc(8)
srat_org@meta.data$idents_for_plot=factor(srat_org@meta.data$idents_for_plot,
                                          levels=c("tumour cluster 1", "tumour cluster 2",
                                                   "tumour cluster 3", "tumour cluster 4",
                                                   "tumour cluster 5", "tumour cluster 6",
                                                   "tumour cluster 7", "tumour cluster 8",
                                                   "tumour cluster 9", "tumour cluster 10",
                                                   "mesenchyme"))


freq_tab=as.data.frame(table(srat_org@meta.data$unique_sample[which(srat_org@meta.data$new_idents=="tumour")]))
freq_tab$type="tumour"

freq_tab2=as.data.frame(table(srat_org@meta.data$unique_sample[which(srat_org@meta.data$new_idents=="mesenchyme")]))
freq_tab2$type="mesenchyme"

freq_tab_all=rbind(freq_tab, freq_tab2)
freq_tab_all$Var1=factor(freq_tab_all$Var1, levels = rev(as.character(freq_tab$Var1)))
freq_tab_all$type=factor(freq_tab_all$type, levels=c("tumour", "mesenchyme"))
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
ggplot(freq_tab_all, aes(x=type, y=Freq, fill=Var1)) +geom_bar(position="fill", stat="identity") + scale_fill_manual(values=tol5qualitative) + theme(panel.background = element_blank())

just_tum_cells=subset(srat_org, idents = c(0,2,4,6,8,11,12,13))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[final_list$gene,]
final_list$cluster=factor(final_list$cluster, levels = c("3_SCPs", "2_bridge", "1_neuroblasts"))
col_fun=col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
names(tol5qualitative)= as.character(freq_tab$Var1[order(freq_tab$Freq)])
just_tum_cells@meta.data$unique_sample=factor(just_tum_cells@meta.data$unique_sample, levels = as.character(freq_tab$Var1[order(freq_tab$Freq)]))
ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = tol5qualitative))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = final_list$cluster, col=col_fun,bottom_annotation = ha)

sum(freq_tab$Freq)
sum(freq_tab2$Freq)

#now do the same for in-house NB samples

all_paths = c("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685341/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685342/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843576/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843577/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843578/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733084/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733085/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733086/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004894/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004902/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004910/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787237/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787238/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787239/outs/filtered_gene_bc_matrices/GRCh38/")
names(all_paths) = c("STDY7685340","STDY7685341","STDY7685342","STDY7843576", "STDY7843577", "STDY7843578","STDY7733084",
                     "STDY7733085", "STDY7733086", "STDY8004894",
                     "STDY8004902", "STDY8004910", "STDY7787237", "STDY7787238", "STDY7787239")
inhouse_srats=Read10X(all_paths)

process_inhouse =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 300 & nCount_RNA > 1000 & percent.mt < mtCut)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
  keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 75)
  srat = FindNeighbors(srat, dims=1:75)
  srat = FindClusters(srat, resolution = clusterRes)
  srat = RunUMAP(srat, dims=1:75, min.dist = finalMinDist, n.neighbors = finalNN)
  return(srat)
}

srat_inhouse=process_inhouse(inhouse_srats)
DimPlot(srat_inhouse, group.by = "unique_samples")

srat_inhouse@meta.data$unique_samples="x"

srat_inhouse@meta.data$unique_samples[which(srat_inhouse@meta.data$orig.ident%in%c("STDY7685340","STDY7685341","STDY7685342"))]="NB1"
srat_inhouse@meta.data$unique_samples[which(srat_inhouse@meta.data$orig.ident%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="NB2"
srat_inhouse@meta.data$unique_samples[which(srat_inhouse@meta.data$orig.ident%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="NB3"
srat_inhouse@meta.data$unique_samples[which(srat_inhouse@meta.data$orig.ident%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="NB4"
srat_inhouse@meta.data$unique_samples[which(srat_inhouse@meta.data$orig.ident%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="GN"

srat_inhouse@meta.data$new_idents="x"
srat_inhouse@meta.data$idents_for_plot="x"

srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))]="tumour"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                         17,15,18,3,6,4,21,22,24))]="immune"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="vascular endothelium"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="mesenchyme"

srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25))]="tumour cluster 1"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(12))]="tumour cluster 2"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(5))]="tumour cluster 3"


srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                    17,15,18,3,6,4,21,22,24))]="immune"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="vascular endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="mesenchyme"
DimPlot(srat_inhouse,group.by = "idents_for_plot" )

srat_inhouse@meta.data$idents_for_plot=factor(srat_inhouse@meta.data$idents_for_plot,
                                              levels=c("tumour cluster 1", "tumour cluster 2",
                                                       "tumour cluster 3","mesenchyme",
                                                       "immune","vascular endothelium", "other"))



freq_tab=as.data.frame(table(srat_inhouse@meta.data$unique_sample[which(srat_inhouse@meta.data$new_idents=="tumour")]))
freq_tab$type="tumour"

freq_tab2=as.data.frame(table(srat_inhouse@meta.data$unique_sample[which(srat_inhouse@meta.data$new_idents=="mesenchyme")]))
freq_tab2$type="mesenchyme"

freq_tab_all=rbind(freq_tab, freq_tab2)
freq_tab_all$Var1=factor(freq_tab_all$Var1, levels = rev(as.character(freq_tab$Var1[order(freq_tab$Freq, decreasing = T)])))
freq_tab_all$type=factor(freq_tab_all$type, levels=c("tumour", "mesenchyme"))
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
ggplot(freq_tab_all, aes(x=type, y=Freq, fill=Var1)) +geom_bar(position="fill", stat="identity") + scale_fill_manual(values=tol5qualitative) + theme(panel.background = element_blank())

just_tum_cells=subset(srat_inhouse, idents = c(1,5,8,12,16,25))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[final_list$gene,]
final_list$cluster=factor(final_list$cluster, levels = c("3_SCPs", "2_bridge", "1_neuroblasts"))
col_fun=col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
names(tol5qualitative)= as.character(freq_tab$Var1[order(freq_tab$Freq)])
just_tum_cells@meta.data$unique_sample=factor(just_tum_cells@meta.data$unique_sample, levels = as.character(freq_tab$Var1[order(freq_tab$Freq)]))
ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = tol5qualitative))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = final_list$cluster, col=col_fun,bottom_annotation = ha)

sum(freq_tab$Freq)
sum(freq_tab2$Freq)
pdf("inhouse_umap.pdf", height = 7, width=7,useDingbats = F)
gg=DimPlot(srat_inhouse, group.by = "idents_for_plot", label=T,
        cols = c(colfunc(3),"#7E481C", "#DDCC77", "#CC6677", "grey" ),
        pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
pdf("org_umap.pdf", height = 7, width=7,useDingbats = F)
gg=DimPlot(srat_org, group.by = "idents_for_plot",
        cols = c(colfunc(10), "#7E481C"), label=T,
        pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank() ,legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()
pdf("dutch_umap.pdf", height = 7, width=7,useDingbats = F)
gg=DimPlot(srat_tumour, group.by = "idents_for_plot",label=T,
        cols = c("#084594", "#4292C6",viridis(4)[1], "#CC6677", "#7E481C","#DDCC77"),
        pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                              axis.ticks = element_blank(), 
                              legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()



saveRDS(srat_inhouse, "srat_inhouse.rds")
saveRDS(srat_org, "srat_org.rds")
saveRDS(srat_tumour, "srat_dutch.rds")
