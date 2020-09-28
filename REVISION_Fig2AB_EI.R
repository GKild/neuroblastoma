library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
source("gene_sets.R")

s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

#all 10x Neuroblastoma samples
#list of paths to gene_bc_matrices: all_paths = c()
names(all_paths) = c("STDY7685340","STDY7685341","STDY7685342","STDY7843576", "STDY7843577", "STDY7843578","STDY7733084",
                     "STDY7733085", "STDY7733086", "STDY8004894",
                     "STDY8004902", "STDY8004910", "STDY7787237", "STDY7787238", "STDY7787239")
inhouse_srats=Read10X(all_paths)
#processing function
process_inhouse =function(mtx){
  srat=CreateSeuratObject(mtx)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 300 & nCount_RNA > 1000 & percent.mt < 30)
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
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:75, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}
# get 10x tumour Seurat object
srat_inhouse=process_inhouse(inhouse_srats)
#literature markers 
srat_inhouse@meta.data$PD_ID=as.character(srat_inhouse@meta.data$orig.ident)

srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="PD42184"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="PD42752-1"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="PD42752-2"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="PD46693"
srat_inhouse@meta.data$PD_ID[which(srat_inhouse@meta.data$PD_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="PD43255"
lit_genes=c("SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH", "DBH", "PHOX2B", 
            "CHGB","GAP43", "BCL2", "NPY", "PHOX2A", "PHOX2B","MYCN", "PNMT", "PDGFRB", "TCF21", "PLVAP", "PECAM1", "KDR",
            "PTPRB","HBG1", "HBG2", "HBB","STAR", "MC2R", "PTPRC")
FeaturePlot(srat_inhouse, lit_genes)

#annotate louvain clusters 

srat_inhouse@meta.data$new_idents="x"
srat_inhouse@meta.data$idents_for_plot="x"

srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))]="tumour"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                    17,15,18,3,6,4,21,22,24))]="leukocytes"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="endothelium"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="mesenchyme"

srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25))]="Tumour cluster 1"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(12))]="Tumour cluster 2"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(5))]="Tumour cluster 3"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                         17,15,18,3,6,4,21,22,24))]="Leukocytes"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="Endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="Mesenchyme"
srat_inhouse@meta.data$idents_for_plot=factor(srat_inhouse@meta.data$idents_for_plot, levels=c("Tumour cluster 1", "Tumour cluster 2",
                                                                                               "Tumour cluster 3","Mesenchyme",
                                                                                               "Leukocytes","Endothelium"))

#make fig2A
colfunc <- colorRampPalette(c("#084594", "#C6DBEF")) 

pdf("inhouse_umap.pdf", height = 3, width=4,useDingbats = F)
gg=DimPlot(srat_inhouse, group.by = "idents_for_plot", label=T,
           cols = c(colfunc(3),"#7E481C", "#DDCC77", "#CC6677"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()


#############################
### process CEL-seq2 data ###
load("scRNAseq_NB_PMC_Full.RData")
Jan_nb = data
genes_10x = read.table("genes.tsv", sep = "\t", header = F)
#main matrix 
mtx = Jan_nb$data$counts
srcs='extra_samples'
srcs = file.path(srcs,list.files(srcs))
srcs = srcs[grep(".txt", srcs)]
source("gene_sets.R")
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

full_obj=read.table(srcs[1], sep = "\t", header = T)
colnames(full_obj)[2:ncol(full_obj)]=paste0(sapply(strsplit(basename(srcs[1]), "\\."), "[", 1),"-", colnames(full_obj)[2:ncol(full_obj)])
for (x in srcs[2:13]){
  tab=read.table(x, sep = "\t", header = T)
  colnames(tab)[2:ncol(tab)]=paste0(sapply(strsplit(basename(x), "\\."), "[", 1),"-", colnames(tab)[2:ncol(tab)])
  full_obj=merge(full_obj, tab, by="GENEID", all=T)
}

full_obj=full_obj[,-grep("UNK", colnames(full_obj))]
full_obj=full_obj[-grep("cluster", full_obj$GENEID),]
full_obj$GENEID=sapply(strsplit(as.character(full_obj$GENEID), "__"),"[", 1)
full_obj$GENEID

old_df=as.data.frame(as.matrix(mtx))
old_df$GENEID=rownames(old_df)

f=merge(old_df, full_obj, by="GENEID", all=T)

rownames(f)=f$GENEID
f=f[,2:ncol(f)]

joined_mtx=as.matrix(f)
joined_mtx[is.na(joined_mtx)] <- 0
joined_mtx=Matrix(joined_mtx, sparse = T)

#remove ERCCs


joined_mtx = joined_mtx[-grep("^ERCC-",rownames(joined_mtx)),]
rownames(joined_mtx)=sapply(strsplit(rownames(joined_mtx), "__"), "[", 1)

mtx_10x =joined_mtx[which(rownames(joined_mtx)%in%genes_10x$V1),]

rownames(mtx_10x)<-genes_10x$V2[match(rownames(mtx_10x),genes_10x$V1)]
rownames(mtx_10x) <- make.unique(rownames(mtx_10x))
dim(mtx_10x)

patient_ids =unique(Jan_nb$col.annot$PatientId)
cells_each_sample =sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})


cells_list = list("NB039_Tumour_organoids_in_20perc_HP_3D"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
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
                  "NB152_Tumour_organoids_P1"=cells_each_sample$NB152$`Tumour organoids P1`, 
                  "000CGH_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK001", "KK002", "TM230", "TM231"))],
                  "000FQM_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK051", "KK052", "KK053", "KK054"))],
                  "000GGU_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK055", "KK056", "KK057", "KK058", "KK059"))])

meta_data=as.data.frame(colnames(mtx_10x))
colnames(meta_data)="CellId"
rownames(meta_data)=meta_data$CellId

meta_data$unique_sample <- "sample"
for (x in 1:length(cells_list)) {
  meta_data$unique_sample[which(meta_data$CellId%in%cells_list[[x]])] <-names(cells_list[x])
}

process_10x =function(mtx, meta_data){
  srat=CreateSeuratObject(mtx, meta.data = meta_data)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  #exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
  #keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  #srat=subset(srat, features=keep_features)
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


fresh_biopsies =unique(meta_data$unique_sample)[grep("Fresh",unique(meta_data$unique_sample))]
tumour_only_mtx =mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%fresh_biopsies)]]
tumour_only_meta =meta_data[which(meta_data$unique_sample%in%fresh_biopsies),]

srat_dutch =process_10x(tumour_only_mtx,tumour_only_meta)
FeaturePlot(srat_dutch, lit_genes)

#annotate louvain clusters

srat_dutch@meta.data$new_idents="x"
srat_dutch@meta.data$idents_for_plot="x"

srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(13))]="Tumour cluster 1"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(24))]="Tumour cluster 2"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(15))]="X1"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(10))]="Tumour cluster 3"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(17))]="X2"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(6))]="Schwannian stroma"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(7,29))]="Endothelium"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(28,18,23,5,26,0,8,20,9,3,22,25,4,19,21,27))]="Leukocytes"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(1,2,11,12,14,16))]="Mesenchyme"

srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(13,24,10))]="Tumour"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(6))]="Schwannian stroma"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(7,29))]="Endothelium"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(28,18,23,5,26,0,8,20,9,3,22,25,4,19,21,27))]="Leukocytes"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(1,2,11,12,14,16))]="Mesenchyme"

srat_dutch@meta.data$idents_for_plot=factor(srat_dutch@meta.data$idents_for_plot, levels=c(
  "Tumour cluster 1","Tumour cluster 2","Tumour cluster 3","X1", "X2",
  "Schwannian stroma","Endothelium", "Mesenchyme", "Leukocytes"
))

pdf("CEL-seq2_umap.pdf", height = 3, width=4,useDingbats = F)
gg=DimPlot(srat_dutch, group.by = "idents_for_plot",label=T,
           cols = c(colfunc(3),"A5A4A2","A5A4A2","#7c5295", "#CC6677", "#7E481C","#DDCC77"),
           pt.size = 0.1) +theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                 axis.ticks = element_blank(), 
                                 legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")


plot(gg)
dev.off()


######Logistic regresssion against fetal adrenal########
adr_all=readRDS("fAdrenal/processedSeurat.RDS")
fetal_pods=readRDS("fetal_podocytes_forGerda.RDS")
genes_10x = read.table("features.tsv.gz", sep = "\t", header = F)

fetal_pds_counts=fetal_pods@assays$RNA@counts
rownames(fetal_pds_counts)=genes_10x$V2
merged_cells=cbind(adr_all@assays$RNA@counts, fetal_pds_counts)
comm_genes=intersect(rownames(merged_cells), rownames(srat_dutch)) 
comm_genes=intersect(comm_genes, rownames(srat_inhouse))

cell_labs=c(as.character(adr_all@meta.data$new_clust), rep("Podocytes", 278))
#train logistic regression model
fit = trainModel(merged_cells[comm_genes, ], cell_labs, workers=NULL)

ps_dutch = predictSimilarity(fit, srat_dutch@assays$RNA@counts[comm_genes,],srat_dutch@meta.data$idents_for_plot)
ps_inhouse=predictSimilarity(fit, srat_inhouse@assays$RNA@counts[comm_genes,], srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.8)
ps_self=predictSimilarity(fit, merged_cells, cell_labs)


lr_col_ord=c("SCPs","Bridge","Chromaffin","Sympathoblastic","Mesenchyme","Cortex","Leukocytes", "Endothelium", "Erythroblasts","Other", "Podocytes")
pdf("dutch_lr.pdf", height = 4, width = 7)
print(similarityHeatmap(ps_dutch, column_order=lr_col_ord,row_order=as.character(c("Tumour cluster 1", "Tumour cluster 2",
                                                                                   "Tumour cluster 3", "X1", "X2", 
                                                                                   "Schwannian stroma",
                                                                                   "Mesenchyme", "Leukocytes", "Endothelium"))))
dev.off()

pdf("inhouse_lr.pdf", height = 4, width = 7)
print(similarityHeatmap(ps_inhouse, column_order=lr_col_ord, row_order=c("Tumour cluster 1", "Tumour cluster 2",
                                                                         "Tumour cluster 3","Mesenchyme",
                                                                         "Leukocytes", "Endothelium")))
dev.off()


####symp, pos and neg ctrl signal######

dutch_tum_only=subset(srat_dutch, cells=rownames(srat_dutch@meta.data)[grep("Tumour", srat_dutch@meta.data$idents_for_plot)])

inhouse_tum_only=subset(srat_inhouse, cells=rownames(srat_inhouse@meta.data)[grep("Tumour", srat_inhouse@meta.data$idents_for_plot)])

dutch_tum_only@meta.data$risk="high risk"
dutch_tum_only@meta.data$risk[which(dutch_tum_only@meta.data$unique_sample%in%c("NB152_Fresh_biopsy"))]="low/intermediate risk"

inhouse_tum_only@meta.data$risk="low/intermediate risk"
inhouse_tum_only@meta.data$risk[which(inhouse_tum_only@meta.data$PD_ID=="PD43255")]="high risk"

tum_risk_labs=c(dutch_tum_only@meta.data$risk, inhouse_tum_only@meta.data$risk)
tum_genes=intersect(rownames(dutch_tum_only), rownames(inhouse_tum_only))
tum_together=cbind(dutch_tum_only@assays$RNA@counts[tum_genes, ], inhouse_tum_only@assays$RNA@counts[tum_genes, ])
tum_together_labs=c(paste0("PMC ", dutch_tum_only@meta.data$idents_for_plot),
                    paste0("GOSH ", inhouse_tum_only@meta.data$idents_for_plot))

dutch_no_15=subset(srat_dutch, idents = c(15,17), invert=T)

dutch_inh_together=cbind(srat_inhouse@assays$RNA@counts[tum_genes, ], dutch_no_15@assays$RNA@counts[tum_genes,])
dutch_inh_together_labs=c(as.character(srat_inhouse@meta.data$idents_for_plot), as.character(dutch_no_15@meta.data$idents_for_plot))

dutch_inh_together_labs[grep("Tumour", dutch_inh_together_labs)]="Tumour"

dutch_inh_together_train=trainModel(dutch_inh_together, dutch_inh_together_labs, workers=NULL)


together_against_tum=predictSimilarity(dutch_inh_together_train, tum_together, logits = F)
adr_against_tum=predictSimilarity(fit, tum_together)
adr_against_tum2=adr_against_tum[,c(2,4,7,5,8,10,1)]
similarityHeatmap(adr_against_tum2,cluster_rows=T, row_split=factor(tum_risk_labs, levels=c("low/intermediate risk", "high risk")))
ps_train_podos_tums=predictSimilarity(fit,tum_together, logits = F)
high_risk_symps=ps_train_podos_tums[rownames(ps_train_podos_tums)[which(tum_risk_labs=="high risk")],7]
high_risk_only=ps_train_podos_tums[rownames(ps_train_podos_tums)[which(tum_risk_labs=="high risk")],]

kmeans_2=kmeans(high_risk_symps,centers = 2)
similarityHeatmap(high_risk_only, row_split=kmeans_2$cluster)

qc_df=data.frame(cells=names(kmeans_2$cluster), kmeans=kmeans_2$cluster)

qc_df$symp_cat="high_symp"
qc_df$symp_cat[which(qc_df$kmeans==2)]="low_symp"
qc_df$cohort="PMC"
qc_df$cohort[grep("STDY", qc_df$cells)]="GOSH"

qc_df$symp_and_cohort=paste0(qc_df$symp_cat, "_", qc_df$cohort)

qc_df$ngenes=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$nFeature_RNA,
               srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$nFeature_RNA)
qc_df$mito=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$percent.mt,
             srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$percent.mt)
qc_df$reads=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$nCount_RNA,
              srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$nCount_RNA)

ggplot(qc_df, aes(x=symp_cat, y=log10(ngenes))) +geom_boxplot()
ggplot(qc_df, aes(x=symp_cat, y=log10(reads))) +geom_boxplot()
ggplot(qc_df, aes(x=symp_cat, y=mito)) +geom_boxplot()

just_high_mtx=tum_together[,rownames(qc_df)]
excludeGenes=c(riboGenes, mtGenes, riboGenes)
keep_genes=which(!rownames(just_high_mtx)%in%excludeGenes)
just_high_srat=CreateSeuratObject(just_high_mtx)
just_high_srat=subset(just_high_srat, features=keep_genes)
just_high_srat=NormalizeData(just_high_srat)
just_high_srat=ScaleData(just_high_srat)
just_high_srat=FindVariableFeatures(just_high_srat)
just_high_srat=RunPCA(just_high_srat, npcs = 50)
just_high_srat = FindNeighbors(just_high_srat, dims = 1:50)
just_high_srat = FindClusters(just_high_srat, resolution = 1)
just_high_srat=RunUMAP(just_high_srat,dims=1:50, min.dist = 0.5, n.neighbors = 50)
symp_sg=qc_df$symp_cat
names(symp_sg)=names(just_high_srat@active.ident)
just_high_srat@active.ident=factor(symp_sg)
high_low_symp_de=FindMarkers(just_high_srat, ident.1 = "high_symp", ident.2 = "low_symp")
high_low_symp_de=high_low_symp_de[which(high_low_symp_de$p_val_adj<0.01),]
de_dediff_annot=read.csv("de_dediff.csv")
rownames(de_dediff_annot)=de_dediff_annot$Gene
high_low_symp_de$gene=rownames(high_low_symp_de)
high_low_symp_de=high_low_symp_de[,c(6,1:5)]
de_dediff_annot=de_dediff_annot[which(de_dediff_annot$Gene%in%rownames(high_low_symp_de)[1:50]),]
bla=merge(high_low_symp_de, de_dediff_annot,by=0, all=TRUE)
x1=bla[order(bla$p_val_adj, decreasing = F),]
rownames(x1)=1:nrow(x1)

x1=x1[,c(2:7,9:13)]
length(which(de_dediff_annot$Gene%in%rownames(high_low_symp_de)[1:50]))

write.table(x1, "high_low_symp_de.txt", sep = "\t", row.names = F, col.names = T, quote = F)

adr_against_tum=predictSimilarity(fit, tum_together, logits = F)
tum_sig_df=data.frame(positive_control=together_against_tum[,3], negative_control=together_against_tum[,5], sympathoblast=adr_against_tum[,7])
tum_sig_df=melt(tum_sig_df)



ggplot(data = tum_sig_df, aes(x=variable,y=value)) +
  geom_quasirandom(shape=19, alpha=0.05, size=0.2)+
  stat_summary(mapping = aes(x = variable, y = value),
               geom="pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, shape=22, fill="black",
               position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))

#make figs H/I 

inhouse_genotyping=read.delim("genotyping_metadata_Sanger.tsv", header = T)
dutch_genotyping=read.delim("genotyping_metadata_PMC.tsv", header = T)

inhouse_genotyping_tum=inhouse_genotyping[which(inhouse_genotyping$new_idents=="tumour"),]
dutch_genotyping_tum=dutch_genotyping[which(dutch_genotyping$idents_for_plot%in%c("Tumour cluster 1", "Tumour cluster 2",
                                                                                  "Tumour cluster 4")),]

sapply(as.character(unique(inhouse_genotyping_tum$scID)), function(x){
  bla=inhouse_genotyping_tum[which(inhouse_genotyping_tum$scID==x),]
  table(bla$call)
})

sapply(as.character(unique(dutch_genotyping_tum$unique_sample)), function(x){
  bla=dutch_genotyping_tum[which(dutch_genotyping_tum$unique_sample==x),]
  table(bla$call)
})

sample_labs=c(dutch_tum_only@meta.data$unique_sample, inhouse_tum_only@meta.data$PD_ID)
tum_ps=predictSimilarity(fit, tum_together, sapply(strsplit(sample_labs, "_"), "[",1))


similarityHeatmap(tum_ps, column_order=lr_col_ord)

#get markers for fig3 A-C
quick_with_podos=quickMarkers(merged_cells, cell_labs, FDR=0.01, N=Inf)
quick_with_podos_rel=quick_with_podos[which(quick_with_podos$cluster%in%c("SCPs","Bridge","Sympathoblastic", "Chromaffin","Podocytes")),]
podo_marks_filt=dplyr::filter(quick_with_podos_rel, tfidf>1 & geneFrequencySecondBest <0.2)




