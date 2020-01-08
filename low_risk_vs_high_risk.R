#define SCP, neuroblast and mature medulla genes

#combine baby adrenal 2 and bilateral adrenal
exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
bilat_srat=readRDS("/lustre/scratch117/casm/team274/gk14/bilateral_adrenal/bilateralAdrenals_Seurat.RDS")
Idents(bilat_srat)
just_medulla = subset(bilat_srat,idents = c(15,16,20))
just_medulla@meta.data$old_ident=just_medulla@meta.data$seurat_clusters
just_scp_traj=subset(adr, idents = c(19, 7,11,17,10))
just_scp_traj@meta.data$old_ident=just_scp_traj@meta.data$seurat_clusters
combined=merge(just_scp_traj,y=just_medulla, project="def_genes")


keep_features=rownames(combined)[which(!rownames(combined)%in%exclude_genes)]
combined=subset(combined, features = keep_features)
combined=NormalizeData(combined)
combined=FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined=ScaleData(combined, features = rownames(combined))
combined = CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
combined = RunPCA(combined, npcs = defNumPCs)
combined = FindNeighbors(combined, dims=1:defNumPCs)
combined = FindClusters(combined, resolution = clusterRes)
combined = RunUMAP(combined, dims=1:defNumPCs, min.dist = finalMinDist, n.neighbors = finalNN)


DimPlot(combined, group.by = "old_ident", label=T, label.size = 6)
combined@meta.data


DimPlot(combined, label=T, label.size = 8)

FeaturePlot(combined, features = c("SLC18A1", "ISL1", "HAND2", "PHOX2A", "PHOX2B", "CHRNA3", "SOX2", "SOX10"))
combined@meta.data
adr_quick_markers=quickMarkers(adr@assays$RNA@counts[,which(adr@meta.data$seurat_clusters%in%c(7,17))], 
                               as.character(adr@meta.data$seurat_clusters[which(adr@meta.data$seurat_clusters%in%c(7,17))]), N=100)
scp_and_neurobl_markers=adr_quick_markers[which(adr_quick_markers$cluster%in%c(7,17)),]
rownames(scp_and_neurobl_markers)= 1:nrow(scp_and_neurobl_markers)
markers_filtered=scp_and_neurobl_markers[which(scp_and_neurobl_markers$geneFrequencyOutsideCluster<0.1),]

just_19_and_7=subset(adr, idents = c(7,17))


dim(just_19_and_7)

mark =FindAllMarkers(just_19_and_7,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mark %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

#bulk NB 
norm_bulk_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/GTEX/GTEX_BulkMetadata.tsv", sep = "\t")

adr_samples=norm_bulk_meta[which(norm_bulk_meta$Tissue=="Adrenal Gland"),]$UniqueSampleID

gtex_bulk=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/GTEX/GTEX_BulkCountsFrag.tsv", sep="\t")


just_adr_counts=gtex_bulk[,adr_samples]
write.table(just_adr_counts, "useful_files/gtex_adrenals.tsv", sep = "\t", row.names = T, col.names = T, quote=F)
gtex_gene_length=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/GTEX/GTEX_BulkGeneLengths.tsv", sep="\t")
gtex_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/GTEX/GTEX_BulkGenes.tsv", sep = "\t", header = T)


colnames(gtex_gene_length)

adr_gene_length=gtex_gene_length[,adr_samples]
dim(adr_gene_length)


#TARGET bulk data

target_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv", header = T, sep= "\t")
target_gene_length= read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGeneLengths.tsv", sep = "\t")
target_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t", header = T)
target_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")

target_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv", sep = "\t", header = T)




colnames(target_counts)=sapply(strsplit(colnames(target_counts), "\\."), "[", 4)
target_clinical_data$TARGET.USI=sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 3)

target_length_norm =target_counts/(target_gene_length/1000)

scale_factor=colSums(target_length_norm)/1000000
target_tpm=target_length_norm/scale_factor
marker_gene_symbols=as.character(unique(markers_filtered$gene))
length(which(marker_gene_symbols%in%target_genes$external_gene_name))
target_genes[which(target_genes$external_gene_name%in%marker_gene_symbols),]

low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event=="Relapse")

scp_markers=as.character(filtered1$gene[which(filtered1$cluster==7)])
scp_markers=c(scp_markers, "PHOX2A", "PTPRC", "SOX10")
neurobl_markers=as.character(filtered1$gene[which(filtered1$cluster==17)])
neurobl_markers=c(neurobl_markers, "PHOX2A", "PTPRC", "SOX10")

rel_patients=c(low_risk_4s$TARGET.USI, high_risk_relapse$TARGET.USI)


scp_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%scp_markers),]$ensembl_gene_id)

neurobl_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%neurobl_markers),]$ensembl_gene_id)

scp_tab=target_tpm[scp_ensembl, rel_patients]
neurobl_tab=target_tpm[neurobl_ensembl, rel_patients]

rownames(scp_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(scp_tab))]
rownames(neurobl_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(neurobl_tab))]

scp_tab=as.data.frame(t(scp_tab))
neurobl_tab=as.data.frame(t(neurobl_tab))

scp_tab$patient=rownames(scp_tab)
neurobl_tab$patient=rownames(neurobl_tab)

scp_tab$risk="x"
scp_tab$risk[1:13]="low_risk_4s"
scp_tab$risk[14:52]="high_risk_relapse"

neurobl_tab$risk="x"
neurobl_tab$risk[1:13]="low_risk_4s"
neurobl_tab$risk[14:52]="high_risk_relapse"
scp_tab_melted$risk=factor(scp_tab_melted$risk, levels = c("low_risk_4s", "high_risk_relapse"))
scp_tab_melted=melt(scp_tab, id.vars=c("patient", "risk"))
neurobl_tab_melted=melt(neurobl_tab, id.vars=c("patient", "risk"))
neurobl_tab_melted$risk=factor(neurobl_tab_melted$risk, levels = c("low_risk_4s", "high_risk_relapse"))
scp_tab_melted$value=log10(scp_tab_melted$value)
neurobl_tab_melted$value=log10(neurobl_tab_melted$value)
ggplot(neurobl_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(scp_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


count2 =read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv", header = T, sep= "\t")
x=colnames(count2)[which(sapply(strsplit(colnames(count2), "\\."), "[", 4)%in%high_risk_relapse$TARGET.USI)]
x1=colnames(count2)[which(sapply(strsplit(colnames(count2), "\\."), "[", 4)%in%low_risk_4s$TARGET.USI)]
df1=data.frame(patient_id=x, risk="high_risk_relapse")
df2=data.frame(patient_id=x1, risk="low_risk_4s")
rel_df=rbind(df1, df2)

write.table(rel_df,"/lustre/scratch117/casm/team274/gk14/target_low_and_high_risk.txt", sep = "\t", col.names = T, row.names = F, quote=F)
