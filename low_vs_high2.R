bilateral=readRDS("/lustre/scratch117/casm/team274/gk14/bilateral_adrenal/bilateralAdrenals_Seurat.RDS")
mature_count_table =read.table("mature_adrenal/matureAdrenalNeuronalCounts.tsv", header = T)
mature_metadata = read.table("mature_adrenal/matureAdrenalNeuronalMetadata.tsv", header = T)
#make it into seurat object
mature_matrix = as.matrix(mature_count_table)
colnames(mature_matrix) = rownames(mature_metadata)
mature_srat=CreateSeuratObject(mature_matrix, meta.data = mature_metadata)

mature_small_colnames=read.table("mature_adrenal/super_filtered/matureAdrenalRef_columnNames.tsv")

mature_small_srat=subset(mature_srat,cells=gsub(".*X","", mature_small_colnames$V1) )

mature_small_srat@meta.data$seurat_clusters=1
DimPlot(bilateral, label=T)
match(rownames(bilateral),rownames(mature_small_counts))
dim(scps_and_neurobl)
FeaturePlot(bilateral, c("PHOX2A", "SOX10", "ERBB3", "MPZ", "PHOX2B", "ISL1"))

FeaturePlot(adr, c("PHOX2A", "SOX10", "ERBB3", "MPZ", "PHOX2B", "SOX2"))

DimPlot(bilateral, group.by = "Phase")

scps_and_neurobl=subset(bilateral, idents = c(12,22))


no_cc_cells=rownames(scps_and_neurobl@meta.data)[which(scps_and_neurobl@meta.data$Phase=="G1")]

scps_and_neurobl=subset(scps_and_neurobl, cells = no_cc_cells)

scps_and_neurobl@meta.data$seurat_clusters


merged_srat=merge(scps_and_neurobl, y=mature_small_srat)

all_markers=quickMarkers(merged_srat@assays$RNA@counts, merged_srat@meta.data$seurat_clusters, N=100)

markers_filtered=dplyr::filter(all_markers, geneFrequency>0.5&geneFrequencyOutsideCluster<0.05)

scp_markers=as.character(markers_filtered$gene[which(markers_filtered$cluster==22)])

neurobl_markers=as.character(markers_filtered$gene[which(markers_filtered$cluster==12)])

mature_markers=as.character(markers_filtered$gene[which(markers_filtered$cluster==1)])


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

low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event=="Relapse")

rel_patients=c(low_risk_4s$TARGET.USI, high_risk_relapse$TARGET.USI)


scp_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%scp_markers),]$ensembl_gene_id)

neurobl_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%neurobl_markers),]$ensembl_gene_id)

mature_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%mature_markers),]$ensembl_gene_id)

scp_tab=target_tpm[scp_ensembl, rel_patients]
neurobl_tab=target_tpm[neurobl_ensembl, rel_patients]
mature_tab=target_tpm[mature_ensembl, rel_patients]

rownames(scp_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(scp_tab))]
rownames(neurobl_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(neurobl_tab))]
rownames(mature_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(mature_tab))]

scp_tab=as.data.frame(t(scp_tab))
neurobl_tab=as.data.frame(t(neurobl_tab))
mature_tab=as.data.frame(t(mature_tab))

scp_tab$patient=rownames(scp_tab)
neurobl_tab$patient=rownames(neurobl_tab)
mature_tab$patient=rownames(mature_tab)


scp_tab$risk="x"
scp_tab$risk[1:13]="low_risk_4s"
scp_tab$risk[14:52]="high_risk_relapse"

neurobl_tab$risk="x"
neurobl_tab$risk[1:13]="low_risk_4s"
neurobl_tab$risk[14:52]="high_risk_relapse"

mature_tab$risk="x"
mature_tab$risk[1:13]="low_risk_4s"
mature_tab$risk[14:52]="high_risk_relapse"


scp_tab_melted=melt(scp_tab, id.vars=c("patient", "risk"))
scp_tab_melted$risk=factor(scp_tab_melted$risk, levels = c("low_risk_4s", "high_risk_relapse"))

neurobl_tab_melted=melt(neurobl_tab, id.vars=c("patient", "risk"))
neurobl_tab_melted$risk=factor(neurobl_tab_melted$risk, levels = c("low_risk_4s", "high_risk_relapse"))

mature_tab_melted=melt(mature_tab, id.vars=c("patient", "risk"))
mature_tab_melted$risk=factor(mature_tab_melted$risk, levels = c("low_risk_4s", "high_risk_relapse"))

scp_tab_melted$value=log10(scp_tab_melted$value)
neurobl_tab_melted$value=log10(neurobl_tab_melted$value)

mature_tab_melted$value=log10(mature_tab_melted$value)

scp_tab_melted$variable=factor(scp_tab_melted$variable, levels =  as.character(unique(scp_tab_melted$variable)[match(scp_markers,unique(scp_tab_melted$variable))]))
neurobl_tab_melted$variable=factor(neurobl_tab_melted$variable, levels =  as.character(unique(neurobl_tab_melted$variable)[match(neurobl_markers,unique(neurobl_tab_melted$variable))]))
mature_tab_melted$variable=factor(mature_tab_melted$variable, levels =  as.character(unique(mature_tab_melted$variable)[match(mature_markers,unique(mature_tab_melted$variable))]))


ggplot(neurobl_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(scp_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mature_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


scp_tab=counts_gs[which(counts_gs$gene_symbol%in%scp_markers),c(relevant_patients$Tumor.ID, "gene_symbol")]
neurobl_tab=counts_gs[which(counts_gs$gene_symbol%in%neurobl_markers),c(relevant_patients$Tumor.ID, "gene_symbol")]
mature_tab=counts_gs[which(counts_gs$gene_symbol%in%mature_markers),c(relevant_patients$Tumor.ID, "gene_symbol")]

scp_tab$gene_symbol = make.unique(as.character(scp_tab$gene_symbol))
neurobl_tab$gene_symbol=make.unique(as.character(neurobl_tab$gene_symbol))
mature_tab$gene_symbol=make.unique(as.character(mature_tab$gene_symbol))

rownames(scp_tab)=scp_tab$gene_symbol
rownames(neurobl_tab)=neurobl_tab$gene_symbol
rownames(mature_tab)=mature_tab$gene_symbol

scp_tab=scp_tab[,1:210]
neurobl_tab=neurobl_tab[,1:210]
mature_tab=mature_tab[,1:210]

scp_transposed=as.data.frame(t(scp_tab))
neurobl_transposed=as.data.frame(t(neurobl_tab))
mature_transposed=as.data.frame(t(mature_tab))

scp_transposed$patient=rownames(scp_transposed)
neurobl_transposed$patient=rownames(neurobl_transposed)
mature_transposed$patient=rownames(mature_transposed)


scp_transposed$risk= relevant_patients$risk_level
neurobl_transposed$risk=relevant_patients$risk_level
mature_transposed$risk=relevant_patients$risk_level

scp_melted=melt(scp_transposed, id.vars=c("patient", "risk"))
neurobl_melted=melt(neurobl_transposed,id.vars=c("patient", "risk"))
mature_melted=melt(mature_transposed,id.vars=c("patient", "risk"))

scp_melted$value=log10(scp_melted$value)
neurobl_melted$value=log10(neurobl_melted$value)
mature_melted$value=log10(mature_melted$value)

scp_melted$risk=factor(scp_melted$risk, levels = c("low_risk", "4s", "high_risk"))
neurobl_melted$risk=factor(neurobl_melted$risk, levels = c("low_risk", "4s", "high_risk"))
mature_melted$risk=factor(mature_melted$risk, levels = c("low_risk", "4s", "high_risk"))

scp_melted$variable=factor(scp_melted$variable, levels =  as.character(unique(scp_melted$variable)[match(scp_markers,unique(scp_melted$variable))]))
neurobl_melted$variable=factor(neurobl_melted$variable, levels =  as.character(unique(neurobl_melted$variable)[match(neurobl_markers,unique(neurobl_melted$variable))]))
neurobl_melted2=dplyr::filter(neurobl_melted, !is.na(variable))
mature_melted$variable=factor(mature_melted$variable, levels =  as.character(unique(mature_melted$variable)[match(mature_markers,unique(mature_melted$variable))]))

ggplot(neurobl_melted2, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(scp_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(mature_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))




