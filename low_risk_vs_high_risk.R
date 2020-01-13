#define SCP, neuroblast and mature medulla genes

#combine baby adrenal 2 and bilateral adrenal
exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
adr=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")




just_19_and_7=subset(adr, idents = c(7,17))


dim(just_19_and_7)

mark =FindAllMarkers(just_19_and_7,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mark =mark %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

markers_filtered =mark[mark$pct.2<0.1,]
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

scp_markers=as.character(markers_filtered$gene[which(markers_filtered$cluster==7)])
scp_markers=c(scp_markers, "PHOX2A", "PTPRC", "SOX10")
neurobl_markers=as.character(markers_filtered$gene[which(markers_filtered$cluster==17)])
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


# German cohort 
german_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_series_matrix.txt", header = T, sep= "\t", skip = 63, nrows =44708)

german_genes=read.delim("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_family_soft.txt", sep = "\t", header = T, skip=471, nrows = 44708, na.strings = "")


german_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/bulk_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)

german_sample_ids=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_series_matrix.txt", header = F, skip=26, nrows = 1)

x=as.data.frame(t(german_sample_ids))

colnames(german_counts)[2:395]=as.character(x$V1[2:395])


no_na_gs=german_genes[!is.na(german_genes$GeneSymbol),]
no_na_gs$ID
counts_gs=german_counts[which(german_counts$ID_REF%in%no_na_gs$ID),] 

counts_gs$gene_symbol=no_na_gs$GeneSymbol 

#get correct patients 
german_clinical_data=german_clinical_data[which(german_clinical_data$Tumor.ID%in%colnames(counts_gs)),]
german_clinical_data$Age.at.Diagnosis..d.=as.numeric(as.character(german_clinical_data$Age.at.Diagnosis..d.))
low_risk=dplyr::filter(german_clinical_data,
                       german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &german_clinical_data$Status!="relapse/progression")
high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified")

risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")

low_risk$risk_level="low_risk"
high_risk$risk_level="high_risk"
risk_4s$risk_level="4s"


relevant_patients=rbind(low_risk, high_risk, risk_4s)


scp_tab=counts_gs[which(counts_gs$gene_symbol%in%scp_markers),c(relevant_patients$Tumor.ID, "gene_symbol")]
neurobl_tab=counts_gs[which(counts_gs$gene_symbol%in%neurobl_mesen),c(relevant_patients$Tumor.ID, "gene_symbol")]

scp_tab$gene_symbol = make.unique(as.character(scp_tab$gene_symbol))
neurobl_tab$gene_symbol=make.unique(as.character(neurobl_tab$gene_symbol))
rownames(scp_tab)=scp_tab$gene_symbol
rownames(neurobl_tab)=neurobl_tab$gene_symbol
scp_tab=scp_tab[,1:210]
neurobl_tab=neurobl_tab[,1:210]
scp_transposed=as.data.frame(t(scp_tab))
neurobl_transposed=as.data.frame(t(neurobl_tab))

scp_transposed$patient=rownames(scp_transposed)
neurobl_transposed$patient=rownames(neurobl_transposed)


scp_transposed$risk= relevant_patients$risk_level
neurobl_transposed$risk=relevant_patients$risk_level
is.factor(scp_transposed$PTPRC)
scp_melted=melt(scp_transposed, id.vars=c("patient", "risk"))
neurobl_melted=melt(neurobl_transposed,id.vars=c("patient", "risk"))

scp_melted$value=log10(scp_melted$value)
neurobl_melted$value=log10(neurobl_melted$value)

scp_melted$risk=factor(scp_melted$risk, levels = c("low_risk", "4s", "high_risk",))
neurobl_melted$risk=factor(neurobl_melted$risk, levels = c("low_risk", "4s", "high_risk"))  
ggplot(neurobl_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(scp_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


