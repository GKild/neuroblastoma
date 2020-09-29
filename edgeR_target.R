#prep the target data
target_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv",
                         header = T, sep= "\t")
target_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t",
                        header = T)
target_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")
target_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv",
                                sep = "\t", header = T)
#sort out clinical data NAs
names(which(table(sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 3))==2))

df_bla=target_clinical_data[grep("PARBAJ|PARHAM|PASGAP|PASNPG|PATYIL|PAUDDK",target_clinical_data$TARGET.USI),]

target_clinical_data[39,2:ncol(target_clinical_data)]=target_clinical_data[40,2:ncol(target_clinical_data)]
target_clinical_data[42,2:ncol(target_clinical_data)]=target_clinical_data[78,2:ncol(target_clinical_data)]
target_clinical_data[76,2:ncol(target_clinical_data)]=target_clinical_data[125,2:ncol(target_clinical_data)]
target_clinical_data[96,2:ncol(target_clinical_data)]=target_clinical_data[41,2:ncol(target_clinical_data)]
target_clinical_data[102,2:ncol(target_clinical_data)]=target_clinical_data[80,2:ncol(target_clinical_data)]
target_clinical_data[142,2:ncol(target_clinical_data)]=target_clinical_data[1,2:ncol(target_clinical_data)]

#unify the IDs in some capacity
colnames(target_counts)=paste0(sapply(strsplit(colnames(target_counts), "\\."), "[", 4), "-", sapply(strsplit(colnames(target_counts),
                                                                                                              "\\."), "[", 5))
target_clinical_data$TARGET_SHORT=gsub("NA", "01A",paste0(sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 3),
                                                          "-", sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 4)))


#define risk groups
low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk_relapse"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"
all=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
extremes=rbind(low_risk_4s, high_risk_relapse)

#get tabs all in the same order 

count_tab=target_counts[,all$TARGET_SHORT]

extreme_count_tab=target_counts[,extremes$TARGET_SHORT]


edg_file=DGEList(extreme_count_tab, samples = extremes, genes = target_genes, group = extremes$risk)
edg_file$common.dispersion=0.4
o <- order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file <- edg_file[o,]
d <- duplicated(edg_file$genes$external_gene_name)
edg_file <- edg_file[!d,]
nrow(edg_file)


keep <- filterByExpr(edg_file)
edg_file <- edg_file[keep, , keep.lib.sizes=FALSE]

rownames(edg_file)=edg_file$genes$external_gene_name[match(rownames(edg_file), edg_file$genes$ensembl_gene_id)]
rownames(edg_file$genes)=edg_file$genes$external_gene_name

edg_file <- calcNormFactors(edg_file)


edg_file$samples$risk=factor(edg_file$samples$risk, levels = c("low_risk_4s", "high_risk_relapse"))
design=model.matrix(~risk+Age.at.Diagnosis.in.Days, data = edg_file$samples)


fit <- glmQLFit(edg_file, design)
qlf <- glmQLFTest(fit,coef=2)

qlf$table[c("ERBB3", "SOX2", "MPZ", "PLP1"),]
et_man_disp <- exactTest(edg_file)

et$table[which(et$genes$external_gene_name%in%c("ERBB3", "SOX2", "MPZ", "PLP1")),]
et_man_disp$table[which(et_man_disp$genes$external_gene_name%in%c("ERBB3", "SOX2", "MPZ", "PLP1")),]
target_genes[which(target_genes$external_gene_name=="PLP1"),]

all_de=et_man_disp$table
topTags(et_man_disp)
all_de$gene=rownames(all_de)
all_de$qvalue=p.adjust(all_de$PValue, method = "BH")
de_filtered=dplyr::filter(all_de,  qvalue<0.01)

de_filtered[which(de_filtered$gene%in%c("ERBB3", "SOX2", "MPZ", "PLP1", "TMEM176A", "TMEM176B", "FOXD3-AS1", "SCRG1")),]

get_target_marker_plots(target_tpm,rownames(topTags(et_man_disp)),rownames(topTags(et_man_disp)),rownames(topTags(et_man_disp)), "bulk_df")


pch <- c(0,1)
colors <- c("darkgreen", "red")
plotMDS(edg_file, col=colors[edg_file$samples$group], pch=pch[edg_file$samples$group])
legend("topleft", legend=levels(edg_file$samples$group), pch=pch, col=colors, ncol=2)

DimPlot(adr_all_med, group.by = "devTime")

FeaturePlot(adr_all_med, c("MYCN", "ALK", "RET", "TERT", "ATRX"))
