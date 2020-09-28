####################prep the target data####################
target_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv",
                         header = T, sep= "\t")
target_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t",
                        header = T)
target_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")
target_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv",
                                sep = "\t", header = T)
target_gene_length= read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGeneLengths.tsv", sep = "\t")
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

target_clinical_data=target_clinical_data[which(!target_clinical_data$TARGET_SHORT%in%c("PATNKP-02A", "PAPVEB-04A")),]
target_counts=target_counts[, target_clinical_data$TARGET_SHORT]
target_length_norm =apply(target_counts, 2, function(x){x/(target_gene_length$TARGET.NBL_TARGET.30.PATHVK.01A.01R/1000)})

scale_factor=colSums(target_length_norm)/1000000
target_tpm=target_length_norm/scale_factor
target_log2tpm=log2(target_tpm+1)
rownames(target_log2tpm)=make.unique(as.character(target_genes$external_gene_name))
rownames(target_tpm)=make.unique(as.character(target_genes$external_gene_name))
rownames(target_counts)=make.unique(as.character(target_genes$external_gene_name))
percent_above_threshold=apply(target_log2tpm, 1, function(x){
  sum(x>5)/length(x)*100
})

#############prep SEQC data ########


seqc_nb=read.table("seqc/SEQC_RNAseq_gene_counts.txt", header = T, sep = "\t")
seqc_meta=read.csv("seqc/seqc_meta.csv")
adr_all$
rownames(seqc_nb)=seqc_nb$X
seqc_nb=seqc_nb[,2:ncol(seqc_nb)]

adr_mtx=adr_all@assays$RNA@counts
rownames(adr_mtx)=genes_10x$V1
merged_cells=cbind(adr_mtx, fetal_pods@assays$RNA@counts)
tg_comm_genes=intersect(rownames(merged_cells), rownames(target_counts))
cell_labs=c(adr_all@meta.data$new_clust, rep("Podocytes", 278))

train_with_podos_tg=trainModel(merged_cells[tg_comm_genes,], cell_labs, workers=NULL)
tg_against_adr=predictSimilarity(train_with_podos_tg, as.matrix(target_counts)[tg_comm_genes,], lengths = target_gene_length[tg_comm_genes,]$TARGET.NBL_TARGET.30.PATHVK.01A.01R, logits = F)

ha = rowAnnotation(Risk = factor(target_clinical_data$COG.Risk.Group, levels = c("Low Risk", "Intermediate Risk", "High Risk")), 
                   Age_at_diagnosis= target_clinical_data$Age.at.Diagnosis.in.Days,
                   Vital_status=target_clinical_data$Vital.Status, 
                   Histology=target_clinical_data$Histology, 
                   col=list(Risk=c("Low Risk"="grey90", "Intermediate Risk"="grey60","High Risk"= "black"), 
                            Age_at_diagnosis=colorRamp2(c(0, 8000), c("#FEE0D2", "#99000D")),
                            Vital_status=c("Alive"="#117733", "Dead"="#332288"), 
                            Histology=c("Favorable"="#AA4499", "Unfavorable"="#DDCC77", "Unknown"="grey")))

similarityHeatmap(tg_against_adr, cluster_rows=F, left_annotation=ha, column_order=c("Mesenchyme","SCPs","Bridge","Sympathoblastic", "Chromaffin","Cortex",
                                                                                     "Leukocytes","Endothelium", "Erythroblasts","Other","Podocytes"),
                  row_order=names(tg_against_adr[,8][order(tg_against_adr[,8], decreasing = T)]))
colnames(target_counts)

genes_10x <- read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)

rownames(seqc_nb)[grep("or", rownames(seqc_nb))]

seqc_strspl=sapply(strsplit(rownames(seqc_nb), "and"), "[", 1)

length(intersect(seqc_strspl, rownames(target_counts)))


test_seqc=seqc_nb
rownames(test_seqc)=make.unique(sapply(strsplit(rownames(seqc_nb), "and"), "[", 1))
test_seqc=test_seqc[intersect(rownames(test_seqc), rownames(target_counts)),]
dim(test_seqc)
sum(rowSums(seqc_nb))
sum(rowSums(test_seqc))

identical(rownames(merged_cells), genes_10x$V2)



######mouse stuff####### 

m13=read.delim("GSE99933_E13.5_counts.txt.gz")
m13=as.matrix(m13)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(m13) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)


rownames(genesV2)=make.unique(genesV2$MGI.symbol)

comm_rows=intersect(rownames(genesV2), rownames(m13))

m13=m13[comm_rows, ]
rownames(m13)=genesV2[comm_rows,]$HGNC.symbol


m13_srat=CreateSeuratObject(m13)
m13_srat = NormalizeData(m13_srat)
m13_srat =FindVariableFeatures(m13_srat, selection.method = "vst", nfeatures = 2000)
m13_srat = ScaleData(m13_srat, features = rownames(m13_srat))
m13_srat = RunPCA(m13_srat, npcs = 25)
m13_srat = FindNeighbors(m13_srat, dims=1:25)
m13_srat = FindClusters(m13_srat, resolution = 0.6)
m13_srat = RunUMAP(m13_srat, dims=1:25, min.dist = 0.5, n.neighbors = 30)

m13_srat@meta.data$new_clust="x"

m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==0)]="Sympathoblast"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==1)]="Chromaffin"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==2)]="SCPs"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==3)]="Bridge"

DimPlot(m13_srat, group.by = "new_clust", label=T)


library(SoupX)


mouse_marks=quickMarkers(m13_srat@assays$RNA@counts, m13_srat@meta.data$new_clust, FDR=0.01, N = Inf)
mouse_marks_filt=dplyr::filter(mouse_marks, tfidf>0.8)

m13_srat@active.ident
mouse_cell_names=names(m13_srat@active.ident)
new_mouse_ident=m13_srat@meta.data$new_clust
names(new_mouse_ident)=mouse_cell_names
m13_srat@active.ident=as.factor(new_mouse_ident)
mouse_marks=FindAllMarkers(m13_srat, logfc.threshold = 1)
mouse_marks_filt=dplyr::filter(mouse_marks, avg_logFC>1)


tum_cell_names=names(inhouse_tum_only@active.ident)
new_tum_ident=inhouse_tum_only@meta.data$risk
names(new_tum_ident)=tum_cell_names
inhouse_tum_only@active.ident=as.factor(new_tum_ident)
inh_tum_marks=FindAllMarkers(inhouse_tum_only, logfc.threshold = 1)
inh_tum_marks_filt=dplyr::filter(inh_tum_marks, avg_logFC>2)
inh_tum_marks=quickMarkers(inhouse_tum_only@assays$RNA@counts, inhouse_tum_only@meta.data$risk, FDR=0.01, N=Inf)
inh_tum_marks_filt=dplyr::filter(inh_tum_marks, tfidf>0.4 & geneFrequencySecondBest <0.2)
table(inh_tum_marks_filt$cluster)


low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

target_risks=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
target_risks$patient=target_risks$TARGET_SHORT
unique(target_risks$risk)

scp_marks=as.character(mouse_marks_filt[which(mouse_marks_filt$cluster=="SCPs"),]$gene)
bridge_marks=as.character(mouse_marks_filt[which(mouse_marks_filt$cluster=="Bridge"),]$gene)
chromaffin_marks=as.character(mouse_marks_filt[which(mouse_marks_filt$cluster=="Chromaffin"),]$gene)
sympathoblastic_marks=as.character(mouse_marks_filt[which(mouse_marks_filt$cluster=="Sympathoblast"),]$gene)


t_scp=as.data.frame(percent_above_threshold[scp_marks[which(scp_marks%in%names(percent_above_threshold))]])
t_bridge=as.data.frame(percent_above_threshold[bridge_marks[which(bridge_marks%in%names(percent_above_threshold))]])
t_sympathoblastic=as.data.frame(percent_above_threshold[sympathoblastic_marks[which(sympathoblastic_marks%in%names(percent_above_threshold))]])
t_chromaffin=as.data.frame(percent_above_threshold[chromaffin_marks[which(chromaffin_marks%in%names(percent_above_threshold))]])

t_scp$marker="SCP"
t_bridge$marker="Bridge"
t_sympathoblastic$marker="Sympathoblastic"
t_chromaffin$marker="Chromaffin"

t_scp$gene=rownames(t_scp)
t_bridge$gene=rownames(t_bridge)
t_sympathoblastic$gene=rownames(t_sympathoblastic)
t_chromaffin$gene=rownames(t_chromaffin)

rownames(t_scp)=1:nrow(t_scp)
rownames(t_bridge)=1:nrow(t_bridge)
rownames(t_sympathoblastic)=1:nrow(t_sympathoblastic)
rownames(t_chromaffin)=1:nrow(t_chromaffin)

colnames(t_scp)=c("value", "marker", "gene")
colnames(t_bridge)=c("value", "marker", "gene")
colnames(t_sympathoblastic)=c("value", "marker", "gene")
colnames(t_chromaffin)=c("value", "marker","gene")

t_all_m=rbind(t_scp, t_bridge,t_sympathoblastic,t_chromaffin)

ggplot(t_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >5 log2(TPM) in TARGET cohort (N=154)")


get_per_risk=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
    })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  scp_markers_in_rg=values_per_risk_group[scp_marks[which(scp_marks%in%rownames(values_per_risk_group))],]
  bridge_markers_in_rg=values_per_risk_group[bridge_marks[which(bridge_marks%in%rownames(values_per_risk_group))],]
  sympathoblastic_markers_in_rg=values_per_risk_group[sympathoblastic_marks[which(sympathoblastic_marks%in%rownames(values_per_risk_group))],]
  chromaffin_markers_in_rg=values_per_risk_group[chromaffin_marks[which(chromaffin_marks%in%rownames(values_per_risk_group))],]
  
  
  scp_markers_in_rg=melt(scp_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  bridge_markers_in_rg=melt(bridge_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  sympathoblastic_markers_in_rg=melt(sympathoblastic_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  chromaffin_markers_in_rg=melt(chromaffin_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))

  
  scp_markers_in_rg$marker="SCP"
  bridge_markers_in_rg$marker="Bridge"
  sympathoblastic_markers_in_rg$marker="Sympathoblastic"
  chromaffin_markers_in_rg$marker="Chromaffin"
 
  
  all_in_groups=rbind(scp_markers_in_rg, bridge_markers_in_rg, sympathoblastic_markers_in_rg, chromaffin_markers_in_rg)
  
  return(all_in_groups)
  
}

target_per_risk=get_per_risk(target_log2tpm, target_risks, 5)
target_per_risk$marker=factor(target_per_risk$marker, levels=c("SCP", "Bridge", "Chromaffin", "Sympathoblastic"))

ggplot(data = target_per_risk, aes(x=marker,y=value)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5)+facet_wrap(~variable)+
  stat_summary(mapping = aes(x = marker, y = value),
                  geom="pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))


get_per_risk_tum=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
    })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  high_markers_in_rg=values_per_risk_group[high_marks[which(high_marks%in%rownames(values_per_risk_group))],]
  low_markers_in_rg=values_per_risk_group[low_marks[which(low_marks%in%rownames(values_per_risk_group))],]
  
  high_markers_in_rg=melt(high_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  low_markers_in_rg=melt(low_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))

  
  high_markers_in_rg$marker="high risk"
  low_markers_in_rg$marker="low/intermediate risk"
  
  
  all_in_groups=rbind(high_markers_in_rg, low_markers_in_rg)
  
  return(all_in_groups)
  
}

high_marks=as.character(inh_tum_marks_filt[which(inh_tum_marks_filt$cluster=="high risk"),]$gene)
low_marks=as.character(inh_tum_marks_filt[which(inh_tum_marks_filt$cluster=="intermediate/low risk"),]$gene)
target_per_risk_tum=get_per_risk_tum(target_log2tpm, target_risks, 5)
target_per_risk_tum$marker=factor(target_per_risk_tum$marker, levels=c("low/intermediate risk","high risk"))
ggplot(data = target_per_risk_tum, aes(x=marker,y=value)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5)+facet_wrap(~variable)+
  stat_summary(mapping = aes(x = marker, y = value),
               geom="pointrange",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))

DimPlot(just_med, group.by = "timepoints")
just_med@meta.data$timepoints=just_med@meta.data$sample_name
just_med@meta.data$timepoints[grep("w21", just_med@meta.data$sample_name)]="w21"
just_med@meta.data$timepoints[grep("w10", just_med@meta.data$sample_name)]="w10d5"
just_med@meta.data$timepoints=factor(just_med@meta.data$timepoints, levels=c("w8", "w8d6", "w10d5","w11", "w21"))



med_cell_names=names(just_med@active.ident)
new_med_ident=just_med@meta.data$time_ct
names(new_med_ident)=med_cell_names
just_med@active.ident=as.factor(new_med_ident)
source("gene_sets.R")
tf_table=read.table("useful_files/Homo_sapiens_TF.txt", sep = "\t", header = T)
tfs=as.character(tf_table$Symbol)
excludeGenes=c(riboGenes, mtGenes, riboGenes)
keep_genes=rownames(just_med)[which(!rownames(just_med)%in%excludeGenes)]
just_med=subset(just_med, features=tfs)
just_med=ScaleData(just_med, features=rownames(just_med))
scp_de = FindMarkers(just_med, ident.1 = "early_SCPs", ident.2 = "late_SCPs", min.pct = 0.25)
bridge_de=FindMarkers(just_med, ident.1 = "early_Bridge", ident.2 = "late_Bridge", min.pct = 0.25)
symp_de=FindMarkers(just_med, ident.1 = "early_Sympathoblastic", ident.2 = "late_Sympathoblastic", min.pct = 0.25)
chromaffin_de=FindMarkers(just_med, ident.1 = "early_Chromaffin", ident.2 = "late_Chromaffin", min.pct = 0.25)


scp_de=scp_de[which(scp_de$p_val_adj<0.01),]

scp_de$lh="x"
scp_de$lh[which(scp_de$avg_logFC>0)]="early"
scp_de$lh[which(scp_de$avg_logFC<0)]="late"
bridge_de=bridge_de[which(bridge_de$p_val_adj<0.01),]
symp_de=symp_de[which(symp_de$p_val_adj<0.01),]
chromaffin_de=chromaffin_de[which(chromaffin_de$p_val_adj<0.01),]

rownames(scp_de[1:20,][order(scp_de[1:20,]$avg_logFC, decreasing = T),])
med_scp=subset(just_med, idents=c("early_SCPs", "late_SCPs"))
med_scp=ScaleData(med_scp, features = rownames(med_scp))
col_fun = colorRamp2(c(-3,-1,0,1,3), rev(brewer.pal(n = 5, name = "RdYlBu")))
Heatmap(as.matrix(med_scp@assays$RNA@scale.data[rownames(scp_de),]),col=col_fun, cluster_rows = T, cluster_columns = F,
        show_column_names = F, column_split=med_scp@meta.data$time_ct, row_split = scp_de$lh, cluster_row_slices = F, show_row_dend = F, cluster_column_slices = F)

DotPlot(just_med, idents=c("early_SCPs", "late_SCPs"), features=c("MPZ", "SOX10", "SOX2", "PLP1", "ERBB3",rownames(scp_de[1:20,][order(scp_de[1:20,]$avg_logFC, decreasing = T),]))) + RotatedAxis() +
    scale_color_gradient2(low="#313695", mid = "#ffffbf", high = "#a50026") 

DotPlot(just_med, idents=c("early_Bridge", "late_Bridge"), features=c("DLL3", "HOXB9", "TUBB3", "ASCL1", "SHD",rownames(bridge_de[1:20,][order(bridge_de[1:20,]$avg_logFC, decreasing = T),]))) + RotatedAxis() +
  scale_color_gradient2(low="#313695", mid = "#ffffbf", high = "#a50026") 

DotPlot(just_med, idents=c("early_Sympathoblastic", "late_Sympathoblastic"), features=c("STMN4", "PRPH", "ISL1", "ELAVL4", "TTC9B",rownames(symp_de[1:20,][order(symp_de[1:20,]$avg_logFC, decreasing = T),]))) + RotatedAxis() +
  scale_color_gradient2(low="#313695", mid = "#ffffbf", high = "#a50026") 

DotPlot(just_med, idents=c("late_Chromaffin", "early_Chromaffin"), features=c("PENK", "CARTPT", "INSM1", "NPAS4", "SCG2",rownames(chromaffin_de[1:20,][order(chromaffin_de[1:20,]$avg_logFC, decreasing = T),]))) + RotatedAxis() +
  scale_color_gradient2(low="#313695", mid = "#ffffbf", high = "#a50026") 





adr_all_med=subset(adr_all, idents=c(25,24,13,19,20))

tum_cells_inhouse=subset(srat_inhouse, idents = c(1,5,8,12,16,25))
tum_cells_dutch=subset(srat_dutch, idents = c(24,13,10,17))

#assign all tumour cells the same identity  
tum_cells_inhouse@meta.data$new_idents="tumour"
tum_cells_dutch@meta.data$new_idents="tumour"
adr_all_med@meta.data$new_idents=as.character(adr_all_med@meta.data$seurat_clusters) 
#combine each dataset with medulla cells and then find markers
combined_inhouse=merge(adr_all_med, tum_cells_inhouse)
inhouse_markers=quickMarkers(combined_inhouse@assays$RNA@counts, combined_inhouse@meta.data$new_idents, FDR = 0.01, N=Inf)
combined_dutch=merge(adr_all_med, tum_cells_dutch)
dutch_markers=quickMarkers(combined_dutch@assays$RNA@counts, combined_dutch@meta.data$new_idents, FDR = 0.01, N=Inf)
#filter by frequency in 2nd best cluster
dutch_tc=dutch_markers[which(dutch_markers$cluster=="tumour"),]
dutch_tc_filt=dutch_tc[which(dutch_tc$geneFrequencySecondBest<0.2),]


inhouse_tc=inhouse_markers[which(inhouse_markers$cluster=="tumour"),]
inhouse_tc_filt=inhouse_tc[which(inhouse_tc$geneFrequencySecondBest<0.2),]

#calculate mean tfidf and sort by it
no_org_lj=full_join(inhouse_tc_filt, dutch_tc_filt, by='gene')
no_org_lj$mean_tfidf=rowMeans(no_org_lj[,c(8,17)], na.rm = T)
no_org_final=no_org_lj[no_org_lj$mean_tfidf>0.85,]
no_org_final_sorted=no_org_final[order(no_org_final$mean_tfidf, decreasing = T),]


no_org_final=no_org_final[which(no_org_final$gene%in%rownames(target_log2tpm)),]
tum_de_genes=as.character(no_org_final$gene)

#exclude genes expressed in more than 25% of leukocyte cells in either sample as these are likely contaminants
inhouse_immune_cells=rownames(srat_inhouse@meta.data)[which(srat_inhouse@meta.data$idents_for_plot=="Leukocytes")]
inhouse_immune_cellexp=apply(srat_inhouse@assays$RNA@counts[,inhouse_immune_cells],1,function(y){
  sum(y>0)/length(y)*100
})

dutch_immune_cells=rownames(srat_dutch@meta.data)[which(srat_dutch@meta.data$idents_for_plot=="Leukocytes")]
dutch_immune_cellexp=apply(srat_dutch@assays$RNA@counts[,dutch_immune_cells],1,function(y){
  sum(y>0)/length(y)*100
})

immune_inhouse_t25=tum_de_genes[which(tum_de_genes%in%names(inhouse_immune_cellexp[inhouse_immune_cellexp>25]))]
dutch_immune_t25=tum_de_genes[which(tum_de_genes%in%names(dutch_immune_cellexp[dutch_immune_cellexp>25]))]

immune_sc=unique(c(immune_inhouse_t25, dutch_immune_t25))
final_de=no_org_final[which(!no_org_final$gene%in%immune_sc),]

a=gsub("\\.x", ".GOSH", colnames(final_de))
b=gsub("\\.y", ".PMC", a)

colnames(final_de)=b
rownames(final_de)=1:nrow(final_de)

final_de_sorted=final_de[order(final_de$mean_tfidf, decreasing = T),]


write.table(final_de_sorted, "final_de_sorted.txt", sep = "\t", row.names = F, col.names = T, quote = F)


edg_file=DGEList(counts=extreme_count_tab, samples = extremes, genes = target_genes, group=extremes$risk)

o = order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file = edg_file[o,]
d = duplicated(edg_file$genes$external_gene_name)
edg_file = edg_file[!d,]

rownames(edg_file)=edg_file$genes$external_gene_name

keep = filterByExpr(edg_file)
table(keep)

edg_file = edg_file[keep, , keep.lib.sizes=FALSE]

edg_file=calcNormFactors(edg_file)

edg_file$samples$group=factor(edg_file$samples$group, levels = c("low_risk_4s", "high_risk"))
#design a model matrix. age at diagnosis and mycn status included as covariates
design=model.matrix(~group+Age.at.Diagnosis.in.Days+MYCN.status,edg_file$samples)
fit=glmQLFit(edg_file, design, robust=TRUE, dispersion = 0.4^2)
qlf=glmQLFTest(fit,coef=2)

de_of_de=qlf$table[as.character(final_de_sorted$gene)[which(as.character(final_de_sorted$gene)%in%rownames(qlf$table))],]
de_of_de_sorted=de_of_de[order(abs(de_of_de$logFC), decreasing = T),]
de_of_de_sorted$qval=p.adjust(de_of_de_sorted$PValue, method = "fdr")


de_of_de_pat=as.data.frame(percent_above_threshold[rownames(de_of_de_sorted)])
de_of_de_pat$gene=rownames(de_of_de_pat)

colnames(de_of_de_pat)=c("value", "gene")
de_of_de_pat=de_of_de_pat[order(de_of_de_pat$value, decreasing = T),]
de_of_de_pat$gene=factor(de_of_de_pat$gene, levels = de_of_de_pat$gene)

#% TARGET samples expressed above threshold 
ggplot(de_of_de_pat, aes(x=reorder(gene, value), y=value))+geom_bar(stat="identity", width = 0.6) +coord_flip() +
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", face = "italic", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 12)) + labs( y ="% of samples above expression threshold", x = "Gene")

#import SEQC data

seqc_counts=read.table("seqc/GSE49711_SEQC_NB_TUC_G_log2.txt", sep = "\t", header = T)
seqc_meta=read.csv("seqc/seqc_meta.csv")

rownames(seqc_counts)=seqc_counts$X00gene_id

seqc_counts=seqc_counts[,2:ncol(seqc_counts)]
colnames(seqc_counts)=sapply(strsplit(colnames(seqc_counts), "_"), "[", 2)


seqc_unlist=apply(seqc_counts, 1, function(x){unlist(x)})

seqc_unlist_t=t(seqc_unlist)


#transform from log2(FPKM+1)
untransf=(2^seqc_unlist_t)-1

#transform to log2(tpm+1)
seqc_tpm=apply(untransf,2,function(x){
  (x/sum(x))*1000000
})


seqc_tpm_log2=log2(seqc_tpm+1)

seqc_percent_above_threshold=apply(seqc_tpm_log2, 1, function(x){
  sum(x>5)/length(x)*100
})
#SEQC risks
seqc_lr=dplyr::filter(seqc_meta, Age<540 & MYCN.status=="1", INSS.Stage!="4S")
seqc_4s=dplyr::filter(seqc_meta, INSS.Stage=="4S")
seqc_hr=dplyr::filter(seqc_meta, Age>540 & MYCN.status=="Amp" &INSS.Stage!="4S")
seqpt1=rbind(seqc_lr, seqc_4s, seqc_hr)
seqc_intermediate=seqc_meta[-which(seqc_meta$NB.ID%in%seqpt1$NB.ID),]

seqc_lr$risk="low_risk"
seqc_4s$risk="low_risk_4s"
seqc_hr$risk="high_risk"
seqc_intermediate$risk="intermediate_risk"
seqc_risks=rbind(seqc_lr, seqc_4s, seqc_hr, seqc_intermediate)

seqc_risks$patient=seqc_risks$NB.ID

total_genes=length(rownames(edg_file))
genes_over_50=length(rownames(de_of_de_pat)[which(de_of_de_pat$value>50)])
de_genes=length(rownames(de_of_de_pat))
not_de_genes=de_genes-genes_over_50
genes_not_assoc=total_genes-de_genes

symp_to_bridge_ratio=tg_against_adr[,7]/tg_against_adr[,4]
symp_sig=tg_against_adr[,7]
target_clinical_data$new_cats=paste0(target_clinical_data$COG.Risk.Group, "_", target_clinical_data$Vital.Status)
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse"))

target_clinical_data$new_cats[which(target_clinical_data$TARGET_SHORT%in%high_risk_relapse$TARGET_SHORT)]="High Risk_Recurrent"

target_clinical_data$new_cats[which(target_clinical_data$new_cats%in%c("Low Risk_Alive", "Intermediate Risk_Dead","Intermediate Risk_Alive"))]="Low/Intermediate Risk"


symp_to_b_df=data.frame(symp_to_bridge_ratio, target_clinical_data$new_cats, symp_sig, target_clinical_data$Percent.Tumor)
symp_to_b_df$target_clinical_data.Percent.Tumor=as.character(symp_to_b_df$target_clinical_data.Percent.Tumor)
ggplot(symp_to_b_df, aes(target_clinical_data.new_cats, symp_sig)) + geom_boxplot()


symp_to_b_df$target_clinical_data.Percent.Tumor[which(symp_to_b_df$target_clinical_data.Percent.Tumor=="")]="Unknown"

symp_to_b_df$target_clinical_data.Percent.Tumor=factor(symp_to_b_df$target_clinical_data.Percent.Tumor,
                                                       levels = c("10", "10-20", "20",   "20-30","30","30-40","35",
                                                                  "40", "50","55","60", "65", "70","70-80","75","75-85",
                                                                  "80", "80-90","85","85-90", "90", "90-95","95","97-98", "100", "Unknown"))

ggplot(symp_to_b_df, aes(target_clinical_data.Percent.Tumor, symp_to_bridge_ratio)) + geom_point()

symp_to_b_df$rg=as.character(target_clinical_data$COG.Risk.Group)
symp_to_b_df$rg[which(target_clinical_data$INSS.Stage=="Stage 4s")]="4S"
symp_to_b_df$rg=factor(symp_to_b_df$rg, levels=c("4S", "Intermediate Risk", "High Risk"))

ggplot(symp_to_b_df, aes(rg, symp_to_bridge_ratio)) + geom_boxplot()
