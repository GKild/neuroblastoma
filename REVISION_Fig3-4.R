library(dplyr)
library(Seurat)
####################prep the target data####################
target_counts=read.table("Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv",
                         header = T, sep= "\t")
target_genes=read.table("Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t",
                        header = T)
target_meta=read.table("Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")
target_clinical_data=read.table("Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv",
                                sep = "\t", header = T)
target_gene_length= read.table("Neuroblastoma/bulk_NB/TARGET-NBL_BulkGeneLengths.tsv", sep = "\t")
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

###### LR of target against fetal adrenal ##############

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


#####define target RGs########

#define risk groups
#all TARGET low risk are stage 4S
low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk")

low_risk_4s$risk="low_risk_4s"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"
target_risks=rbind(low_risk_4s,intermediate_risk, high_risk)

#######seqc data##########

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
target_risks$patient=target_risks$TARGET_SHORT

scp_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="SCPs"),]$gene)
bridge_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="Bridge"),]$gene)
chromaffin_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="Chromaffin"),]$gene)
sympathoblastic_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="Sympathoblastic"),]$gene)
pod_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="Podocytes"),]$gene)

t_scp=as.data.frame(percent_above_threshold[scp_marks[which(scp_marks%in%names(percent_above_threshold))]])
t_bridge=as.data.frame(percent_above_threshold[bridge_marks[which(bridge_marks%in%names(percent_above_threshold))]])
t_sympathoblastic=as.data.frame(percent_above_threshold[sympathoblastic_marks[which(sympathoblastic_marks%in%names(percent_above_threshold))]])
t_chromaffin=as.data.frame(percent_above_threshold[chromaffin_marks[which(chromaffin_marks%in%names(percent_above_threshold))]])
t_pods=as.data.frame(percent_above_threshold[pod_marks[which(pod_marks%in%names(percent_above_threshold))]])

t_scp$marker="scp"
t_bridge$marker="bridge"
t_sympathoblastic$marker="sympathoblastic"
t_chromaffin$marker="chromaffin"
t_pods$marker="podocyte"

t_scp$gene=rownames(t_scp)
t_bridge$gene=rownames(t_bridge)
t_sympathoblastic$gene=rownames(t_sympathoblastic)
t_chromaffin$gene=rownames(t_chromaffin)
t_pods$gene=rownames(t_pods)

rownames(t_scp)=1:nrow(t_scp)
rownames(t_bridge)=1:nrow(t_bridge)
rownames(t_sympathoblastic)=1:nrow(t_sympathoblastic)
rownames(t_chromaffin)=1:nrow(t_chromaffin)
rownames(t_pods)=1:nrow(t_pods)

colnames(t_scp)=c("value", "marker", "gene")
colnames(t_bridge)=c("value", "marker", "gene")
colnames(t_sympathoblastic)=c("value", "marker", "gene")
colnames(t_chromaffin)=c("value", "marker","gene")
colnames(t_pods)=c("value", "marker","gene")

t_all_m=rbind(t_scp, t_bridge,t_sympathoblastic,t_chromaffin, t_pods)

seqc_scp=as.data.frame(seqc_percent_above_threshold[scp_marks[which(scp_marks%in%names(seqc_percent_above_threshold))]])
seqc_bridge=as.data.frame(seqc_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(seqc_percent_above_threshold))]])
seqc_sympathoblastic=as.data.frame(seqc_percent_above_threshold[sympathoblastic_marks[which(sympathoblastic_marks%in%names(seqc_percent_above_threshold))]])
seqc_chromaffin=as.data.frame(seqc_percent_above_threshold[chromaffin_marks[which(chromaffin_marks%in%names(seqc_percent_above_threshold))]])
seqc_pods=as.data.frame(seqc_percent_above_threshold[pod_marks[which(pod_marks%in%names(seqc_percent_above_threshold))]])

seqc_scp$marker="scp"
seqc_bridge$marker="bridge"
seqc_sympathoblastic$marker="sympathoblastic"
seqc_chromaffin$marker="chromaffin"
seqc_pods$marker="podocyte"

seqc_scp$gene=rownames(seqc_scp)
seqc_bridge$gene=rownames(seqc_bridge)
seqc_sympathoblastic$gene=rownames(seqc_sympathoblastic)
seqc_chromaffin$gene=rownames(seqc_chromaffin)
seqc_pods$gene=rownames(seqc_pods)

rownames(seqc_scp)=1:nrow(seqc_scp)
rownames(seqc_bridge)=1:nrow(seqc_bridge)
rownames(seqc_sympathoblastic)=1:nrow(seqc_sympathoblastic)
rownames(seqc_chromaffin)=1:nrow(seqc_chromaffin)
rownames(seqc_pods)=1:nrow(seqc_pods)

colnames(seqc_scp)=c("value", "marker", "gene")
colnames(seqc_bridge)=c("value", "marker", "gene")
colnames(seqc_sympathoblastic)=c("value", "marker", "gene")
colnames(seqc_chromaffin)=c("value", "marker","gene")
colnames(seqc_pods)=c("value", "marker","gene")

seqc_all_m=rbind(seqc_scp, seqc_bridge,seqc_sympathoblastic,seqc_chromaffin, seqc_pods)
t_all_m$dataset="TARGET"
seqc_all_m$dataset="SEQC"
joined_3a=rbind(t_all_m,seqc_all_m)


joined_3a$marker=factor(joined_3a$marker, levels = c("podocyte","scp", "bridge", "sympathoblastic", "chromaffin"))
joined_3a$dataset=factor(joined_3a$dataset, levels = c("TARGET", "SEQC"))
#Figure 3A
pdf("marks_all.pdf", width = 4, height = 3, useDingbats = F)
gg=ggplot(data = joined_3a, aes(x=marker,y=value, group=dataset)) +
  geom_quasirandom(mapping = aes(x=marker,y=value, group=dataset), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  stat_summary(mapping = aes(x = marker, y = value, shape=dataset),
                  geom = "pointrange",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black", 
                  position=position_dodge(width=0.8)) +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Markers", y = "% samples above expression threshold")
plot(gg)
dev.off()



#non-medullary samples (target)
#exclude adrenal or putatively adrrnal samples
leftover_ids=as.character(target_clinical_data$ICDO.Description)[-grep("adrenal|abdominal|kidney|abdomen|unknown|retroperitoneum|other", tolower(as.character(target_clinical_data$ICDO.Description)))]
#exclude empty strings
leftover_ids=leftover_ids[c(1:13,16:23)]

leftover_samples=target_clinical_data$TARGET_SHORT[which(target_clinical_data$ICDO.Description%in%leftover_ids)]

os_percent_above_threshold=apply(target_log2tpm[,leftover_samples], 1, function(y){
  sum(y>5)/length(y)*100
})





os_scp=as.data.frame(os_percent_above_threshold[scp_marks[which(scp_marks%in%names(os_percent_above_threshold))]])
os_bridge=as.data.frame(os_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(os_percent_above_threshold))]])
os_sympathoblastic=as.data.frame(os_percent_above_threshold[sympathoblastic_marks[which(sympathoblastic_marks%in%names(os_percent_above_threshold))]])
os_chromaffin=as.data.frame(os_percent_above_threshold[chromaffin_marks[which(chromaffin_marks%in%names(os_percent_above_threshold))]])
os_pods=as.data.frame(os_percent_above_threshold[pod_marks[which(pod_marks%in%names(os_percent_above_threshold))]])

os_scp$marker="scp"
os_bridge$marker="bridge"
os_sympathoblastic$marker="sympathoblastic"
os_chromaffin$marker="chromaffin"
os_pods$marker="podocytes"

os_scp$gene=rownames(os_scp)
os_bridge$gene=rownames(os_bridge)
os_sympathoblastic$gene=rownames(os_sympathoblastic)
os_chromaffin$gene=rownames(os_chromaffin)
os_pods$gene=rownames(os_pods)

rownames(os_scp)=1:nrow(os_scp)
rownames(os_bridge)=1:nrow(os_bridge)
rownames(os_sympathoblastic)=1:nrow(os_sympathoblastic)
rownames(os_chromaffin)=1:nrow(os_chromaffin)
rownames(os_pods)=1:nrow(os_pods)

colnames(os_scp)=c("value", "marker", "gene")
colnames(os_bridge)=c("value", "marker", "gene")
colnames(os_sympathoblastic)=c("value", "marker", "gene")
colnames(os_chromaffin)=c("value", "marker","gene")
colnames(os_pods)=c("value", "marker","gene")

os_all_m=rbind(os_scp, os_bridge,os_sympathoblastic,os_chromaffin, os_pods)

os_all_m$marker=factor(os_all_m$marker, levels=c("podocytes" , "scp","bridge","sympathoblastic","chromaffin"))
#Figure 3B
pdf("marks_om.pdf", width = 3, height = 3, useDingbats = F)
ggplot(data = os_all_m, aes(x=marker,y=value)) +
  geom_quasirandom(mapping = aes(x=marker,y=value), alpha=0.2, cex=0.5, shape=19) +
  stat_summary(mapping = aes(x = marker, y = value),
                  geom="pointrange", shape=22,
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#per risk group function

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
  podo_markers_in_rg=values_per_risk_group[pod_marks[which(pod_marks%in%rownames(values_per_risk_group))],]
  
  scp_markers_in_rg=melt(scp_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  bridge_markers_in_rg=melt(bridge_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  sympathoblastic_markers_in_rg=melt(sympathoblastic_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  chromaffin_markers_in_rg=melt(chromaffin_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  podo_markers_in_rg=melt(podo_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  
  scp_markers_in_rg$marker="scp"
  bridge_markers_in_rg$marker="bridge"
  sympathoblastic_markers_in_rg$marker="sympathoblastic"
  chromaffin_markers_in_rg$marker="chromaffin"
  podo_markers_in_rg$marker="podocyte"
  
  all_in_groups=rbind(scp_markers_in_rg, bridge_markers_in_rg, sympathoblastic_markers_in_rg, chromaffin_markers_in_rg, podo_markers_in_rg)
  
  return(all_in_groups)
  
}

target_per_risk=get_per_risk(target_log2tpm, target_risks, 5)
seqc_per_risk=get_per_risk(seqc_tpm_log2, seqc_risks, 5)

target_per_risk$cohort="TARGET"
seqc_per_risk$cohort="SEQC"

combined_per_risk=rbind(target_per_risk, seqc_per_risk)

colnames(combined_per_risk)=c("gene", "risk", "value", "marker", "cohort")

combined_per_risk$risk=factor(combined_per_risk$risk, levels=c("low_risk","low_risk_4s", "intermediate_risk", "high_risk"))

combined_per_risk$marker=factor(combined_per_risk$marker, levels=c("podocyte" , "scp","bridge","sympathoblastic","chromaffin"))

combined_per_risk$cohort=factor(combined_per_risk$cohort, levels = c("TARGET", "SEQC"))
#Figure 3C
ggplot(data = combined_per_risk, aes(x=marker,y=value)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5)+facet_grid(vars(risk), vars(cohort))+
  stat_summary(mapping = aes(x = marker, y = value),
                  geom = "pointrange",
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

########difference in sympathoblast signal in  dead/alive/relapse#############

target_clinical_data$new_cats=paste0(target_clinical_data$COG.Risk.Group, "_", target_clinical_data$Vital.Status)
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse"))

target_clinical_data$new_cats[which(target_clinical_data$TARGET_SHORT%in%high_risk_relapse$TARGET_SHORT)]="High Risk_Recurrent"

target_clinical_data$new_cats[which(target_clinical_data$new_cats%in%c("Low Risk_Alive", "Intermediate Risk_Dead","Intermediate Risk_Alive"))]="Low/Intermediate Risk"


values_per_risk_group=sapply(unique(target_clinical_data$new_cats), function(x){
  pts=target_clinical_data$TARGET_SHORT[which(target_clinical_data$new_cats==x)]
  tab=target_log2tpm[,pts]
  apply(tab, 1, function(y){
    sum(y>5)/length(y)*100
  })
})

sympathoblastic_markers_in_rg=values_per_risk_group[sympathoblastic_marks[which(sympathoblastic_marks%in%rownames(values_per_risk_group))],]


sympatho_melted=melt(sympathoblastic_markers_in_rg)


ggplot(data = sympatho_melted, aes(x=Var2,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+
  stat_summary(mapping = aes(x = Var2, y = value),
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


##########differential expression for fig4###############
adr_all=readRDS("fAdrenal/processedSeurat.RDS")
srat_inhouse=readRDS("srat_inhouse.rds")
srat_dutch=readRDS("srat_dutch.rds")

adr_all_med=subset(adr_all, idents=c(25,24,13,19,20))

tum_cells_inhouse=subset(srat_inhouse, idents = c(1,5,8,12,16,25))
tum_cells_dutch=subset(srat_dutch, idents = c(24,13,10))

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
sortedCols=c("gene","mean_tfidf","tfidf.PMC","tfidf.GOSH","qval.PMC","qval.GOSH","geneFrequency.PMC","geneFrequency.GOSH","geneFrequencyOutsideCluster.PMC",
             "geneFrequencyOutsideCluster.GOSH",	"geneFrequencySecondBest.PMC","geneFrequencySecondBest.GOSH")
final_de_sorted[,sortedCols]

write.table(final_de_sorted[,sortedCols], "final_de_sorted.txt", sep = "\t", row.names = F, col.names = T, quote = F)

extremes=rbind(low_risk_4s, high_risk)

#subset to only include low and high risk samples
extreme_count_tab=target_counts[,extremes$TARGET_SHORT]
edg_file=DGEList(counts=extreme_count_tab, samples = extremes, genes = target_genes, group=extremes$risk)

o = order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file = edg_file[o,]
d = duplicated(edg_file$genes$external_gene_name)
edg_file = edg_file[!d,]

rownames(edg_file)=edg_file$genes$external_gene_name

keep = filterByExpr(edg_file)
table(keep)

edg_file = edg_file[keep,keep.lib.sizes=FALSE]

edg_file=calcNormFactors(edg_file)

edg_file$samples$group=factor(edg_file$samples$group, levels = c("low_risk_4s", "high_risk"))
which(edg_file$samples$MYCN.status=="Unknown")
edg_file$samples$MYCN.status[130]="Not Amplified"
edg_file$samples$MYCN.status=factor(edg_file$samples$MYCN.status,levels=c("Not Amplified", "Amplified"))

design=model.matrix(~group+Age.at.Diagnosis.in.Days+MYCN.status,edg_file$samples)
fit=glmQLFit(edg_file, design, robust=TRUE, dispersion = 0.4^2)
qlf=glmQLFTest(fit,coef=2)

de_of_de=qlf$table[as.character(final_de_sorted$gene)[which(as.character(final_de_sorted$gene)%in%rownames(qlf$table))],]
de_of_de_sorted=de_of_de[order(abs(de_of_de$logFC), decreasing = T),]
de_of_de_sorted$qval=p.adjust(de_of_de_sorted$PValue, method = "fdr")
de_of_de_sorted$gene=rownames(de_of_de_sorted)
de_of_de_sorted=de_of_de_sorted[,c(6,1:5)]
write.table(de_of_de_sorted, "de_of_de_sorted.txt", row.names = F, col.names = T, quote = F, sep = "\t")


de_of_de_pat=as.data.frame(percent_above_threshold[rownames(de_of_de_sorted)])
de_of_de_pat$gene=rownames(de_of_de_pat)

colnames(de_of_de_pat)=c("value", "gene")
de_of_de_pat=de_of_de_pat[order(de_of_de_pat$value, decreasing = T),]
de_of_de_pat$gene=factor(de_of_de_pat$gene, levels = de_of_de_pat$gene)
length(which(de_of_de_pat$value>=50))

#% TARGET samples expressed above threshold 
ggplot(de_of_de_pat, aes(x=reorder(gene, value), y=value))+geom_bar(stat="identity", width = 0.6) +coord_flip() +
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", face = "italic", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 12)) + labs( y ="% of samples above expression threshold", x = "Gene")

shared_genes=percent_above_threshold[rownames(edg_file)]
length(which(shared_genes>50))
phyper(37, 5923, 23924, 89, lower.tail = F)
above_log=de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),][which(abs(de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),]$logFC)>1),]
above_log=above_log[which(above_log$qval<0.05),]
genes_for_kap=rownames(above_log)

de_of_de_facet_tab=target_log2tpm[genes_for_kap, extremes$TARGET_SHORT]

de_tab=as.data.frame(t(de_of_de_facet_tab))  
de_tab$risk=extremes$risk
de_tab$patient=rownames(de_tab)

de_tab_melted=melt(de_tab, id.vars=c("patient", "risk"))
de_tab_melted$risk=factor(de_tab_melted$risk, levels = c("low_risk_4s","high_risk"))


seqc_risk2=seqc_risks

seqc_risk2$risk[which(seqc_risk2$risk=="low_risk")]="low_risk_4s"
seqc_risk2=seqc_risk2[which(seqc_risk2$risk%in%c("low_risk_4s", "high_risk")), ]

genes_for_seqc=genes_for_kap[which(genes_for_kap%in%rownames(seqc_tpm_log2))]


seqc_facet_tab=seqc_tpm_log2[genes_for_seqc, seqc_risk2$patient]
seqc_facet_tab=as.data.frame(t(seqc_facet_tab))

seqc_facet_tab$risk=seqc_risk2$risk
seqc_facet_tab$patient=rownames(seqc_facet_tab)

seqc_facet_tab_melted=melt(seqc_facet_tab, id.vars=c("patient", "risk"))

seqc_facet_tab_melted$risk=factor(seqc_facet_tab_melted$risk, levels = c("low_risk_4s","high_risk"))

de_tab_melted$dataset="TARGET"
seqc_facet_tab_melted$dataset="SEQC"

de_of_de_two_datasets=rbind(de_tab_melted, seqc_facet_tab_melted)
de_of_de_two_datasets$risk=as.character(de_of_de_two_datasets$risk)
de_of_de_two_datasets$risk[which(de_of_de_two_datasets$risk=="low_risk_4s")]="low_risk"
de_of_de_two_datasets$risk=factor(de_of_de_two_datasets$risk, levels = c("low_risk", "high_risk"))


ggplot(data=de_of_de_two_datasets, aes(x=risk, y=value, group=dataset)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5,dodge.width=.8)+facet_wrap(~variable,scales="free_y")+
  stat_summary(mapping = aes(x = risk, y = value, shape=dataset, fill=dataset),
               geom="pointrange", colour="#1c2e4a",
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median, fill="black",
               position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))+ 
  labs(y="Expression (log2(TPM))")


gen <- read.delim( "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm (1).gct", skip = 2, header = TRUE )

## First two columns corresponds to "Transcript name" and "Gene name"
gen_ann <- gen[ , c(1,2) ]
gen_clean <- as.matrix( gen[ , -seq( 2 ) ] )
rm(gen)
# Load map file
sample_map <- read.delim( "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" )
sample_map$SAMPID_new <- gsub( "-", "\\.", sample_map$SAMPID )
rownames( sample_map ) <- sample_map$SAMPID

sample_map=sample_map[which(sample_map$SAMPID_new%in%colnames(gen_clean)),]
srat_dutch=subset(srat_dutch, idents=17, invert=T)

dutch_tum_marks=quickMarkers(srat_dutch@assays$RNA@counts, srat_dutch@meta.data$new_idents, N = Inf, FDR=0.01)
inh_tum_marks=quickMarkers(srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$new_idents, N=Inf, FDR = 0.01)
adr_all_marks=quickMarkers(adr_all@assays$RNA@counts, adr_all@meta.data$new_clust, N=Inf, FDR = 0.01)

dutch_marks_filt=filter(dutch_tum_marks, cluster=="Tumour", geneFrequency>0.5, geneFrequencySecondBest<0.2)
inh_marks_filt=filter(inh_tum_marks, cluster=="tumour", geneFrequency>0.5, geneFrequencySecondBest<0.2)


tum_genes_high=intersect(inh_marks_filt$gene, dutch_marks_filt$gene)
just_med=subset(adr_all,
                cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$new_clust%in%c("SCPs", "Bridge",
                                                                                         "Sympathoblastic", "Chromaffin"))])
med_gene_frac=apply(just_med@assays$RNA@counts, 1, function(x){sum(x>0)/6451*100})
med_genes_high=names(med_gene_frac[which(med_gene_frac>10)])
med_and_tum=intersect(tum_genes_high, med_genes_high)

genes_of_i=c(med_and_tum, genes_for_kap, "B4GALNT1", "ST3GAL5", "ST8SIA1") 
genes_of_i=intersect(genes_of_i, rownames(target_tpm))

length(which(genes_of_i%in%gen_ann$Description))
rownames(gen_clean)=gen_ann$Description

lr_target_tab=target_tpm[genes_of_i,target_clinical_data$TARGET_SHORT[which(target_clinical_data$COG.Risk.Group=="Low Risk")]]
lr_target_means=rowMeans(lr_target_tab)
hr_target_tab=target_tpm[genes_of_i,target_clinical_data$TARGET_SHORT[which(target_clinical_data$COG.Risk.Group=="High Risk")]]
hr_target_means=rowMeans(hr_target_tab)

gtex_goi=gen_clean[genes_of_i,]

bla=sapply(as.character(unique(sample_map$SMTSD)), function(x){
  rowMeans(gtex_goi[,sample_map$SAMPID_new[which(sample_map$SMTSD==x)]])
})

bla=as.data.frame(bla)

bla=cbind(lr_target_means, hr_target_means, bla)
bla_log2=log2(bla+1)

labels_for_gtex=read.csv("labels_for_Gtex.csv")
rownames(labels_for_gtex)=labels_for_gtex$OLD.LABEL
bla_log2=bla_log2[,rownames(labels_for_gtex)]
colnames(bla_log2)=labels_for_gtex$NEW.LABEL

no_brain=bla_log2[,17:ncol(bla_log2)]
row_ord=c(rep("tum-med shared", 95), rep("tum-med different", 10), rep("GD2", 3))
gene_annot=data.frame(no_brain_mean=rowMeans(no_brain), annot=row_ord)


shared_gene_rep=read.csv('shared_gene_rep.csv')
rownames(shared_gene_rep)=shared_gene_rep$Gene
shared_gene_rep=shared_gene_rep[which(rownames(shared_gene_rep)%in%med_and_tum),]
length(which(shared_gene_rep$Gene%in%rownames(bla_log2)))
gene_annot$surface="x"
gene_annot$protein_type="x"
gene_annot$protein_coding="x"
gene_annot[rownames(shared_gene_rep),]$surface=as.character(shared_gene_rep$Cell.surface.)
gene_annot[rownames(shared_gene_rep),]$protein_type=as.character(shared_gene_rep$Protein.type)
gene_annot[rownames(shared_gene_rep),]$protein_coding=as.character(shared_gene_rep$Protein.coding.)
gene_annot$surface[which(gene_annot$surface=="x")]=NA
gene_annot$surface[which(gene_annot$surface=="")]=NA
gene_annot$protein_type[which(gene_annot$protein_type=="x")]=NA
gene_annot$protein_type[which(gene_annot$protein_type=="")]=NA
gene_annot$protein_coding[which(gene_annot$protein_coding=="x")]=NA
gene_annot$protein_coding[which(gene_annot$protein_coding=="")]=NA

gene_annot_ord=gene_annot[with(gene_annot, order(annot, no_brain_mean)), ]
as.matrix(bla_log2)[rownames(gene_annot_ord),]

#pick the top 20 tum-med shared genes. 
top_25=gene_annot_ord[c(1:38, 104:108),]
top_25$avg_exp=rowMeans(bla_log2[rownames(top_20),])
top_25=top_20[with(top_25, order(annot, no_brain_mean, decreasing = F)), ]

top_25$avg_exp=rowMeans(bla_log2[rownames(top_25),])

labels_for_gtex$GROUP=factor(labels_for_gtex$GROUP, levels=c("Neuroblastoma", "CNS", "PNS","Skin", "Digestive_system",
                                                             "Circulation", "Blood_immune", "Musculoskeletal", "Respiratory", 
                                                             "Urinary_system", "Female reproduction","Male_reproduction", "Fat", "Cell_line", 
                                                             "Endocrine", "Gonads"))
col_fun=colorRamp2(c(0, 5, 10), c("#edf8b1", "#7fcdbb", "#2c7fb8"))

ha = rowAnnotation(Surface = top_25$surface, Protein_coding=top_25$protein_coding, col=list(Surface=c("Yes"="red","No"="black"), Protein_coding=c("Yes"="purple", "No"="yellow")))
Heatmap(as.matrix(bla_log2)[rownames(top_25),],col=col_fun, cluster_rows = T, cluster_columns = F,
        cluster_row_slices = F, show_row_dend = F, column_split = labels_for_gtex$GROUP, row_split = top_25$annot, right_annotation = ha )


