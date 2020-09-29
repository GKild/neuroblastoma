#process the bulks

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


target_length_norm =target_counts/(target_gene_length/1000)

scale_factor=colSums(target_length_norm)/1000000
target_tpm=target_length_norm/scale_factor
target_log2tpm=log2(target_tpm+1)

rownames(target_log2tpm)=make.unique(as.character(target_genes$external_gene_name))
rownames(target_tpm)=make.unique(as.character(target_genes$external_gene_name))
genes_above_threshold=names(rowMeans(target_log2tpm)[rowMeans(target_log2tpm)>5])
percent_above_threshold=apply(target_log2tpm, 1, function(x){
  sum(x>5)/length(x)*100
})

#prep seqc data

seqc_counts=read.table("seqc/GSE49711_SEQC_NB_TUC_G_log2.txt", sep = "\t", header = T, row.names = T)
seqc_meta=read.csv("seqc/seqc_meta.csv")

rownames(seqc_counts)=seqc_counts$X00gene_id

seqc_counts=seqc_counts[,2:ncol(seqc_counts)]
colnames(seqc_counts)=sapply(strsplit(colnames(seqc_counts), "_"), "[", 2)


seqc_unlist=apply(seqc_counts, 1, function(x){unlist(x)})

seqc_unlist_t=t(seqc_unlist)



untransf=(2^seqc_unlist_t)-1

threshold_10 <-apply(untransf, 1,
                     function(x){
                       sum(x >0)>=50 })
filtered_seqc=untransf[threshold_10,]

seqc_tpm=apply(filtered_seqc,2,function(x){
  (x/sum(x))*1000000
})


seqc_tpm_log2=log2(seqc_tpm+1)

seqc_percent_above_threshold=apply(seqc_tpm_log2, 1, function(x){
  sum(x>5)/length(x)*100
})

#prep german data

#array data DE 

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
german_clinical_data=german_clinical_data[which(german_clinical_data$Tumor.ID%in%colnames(counts_gs)),]
german_clinical_data$Age.at.Diagnosis..d.=as.numeric(as.character(german_clinical_data$Age.at.Diagnosis..d.))
counts_gs$gene_symbol=make.unique(as.character(counts_gs$gene_symbol))
rownames(counts_gs)=counts_gs$gene_symbol
counts_gs=counts_gs[,2:(ncol(counts_gs)-1)]

extreme_counts_german=counts_gs[,extremes$Tumor.ID]
bla=apply(counts_gs, 1, function(x){unlist(x)})
bla_t=t(bla)
bla_t=bla_t[,extremes$Tumor.ID]
german_counts_log2=log2(bla_t)

german_percent_above_threshold=apply(german_counts_log2, 1, function(x){
  sum(x>10)/length(x)*100
})


#get the markers 

quick_with_podos=quickMarkers(merged_cells, cell_labs, N=Inf)
quick_with_podos_rel=quick_with_podos[which(quick_with_podos$cluster%in%c("SCPs","bridge","sympathoblastic", "chromaffin","podocytes")),]
podo_marks_filt=dplyr::filter(quick_with_podos_rel, tfidf>1 & geneFrequencySecondBest <0.2)

#stratify by risk groups - target

low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

target_risks=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)

unique(target_risks$risk)

#stratify by risk group - german 

low_risk=dplyr::filter(german_clinical_data,
                       german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &!german_clinical_data$Status%in%c("relapse/progression","death of disease"))
high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified", german_clinical_data$Status%in%c("relapse/progression","death of disease"))
risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")
pt1=rbind(low_risk, high_risk, risk_4s)
intermediate_risk=german_clinical_data[-which(german_clinical_data$Tumor.ID%in%pt1$Tumor.ID),]

low_risk$risk="low_risk"
risk_4s$risk="low_risk_4s"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

german_risks=rbind(low_risk, risk_4s, intermediate_risk, high_risk)

# stratify by risk group - seqc 


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
dim(german_risks)


german_risks$patient=german_risks$Tumor.ID
seqc_risks$patient=seqc_risks$NB.ID
target_risks$patient=target_risks$TARGET_SHORT


#all samples plot

scp_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="SCPs"),]$gene)
bridge_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="bridge"),]$gene)
chromaffin_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="chromaffin"),]$gene)
sympathoblastic_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="sympathoblastic"),]$gene)
pod_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="podocytes"),]$gene)

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


ggplot(t_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >5 log2(TPM) in TARGET cohort (N=154)")


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


ggplot(seqc_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() + ggtitle("% samples a marker is expressed at >5 log2(TPM) in SEQC cohort (N=498)")




german_scp=as.data.frame(german_percent_above_threshold[scp_marks[which(scp_marks%in%names(german_percent_above_threshold))]])
german_bridge=as.data.frame(german_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(german_percent_above_threshold))]])
german_sympathoblastic=as.data.frame(german_percent_above_threshold[sympathoblastic_marks[which(sympathoblastic_marks%in%names(german_percent_above_threshold))]])
german_chromaffin=as.data.frame(german_percent_above_threshold[chromaffin_marks[which(chromaffin_marks%in%names(german_percent_above_threshold))]])
german_pods=as.data.frame(german_percent_above_threshold[pod_marks[which(pod_marks%in%names(german_percent_above_threshold))]])

german_scp$marker="scp"
german_bridge$marker="bridge"
german_sympathoblastic$marker="sympathoblastic"
german_chromaffin$marker="chromaffin"
german_pods$marker="podocyte"

german_scp$gene=rownames(german_scp)
german_bridge$gene=rownames(german_bridge)
german_sympathoblastic$gene=rownames(german_sympathoblastic)
german_chromaffin$gene=rownames(german_chromaffin)
german_pods$gene=rownames(german_pods)

rownames(german_scp)=1:nrow(german_scp)
rownames(german_bridge)=1:nrow(german_bridge)
rownames(german_sympathoblastic)=1:nrow(german_sympathoblastic)
rownames(german_chromaffin)=1:nrow(german_chromaffin)
rownames(german_pods)=1:nrow(german_pods)

colnames(german_scp)=c("value", "marker", "gene")
colnames(german_bridge)=c("value", "marker", "gene")
colnames(german_sympathoblastic)=c("value", "marker", "gene")
colnames(german_chromaffin)=c("value", "marker","gene")
colnames(german_pods)=c("value", "marker","gene")

german_all_m=rbind(german_scp, german_bridge,german_sympathoblastic,german_chromaffin, german_pods)


ggplot(german_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >10 log2(intensity) in GERMAN cohort (N=394)")

which(rownames(german_counts_log2)%in%"ERBB3.4")


t_all_m$dataset="TARGET"
german_all_m$dataset="GERMAN"
seqc_all_m$dataset="SEQC"


joined_dat=rbind(t_all_m, german_all_m, seqc_all_m)


joined_dat$marker=factor(joined_dat$marker, levels = c("podocyte","scp", "bridge", "sympathoblastic", "chromaffin"))
joined_dat$dataset=factor(joined_dat$dataset, levels = c("TARGET", "SEQC", "GERMAN"))




f1=function(dat){
  l1=sapply(as.character(unique(dat$marker)), function(x){
  dat$value[which(dat$marker==x)]
})
  
  combs=t(combn(names(l1), 2))
  comp_value <- apply(combs,1,function(y){
    wilcox.test(l1[[y[1]]], l1[[y[2]]])$p.value
  })
  
  df_out <- data.frame(combs, comp_value)
  names(df_out) <- c("mod_1", "mod_2", "comp_value")
  
  return(df_out)
}

pairwise_tg=f1(t_all_m)
pairwise_seqc=f1(seqc_all_m)
pairwise_german=f1(german_all_m)


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
german_per_risk=get_per_risk(german_counts_log2, german_risks, 10)
seqc_per_risk=get_per_risk(seqc_tpm_log2, seqc_risks, 5)

target_per_risk$cohort="target"
german_per_risk$cohort="german"
seqc_per_risk$cohort="seqc"


combined_per_risk=rbind(target_per_risk, german_per_risk, seqc_per_risk)


colnames(combined_per_risk)=c("gene", "risk", "value", "marker", "cohort")

combined_per_risk$risk=factor(combined_per_risk$risk, levels=c("low_risk","low_risk_4s", "intermediate_risk", "high_risk", "high_risk_relapse"))

combined_per_risk$marker=factor(combined_per_risk$marker, levels=c("podocyte" , "scp","bridge","sympathoblastic","chromaffin"))

combined_per_risk$cohort=factor(combined_per_risk$cohort, levels = c("target", "seqc", "german"))


#non-medullary samples (target)

leftover_ids=as.character(target_clinical_data$ICDO.Description)[-grep("adrenal|abdominal|kidney|abdomen|unknown|retroperitoneum|other", tolower(as.character(target_clinical_data$ICDO.Description)))]
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
pdf("marks_all.pdf", width = 5, height = 3, useDingbats = F)
gg=ggplot(data = joined_dat, aes(x=marker,y=value, group=dataset)) +
  geom_quasirandom(mapping = aes(x=marker,y=value, group=dataset), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  geom_pointrange(mapping = aes(x = marker, y = value, shape=dataset),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black", 
                  position=position_dodge(width=0.8)) +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Markers", y = "% samples above expression threshold")
plot(gg)
dev.off()

ggplot(data = combined_per_risk, aes(x=marker,y=value)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5)+facet_grid(vars(risk), vars(cohort))+
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))

pdf("marks_om.pdf", width = 3, height = 3, useDingbats = F)
ggplot(data = os_all_m, aes(x=marker,y=value)) +
  geom_quasirandom(mapping = aes(x=marker,y=value), alpha=0.2, cex=0.5, shape=19) +
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary", shape=22,
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black") +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


joined_no_arrays=rbind(t_all_m,seqc_all_m)


joined_no_arrays$marker=factor(joined_no_arrays$marker, levels = c("podocyte","scp", "bridge", "sympathoblastic", "chromaffin"))
joined_no_arrays$dataset=factor(joined_no_arrays$dataset, levels = c("TARGET", "SEQC"))

pdf("marks_no_german.pdf", width = 4, height = 3, useDingbats = F)
gg=ggplot(data = joined_no_arrays, aes(x=marker,y=value, group=dataset)) +
  geom_quasirandom(mapping = aes(x=marker,y=value, group=dataset), dodge.width=.8, shape=19, cex=0.5, alpha=0.2) +
  geom_pointrange(mapping = aes(x = marker, y = value, shape=dataset),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black", 
                  position=position_dodge(width=0.8)) +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x="Markers", y = "% samples above expression threshold")
plot(gg)
dev.off()

