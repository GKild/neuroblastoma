#prep the target data
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

# prep the german data 

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
german_clinical_data=german_clinical_data[which(german_clinical_data$Tumor.ID%in%colnames(counts_gs)),]
german_clinical_data$Age.at.Diagnosis..d.=as.numeric(as.character(german_clinical_data$Age.at.Diagnosis..d.))


#function to produce the marker expression level plots for target data

get_target_marker_plots=function(exp_dat, scp_markers, neurobl_markers, mature_markers, mark_list){
  #relevant ensembl IDs for the genes of interest
  scp_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%scp_markers),]$ensembl_gene_id)
  neurobl_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%neurobl_markers),]$ensembl_gene_id)
  mature_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%mature_markers),]$ensembl_gene_id)
  
  low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
  high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
  intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
  high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))
  
  low_risk_4s$risk="low_risk_4s"
  high_risk_relapse$risk="high_risk_relapse"
  intermediate_risk$risk="intermediate_risk"
  high_risk$risk="high_risk"
  
  extremes=rbind(low_risk_4s, high_risk_relapse)
  all=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
  
  scp_tab=exp_dat[scp_ensembl, all$TARGET.USI]
  neurobl_tab=exp_dat[neurobl_ensembl, all$TARGET.USI]
  mature_tab=exp_dat[mature_ensembl, all$TARGET.USI]

  
  rownames(scp_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(scp_tab))]
  rownames(neurobl_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(neurobl_tab))]
  rownames(mature_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(mature_tab))]
  
  scp_tab=as.data.frame(t(scp_tab))
  neurobl_tab=as.data.frame(t(neurobl_tab))
  mature_tab=as.data.frame(t(mature_tab))
  
  scp_tab$risk=all$risk
  neurobl_tab$risk=all$risk
  mature_tab$risk=all$risk
  
  scp_tab$patient=rownames(scp_tab)
  neurobl_tab$patient=rownames(neurobl_tab)
  mature_tab$patient=rownames(mature_tab)
  
  
  
  scp_tab_melted=melt(scp_tab, id.vars=c("patient", "risk"))
  scp_tab_melted$risk=factor(scp_tab_melted$risk, levels = c("low_risk_4s","intermediate_risk", "high_risk", "high_risk_relapse"))
  
  neurobl_tab_melted=melt(neurobl_tab, id.vars=c("patient", "risk"))
  neurobl_tab_melted$risk=factor(neurobl_tab_melted$risk, levels = c("low_risk_4s","intermediate_risk", "high_risk" ,"high_risk_relapse"))
  
  mature_tab_melted=melt(mature_tab, id.vars=c("patient", "risk"))
  mature_tab_melted$risk=factor(mature_tab_melted$risk, levels = c("low_risk_4s","intermediate_risk", "high_risk", "high_risk_relapse"))
  
  scp_tab_melted$value=log10(scp_tab_melted$value)
  neurobl_tab_melted$value=log10(neurobl_tab_melted$value)
  mature_tab_melted$value=log10(mature_tab_melted$value)
  
  scp_tab_melted$variable=factor(scp_tab_melted$variable, levels =  as.character(scp_markers[which(scp_markers%in%scp_tab_melted$variable)]))
  neurobl_tab_melted$variable=factor(neurobl_tab_melted$variable, levels =  as.character(unique(neurobl_tab_melted$variable)[match(neurobl_markers,unique(neurobl_tab_melted$variable))]))

  pdf(paste0("plts/new_plts/target_all_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  pdf(paste0("plts/new_plts/target_extreme_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted[which(neurobl_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted[which(scp_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted[which(mature_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  
}

# function to produce the plots for german data

get_german_marker_plots= function(exp_dat, scp_markers, neurobl_markers, mature_markers, mark_list){
  low_risk=dplyr::filter(german_clinical_data,
                         german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &!german_clinical_data$Status%in%c("relapse/progression","death of disease"))
  high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified", german_clinical_data$Status%in%c("relapse/progression","death of disease"))
  risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")
  pt1=rbind(low_risk, high_risk, risk_4s)
  pt2=german_clinical_data[-which(german_clinical_data$Tumor.ID%in%pt1$Tumor.ID),]
  inter_death_or_relapse=dplyr::filter(pt2,Status%in%c("relapse/progression","death of disease"))
  inter_no_event=pt2[-which(pt2$Tumor.ID%in%inter_death_or_relapse$Tumor.ID),]
  low_risk$risk_level="low_risk"
  inter_no_event$risk_level="inter_no_event"
  inter_death_or_relapse$risk_level="inter_death_or_relapse"
  high_risk$risk_level="high_risk"
  risk_4s$risk_level="4s"
  
  extremes=rbind(low_risk, risk_4s, high_risk)
  all=rbind(low_risk, risk_4s, inter_no_event, inter_death_or_relapse, high_risk)
  
  scp_tab=exp_dat[which(counts_gs$gene_symbol%in%scp_markers),c(all$Tumor.ID, "gene_symbol")]
  neurobl_tab=exp_dat[which(counts_gs$gene_symbol%in%neurobl_markers),c(all$Tumor.ID, "gene_symbol")]
  mature_tab=exp_dat[which(counts_gs$gene_symbol%in%mature_markers),c(all$Tumor.ID, "gene_symbol")]
  
  scp_tab$gene_symbol = make.unique(as.character(scp_tab$gene_symbol))
  neurobl_tab$gene_symbol=make.unique(as.character(neurobl_tab$gene_symbol))
  mature_tab$gene_symbol=make.unique(as.character(mature_tab$gene_symbol))
  
  rownames(scp_tab)=scp_tab$gene_symbol
  rownames(neurobl_tab)=neurobl_tab$gene_symbol
  rownames(mature_tab)=mature_tab$gene_symbol
  scp_tab=scp_tab[,1:ncol(scp_tab)-1]
  neurobl_tab=neurobl_tab[,1:ncol(neurobl_tab)-1]
  mature_tab=mature_tab[,1:ncol(mature_tab)-1]
  
  scp_transposed=as.data.frame(t(scp_tab))
  neurobl_transposed=as.data.frame(t(neurobl_tab))
  mature_transposed=as.data.frame(t(mature_tab))
  
  scp_transposed$patient=rownames(scp_transposed)
  neurobl_transposed$patient=rownames(neurobl_transposed)
  mature_transposed$patient=rownames(mature_transposed)
  
  
  scp_transposed$risk= all$risk_level
  neurobl_transposed$risk=all$risk_level
  mature_transposed$risk=all$risk_level
  
  scp_melted=melt(scp_transposed, id.vars=c("patient", "risk"))
  neurobl_melted=melt(neurobl_transposed,id.vars=c("patient", "risk"))
  mature_melted=melt(mature_transposed,id.vars=c("patient", "risk"))
  
  scp_melted$value=log10(scp_melted$value)
  neurobl_melted$value=log10(neurobl_melted$value)
  mature_melted$value=log10(mature_melted$value)
  
  scp_melted$risk=factor(scp_melted$risk, levels = c("low_risk", "4s","inter_no_event","inter_death_or_relapse","high_risk"))
  neurobl_melted$risk=factor(neurobl_melted$risk, levels = c("low_risk", "4s","inter_no_event","inter_death_or_relapse","high_risk"))
  mature_melted$risk=factor(mature_melted$risk, levels = c("low_risk", "4s","inter_no_event","inter_death_or_relapse","high_risk"))
  
  scp_melted$variable=sapply(strsplit(as.character(scp_melted$variable), "\\."),"[",1)
  neurobl_melted$variable=sapply(strsplit(as.character(neurobl_melted$variable), "\\."),"[",1)
  mature_melted$variable=sapply(strsplit(as.character(mature_melted$variable), "\\."),"[",1)
  
  scp_melted$variable=factor(scp_melted$variable, levels =  as.character(scp_markers[which(scp_markers%in%scp_melted$variable)]))
  neurobl_melted$variable=factor(neurobl_melted$variable, levels =  as.character(neurobl_markers[which(neurobl_markers%in%neurobl_melted$variable)]))
  
  pdf(paste0("plts/new_plts/german_all_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  pdf(paste0("plts/new_plts/german_extreme_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_melted[which(neurobl_melted$patient%in%extremes$Tumor.ID),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_melted[which(scp_melted$patient%in%extremes$Tumor.ID),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_melted[which(mature_melted$patient%in%extremes$Tumor.ID),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  
}

get_target_marker_plots(target_tpm, scp_markers, neurobl_markers, mature_markers, "no_mesenchyme")
get_german_marker_plots(counts_gs, scp_markers, neurobl_markers, mature_markers, "no_mesenchyme")
#get the marker list with mesenchyme
x=subset(bilateral, idents = c(8,12,22))
x_nocc=rownames(x@meta.data)[which(x@meta.data$Phase=="G1")]
x=subset(x, cells = x_nocc)
z=mature_small_srat
z@meta.data$seurat_clusters=1

u=merge(x,z)

q=quickMarkers(u@assays$RNA@counts, u@meta.data$seurat_clusters, N=100)
f=dplyr::filter(q, geneFrequency>0.5, geneFrequencyOutsideCluster<0.05)

scp_mesen=as.character(f$gene[which(f$cluster==22)])
neurobl_mesen=as.character(f$gene[which(f$cluster==12)])
mature_mesen=as.character(f$gene[which(f$cluster==1)])

get_target_marker_plots(target_tpm, scp_mesen, neurobl_mesen, mature_mesen, "mesenchyme")
get_german_marker_plots(counts_gs, scp_mesen, neurobl_mesen, mature_mesen, "mesenchyme")

#get the marker list with everything

bilat_nocc=rownames(bilateral@meta.data)[which(bilateral@meta.data$Phase=="G1")]
bilat2=subset(bilateral, cells = bilat_nocc)
mature_small=mature_small_srat
mature_small@meta.data$seurat_clusters=25

bilat_and_mature=merge(bilat2,mature_small)

d=quickMarkers(bilat_and_mature@assays$RNA@counts, bilat_and_mature@meta.data$seurat_clusters, N=100)
e=dplyr::filter(d, geneFrequency>0.5, geneFrequencyOutsideCluster<0.05)

scp_bilat=as.character(e$gene[which(e$cluster==22)])
neurobl_bilat=as.character(e$gene[which(e$cluster==12)])
mature_bilat=as.character(e$gene[which(e$cluster==25)])

get_target_marker_plots(target_tpm, scp_bilat, neurobl_bilat, mature_bilat, "bilat_and_mature")
get_german_marker_plots(counts_gs, scp_bilat, neurobl_bilat, mature_bilat, "bilat_and_mature")


