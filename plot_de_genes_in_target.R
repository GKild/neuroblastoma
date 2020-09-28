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
  mature_tab_melted$variable=factor(mature_tab_melted$variable, levels =  as.character(unique(mature_tab_melted$variable)[match(mature_markers,unique(mature_tab_melted$variable))]))
  
  pdf(paste0("plts/new_plts/target_all_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Dutch DE genes")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Inhouse DE genes") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Dutch organoid DE genes") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  pdf(paste0("plts/new_plts/target_extreme_groups_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted[which(neurobl_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Dutch DE genes")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted[which(scp_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Inhouse DE genes") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted[which(mature_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Dutch organoid DE genes") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  
}

scp_markers=as.character(inhouse_de$Var1)
neurobl_markers=as.character(dutch_de$Var1)
mature_markers=as.character(org_de$Var1)

get_target_marker_plots(target_tpm, scp_markers, neurobl_markers, mature_markers ,"diff_exp")
