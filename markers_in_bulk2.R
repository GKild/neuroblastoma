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

#prep german data

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



get_target_marker_plots=function(exp_dat, gene_tab){
  #relevant ensembl IDs for the genes of interest
  gene_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%gene_tab$gene),]$ensembl_gene_id)
  
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
  extremes$risk=factor(extremes$risk, levels=c("low_risk_4s", "high_risk_relapse"))
  
  gene_counts=exp_dat[gene_ensembl, all$TARGET.USI]
  gene_counts=as.matrix(gene_counts)
  gene_counts=log10(gene_counts+1)
  rownames(gene_counts)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(gene_counts))]
  gene_split=gene_tab[which(gene_tab$gene%in%rownames(gene_counts)),]
 
  gene_counts=t(scale(t(gene_counts)))
  
  low_risk_4s_mean=rowMeans(gene_counts[,all$TARGET.USI[which(all$risk=="low_risk_4s")]])
  intermediate_risk_mean=rowMeans(gene_counts[,all$TARGET.USI[which(all$risk=="intermediate_risk")]])
  high_risk_mean=rowMeans(gene_counts[,all$TARGET.USI[which(all$risk=="high_risk")]])
  high_risk_relapse_mean=rowMeans(gene_counts[,all$TARGET.USI[which(all$risk=="high_risk_relapse")]])
  
  all_means=cbind(low_risk_4s_mean,intermediate_risk_mean,high_risk_mean, high_risk_relapse_mean)
  rownames(all_means)=rownames(gene_counts)
  
  all_means=as.matrix(all_means)

  
  #col_fun = colorRamp2(c(-3, 0, 3), c("purple", "black", "yellow"))
  hm=Heatmap(gene_counts, cluster_columns = F, show_row_dend = F, column_split = as.character(all$risk[match(colnames(gene_counts), all$TARGET.USI)]))
  return(hm)                                                                                            
  
  
}

x=get_target_marker_plots(target_tpm, a)
draw(x)


get_german_marker_plots=function(exp_dat, gene_tab){
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
  
  gene_counts=exp_dat[which(counts_gs$gene_symbol%in%gene_tab$gene),c(all$Tumor.ID, "gene_symbol")]
  all=all[which(all$Tumor.ID%in%colnames(gene_counts)),]
  gene_counts$gene_symbol=make.unique(as.character(gene_counts$gene_symbol))
  rownames(gene_counts)=gene_counts$gene_symbol
  gene_counts=gene_counts[,1:ncol(gene_counts)-1]
  gene_counts=gene_counts[which(rownames(gene_counts)%in%gene_tab$gene),]
  gene_counts=log10(gene_counts)
  gene_split=gene_tab[which(gene_tab$gene%in%rownames(gene_counts)),]
  gene_counts=t(scale(t(gene_counts)))
  
  low_risk_mean=rowMeans(gene_counts[,all$Tumor.ID[which(all$risk=="low_risk")]])
  inter_no_event_mean=rowMeans(gene_counts[,all$Tumor.ID[which(all$risk=="inter_no_event")]])
  inter_death_or_relapse_mean=rowMeans(gene_counts[,all$Tumor.ID[which(all$risk=="inter_death_or_relapse")]])
  high_risk_mean=rowMeans(gene_counts[,all$Tumor.ID[which(all$risk=="high_risk")]])
  risk_4s_mean=rowMeans(gene_counts[,all$Tumor.ID[which(all$risk=="4s")]])
  
  all=cbind(low_risk_mean,inter_no_event_mean, inter_death_or_relapse_mean,high_risk_mean,risk_4s_mean)
  rownames(all)=rownames(gene_counts)
  
  all=as.matrix(all)
  hm=Heatmap(all, cluster_columns = F, cluster_rows = T, show_row_dend = F, column_s = c("low_risk_mean", "risk_4s_mean", 
                                                                           "inter_no_event_mean", "inter_death_or_relapse_mean", "high_risk_mean"))
  return(hm)
}

gene=c("PHOX2A", "PHOX2B", "SOX2", "SOX10")
split=rep("split", 4)

df=data.frame(gene=gene, split=split)

x=get_german_marker_plots(counts_gs, a)
draw(x)


tabdat=read.table("de_genes_tfidf_na.txt", sep = "\t", header = T)

a=data.frame(gene=tabdat[c(1:100),1], split="split")

