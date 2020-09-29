adr_all_new =adr_all
adr_all_new@meta.data$new_clust= as.character(adr_all_new@meta.data$seurat_clusters)

adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("13"))]="left_clust"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("24"))]="bridge"
adr_all_new@meta.data$new_clust[which(adr_all_new@meta.data$new_clust%in%c("19", "20"))]="right_clust"

tr2mark=quickMarkers(adr_all_new@assays$RNA@counts, adr_all_new@meta.data$new_clust, N=30)
tr2mark_rel=tr2mark[which(tr2mark$cluster%in%c("SCPs", "left_clust", "bridge", "right_clust")),]
length(unique(tr2mark_rel$gene))

tr2mark_rel$gene=make.unique(as.character(tr2mark_rel$gene))

tr2mark_final=tr2mark_rel[(which(tr2mark_rel$gene%in%rownames(srat_org))),]

tr2mark_final$cluster=factor(tr2mark_final$cluster, levels=c("1_SCPs","3_bridge", "2_neuroblasts", "4_ganglia"))


just_tum_cells=subset(srat_tumour, idents = c(9,12,15,23))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[tr2mark_final$gene,]

col_fun=col_fun = colorRamp2(c(-3,0,3), c("blue", "black", "red") )
tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
names(tol14rainbow)= unique(srat_tumour@meta.data$unique_sample)
ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = tol14rainbow))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = tr2mark_final$cluster, col=col_fun,bottom_annotation = ha)


just_tum_cells=subset(srat_org, idents = c(0,1,2,3,4,6,7,8,9,10,12))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[tr2mark_final$gene,]
cols_4=c("#882E72","#1965B0","#4EB265","#F7EE55")
names(cols_4)=unique(just_tum_cells@meta.data$unique_sample)
just_tum_cells@meta.data$unique_sample=factor(just_tum_cells@meta.data$unique_sample, levels = as.character(freq_tab$Var1[order(freq_tab$Freq)]))
ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = cols_4))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = tr2mark_final$cluster, col=col_fun,bottom_annotation = ha)


just_tum_cells=subset(srat_inhouse, idents = c(1,5,8,12,16,25))
just_tum_cells=NormalizeData(just_tum_cells)
just_tum_cells=ScaleData(just_tum_cells, features = rownames(just_tum_cells))
dutch_mat_rel_genes=just_tum_cells@assays$RNA@scale.data[tr2mark_final$gene,]

cols_5=c("#882E72","#1965B0","#4EB265","#F7EE55","#F1932D")
names(cols_5)=unique(just_tum_cells@meta.data$unique_sample)


ha = HeatmapAnnotation(nb_sample = just_tum_cells@meta.data$unique_sample,
                       col = list(nb_sample = cols_5))
Heatmap(dutch_mat_rel_genes, cluster_rows = F, show_column_names = F, show_column_dend = F, row_split = tr2mark_final$cluster, col=col_fun,bottom_annotation = ha)

neurobl=as.character(neurobl_markers$gene)

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
  
  pdf(paste0("plts/new_plts/neurobl_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  pdf(paste0("plts/new_plts/neurobl_",mark_list,".pdf"), width = 15, height = 15)
  gg = ggplot(neurobl_tab_melted[which(neurobl_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Neuroblast markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(scp_tab_melted[which(scp_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("SCP markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  gg = ggplot(mature_tab_melted[which(mature_tab_melted$patient%in%extremes$TARGET.USI),], aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("Mature markers") +facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(gg)
  dev.off()
  
  
}

get_target_marker_plots(target_tpm, neurobl, neurobl, neurobl, "test_neurobl")

srats=list()
for (x in 1:length(cells_list)){
  message(names(cells_list[x]))
  srats[[names(cells_list[x])]]=list()
  each_sample = mtx_10x[,cells_list[[x]]]
  each_meta = meta_data[cells_list[[x]],]
  #avoid NAs in MTfract column - if NA, then total reads is 0 and will fail QC anyways
  each_meta$mt_fract[is.na(each_meta$mt_fract)]=0
  #which cells are passing QC
  passed_cells = dplyr::filter(each_meta, mt_fract<0.2, umi_count>500, nFeatures>200)
  w =nrow(passed_cells)
  message(sprintf("After filtering, keeping %d cells of %d",w,nrow(each_meta)))
  if(w<100){
    message(sprintf("Sample has only %d cells passing QC.  Too shit to continue.",w))
    next
  }
  #create seurat object and add it into a list of Seurat objects, run LR, make a umap, some featureplots, and LR heatmap
  seurat_obj =CreateSeuratObject(each_sample[,passed_cells$CellId],meta.data=each_meta[passed_cells$CellId,], min.cells = )
  srats[[names(cells_list[x])]] = process_dutch(seurat_obj)
  # ps = predictSimilarity(fit,seurat_obj@assays$RNA@counts[rownames(seurat_obj@assays$RNA@counts)%in%rownames(dat),],
  #                        classes=as.character(srats[[names(cells_list[x])]]@active.ident))
  # genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
  # col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
  # pdf(paste0("all_tum/",names(cells_list[x]),'.pdf'),width=14,height=14)
  # plot(DimPlot(srats[[names(cells_list[x])]],label=TRUE, label.size = 10, pt.size = 4)+guides(colour=FALSE))
  # plot(FeaturePlot(srats[[names(cells_list[x])]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB')))
  # plot(FeaturePlot(srats[[names(cells_list[x])]],c('MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')))
  # print(Heatmap(srats[[names(cells_list[x])]]@assays$RNA@scale.data[genes_plot,],col = col_fun, cluster_columns = F, show_column_names = F,
  #               column_split=srats[[names(cells_list[x])]]@active.ident))
  # print(similarityHeatmap(ps,
  #                         row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
  #                         column_order = c("10","17","11","7","19","0","8","1","2","4","6","12","13","5","16","14","18","3","9","15")))
  # dev.off()
}

inhouse_srats=list()
for(x in c(1:length(nb_names))){
  inhouse_srats[[nb_names[x]]] = process_inhouse(nb_list[[nb_names[x]]])
  # ps = predictSimilarity(fit,inhouse_srats[[nb_names[x]]]@assays$RNA@counts[rownames(inhouse_srats[[nb_names[x]]]@assays$RNA@counts)%in%rownames(dat),],
  #                        classes=as.character(inhouse_srats[[nb_names[x]]]@active.ident))
  # genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
  # col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
  # pdf(paste0("all_tum/",nb_names[x],'.pdf'),width=14,height=14)
  # plot(DimPlot(inhouse_srats[[nb_names[x]]],label=TRUE, label.size = 10, pt.size = 4)+guides(colour=FALSE))
  # plot(FeaturePlot(inhouse_srats[[nb_names[x]]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB')))
  # plot(FeaturePlot(inhouse_srats[[nb_names[x]]],c('MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')))
  # print(Heatmap(inhouse_srats[[nb_names[x]]]@assays$RNA@scale.data[genes_plot,],col = col_fun, cluster_columns = F, show_column_names = F,
  #               column_split=inhouse_srats[[nb_names[x]]]@active.ident))
  # print(similarityHeatmap(ps,
  #                         row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
  #                         column_order = c("10","17","11","7","19","0","8","1","2","4","6","12","13","5","16","14","18","3","9","15")))
  # dev.off()
}

length(names(srats))

levels(adr_all)=c(25,24,13,19,20,11,18,23,3,15,14,17,5,22,28,0,1,2,10,6,8,9,21,27,30,7,4,12,16,26,29,31)
adr_avg=AverageExpression(adr_all, features = rownames(adr_all), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)


genes_plot=c("SOX2", "SOX10", "PHOX2A", "PHOX2B", "PDGFRB","PECAM1","STAR", "HBB", "HBA1", "PTPRC")

WhichCells(adr_all, idents = c(25,24,13,19,20))
DimPlot(adr_all, cells.highlight = WhichCells(adr_all, idents = c(25,24,13,19,20)))

adr_med=subset(adr_all, idents=c(25,24,13,19,20))

tum_cells_inhouse=subset(srat_inhouse, idents = c(1,5,8,12,16,25))
tum_cells_dutch=subset(srat_tumour, idents = c(22,8))
tum_cells_org=subset(srat_org, idents = c(0,1,2,3,4,6,7,8,9,10,11,13))
tum_cells_inhouse$new_idents="tumour"
tum_cells_dutch$new_idents="tumour"
tum_cells_org$new_idents="tumour"
adr_med@meta.data$new_idents=as.character(adr_med@meta.data$seurat_clusters) 
combined_inhouse=merge(adr_med, tum_cells_inhouse)
inhouse_markers=quickMarkers(combined_inhouse@assays$RNA@counts, combined_inhouse@meta.data$new_idents, FDR = 0.01, N=Inf)
combined_dutch=merge(adr_med, tum_cells_dutch)
dutch_markers=quickMarkers(combined_dutch@assays$RNA@counts, combined_dutch@meta.data$new_idents, FDR = 0.01, N=Inf)
combined_org=merge(adr_med, tum_cells_org)
org_markers=quickMarkers(combined_org@assays$RNA@counts, combined_org@meta.data$new_idents, FDR = 0.01, N=Inf)

org_tc=org_markers[which(org_markers$cluster=="tumour"),]
org_tc_filt=org_tc[which(org_tc$geneFrequencySecondBest<0.2),]


inhouse_tc=inhouse_markers[which(inhouse_markers$cluster=="tumour"),]
inhouse_tc_filt=inhouse_tc[which(inhouse_tc$geneFrequencySecondBest<0.2),]

dutch_tc=dutch_markers[which(dutch_markers$cluster=="tumour"),]
dutch_tc_filt=dutch_tc[which(dutch_tc$geneFrequencySecondBest<0.2),]


lj_genes=full_join(inhouse_tc_filt, dutch_tc_filt, by='gene') %>% full_join(.,org_tc_filt , by='gene') 


joined_rel_cols=lj_genes[,c(1,3,12,21,8,17,26)]

rm_org_only_genes=as.character(dplyr::filter(joined_rel_cols, is.na(geneFrequency.x) & is.na(geneFrequency.y))$gene)

joined_rel_cols2=joined_rel_cols[which(!joined_rel_cols$gene%in%rm_org_only_genes),]

joined_rel_cols2$mean_tfidf=rowMeans(joined_rel_cols2[,c(5:7)], na.rm = T)

tum_n_de_final=joined_rel_cols2[joined_rel_cols2$mean_tfidf>0.75,]

dutch_tc_filt
org_tc_filt
inhouse_tc_filt

tum_n_de_final
colnames(tum_n_de_final)=c("gene", "gene_freq_inhouse", "gene_freq_dutch", "gene_freq_org", "tfidf_inhouse", "tfidf_dutch","tfidf_org", "mean_tfidf")


tum_n_de_final_sorted=tum_n_de_final[order(tum_n_de_final$mean_tfidf, decreasing = T),]


de_in_tg=as.data.frame(percent_above_threshold[as.character(tum_n_de_final_sorted$gene)[which(as.character(tum_n_de_final_sorted$gene)%in%names(percent_above_threshold))]])
de_in_tg$gene=rownames(de_in_tg)

colnames(de_in_tg)=c("value", "gene")

de_in_tg=de_in_tg[order(de_in_tg$value, decreasing = T), ]
de_in_tg$gene=factor(de_in_tg$gene, levels = de_in_tg$gene)

qlf3$table[as.character(de_in_tg$gene)[which(as.character(de_in_tg$gene)%in%rownames(qlf3$table))],]

ggplot(de_in_tg,aes(x=reorder(gene, value), y=value)) + geom_bar(stat="identity") + coord_flip()

no_org_lj=full_join(inhouse_tc_filt, dutch_tc_filt, by='gene')

no_org_lj$mean_tfidf=rowMeans(no_org_lj[,c(8,17)], na.rm = T)

no_org_final=no_org_lj[no_org_lj$mean_tfidf>0.85,]
no_org_final_sorted=no_org_final[order(no_org_final$mean_tfidf, decreasing = T),]



no_org_final=no_org_final[which(no_org_final$gene%in%rownames(target_log2tpm)),]

final_de=no_org_final[which(!no_org_final$gene%in%immune_sc),]

a=gsub("\\.x", ".GOSH", colnames(final_de))
b=gsub("\\.y", ".PMC", a)

colnames(final_de)=b
rownames(final_de)=1:nrow(final_de)

final_de_sorted=final_de[order(final_de$mean_tfidf, decreasing = T),]


write.table(final_de_sorted, "final_de_sorted.txt", sep = "\t", row.names = F, col.names = T, quote = F)
