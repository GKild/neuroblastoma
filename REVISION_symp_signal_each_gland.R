source("/home/jovyan/Dediff/logisticRegression.R")
workers=NULL
adr_all=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")
adr_all@meta.data$sample_name="x"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]="w21_1"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(adr_all@active.ident), "_"), "[", 6)%in%c("Adr8710634", "Adr8710635"))]="w21_2"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident=="babyAdrenal1")]="w8"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident=="babyAdrenal2")]="w8d6"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5388STDY7717452","5388STDY7717453",
                                                                      "5388STDY7717454","5388STDY7717455"))]="w10d5_1"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5388STDY7717456","5388STDY7717458",
                                                                      "5388STDY7717459"))]="w10d5_2"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident%in%c("5698STDY7839907","5698STDY7839909",
                                                                      "5698STDY7839917"))]="w11"

adr_all@meta.data$new_clust=as.character(adr_all@meta.data$seurat_clusters)
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("13"))]="Chromaffin"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("24"))]="Bridge"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("19", "20"))]="Sympathoblastic"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("3", "15", "14","17"))]="Endothelium"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("23", "11", "18"))]="Mesenchyme"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("22", "5", "28", "0", "1", "2", "10", "6", "8", "9","21"))]="Cortex"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("16", "26"))]="Leukocytes"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("30", "7", "4", "12"))]="Erythroblasts"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("27", "31", "29"))]="Other"

srat_inhouse=readRDS("/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_inhouse.rds")
srat_dutch=readRDS('/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_dutch.rds')


srat_dutch@meta.data$idents_for_plot="x"

srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(13))]="Tumour cluster 1"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(24))]="Tumour cluster 2"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(15))]="Tumour cluster 3"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(10))]="Tumour cluster 4"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(17))]="Tumour cluster 5"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(6))]="Schwannian stroma"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(7,29))]="Endothelium"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(28,18,23,5,26,0,8,20,9,3,22,25,4,19,21,27))]="Leukocytes"
srat_dutch@meta.data$idents_for_plot[which(srat_dutch@meta.data$seurat_clusters%in%c(1,2,11,12,14,16))]="Mesenchyme"

srat_inhouse@meta.data$new_idents="x"
srat_inhouse@meta.data$idents_for_plot="x"

srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))]="tumour"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                    17,15,18,3,6,4,21,22,24))]="leukocytes"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="endothelium"
srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="mesenchyme"

srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25))]="Tumour cluster 1"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(12))]="Tumour cluster 2"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(5))]="Tumour cluster 3"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(20,14,29,19,0,7,2,9,23,13,
                                                                                         17,15,18,3,6,4,21,22,24))]="Leukocytes"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(10))]="Endothelium"
srat_inhouse@meta.data$idents_for_plot[which(srat_inhouse@meta.data$seurat_clusters%in%c(11))]="Mesenchyme"

srat_inhouse@meta.data$GOSH_ID=as.character(srat_inhouse@meta.data$orig.ident)
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="GOSH014"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="GOSH023"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="GOSH019"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="GOSH025"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="GOSH021"



dutch_tum_only=subset(srat_dutch, cells=rownames(srat_dutch@meta.data)[grep("Tumour", srat_dutch@meta.data$idents_for_plot)])
dutch_tum_only=subset(dutch_tum_only, idents = 15, invert=T)
inhouse_tum_only=subset(srat_inhouse, cells=rownames(srat_inhouse@meta.data)[grep("Tumour", srat_inhouse@meta.data$idents_for_plot)])

dutch_tum_only@meta.data$risk="high risk"
dutch_tum_only@meta.data$risk[which(dutch_tum_only@meta.data$unique_sample%in%c("NB152_Fresh_biopsy"))]="low/intermediate risk"

inhouse_tum_only@meta.data$risk="low/intermediate risk"
inhouse_tum_only@meta.data$risk[which(inhouse_tum_only@meta.data$GOSH_ID=="GOSH021")]="high risk"

tum_risk_labs=c(dutch_tum_only@meta.data$risk, inhouse_tum_only@meta.data$risk)
tum_genes=intersect(rownames(dutch_tum_only), rownames(inhouse_tum_only))
tum_together=cbind(dutch_tum_only@assays$RNA@counts[tum_genes, ], inhouse_tum_only@assays$RNA@counts[tum_genes, ])
tum_together_labs=c(paste0("PMC ", dutch_tum_only@meta.data$idents_for_plot),
                    paste0("GOSH ", inhouse_tum_only@meta.data$idents_for_plot))

dutch_no_15=subset(srat_dutch, idents = 15, invert=T)

dutch_inh_together=cbind(srat_inhouse@assays$RNA@counts[tum_genes, ], dutch_no_15@assays$RNA@counts[tum_genes,])
dutch_inh_together_labs=c(srat_inhouse@meta.data$idents_for_plot, dutch_no_15@meta.data$idents_for_plot)

dutch_inh_together_labs[grep("Tumour", dutch_inh_together_labs)]="Tumour"

dutch_inh_together_train=trainModel(dutch_inh_together, dutch_inh_together_labs, workers=NULL)


together_against_tum=predictSimilarity(dutch_inh_together_train, tum_together, logits = F)

similarityHeatmap(together_against_tum)

adr_tum_comm_genes=intersect(rownames(adr_all), rownames(tum_together))

week8=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w8")])
week8d6=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w8d6")])
week8d6=subset(week8d6, cells=rownames(week8d6@meta.data)[which(week8d6@meta.data$new_clust!="Erythroblasts")])
week10d5=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name%in%c("w10d5_1","w10d5_2"))])
week11=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w11")])
week21=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name%in%c("w21_1", "w21_2"))])
week21=subset(week21, cells=rownames(week21@meta.data)[which(week21@meta.data$new_clust!="Bridge")])


week8_train=trainModel(week8@assays$RNA@counts[adr_tum_comm_genes,], week8@meta.data$new_clust, workers=NULL)
week8d6_train=trainModel(week8d6@assays$RNA@counts[adr_tum_comm_genes,], week8d6@meta.data$new_clust, workers=NULL)
week10d5_train=trainModel(week10d5@assays$RNA@counts[adr_tum_comm_genes,], week10d5@meta.data$new_clust, workers=NULL)
week11_train=trainModel(week11@assays$RNA@counts[adr_tum_comm_genes,], week11@meta.data$new_clust, workers=NULL)
week21_train=trainModel(week21@assays$RNA@counts[adr_tum_comm_genes,], week21@meta.data$new_clust, workers=NULL)

week8_pred=predictSimilarity(week8_train, tum_together, logits = F)
week8d6_pred=predictSimilarity(week8d6_train, tum_together, logits = F)
week10d5_pred=predictSimilarity(week10d5_train, tum_together, logits = F)
week11_pred=predictSimilarity(week11_train, tum_together, logits = F)
week21_pred=predictSimilarity(week21_train, tum_together, logits = F)


cells=rownames(week8_pred)
pos_ctrl=unname(together_against_tum[,3])
neg_ctrl=unname(together_against_tum[,5])
week8d6_symp=unname(week8d6_pred[,3])
week10d5_symp=unname(week10d5_pred[,6])
week11_symp=unname(week11_pred[,7])
week21_symp=unname(week21_pred[,4])

all_symp_df=data.frame(cells, pos_ctrl, neg_ctrl, week8d6_symp, week10d5_symp, week11_symp, week21_symp)

all_symp_df_melted=melt(all_symp_df, id.vars = "cells")

ggplot(data = all_symp_df_melted, aes(x=variable,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+
  stat_summary(mapping = aes(x = variable, y = value),
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

adr_clust_train=trainModel(adr_all@assays$RNA@counts[tum_genes,], adr_all@meta.data$seurat_clusters, workers=NULL)
bilat_med=subset(week10d5, idents = c(25,24,13,19,20))
tech_med=subset(week11, idents = c(25,24,13,19,20))
tech_clust_train=trainModel(tech_med@assays$RNA@counts[tum_genes,], tech_med@meta.data$seurat_clusters, workers=NULL)
bilat_clust_train=trainModel(bilat_med@assays$RNA@counts[tum_genes,], bilat_med@meta.data$seurat_clusters, workers=NULL)
clust_tum_pred=predictSimilarity(adr_clust_train, tum_together)
bilat_clust_pred=predictSimilarity(bilat_clust_train, tum_together,logits =F)
tech_clust_pred=predictSimilarity(tech_clust_train, tum_together,logits =F)
similarityHeatmap(clust_tum_pred)
similarityHeatmap(bilat_clust_pred, row_split=tum_risk_labs)
similarityHeatmap(tech_clust_pred,row_split=tum_risk_labs)
row_split=tum_risk_labs


table(adr_all@meta.data$sample_name[which(adr_all@meta.data$seurat_clusters==20)])

get_symp_signal=function(srat){
  remove_small_clusts=function(x){
    keep_idents=names(which(table(x@meta.data$new_clust)>=10))
    x=subset(x, cells=rownames(x@meta.data)[which(x@meta.data$new_clust%in%keep_idents)])
    return(x)
  }
  week8=subset(srat, cells=rownames(srat@meta.data)[which(srat@meta.data$sample_name=="w8")])
  week8=remove_small_clusts(week8)
  week8d6=subset(srat, cells=rownames(srat@meta.data)[which(srat@meta.data$sample_name=="w8d6")])
  week8d6=remove_small_clusts(week8d6)
  week10d5=subset(srat, cells=rownames(srat@meta.data)[which(srat@meta.data$sample_name%in%c("w10d5_1","w10d5_2"))])
  week10d5=remove_small_clusts(week10d5)
  week11=subset(srat, cells=rownames(srat@meta.data)[which(srat@meta.data$sample_name=="w11")])
  week11=remove_small_clusts(week11)
  week21=subset(srat, cells=rownames(srat@meta.data)[which(srat@meta.data$sample_name%in%c("w21_1", "w21_2"))])
  week21=remove_small_clusts(week21)
  
  week8_train=trainModel(week8@assays$RNA@counts[adr_tum_comm_genes,], week8@meta.data$new_clust, workers=NULL)
  week8d6_train=trainModel(week8d6@assays$RNA@counts[adr_tum_comm_genes,], week8d6@meta.data$new_clust, workers=NULL)
  week10d5_train=trainModel(week10d5@assays$RNA@counts[adr_tum_comm_genes,], week10d5@meta.data$new_clust, workers=NULL)
  week11_train=trainModel(week11@assays$RNA@counts[adr_tum_comm_genes,], week11@meta.data$new_clust, workers=NULL)
  week21_train=trainModel(week21@assays$RNA@counts[adr_tum_comm_genes,], week21@meta.data$new_clust, workers=NULL)
  
  week8_pred=predictSimilarity(week8_train, tum_together, logits = F)
  week8_pred=as.data.frame(week8_pred)
  colnames(week8_pred)=paste0("week8_",colnames(week8_pred))
  week8d6_pred=predictSimilarity(week8d6_train, tum_together, logits = F)
  week8d6_pred=as.data.frame(week8d6_pred)
  colnames(week8d6_pred)=paste0("week8d6_",colnames(week8d6_pred))
  week10d5_pred=predictSimilarity(week10d5_train, tum_together, logits = F)
  week10d5_pred=as.data.frame(week10d5_pred)
  colnames(week10d5_pred)=paste0("week10d5_",colnames(week10d5_pred))
  week11_pred=predictSimilarity(week11_train, tum_together, logits = F)
  week11_pred=as.data.frame(week11_pred)
  colnames(week11_pred)=paste0("week11_",colnames(week11_pred))
  week21_pred=predictSimilarity(week21_train, tum_together, logits = F)
  week21_pred=as.data.frame(week21_pred)
  colnames(week21_pred)=paste0("week21_",colnames(week21_pred))
  
  all_df=cbind(week8_pred, week8d6_pred, week10d5_pred, week11_pred, week21_pred)
  just_symp=all_df[,grep("Sympathoblastic", colnames(all_df))]
  return(just_symp)
}

symp_adr_all=get_symp_signal(adr_all)
symp_adr_all_pos_neg=cbind(symp_adr_all, pos_ctrl, neg_ctrl)
symp_adr_all_pos_neg_melted=melt(symp_adr_all_pos_neg)
ggplot(data = symp_adr_all_pos_neg_melted, aes(x=variable,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+
  stat_summary(mapping = aes(x = variable, y = value),
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


adr_all_19=subset(adr_all, idents = 20, invert=T)
adr_all_20=subset(adr_all, idents = 19, invert=T)

symp_adr_19=get_symp_signal(adr_all_19)
symp_adr_19_pos_neg=cbind(symp_adr_19, pos_ctrl, neg_ctrl)
symp_adr_19_pos_neg_melted=melt(symp_adr_19_pos_neg)

ggplot(data = symp_adr_19_pos_neg_melted, aes(x=variable,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+
  stat_summary(mapping = aes(x = variable, y = value),
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


symp_adr_20=get_symp_signal(adr_all_20)
symp_adr_20_pos_neg=cbind(symp_adr_20, pos_ctrl, neg_ctrl)
symp_adr_20_pos_neg_melted=melt(symp_adr_20_pos_neg)

ggplot(data = symp_adr_20_pos_neg_melted, aes(x=variable,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+
  stat_summary(mapping = aes(x = variable, y = value),
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

med_just_tfs=subset(just_med, features = as.character(tf_table$Symbol[which(tf_table$Symbol%in%rownames(just_med))]))


symp_de=FindMarkers(med_just_tfs, ident.1 = 19, ident.2 = 20)
symp_de_filt=symp_de[symp_de$p_val_adj<0.05,]

rownames(symp_de_filt[order(symp_de_filt$avg_logFC, decreasing = T),])
symp_de_tfs=subset(just_med, features = rownames(symp_de_filt[order(symp_de_filt$avg_logFC, decreasing = T),]), idents = c(19,20))

symp_de_tfs=ScaleData(symp_de_tfs, features = rownames(symp_de_tfs))

col_fun = colorRamp2(c(-3,-1,0,1,3), rev(brewer.pal(n = 5, name = "RdYlBu")))
Heatmap(symp_de_tfs@assays$RNA@scale.data, col=col_fun, show_column_names = F, cluster_column_slices = F,
        cluster_columns = T, column_split = symp_de_tfs@meta.data$seurat_clusters, show_column_dend = F,
        show_row_dend = F, row_order =  rownames(symp_de_filt[order(symp_de_filt$avg_logFC, decreasing = T),]),
        row_split = c(rep("pos_19", 13), rep("pos_20", 3)))

rowMeans(symp_de_tfs@assays$RNA@counts)
