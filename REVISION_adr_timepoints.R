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

srat_inhouse=readRDS("/home/jovyan/Neuroblastoma/srat_inhouse.rds")
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


just_bilat=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name%in%c("w10d5_1", "w10d5_2"))])
bilat_med=subset(just_bilat,
                     cells=rownames(just_bilat@meta.data)[which(just_bilat@meta.data$new_clust%in%c("SCPs", "Bridge",
                                                                                              "Sympathoblastic", "Chromaffin", 
                                                                                              "Endothelium", "Mesenchyme"))])
bone_scps=subset(scp_seurat,
                 cells=rownames(scp_seurat@meta.data)[which(scp_seurat@meta.data$new_labs=="bone SCPs")])


just_med=subset(adr_all,
                cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$new_clust%in%c("SCPs", "Bridge",
                                                                                              "Sympathoblastic", "Chromaffin"))])
bone_med_genes=intersect(rownames(bone_scps), rownames(just_med))

med_and_bone=cbind(just_med@assays$RNA@counts[bone_med_genes,], bone_scps@assays$RNA@counts[bone_med_genes,])

med_and_bone_labs=c(just_med@meta.data$new_clust, bone_scps@meta.data$new_labs)

bilat_train=trainModel(bilat_med@assays$RNA@counts, bilat_med@meta.data$new_clust, workers=NULL)

med_ag_bilat=predictSimilarity(bilat_train, med_and_bone, minGeneMatch = 0.8)
sample_ha=c(just_med@meta.data$sample_name, rep("bone SCPs", 299))
logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')

tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
names(tol8qualitative)=as.character(unique(sample_ha))

ha = rowAnnotation(Sample = factor(sample_ha, levels = c("w8", "w8d6", "w10d5_1", "w10d5_2", "w11", "w21_1", "w21_2", "bone SCPs")),
                   col = list(Sample = tol8qualitative))

Heatmap(med_ag_bilat, col=circlize::colorRamp2(seq(-5,5,length.out=length(logitCols)),logitCols),
        show_row_names = F, cluster_rows=F, cluster_columns = F, name = "Predicted\nSimilarity\n(Logit)",
        row_split=factor(med_and_bone_labs, levels=c("bone SCPs", "SCPs", "Bridge", "Sympathoblastic", "Chromaffin")), 
        column_order = c("SCPs", "Bridge", "Sympathoblastic", "Chromaffin", "Mesenchyme", "Endothelium"), left_annotation = ha)

just_med@meta.data$broad_time="x"

just_med@meta.data$broad_time[which(just_med@meta.data$sample_name%in%c("w8", "w8d6","w10d5_1", "w10d5_2", "w11"))]="early"
just_med@meta.data$broad_time[which(just_med@meta.data$sample_name%in%c("w21_1", "w21_2"))]="late"

just_med@meta.data$time_ct=paste0(just_med@meta.data$broad_time, "_", just_med@meta.data$new_clust)
table(just_med@meta.data$time_ct)

time_ct_train=trainModel(just_med@assays$RNA@counts, just_med@meta.data$time_ct, workers=NULL)

dutch_tum_only=subset(srat_dutch, cells=rownames(srat_dutch@meta.data)[grep("Tumour", srat_dutch@meta.data$idents_for_plot)])
inhouse_tum_only=subset(srat_inhouse, cells=rownames(srat_inhouse@meta.data)[grep("Tumour", srat_inhouse@meta.data$idents_for_plot)])

tum_genes=intersect(rownames(dutch_tum_only), rownames(inhouse_tum_only))
tum_together=cbind(dutch_tum_only@assays$RNA@counts[tum_genes, ], inhouse_tum_only@assays$RNA@counts[tum_genes, ])
tum_together_labs=c(paste0("PMC ", dutch_tum_only@meta.data$idents_for_plot),
                    paste0("GOSH ", inhouse_tum_only@meta.data$idents_for_plot))

time_ct_dutch=predictSimilarity(time_ct_train, dutch_tum_only@assays$RNA@counts, minGeneMatch = 0.8)
time_ct_inhouse=predictSimilarity(time_ct_train, inhouse_tum_only@assays$RNA@counts,minGeneMatch = 0.8)

time_ct_all=predictSimilarity(time_ct_train, tum_together, minGeneMatch = 0.8)
similarityHeatmap(time_ct_dutch, column_order=c("early SCPs","late SCPs", "early Bridge", "late Bridge",
                                                "early Sympathoblastic", "late Sympathoblastic",
                                                "early Chromaffin", "late Chromaffin"), row_split=dutch_tum_only@meta.data$idents_for_plot, cluster_rows=T, cluster_row_slices=F)
 
similarityHeatmap(time_ct_inhouse, column_order=c("early SCPs","late SCPs", "early Bridge", "late Bridge",
                                                  "early Sympathoblastic", "late Sympathoblastic",
                                                  "early Chromaffin", "late Chromaffin"), row_split=inhouse_tum_only@meta.data$idents_for_plot, cluster_rows=T, cluster_row_slices=F)

similarityHeatmap(time_ct_all, column_order=c("early_SCPs","late_SCPs", "early_Bridge", "late_Bridge",
                                                  "early_Sympathoblastic", "late_Sympathoblastic",
                                                  "early_Chromaffin", "late_Chromaffin"), row_split=tum_together_labs, cluster_rows=T, cluster_row_slices=F)


just_early=subset(just_med, cells=rownames(just_med@meta.data)[which(just_med@meta.data$broad_time=="early")])
just_late=subset(just_med, cells=rownames(just_med@meta.data)[which(just_med@meta.data$broad_time=="late")])


train_early=trainModel(just_early@assays$RNA@counts, just_early@meta.data$time_ct, workers=NULL)
train_late=trainModel(just_late@assays$RNA@counts, just_late@meta.data$time_ct, workers=NULL)

tum_early_ps=predictSimilarity(train_early, tum_together, minGeneMatch = 0.8)
tum_late_ps=predictSimilarity(train_late, tum_together, minGeneMatch = 0.8)

similarityHeatmap(tum_early_ps, column_order=c("early_SCPs","early_Bridge","early_Sympathoblastic",  "early_Chromaffin"), row_split=tum_together_labs)
similarityHeatmap(tum_late_ps, column_order=c("late_SCPs", "late_Bridge","late_Sympathoblastic", "late_Chromaffin"), row_split=tum_together_labs)


time_ct_marks=getMarkers(time_ct_train)

FeaturePlot(just_late,as.character(time_ct_marks[which(time_ct_marks$class=="late_Sympathoblastic"),]$gene[70:79]))
