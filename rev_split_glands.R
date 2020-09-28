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

chinese_fa2=read.table("chinese_fa/Fetal-Adrenal-Gland2_dge.txt", sep = ",", header = T)
chinese_fa2_meta=read.csv("chinese_fa/Fetal-Adrenal-Gland2_Anno.csv")
chinese_fa3=read.table("chinese_fa/Fetal-Adrenal-Gland3_dge.txt", sep = ",", header = T)
chinese_fa3_meta=read.csv("chinese_fa/Fetal-Adrenal-Gland3_Anno.csv")
chinese_fa4=read.table("chinese_fa/Fetal-Adrenal-Gland4_dge.txt", sep = ",", header = T)
chinese_fa4_meta=read.csv("chinese_fa/Fetal-Adrenal-Gland4_Anno.csv")

make_seurat=function(x,y){
  rownames(x)=x$X
  x=x[,c(2:ncol(x))]
  x=as.matrix(x)
  rownames(y)=y$X
  srat=CreateSeuratObject(x, meta.data = y)
  return(srat)
  
}

ch_fa2_srat=make_seurat(chinese_fa2, chinese_fa2_meta)
ch_fa3_srat=make_seurat(chinese_fa3, chinese_fa3_meta)
ch_fa4_srat=make_seurat(chinese_fa4, chinese_fa4_meta)
srat_inhouse=readRDS("/home/jovyan/Neuroblastoma/srat_inhouse.rds")
srat_dutch=readRDS('/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_dutch.rds')

#annotate the dutch samples

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

DimPlot(srat_dutch, group.by = "idents_for_plot", label=T)

baby1=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w8")])
baby2=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w8d6")])
bilat1=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w10d5_1")])
bilat2=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w10d5_2")])
tech=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w11")])
late1=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w21_1")])
late2=subset(adr_all, cells=rownames(adr_all@meta.data)[which(adr_all@meta.data$sample_name=="w21_2")])

all_train=trainModel(adr_all@assays$RNA@counts, adr_all@meta.data$new_clust)
baby1_train=trainModel(baby1@assays$RNA@counts, baby1@meta.data$new_clust, workers=NULL)
baby2_train=trainModel(baby2@assays$RNA@counts, baby2@meta.data$new_clust, workers=NULL)
bilat1_train=trainModel(bilat1@assays$RNA@counts, bilat1@meta.data$new_clust, workers=NULL)
bilat2_train=trainModel(bilat2@assays$RNA@counts, bilat2@meta.data$new_clust, workers=NULL)
tech_train=trainModel(tech@assays$RNA@counts, tech@meta.data$new_clust, workers=NULL)
late1_train=trainModel(late1@assays$RNA@counts, late1_subs@meta.data$new_clust, workers=NULL)
late2_train=trainModel(late2_subs@assays$RNA@counts, late2_subs@meta.data$new_clust, workers=NULL)
try1_train=trainModel(try1@assays$RNA@counts, try1@meta.data$new_clust, workers=NULL)

chinese2_train=trainModel(ch_fa2_srat@assays$RNA@counts, ch_fa2_srat@meta.data$CT, workers=NULL)
chinese3_train=trainModel(ch_fa3_srat@assays$RNA@counts, ch_fa3_srat@meta.data$CT, workers=NULL)
chinese4_train=trainModel(ch_fa4_srat@assays$RNA@counts, ch_fa4_srat@meta.data$CT, workers=NULL)

all_ps_dutch=predictSimilarity(all_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
all_ps_inhouse=predictSimilarity(all_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot)

baby1_ps_dutch=predictSimilarity(baby1_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
baby1_ps_inhouse=predictSimilarity(baby1_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot)

baby2_ps_dutch=predictSimilarity(try1_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
baby2_ps_inhouse=predictSimilarity(try1_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)

bilat1_ps_dutch=predictSimilarity(bilat1_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
bilat1_ps_inhouse=predictSimilarity(bilat1_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)

bilat2_ps_dutch=predictSimilarity(bilat2_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
bilat2_ps_inhouse=predictSimilarity(bilat2_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)

tech_ps_dutch=predictSimilarity(tech_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
tech_ps_inhouse=predictSimilarity(tech_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)

late1_ps_dutch=predictSimilarity(late1_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
late1_ps_inhouse=predictSimilarity(late1_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)

late2_ps_dutch=predictSimilarity(late2_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.5)
late2_ps_inhouse=predictSimilarity(late2_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.5)


ch2_ps_dutch=predictSimilarity(chinese2_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.4)
ch2_ps_inhouse=predictSimilarity(chinese2_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.4)

ch3_ps_dutch=predictSimilarity(chinese3_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.4)
ch3_ps_inhouse=predictSimilarity(chinese3_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.4)

ch4_ps_dutch=predictSimilarity(chinese4_train, srat_dutch@assays$RNA@counts, srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.4)
ch4_ps_inhouse=predictSimilarity(chinese4_train, srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot, minGeneMatch = 0.4)


pdf("sim_heatmaps.pdf" , width=8, height = 6)
print(similarityHeatmap(baby1_ps_dutch, column_title="Baby1 dutch" ))
print(similarityHeatmap(baby1_ps_inhouse,column_title="Baby1 inhouse"))


print(similarityHeatmap(baby2_ps_dutch, column_title="Baby2 dutch"))
print(similarityHeatmap(baby2_ps_inhouse, column_title="Baby2 inhouse"))

print(similarityHeatmap(bilat1_ps_dutch, column_title="Bilat1 dutch"))
print(similarityHeatmap(bilat1_ps_inhouse, column_title="Bilat1 inhouse"))

print(similarityHeatmap(bilat2_ps_dutch,column_title="Bilat2 dutch"))
print(similarityHeatmap(bilat2_ps_inhouse, column_title="Bilat2 inhouse"))

print(similarityHeatmap(tech_ps_dutch, column_title="Tech dutch"))
print(similarityHeatmap(tech_ps_inhouse, column_title="Tech inhouse"))

print(similarityHeatmap(late1_ps_dutch, column_title="Late1 dutch"))
print(similarityHeatmap(late1_ps_inhouse, column_title="Late1 inhouse"))

print(similarityHeatmap(late2_ps_dutch,column_title="Late2 dutch"))
print(similarityHeatmap(late2_ps_inhouse, column_title="Late2 inhouse"))

print(similarityHeatmap(ch2_ps_dutch,column_title="Chinese2 (12w) dutch"))
print(similarityHeatmap(ch2_ps_inhouse, column_title="Chinese2 (12w) inhouse"))

print(similarityHeatmap(ch3_ps_dutch,column_title="Chinese3 (14w) dutch"))
print(similarityHeatmap(ch3_ps_inhouse, column_title="Chinese3 (14w) inhouse"))

print(similarityHeatmap(ch4_ps_dutch,column_title="Chinese4 (12w) dutch"))
print(similarityHeatmap(ch4_ps_inhouse, column_title="Chinese4 (12w) inhouse"))


dev.off()


#add podocytes and RCC to the ref

crc_ref=readRDS("/lustre/scratch117/casm/team274/processedSeurat_RCC.RDS")
fetal_pods=readRDS("/lustre/scratch117/casm/team274/gk14/fetal_podocytes_forGerda.RDS")
genes_10x <- read.table("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t", header = F)
load("common_genes.rdata")

fetal_pds_counts=fetal_pods@assays$RNA@counts
fetal_pds_counts=fetal_pds_counts[which(rownames(fetal_pds_counts)%in%genes_10x$V1),]

rownames(fetal_pds_counts)=as.character(genes_10x$V2)[match(rownames(fetal_pds_counts), genes_10x$V1)]
rownames(fetal_pds_counts)=make.unique(rownames(fetal_pds_counts))
comm_genes=intersect(rownames(fetal_pds_counts), rownames(adr_all))
pods_comm_genes=fetal_pds_counts[comm_genes, ]
adr_comm_genes=adr_all@assays$RNA@counts[comm_genes,]

merged_cells=cbind(adr_comm_genes, pods_comm_genes)
cell_labs=c(adr_all@meta.data$new_clust, rep("Podocytes", 278))

train_with_podos=trainModel(merged_cells, cell_labs, workers=NULL)
scp_seurat@meta.data$new_labs="x"
scp_seurat@meta.data$new_labs[grep("bone",scp_seurat@meta.data$sample)]="bone SCPs"
scp_seurat@meta.data$new_labs[grep("skin",scp_seurat@meta.data$sample)]="skin SCPs"
scp_seurat@meta.data$new_labs[grep("gut",scp_seurat@meta.data$sample)]="gut SCPs"
scp_seurat@meta.data$new_labs[grep("adrenal",scp_seurat@meta.data$sample)]="adrenal SCPs"

scp_against_adr=predictSimilarity(train_with_podos, scp_seurat@assays$RNA@counts, scp_seurat@meta.data$new_labs, minGeneMatch = 0.8)

similarityHeatmap(scp_against_adr, column_order=lr_col_ord)
lr_col_ord=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme","Cortex","Leukocytes", "Endothelium", "Erythroblasts","Other", "Podocytes")


no_medulla=subset(adr_all, cells=rownames(adr_all@meta.data)[which(!adr_all@meta.data$new_clust%in%c("SCPs", "Bridge", "Chromaffin", "Sympathoblastic"))])

no_medulla_comm_genes=no_medulla@assays$RNA@counts[comm_genes,]
no_medulla_with_pods=cbind(no_medulla_comm_genes, pods_comm_genes)
no_medulla_cell_labs=cell_labs=c(no_medulla@meta.data$new_clust, rep("Podocytes", 278))

train_no_medulla=trainModel(no_medulla_with_pods, no_medulla_cell_labs, workers=NULL)

ps_no_med_dutch=predictSimilarity(train_no_medulla, as.matrix(srat_dutch@assays$RNA@counts), srat_dutch@meta.data$idents_for_plot, minGeneMatch = 0.8)
similarityHeatmap(ps_no_med_dutch)

