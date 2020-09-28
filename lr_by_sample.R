srat_inhouse=readRDS("srat_inhouse.rds")
srat_tumour=readRDS("srat_dutch.rds")
load("adr_all_podo_train.RData")

com_rows=rownames(fit$podocytes$glmnet.fit$beta)


srat_inhouse@meta.data$GOSH_ID=as.character(srat_inhouse@meta.data$orig.ident)

srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="GOSH014"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="GOSH023"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="GOSH019"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="GOSH025"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="GOSH021"

sample_labs=c(srat_tumour@meta.data$unique_sample[which(srat_tumour@meta.data$seurat_clusters%in%c(22,8))],
              srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))])
dataset_labs=c(rep("GOSH", 4),rep("PMC", 5))

dutch_mtx=srat_tumour@assays$RNA@counts[com_rows,which(srat_tumour@meta.data$seurat_clusters%in%c(22,8))]
inhouse_mtx=srat_inhouse@assays$RNA@counts[com_rows, which(srat_inhouse@meta.data$seurat_clusters%in%c(1,16,8,25,12,5))]

all_tum_cells=cbind(dutch_mtx,inhouse_mtx)

tum_ps=predictSimilarity(fit, all_tum_cells, sapply(strsplit(sample_labs, "_"), "[",1))

similarityHeatmap(tum_ps)

logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
cols = circlize::colorRamp2(seq(-5,5,length.out=length(logitCols)),logitCols)
colnames(tum_ps)=c("Renal podocytes", "SCPs", "Other", "Bridge","Chromaffin", "Leukocytes","Sympathoblast", "Mesenchyme",
                   "Erythroblasts","Edothelium", "Cortex")
rownames(tum_ps)=c("GOSH014 (N=1)", "GOSH019 (N=24)", "GOSH021 (N=425)", "GOSH023 (N=21)", "GOSH025 (N=1295)",
                   "NB060 (N=67)" ,  "NB086 (N=5)",   "NB124 (N=8)",   "NB125 (N=425)",   "NB130 (N=11)",
                   "NB132 (N=3)",   "NB152 (N=11)")
pdf("blah.pdf", useDingbats = F, width = 7, height = 10)
print(Heatmap(tum_ps, col=cols, cluster_rows=F, show_row_names = T, cluster_columns = F, row_split = dataset_labs)
dev.off()

rownames(tum_ps)=c("GOSH014 (N=1)", "GOSH019 (N=24)", "GOSH021 (N=425)", "GOSH023 (N=21)", "GOSH025 (N=1295)",
                   "NB060 (N=67)" ,  "NB086 (N=5)",   "NB124 (N=8)",   "NB125 (N=425)",   "NB130 (N=11)",
                   "NB132 (N=3)",   "NB152 (N=11)")

Heatmap(tum_ps, col=cols, cluster_rows=F, show_row_names = T, cluster_columns = F, row_split = dataset_labs, name="Predicted similarity (logit)")



tum_ps2=tum_ps[c(2,3,4,5,6,8,9,10,12),]

rownames(tum_ps2)= c("GOSH019 (N=22/24)", "GOSH021 (N=404/425)", "GOSH023 (N=14/21)", "GOSH025 (N=1294/1295)",
                     "NB060 (N=5/67)" , "NB124 (N=3/8)",   "NB125 (N=85/425)",   "NB130 (N=0/11)","NB152 (N=3/11)")
colnames(tum_ps2)=c("Renal podocytes", "SCPs", "Other", "Bridge","Chromaffin", "Leukocytes","Sympathoblast", "Mesenchyme",
                    "Erythroblasts","Endothelium", "Cortex")
col_ord=c("SCPs", "Bridge", "Chromaffin" , "Sympathoblast", "Mesenchyme",  "Cortex" , "Leukocytes", "Endothelium",
          "Erythroblasts" , "Other", "Renal podocytes")
Heatmap(tum_ps2, col=cols, cluster_rows=F, show_row_names = T, cluster_columns = F,
        row_split = dataset_labs, name="Predicted similarity (logit)", column_order = col_ord)


colnames(tum_ps2)
length(col_ord)
dim(tum_ps2)
