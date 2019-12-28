source("Logistic_regression.R")
keep_features=rownames(nb_srat)[which(!rownames(nb_srat)%in%exclude_genes)]
adr@meta.data$ident2=as.character(adr@meta.data$seurat_clusters)
adr@meta.data$ident2[which(adr@meta.data$ident2%in%c("0", "8"))]="cortex_cycling"
adr@meta.data$ident2[which(adr@meta.data$ident2%in%c("1", "2", "4", "6", "12","13"))]="cortex_non_cyling"
adr@meta.data$ident2[which(adr@meta.data$ident2%in%c("3", "9"))]="mesenchyme_non_cycling"
adr@meta.data$ident2[which(adr@meta.data$ident2%in%c("15"))]="mesenchyme_cycling"
mature_srat@meta.data$ident2=0
mature_srat@meta.data$ident2[which(rownames(mature_srat@meta.data)%in%mature_non_cycling_cells)]='mature_non_cycling'
mature_srat@meta.data$ident2[which(rownames(mature_srat@meta.data)%in%mature_cycling_cells)]="mature_cycling"

srat_train_inhouse=merge(adr, mature_srat)


srat_train_dutch=merge(adr_dutch_subs, mature_dutch_subs)
srat_train_dutch@meta.data$ident2=srat_train_inhouse@meta.data$ident2
srat_train_inhouse=subset(srat_train_inhouse, features = keep_features)
srat_train_dutch=subset(srat_train_dutch, features=keep_features)

inhouse_train_mtx=srat_train_inhouse@assays$RNA@counts
dutch_train_mtx=srat_train_dutch@assays$RNA@counts

inhouse_fit=trainModel(inhouse_train_mtx, as.character(srat_train_inhouse@meta.data$ident2), maxCells = 3000)
dutch_fit=trainModel(dutch_train_mtx, as.character(srat_train_dutch@meta.data$ident2), maxCells = 3000)

inhouse_mtx = nb_srat@assays$RNA@counts
dutch_mtx=srat_tumour@assays$RNA@counts
ps=predictSimilarity(inhouse_fit, inhouse_mtx, classes=as.character(nb_srat@meta.data$seurat_clusters))
ps_dutch=predictSimilarity(dutch_fit, dutch_mtx, classes=as.character(srat_tumour@meta.data$seurat_clusters))
col_order= c("19","7", "11" ,"17", "10" ,
             "mature_non_cycling" , "mature_cycling" ,  "mesenchyme_non_cycling", 
             "mesenchyme_cycling" ,"cortex_non_cyling", "cortex_cycling" ,"16", "5", "18", "14")
pdf("inhouse_nb_lr.pdf", width = 14, height = 14)
print(similarityHeatmap(ps,
                        row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
                        column_order=col_order))
dev.off()

pdf("dutch_nb_lr.pdf", width = 14, height = 14)
print(similarityHeatmap(ps_dutch,
                        row_order = rownames(ps_dutch)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps_dutch))))],
                        column_order=col_order))
dev.off()

DimPlot(nb_srat, label = T, label.size = 6)
DimPlot(srat_tumour, label = T, label.size = 4)



