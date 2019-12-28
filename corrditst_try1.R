
dim(rms_1_srat_no_cc@assays$RNA@scale.data)
rms1_srat=CreateSeuratObject(rms_1_mtx, min.features = 300)
dim(rms1_srat)
rms1_srat=NormalizeData(rms1_srat)
rms1_srat=ScaleData(rms1_srat, features = rownames(rms1_srat))
rms1_srat=FindVariableFeatures(rms1_srat)
rms1_srat=RunPCA(rms1_srat)
rms1_srat=RunUMAP(rms1_srat, dims = 1:20)
DimPlot(rms1_srat, group.by = "new_clust")
rms1_srat@meta.data$new_clust = rms_1_srat_no_cc@meta.data$seurat_clusters
rms1_srat
rms1_srat@active.ident= rms1_srat@meta.data$new_clust
names(rms1_srat@active.ident)= rownames(rms1_srat@meta.data)
clust_avg=AverageExpression(rms1_srat, return.seurat = T)

clust0avg=rowMeans(rms1_srat@assays$RNA@scale.data[,rownames(rms1_srat@meta.data)[which(rms1_srat@meta.data$new_clust==0)]])
femur_mtx =readMM("/home/jovyan/RMS/femur.mtx")
femur_cells=read.table("/home/jovyan/RMS/femur_cells.txt", sep = "\t", header = T)
femur_genes=read.table("/home/jovyan/RMS/femur_genes.txt", sep = "\t", header = T)
pseudotime_cells <- read.table("/home/jovyan/RMS/pseudo_nocc.csv")
rownames(femur_mtx)=femur_genes$Feature
colnames(femur_mtx)=rownames(femur_cells)
correct_order =rownames(pseudotime_cells)[order(pseudotime_cells, decreasing = F)]
femur_mtx

ref_cells =CreateSeuratObject(femur_mtx, meta.data = femur_cells)
ref_cells=NormalizeData(ref_cells)
ref_cells=ScaleData(ref_cells, features = rownames(ref_cells))

myo_nocc =subset(ref_cells, cells = rownames(pseudotime_cells))
reference_mtx =myo_nocc@assays$RNA@scale.data[,correct_order]
dim(myo_nocc)
reference_mtx
query_mtx =clust_avg@assays$RNA@scale.data
x=rbind(reference_mtx[,2], clust0avg)
x[1,]
rdist(x, metric = "correlation")
cdist(reference_mtx[,1], query_mtx[,1])

clust0avg = unname(clust0avg)
cell_1_ref=reference_mtx[,1]
cell_1_ref=unname(cell_1_ref)
cell_2_ref=reference_mtx[,2]
cell_2_ref=unname(cell_2_ref)
cell_3_ref=reference_mtx[,3]
cell_3_ref=unname(cell_3_ref)
cell_4_ref=reference_mtx[,4]
cell_4_ref=unname(cell_4_ref)
x=rbind(cell_1_ref,cell_2_ref ,clust0avg)
rdist(x, metric = "correlation")
dist(x)
sqrt((1-cor(cell_2_ref, clust0avg))/2)
cdist(clust0avg, cell_1_ref)
cor(clust0avg, cell_2_ref)
plot(clust0avg, cell_2_ref)
obs=cbind(clust0avg, cell_1_ref)
dist(obs)
rdist(x, metric= "correlation")
