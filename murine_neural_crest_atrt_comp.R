#atrt mouse NC comp

#trunk E9.5
E10.5_post=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E10.5_post-otic_Sox10_counts.txt", header = T, row.names = 1)
E10.5_tail=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E10.5_tail_Wnt1_counts.txt", header = T, row.names = 1)
E9.5_anterior=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E9.5_anterior_Sox10_counts.txt", header = T, row.names = 1)
E9.5_trunk=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E9.5_trunk_Wnt1_counts.txt", header = T, row.names = 1)
E9.5_posterior=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E9.5_posterior_Sox10_counts.txt", header = T, row.names = 1)
E8.5=read.table("/home/jovyan/Neuroblastoma/mouse_neural_crest/GSE129114_E8.5_whole_embryo_Wnt1_counts.txt", header = T, row.names = 1)
mouse_all_g=intersect(rownames(E8.5), rownames(E9.5_anterior))


E10.5_post=E10.5_post[mouse_all_g,]
E10.5_tail=E10.5_tail[mouse_all_g,]
E9.5_anterior=E9.5_anterior[mouse_all_g,]
E9.5_trunk=E9.5_trunk[mouse_all_g,]
E9.5_posterior=E9.5_posterior[mouse_all_g,]
E8.5=E8.5[mouse_all_g,]

mouse_trunk_mtx=as.matrix(cbind(E8.5, E9.5_anterior, E9.5_trunk, E9.5_posterior, E10.5_post, E10.5_tail))

mouse_annot=read.csv("/home/jovyan/Neuroblastoma/mouse_neural_crest/aas9536_table_S9.csv")
mouse_trunk_mtx=mouse_trunk_mtx[,which(colnames(mouse_trunk_mtx)%in%mouse_annot$X)]

process_mouse_10x =function(mtx){
  srat=CreateSeuratObject(mtx)
  #srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  #srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  #srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  #exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
  #keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  #srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  #srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}
mt_srat=process_mouse_10x(mouse_trunk_mtx)

DimPlot(mt_srat, group.by = "annot")

mt_srat@meta.data$annot=as.character(mouse_annot$label)[match(colnames(mt_srat), mouse_annot$X)]
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.character(rownames(mt_srat)) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

genesV2$MGI.symbol=make.unique(genesV2$MGI.symbol)          

mouse_mtx_2=mouse_trunk_mtx[which(rownames(mouse_trunk_mtx)%in%genesV2$MGI.symbol),]

rownames(mouse_mtx_2)=as.character(genesV2$HGNC.symbol)[match(rownames(mouse_mtx_2), genesV2$MGI.symbol)]
rownames(mouse_mtx_2)=make.unique(rownames(mouse_mtx_2))
comm_genes=intersect(intersect(rownames(atrt1_srat), rownames(atrt_inh_srat)), rownames(mouse_mtx_2))

fit_mouse = trainModel(mouse_mtx_2, as.character(mt_srat@meta.data$annot))

ps_atrt1 = predictSimilarity(fit_mouse, atrt1_srat@assays$RNA@counts[comm_genes,] ,atrt1_srat@meta.data$annot)

ps_atrt2 = predictSimilarity(fit_mouse, atrt2_srat@assays$RNA@counts[comm_genes,] ,atrt2_srat@meta.data$annot)

ps_atrt3 = predictSimilarity(fit_mouse, atrt3_srat@assays$RNA@counts[comm_genes,] ,atrt3_srat@meta.data$annot)

ps_inhouse= predictSimilarity(fit_mouse, atrt_inh_srat@assays$RNA@counts[comm_genes,] ,atrt_inh_srat@meta.data$seurat_clusters)

similarityHeatmap(ps_atrt1)
similarityHeatmap(ps_atrt2)
similarityHeatmap(ps_atrt3)
similarityHeatmap(ps_inhouse)


ellie_mouse=readRDS("/lustre/scratch119/realdata/mdt1/team274/ek12/MRT/mouse_neural_crest/mouseWnt1_ALL_trained_NoREF.RDS")
ellie_mouse$delaminatory$glmnet.fit$beta@Dimnames[[1]]

genes_10x=read.table("ATRT_data/scRNAseq/ATRT1/genes.tsv")

length(which(ellie_mouse$delaminatory$glmnet.fit$beta@Dimnames[[1]]%in%genes_10x$V1))

mtx1=atrt1_srat@assays$RNA@counts
mtx1=mtx1[which(rownames(mtx1)%in%genes_10x$V2),]
rownames(mtx1)=genes_10x$V1[match(rownames(mtx1), genes_10x$V2)]
mtx1=mtx1[which(rownames(mtx1)%in%ellie_mouse$delaminatory$glmnet.fit$beta@Dimnames[[1]]),]

pr_1=predictSimilarity(ellie_mouse, mtx1, atrt1_srat@meta.data$annot)

