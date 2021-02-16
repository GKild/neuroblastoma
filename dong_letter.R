library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(ggplot2)
source("gene_sets.R")

s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

adr_all=readRDS("/lustre/scratch117/casm/team274/gk14/adr_all_annot.rds")
dong1=read.table("Dong_adrenals/GSM4088785_F2_gene_cell_exprs_table.xls", sep = "\t", header = T)
dong2=read.table("Dong_adrenals/GSM4088786_F7_gene_cell_exprs_table.xls", sep = "\t", header = T)
dong3=read.table("Dong_adrenals/GSM4088787_F106_gene_cell_exprs_table.xls", sep = "\t", header = T)
dong4=read.table("Dong_adrenals/GSM4088788_F107_gene_cell_exprs_table.xls", sep = "\t", header = T)


dong_annot=read.csv("Dong_adrenals/Adrenal_gland_annotation.csv")
dong_annot$cell_id_und=sapply(strsplit(as.character(dong_annot$cell_id), "_"), "[", 1)

dong_annot1=dong_annot[which(dong_annot$sample_id=="F2"),]
dong_annot2=dong_annot[which(dong_annot$sample_id=="F7"),]
dong_annot3=dong_annot[which(dong_annot$sample_id=="F106"),]
dong_annot4=dong_annot[which(dong_annot$sample_id=="F107"),]

rownames(dong_annot1)=dong_annot1$cell_id_und
rownames(dong_annot2)=dong_annot2$cell_id_und
rownames(dong_annot3)=dong_annot3$cell_id_und
rownames(dong_annot4)=dong_annot4$cell_id_und

rownames(dong1)=make.unique(as.character(dong1$Symbol))
dong1=dong1[,rownames(dong_annot1)]
rownames(dong2)=make.unique(as.character(dong2$Symbol))
dong2=dong2[,rownames(dong_annot2)]
rownames(dong3)=make.unique(as.character(dong3$Symbol))
dong3=dong3[,rownames(dong_annot3)]
rownames(dong4)=make.unique(as.character(dong4$Symbol))
dong4=dong4[,rownames(dong_annot4)]

dong1=as.matrix(dong1)
dong2=as.matrix(dong2)
dong3=as.matrix(dong3)
dong4=as.matrix(dong4)


dong1_srat=CreateSeuratObject(dong1, meta.data = dong_annot1)
dong2_srat=CreateSeuratObject(dong2, meta.data = dong_annot2)
dong3_srat=CreateSeuratObject(dong3, meta.data = dong_annot3)
dong4_srat=CreateSeuratObject(dong4, meta.data = dong_annot4)


dong_all=merge(dong1_srat, y=list(dong2_srat, dong3_srat, dong4_srat))
process_dong =function(srat){
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
  srat = NormalizeData(srat)
  srat = FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


dong_all=process_dong(dong_all)



dong_medulla=subset(dong_all,
                    cells = rownames(dong_all@meta.data)[which(dong_all@meta.data$cell_type%in%c("SCPs",
                                                                                                 "Chromaffin cells", 
                                                                                                 "Sympathoblasts"))])

redo_adr =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 20)
  srat = FindNeighbors(srat, dims=1:20)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:20, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

dong_medulla_reclust=redo_adr(dong_medulla)

#dong signature genes used to define fetal adrenal cell types
pdf("dong_reclust.pdf", height = 10, width = 10)
plot(DimPlot(dong_medulla_reclust, label = T, pt.size = 3, label.size = 3) +theme(legend.position = "none"))
plot(DimPlot(dong_medulla_reclust, group.by = "cell_type", label = T, pt.size = 3, label.size = 3) +theme(legend.position = "none"))
#Dong SCP
plot(FeaturePlot(dong_medulla_reclust, c("SOX10", "PLP1", "MPZ", "FABP7", "FOXD3", "ERBB3")))
#Dong sympathoblast
plot(FeaturePlot(dong_medulla_reclust, c("CARTPT", "INSM1", "DBH", "CHGA", "TH")))
#Dong chromaffin
plot(FeaturePlot(dong_medulla_reclust, c("TH", "CHGA","CHGB", "DBH","HMX1", "NTRK1", "PHOX2B", "ISL1", "PRPH", "NPY")))
#De Preter sympathoblast
plot(FeaturePlot(dong_medulla_reclust, c("BCL2", "GAP43", "NPY")))
#De preter Chromaffin
plot(FeaturePlot(dong_medulla_reclust, c("CHGA", "CHGB", "DBH", "DDC", "TH", "PNMT")))
plot(FeaturePlot(dong_medulla_reclust, c("STAR", "NOV")))
dev.off()


#############
# Libraries #
#############
library(Seurat)
library(biomaRt)
library(edgeR)
library(SoupX)
import('scQC.R',as='scQC')
import('/home/jovyan/Dediff/logisticRegression.R',as='lr')

#############
# Load data #
#############
#Load adrenal reference
#Load mouse data
m12.5 = read.table('GSE99933_E12.5_counts.txt.gz',sep='\t')
m13.5 = read.table('GSE99933_E13.5_counts.txt.gz',sep='\t')
#Get shared genes
#Get orthologs
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl")
gnMap = getLDS(attributes=c("ensembl_gene_id",'mgi_symbol'),
               filters="mgi_symbol", 
               values=rownames(m12.5), 
               mart=mart2,
               attributesL=c("ensembl_gene_id",'hgnc_symbol'), 
               martL=mart1)
gnMap$ambiguous = table(gnMap$MGI.symbol)[gnMap$MGI.symbol]>1 | table(gnMap$HGNC.symbol)[gnMap$HGNC.symbol]>1
gnMap = gnMap[!gnMap$ambiguous,]

#####################
# Process and label #
#####################
srat12.5 = quickCluster(m12.5, clusteringRes = 1)
srat12.5@meta.data$annot = c('4'='Sympathoblasts', #Based on Cartpt
                             '6' = 'SCP', #Based on c('Sox10','Erbb3','Foxd3')
                             '3' = 'SCP', 
                             '8' = 'B2', #Bridge from SCP towards sympathoblasts made up of cycling cells.  Yellow in Furlan et al.
                             '2' = 'Bridge', #Based on c('Thbd','Htr3a','Dll3')
                             '0' = 'Bridge', 
                             '7' = 'Chromaffin', #Based on c('Th','Chgb','Foxq1')
                             '5' = 'Chromaffin',
                             '1' = 'Chromaffin'
)[as.character(srat12.5@meta.data$seurat_clusters)]
srat13.5 = quickCluster(m13.5, clusteringRes = 1)
srat13.5@meta.data$annot = c('1'='Sympathoblasts', #Based on Cartpt
                             '5' = 'Sympathoblasts',
                             '0' = 'SCP', #Based on c('Sox10','Erbb3','Foxd3')
                             '3' = 'Bridge', #Based on c('Thbd','Htr3a','Dll3')
                             '4' = 'Chromaffin', #Based on c('Th','Chgb','Foxq1')
                             '2' = 'Chromaffin'
)[as.character(srat13.5@meta.data$seurat_clusters)]

###############
# Train model #
###############
######################
# Logistic regression
tgtDat = srat12.5@assays$RNA@counts
#Drop genes that we have no mapping for
tgtDat = tgtDat[rownames(tgtDat) %in% gnMap$MGI.symbol,]
#Or that are not in dong_medulla_reclust
tgtDat = tgtDat[gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)] %in% rownames(dong_medulla_reclust),]
#Rename to human symbol
rownames(tgtDat) = gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)]
tgtDat12.5 = tgtDat
#13.5
tgtDat = srat13.5@assays$RNA@counts
#Drop genes that we have no mapping for
tgtDat = tgtDat[rownames(tgtDat) %in% gnMap$MGI.symbol,]
#Or that are not in dong_medulla_reclust
tgtDat = tgtDat[gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)] %in% rownames(dong_medulla_reclust),]
#Rename to human symbol
rownames(tgtDat) = gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)]
tgtDat13.5 = tgtDat
tgtDatAd = dong_medulla_reclust@assays$RNA@counts
#What genes are common
gns = rownames(tgtDat12.5)
gns = gns[gns%in%rownames(tgtDat13.5)]
gns = gns[gns%in%rownames(dong_medulla_reclust)]
#Make basic seurat objects from them
intDat12.5 = CreateSeuratObject(tgtDat12.5[gns,])
intDat12.5 = NormalizeData(intDat12.5)
intDat12.5 = FindVariableFeatures(intDat12.5)
intDat12.5@meta.data$annot = srat12.5@meta.data$annot
intDat13.5 = CreateSeuratObject(tgtDat13.5[gns,])
intDat13.5 = NormalizeData(intDat13.5)
intDat13.5 = FindVariableFeatures(intDat13.5)
intDat13.5@meta.data$annot = srat13.5@meta.data$annot
intDatAd = dong_medulla_reclust@assays$RNA@counts[gns,]
intDatAd = CreateSeuratObject(intDatAd)
intDatAd = NormalizeData(intDatAd)
intDatAd = FindVariableFeatures(intDatAd)
#Just medulla
intDatMed = dong_medulla_reclust@assays$RNA@counts[gns,dong_medulla_reclust@meta.data$cell_type %in% unique(dong_medulla_reclust@meta.data$cell_type)]
intDatMed = intDatMed[rownames(intDatMed) %in% rownames(tgtDat12.5),]
intDatMed = CreateSeuratObject(intDatMed)
intDatMed@meta.data$new_clust = dong_medulla_reclust@meta.data$cell_type[dong_medulla_reclust@meta.data$cell_type %in% unique(dong_medulla_reclust@meta.data$cell_type)]
intDatMed = NormalizeData(intDatMed)
intDatMed = FindVariableFeatures(intDatMed)
#Train model at last
fit12.5 = lr$trainModel(intDat12.5@assays$RNA@counts,srat12.5@meta.data$annot, workers=NULL)
fit13.5 = lr$trainModel(intDat13.5@assays$RNA@counts,srat13.5@meta.data$annot, workers=NULL)
labAnch12.5 = FindTransferAnchors(reference = intDat12.5,query = intDatMed,dims=seq(30))
labAnch13.5 = FindTransferAnchors(reference = intDat13.5,query = intDatMed,dims=seq(30))








#######
# Map #
#######
out = lr$predictSimilarity(fit = fit12.5,tgtData = dong_medulla_reclust@assays$RNA@counts,minGeneMatch = 0,logits = F)
w = abs(out)>4
out[w] = 4*sign(out)[w]
w = dong_medulla_reclust@meta.data$cell_type %in% unique(dong_medulla_reclust@meta.data$cell_type)
tmp = dong_medulla_reclust[,w]
tmp@meta.data = cbind(tmp@meta.data,out[w,])
pdf('dong_medulla_reclustFurlanLR12.5.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('Sympathoblasts','SCP','Bridge','Chromaffin'),colScheme='default') + scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5)
dev.off()
#13.5
out = lr$predictSimilarity(fit = fit13.5,tgtData = dong_medulla_reclust@assays$RNA@counts,minGeneMatch = 0, logits = F)
w = abs(out)>4
out[w] = 4*sign(out)[w]
w = dong_medulla_reclust@meta.data$cell_type %in% unique(dong_medulla_reclust@meta.data$cell_type)
tmp = dong_medulla_reclust[,w]
tmp@meta.data = cbind(tmp@meta.data,out[w,])
pdf('dong_medulla_reclustFurlanLR13.5.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('Sympathoblasts','SCP','Bridge','Chromaffin'),colScheme='default') + scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5)
dev.off()
#Seurat label transfer
cc = c(Sympathoblasts='#E30613',
       Chromaffin = '#1D71B8',
       Bridge='black',
       SCP='#3AAA35')
#12.5
pred = TransferData(labAnch12.5,refdata=intDat12.5$annot,dims=seq(30))
tmp = dong_medulla_reclust
tmp = AddMetaData(tmp,metadata=pred)
dd = cbind(tmp@meta.data,tmp@reductions$umap@cell.embeddings)
dd = dd[dd$predicted.id %in% names(cc),]
pdf('dong_medulla_reclustFurlanSeurat12.5.pdf',width=4,height=4,useDingbats=FALSE)
plot(dd$UMAP_1,dd$UMAP_2,
     xlab='UMAP1',
     ylab='UMAP2',
     main='Seurat 12.5',
     col = cc[dd$predicted.id],
     pch=19,
     cex=0.1)
dev.off()
#13.5
pred = TransferData(labAnch13.5,refdata=intDat13.5$annot,dims=seq(30))
tmp = dong_medulla_reclust
tmp = AddMetaData(tmp,metadata=pred)
dd = cbind(tmp@meta.data,tmp@reductions$umap@cell.embeddings)
dd = dd[dd$predicted.id %in% names(cc),]
pdf('dong_medulla_reclustFurlanSeurat13.5.pdf',width=4,height=4,useDingbats=FALSE)
plot(dd$UMAP_1,dd$UMAP_2,
     xlab='UMAP1',
     ylab='UMAP2',
     main='Seurat 13.5',
     col = cc[dd$predicted.id],
     pch=19,
     cex=0.1)
dev.off()


#12.5 similarity scores
pred = TransferData(labAnch12.5,refdata=intDat12.5$annot,dims=seq(30))
tmp = dong_medulla_reclust
tmp = AddMetaData(tmp,metadata=pred)
pdf('dong_medulla_reclustFurlan_sim_Seurat12.5.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('prediction.score.Sympathoblasts','prediction.score.SCP','prediction.score.Bridge',
                 'prediction.score.Chromaffin'),colScheme='default') + 
  scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5)
dev.off()

#13.5 similarity scores
pred = TransferData(labAnch13.5,refdata=intDat13.5$annot,dims=seq(30))
tmp = dong_medulla_reclust
tmp = AddMetaData(tmp,metadata=pred)
dd = cbind(tmp@meta.data,tmp@reductions$umap@cell.embeddings)
dd = dd[dd$predicted.id %in% names(cc),]
pdf('dong_medulla_reclustFurlan_sim_Seurat13.5.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('prediction.score.Sympathoblasts','prediction.score.SCP','prediction.score.Bridge',
                 'prediction.score.Chromaffin'),colScheme='default') +
  scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5)
dev.off()


mouse_ref_comb=merge(srat12.5, srat13.5)
tgtDat = mouse_ref_comb@assays$RNA@counts
#Drop genes that we have no mapping for
tgtDat = tgtDat[rownames(tgtDat) %in% gnMap$MGI.symbol,]
#Or that are not in dong_medulla_reclust
tgtDat = tgtDat[gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)] %in% rownames(dong_medulla_reclust),]
#Rename to human symbol
rownames(tgtDat) = gnMap$HGNC.symbol[match(rownames(tgtDat),gnMap$MGI.symbol)]
tgtDat_mouse_ref_comb = tgtDat
intDat_mouse_ref_comb= CreateSeuratObject(tgtDat_mouse_ref_comb[gns,])
intDat_mouse_ref_comb = NormalizeData(intDat_mouse_ref_comb)
intDat_mouse_ref_comb = FindVariableFeatures(intDat_mouse_ref_comb)
intDat_mouse_ref_comb@meta.data$annot = mouse_ref_comb@meta.data$annot

labAnch_mouse_ref_comb = FindTransferAnchors(reference = intDat_mouse_ref_comb,query = intDatMed,dims=seq(30))

pred = TransferData(labAnch_mouse_ref_comb,refdata=intDat_mouse_ref_comb$annot,dims=seq(30))
tmp = dong_medulla_reclust
tmp = AddMetaData(tmp,metadata=pred)
pdf('dong_medulla_reclustFurlanSeuratMouse.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('prediction.score.Sympathoblasts','prediction.score.SCP','prediction.score.Bridge',
                 'prediction.score.Chromaffin'),colScheme='default') +
  scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5)
dev.off()

train_mouse_comb=lr$trainModel(intDat_mouse_ref_comb@assays$RNA@counts, intDat_mouse_ref_comb@meta.data$annot, workers=NULL)
out = lr$predictSimilarity(fit = tratrain_mouse_comb,tgtData = dong_medulla_reclust@assays$RNA@counts,minGeneMatch = 0,logits = F)
w = abs(out)>4
out[w] = 4*sign(out)[w]
w = dong_medulla_reclust@meta.data$cell_type %in% unique(dong_medulla_reclust@meta.data$cell_type)
tmp = dong_medulla_reclust[,w]
tmp@meta.data = cbind(tmp@meta.data,out[w,])
pdf('dong_medulla_reclustFurlanLR_comb.pdf',width=4,height=4,useDingbats=FALSE)
scQC$fPlot(tmp,c('Sympathoblasts','SCP','Bridge','Chromaffin'),colScheme='default') + scale_colour_gradient2(low = '#edf8b1', mid = '#7fcdbb', high = '#2c7fb8', midpoint = 0.5, limits=c(0,1))
dev.off()

dong_medulla_reclust@meta.data$cell_type_clust="x"
dong_medulla_reclust@meta.data$cell_type_clust[which(dong_medulla_reclust@meta.data$cell_type=="SCPs")]=1
dong_medulla_reclust@meta.data$cell_type_clust[which(dong_medulla_reclust@meta.data$cell_type=="Sympathoblasts")]=2
dong_medulla_reclust@meta.data$cell_type_clust[which(dong_medulla_reclust@meta.data$cell_type=="Chromaffin cells")]=3

#make figure 1A
pdf("dong_just_medulla.pdf", height = 3, width=3,useDingbats = F)
gg=DimPlot(dong_medulla_reclust, group.by = "cell_type_clust", label=T,label.size = 3,cols=c("grey8", "grey35", "grey65")) +
  ggtitle("") + theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
                      axis.ticks = element_blank(), legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1")
plot(gg)
dev.off()

new_idents=as.factor(dong_medulla_reclust@meta.data$cell_type_clust)
names(new_idents)=names(dong_medulla_reclust@active.ident)
dong_medulla_reclust@active.ident=new_idents
pdf("dong_dotplot.pdf", height = 5, width = 5, useDingbats = F)
plot(DotPlot(dong_medulla_reclust,
        features = rev(c("INSM1", "CARTPT", "SLC18A3", "TH", "CHGA", "CHGB", "DBH","PHOX2B", "BCL2", "GAP43", "NPY", "DDC", "PNMT"))) +
  coord_flip() +
  scale_colour_gradient2(low = "#4575B4", mid = "#ffffbf", high = "#D73027"))
dev.off()


