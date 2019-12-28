library(Matrix)
library(MASS)
library(dplyr)
library(Seurat)
library(ggplot2)
library(glmnet)
#function to do standard 10X processing in Seurat
process_10x <- function(seurat_obj, npcs=50, res=1.0){
  seurat_obj = NormalizeData(seurat_obj)
  seurat_obj = ScaleData(seurat_obj)
  seurat_obj = FindVariableFeatures(seurat_obj)
  seurat_obj = RunPCA(seurat_obj,approx=FALSE)
  seurat_obj = RunUMAP(seurat_obj,dims=seq(npcs))
  seurat_obj = FindNeighbors(seurat_obj,dims=seq(npcs))
  seurat_obj = FindClusters(seurat_obj,res=res)
  return(seurat_obj)
}

#define params
outdir = 'Jan_NB/split_objects/'
npcs=50
mtMax=0.2
nGeneMin=300
nCntsMin=500
minCells=100
#normal reference for LR
adr <- readRDS("babyAdrenalRef/babyAdrenal2.RDS")
#Processing Jan's data
#Loading data
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb <- data
rm(data)
genes_10x <- read.table("Neuroblastoma1/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
#main matrix 
mtx <- Jan_nb$data$counts
#column metadata 
meta_data <- Jan_nb$col.annot
#calculate fraction of spikeins
meta_data$spikeIn_fraction<- colSums(mtx[grep("^ERCC-",rownames(mtx)),])/colSums(mtx)
#remove spikeins
mtx <- mtx[-grep("^ERCC-",rownames(mtx)),]
#calculate fraction of reads lost by using 10X ensembl IDs, and fraction of reads lost from  MT genes that are NOT in 10X. 
mtx_10x <-mtx[which(rownames(mtx)%in%genes_10x$V1),]
not_10x<-mtx[which(!rownames(mtx)%in%genes_10x$V1),]
meta_data$reads_lost <-(colSums(mtx)-colSums(mtx_10x))/colSums(mtx)
rownames(not_10x)<-Jan_nb$row.annot$gene[match(rownames(not_10x),Jan_nb$row.annot$ensemblId)]
meta_data$mt_reads_lost <-colSums(not_10x[grep("^MT-",rownames(not_10x)),])/colSums(mtx)

#get cols for each patient-sample. 
patient_ids <-unique(Jan_nb$col.annot$PatientId)
cells_each_sample <-sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})

cells_list <- list("NB039_Tumour_organoids_in_20perc_HP_3D"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO_TGFb_inh"=cells_each_sample$NB039$`Tumour organoids in OM_WNT_RSPO_TGFb inh`,
                   "NB039_Tumour_organoids_in_OM_WNT_RSPO"=cells_each_sample$NB039$`Tumour organoids in OM _WNT_RSPO`, 
                   "NB060_Fresh_biopsy"=as.character(cells_each_sample$NB060), 
                   "NB067_Fresh_biopsy"=as.character(cells_each_sample$NB067),
                   "NB086_Fresh_biopsy"=as.character(cells_each_sample$NB086),
                   "NB098_Fresh_biopsy"=cells_each_sample$NB098$`Fresh biopsy`,
                   "NB098_peripheral_blood_from_fresh_biopsy"=cells_each_sample$NB098$`peripheral blood from fresh biopsy`,
                   "NB098_TILs_from_fresh_biopsy"=cells_each_sample$NB098$`TILs from fresh biopsy`,
                   "NB098_Tumour_organoids"=cells_each_sample$NB098$`Tumour organoids`, 
                   "NB106_Fresh_biopsy"=as.character(cells_each_sample$NB106),
                   "NB107_Fresh_biopsy"=as.character(cells_each_sample$NB107),
                   "NB123_Fresh_biopsy"=as.character(cells_each_sample$NB123),
                   "NB124_Fresh_biopsy"=as.character(cells_each_sample$NB124),
                   "NB125_Fresh_biopsy"=as.character(cells_each_sample$NB125),
                   "NB130_Fresh_biopsy"=as.character(cells_each_sample$NB130),
                   "NB132_Fresh_biopsy"=as.character(cells_each_sample$NB132),
                   "NB138_Fresh_biopsy"=as.character(cells_each_sample$NB138),
                   "NB151_Fresh_biopsy"=as.character(cells_each_sample$NB151),
                   "NB152_Fresh_biopsy"=cells_each_sample$NB152$`Fresh biopsy`,
                   "NB152_Tumour_organoids_P1"=cells_each_sample$NB152$`Tumour organoids P1`)


#add the unique sample to metadata 

meta_data$unique_sample <- "sample"
for (x in 1:length(cells_list)) {
  meta_data$unique_sample[which(meta_data$CellId%in%cells_list[[x]])] <-names(cells_list[x])
}

# going to be working with only 10X from now on, so use mtx_10x,
# convert to gene names, and force uniqueness on each rowname. 

rownames(mtx_10x)<-genes_10x$V2[match(rownames(mtx_10x),genes_10x$V1)]
rownames(mtx_10x) <- make.unique(rownames(mtx_10x))

#add other QC metrics to metadata - mt frac, UMIs, nGenes with more than 0 reads
meta_data$mt_fract <- colSums(mtx_10x[grep("^MT-",rownames(mtx_10x)),])/colSums(mtx_10x)
meta_data$umi_count <- colSums(mtx_10x)
meta_data$nFeatures <- colSums(mtx_10x>0)

#train the LR model

dat = adr@assays$RNA@counts
dat = dat[rownames(dat) %in% rownames(mtx_10x),]
fit = trainModel(dat,paste0('fAd',as.character(adr@active.ident)),maxCells=5000)

#now process Seurat and run LR 

srats=list()

for (x in 1:length(cells_list)) {
  message(names(cells_list[x]))
  srats[[names(cells_list[x])]]=list()
  each_sample = mtx_10x[,cells_list[[x]]]
  each_meta = meta_data[cells_list[[x]],]
  #avoid NAs in MTfract column - if NA, then total reads is 0 and will fail QC anyways
  each_meta$mt_fract[is.na(each_meta$mt_fract)]=0
  #which cells are passing QC
  passed_cells = dplyr::filter(each_meta, mt_fract<mtMax, umi_count>nCntsMin, nFeatures>nGeneMin)
  w =nrow(passed_cells)
  message(sprintf("After filtering, keeping %d cells of %d",w,nrow(each_meta)))
  if(w<minCells){
    message(sprintf("Sample has only %d cells passing QC.  Too shit to continue.",w))
    next
  }
  #create seurat object and add it into a list of Seurat objects, run LR, make a umap, some featureplots, and LR heatmap
  seurat_obj =CreateSeuratObject(each_sample[,passed_cells$CellId],meta.data=each_meta[passed_cells$CellId,])
  srats[[names(cells_list[x])]] = process_10x(seurat_obj,npcs=npcs)
  ps = predictSimilarity(fit,seurat_obj@assays$RNA@counts[rownames(seurat_obj@assays$RNA@counts)%in%rownames(dat),],
                         classes=as.character(srats[[names(cells_list[x])]]@active.ident))
  pdf(file.path(outdir,paste0(names(cells_list[x]),'.pdf')),width=14,height=14)
  plot(DimPlot(srats[[names(cells_list[x])]],label=TRUE)+guides(colour=FALSE))
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB')))
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('MYCN','SOX10','SOX2','PHOX2A','CHGB','PHOX2B')))
  print(similarityHeatmap(ps,
                          row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
                          column_order = colnames(ps)[order(as.numeric(gsub('[A-Za-z]','',colnames(ps))))]
  ))
  dev.off()
}
#make metadata plot
rl_melted <- melt(meta_data, id.vars = c("unique_sample", "CellId"), measure.vars = c("reads_lost", "mt_reads_lost"))
pdf('qc_metrics.pdf',width=21,height=14)
gg = ggplot(rl_melted,aes(x=value, fill=variable)) +
  geom_density(alpha=0.5) + 
  ggtitle('Reads lost') +
  facet_wrap(~unique_sample)
plot(gg)
gg=ggplot(meta_data, aes(x=mt_fract)) +
  geom_density()+
  ggtitle('Fraction of mitochondrial reads') +
  facet_wrap(~unique_sample)
plot(gg)
gg=ggplot(meta_data, aes(x=log10(umi_count))) +
  geom_density()+
  ggtitle('Log10 UMI count') +
  facet_wrap(~unique_sample)
plot(gg)
gg=ggplot(meta_data, aes(x=log10(nFeatures))) +
  geom_density()+
  ggtitle('Log10 nFeatures') +
  facet_wrap(~unique_sample)
plot(gg)
dev.off()



