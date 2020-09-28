source("Logistic_regression.R")
source("gene_sets.R")
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb <- data
genes_10x <- read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
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
meta_data$mt_fract <- colSums(mtx_10x[grep("^MT-",rownames(mtx_10x)),])/colSums(mtx_10x)
meta_data$umi_count <- colSums(mtx_10x)
meta_data$nFeatures <- colSums(mtx_10x>0)

process_dutch =function(srat){
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 500)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  srat=subset(srat, features=keep_features)
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
  #ps = predictSimilarity(fit,seurat_obj@assays$RNA@counts[rownames(seurat_obj@assays$RNA@counts)%in%rownames(dat),],
  # classes=as.character(srats[[names(cells_list[x])]]@active.ident))
  pdf(paste0("all_tum/",names(cells_list[x]),'.pdf'),width=14,height=14)
  plot(DimPlot(srats[[names(cells_list[x])]],label=TRUE, label.size = 5, pt.size = 2)+guides(colour=FALSE))
  genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
  ComplexHeatmap::Heatmap(srats[[names(cells_list[x])]]@assays$RNA@scale.data[genes_plot,], 
                          )
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB')))
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')))
  #print(similarityHeatmap(ps,
  #row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
  #column_order = c("10","17","11","7","19","0","8","1","2","4","6","12","13","5","16","14","18","3","9","15")))
  dev.off()
}
genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
Heatmap(adr@assays$RNA@scale.data[genes_plot,], cluster_columns = T, show_column_names = F, show_row_names = T)
