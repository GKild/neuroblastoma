source("Logistic_regression.R")
source("gene_sets.R")
adr <- readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
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

exclude_genes=as.character(read.table("useful_files/excludeGenes.tsv", header = F, sep = "\t")$V1)
keep_features=rownames(adr)[which(!rownames(adr)%in%exclude_genes)]
adr_subset=subset(adr,features=keep_features)
dat = adr_subset@assays$RNA@counts
dat = dat[rownames(dat) %in% rownames(mtx_10x),]
fit = trainModel(dat,as.character(adr@active.ident),maxCells=5000)


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

process_inhouse=function(mtx){
  srat=CreateSeuratObject(mtx)
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
  ps = predictSimilarity(fit,seurat_obj@assays$RNA@counts[rownames(seurat_obj@assays$RNA@counts)%in%rownames(dat),],
  classes=as.character(srats[[names(cells_list[x])]]@active.ident))
  genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
  col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
  pdf(paste0("all_tum/",names(cells_list[x]),'.pdf'),width=14,height=14)
  plot(DimPlot(srats[[names(cells_list[x])]],label=TRUE, label.size = 10, pt.size = 4)+guides(colour=FALSE))
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB')))
  plot(FeaturePlot(srats[[names(cells_list[x])]],c('MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')))
  print(Heatmap(srats[[names(cells_list[x])]]@assays$RNA@scale.data[genes_plot,],col = col_fun, cluster_columns = F, show_column_names = F,
           column_split=srats[[names(cells_list[x])]]@active.ident))
  print(similarityHeatmap(ps,
  row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
  column_order = c("10","17","11","7","19","0","8","1","2","4","6","12","13","5","16","14","18","3","9","15")))
  dev.off()
}



all_paths = c("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843576/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843577/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7843578/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733084/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733085/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7733086/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004894/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004902/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY8004910/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787237/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787238/outs/filtered_gene_bc_matrices/GRCh38/",
              "/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7787239/outs/filtered_gene_bc_matrices/GRCh38/")
names(all_paths) = c("STDY7843576", "STDY7843577", "STDY7843578","STDY7733084","STDY7733085", "STDY7733086", "STDY8004894",
                     "STDY8004902", "STDY8004910", "STDY7787237", "STDY7787238", "STDY7787239")

nb1 = Read10X(data.dir = all_paths[1:3])
nb2 = Read10X(data.dir = all_paths[4:6])
nb3 = Read10X(data.dir = all_paths[7:9])
gn = Read10X(data.di=all_paths[10:12])

nb_list=list(nb1, nb2, nb3, gn)
nb_names=c("nb1", "nb2", "nb3", "gn")
names(nb_list)=nb_names
inhouse_srats=list()
for(x in c(1:length(nb_names))){
  inhouse_srats[[nb_names[x]]] = process_inhouse(nb_list[[nb_names[x]]])
  ps = predictSimilarity(fit,inhouse_srats[[nb_names[x]]]@assays$RNA@counts[rownames(inhouse_srats[[nb_names[x]]]@assays$RNA@counts)%in%rownames(dat),],
                         classes=as.character(inhouse_srats[[nb_names[x]]]@active.ident))
  genes_plot=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB', 'MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')
  col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdYlBu")))
  pdf(paste0("all_tum/",nb_names[x],'.pdf'),width=14,height=14)
  plot(DimPlot(inhouse_srats[[nb_names[x]]],label=TRUE, label.size = 10, pt.size = 4)+guides(colour=FALSE))
  plot(FeaturePlot(inhouse_srats[[nb_names[x]]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB', 'PTPRB')))
  plot(FeaturePlot(inhouse_srats[[nb_names[x]]],c('MYCN','CHGB','SOX2','SOX10','PHOX2A','PHOX2B')))
  print(Heatmap(inhouse_srats[[nb_names[x]]]@assays$RNA@scale.data[genes_plot,],col = col_fun, cluster_columns = F, show_column_names = F,
                column_split=inhouse_srats[[nb_names[x]]]@active.ident))
  print(similarityHeatmap(ps,
                          row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
                          column_order = c("10","17","11","7","19","0","8","1","2","4","6","12","13","5","16","14","18","3","9","15")))
  dev.off()
}

inhouse_has_tum=inhouse_srats[2:4]
dutch_has_tum=c(srats[c(4,5,11,14,15,19,20)])
org_has_tum=c(srats[c(1,2,3,10,21)])

inhouse_clust_numbers=list(nb2=2,nb3=c(0,1,2,3,6,7), gn=c(0,2,6)) 
dutch_clust_numbers=list(NB060_Fresh_biopsy=1, NB067_Fresh_biopsy=0, NB106_Fresh_biopsy=8, 
                   NB124_Fresh_biopsy=11, NB125_Fresh_biopsy=c(0,3,10), NB151_Fresh_biopsy=4, NB152_Fresh_biopsy=c(7,8))
org_clust_numbers=list(NB039_Tumour_organoids_in_20perc_HP_3D=c(0,1,2,3,4,5),
                       NB039_Tumour_organoids_in_OM_WNT_RSPO_TGFb_inh=c(0,1,2,3),
                       NB039_Tumour_organoids_in_OM_WNT_RSPO=c(0,1,2,3,4),
                       NB098_Tumour_organoids=2, NB152_Tumour_organoids_P1=8)

adr_scp=subset(adr, idents=c(10,17,11,7,19))
adr_scp@meta.data$seurat_clusters=as.character(adr_scp@meta.data$seurat_clusters)
adr_scp@meta.data$seurat_clusters[which(adr_scp@meta.data$seurat_clusters==10)]="fa_10"
adr_scp@meta.data$seurat_clusters[which(adr_scp@meta.data$seurat_clusters==17)]="fa_17"
adr_scp@meta.data$seurat_clusters[which(adr_scp@meta.data$seurat_clusters==11)]="fa_11"

adr_scp@meta.data$seurat_clusters[which(adr_scp@meta.data$seurat_clusters==7)]="fa_7"

adr_scp@meta.data$seurat_clusters[which(adr_scp@meta.data$seurat_clusters==19)]="fa_19"

inhouse_top_de_genes=list()
all_inhouse_genes=vector()
for(x in c(1:length(inhouse_has_tum))){
  s1=subset(inhouse_has_tum[[x]], idents = inhouse_clust_numbers[[x]])
  s1@meta.data$seurat_clusters="tum_clust"
  m=merge(adr_scp, s1)
  all_markers=quickMarkers(m@assays$RNA@counts, m@meta.data$seurat_clusters, FDR = 0.01)
  all_markers_tum=all_markers[-grep("fa",all_markers$cluster),]
  rel_x=dplyr::filter(all_markers_tum, geneFrequency>0.5, geneFrequencyOutsideCluster<0.05)
  inhouse_top_de_genes[[names(has_tum)[x]]]=rel_x
  all_inhouse_genes=c(all_inhouse_genes, as.character(rel_x$gene))
}


inhouse_genes_df=data.frame(genes=all_inhouse_genes)
inhouse_freq=as.data.frame(table(inhouse_genes_df$genes))
inhouse_de=inhouse_freq[order(inhouse_freq$Freq, decreasing = T),]

dutch_top_de_genes=list()
all_dutch_genes=vector()
for(x in c(1:length(dutch_has_tum))){
  s1=subset(dutch_has_tum[[x]], idents = dutch_clust_numbers[[x]])
  s1@meta.data$seurat_clusters="tum_clust"
  m=merge(adr_scp, s1)
  all_markers=quickMarkers(m@assays$RNA@counts, m@meta.data$seurat_clusters, FDR = 0.01)
  all_markers_tum=all_markers[-grep("fa",all_markers$cluster),]
  rel_x=dplyr::filter(all_markers_tum, geneFrequency>0.5, geneFrequencyOutsideCluster<0.05)
  dutch_top_de_genes[[names(has_tum)[x]]]=rel_x
  all_dutch_genes=c(all_dutch_genes, as.character(rel_x$gene))
}


dutch_genes_df=data.frame(genes=all_dutch_genes)
dutch_freq=as.data.frame(table(dutch_genes_df$genes))
dutch_de=dutch_freq[order(dutch_freq$Freq, decreasing = T),]


org_top_de_genes=list()
all_org_genes=vector()
for(x in c(1:length(org_has_tum))){
  s1=subset(org_has_tum[[x]], idents = org_clust_numbers[[x]])
  s1@meta.data$seurat_clusters="tum_clust"
  m=merge(adr_scp, s1)
  all_markers=quickMarkers(m@assays$RNA@counts, m@meta.data$seurat_clusters, FDR = 0.01)
  all_markers_tum=all_markers[-grep("fa",all_markers$cluster),]
  rel_x=dplyr::filter(all_markers_tum, geneFrequency>0.5, geneFrequencyOutsideCluster<0.05)
  org_top_de_genes[[names(has_tum)[x]]]=rel_x
  all_org_genes=c(all_org_genes, as.character(rel_x$gene))
}


org_genes_df=data.frame(genes=all_org_genes)
org_freq=as.data.frame(table(org_genes_df$genes))
org_de=org_freq[order(org_freq$Freq, decreasing = T),]



