
load("/home/jovyan/Neuroblastoma/Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb = data
genes_10x = read.table("/lustre/scratch119/realdata/mdt1/team274/ek12/Neuroblastoma/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
#main matrix 
mtx = Jan_nb$data$counts
srcs='/home/jovyan/Neuroblastoma/Jan_NB/extra_samples'
srcs = file.path(srcs,list.files(srcs))
srcs = srcs[grep(".txt", srcs)]
source("/home/jovyan/Neuroblastoma/gene_sets.R")
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

full_obj=read.table(srcs[1], sep = "\t", header = T)
colnames(full_obj)[2:ncol(full_obj)]=paste0(sapply(strsplit(basename(srcs[1]), "\\."), "[", 1),"-", colnames(full_obj)[2:ncol(full_obj)])
for (x in srcs[2:13]){
  tab=read.table(x, sep = "\t", header = T)
  colnames(tab)[2:ncol(tab)]=paste0(sapply(strsplit(basename(x), "\\."), "[", 1),"-", colnames(tab)[2:ncol(tab)])
  full_obj=merge(full_obj, tab, by="GENEID", all=T)
}

full_obj=full_obj[,-grep("UNK", colnames(full_obj))]
full_obj=full_obj[-grep("cluster", full_obj$GENEID),]
full_obj$GENEID=sapply(strsplit(as.character(full_obj$GENEID), "__"),"[", 1)
full_obj$GENEID

old_df=as.data.frame(as.matrix(mtx))
old_df$GENEID=rownames(old_df)

f=merge(old_df, full_obj, by="GENEID", all=T)

rownames(f)=f$GENEID
f=f[,2:ncol(f)]

joined_mtx=as.matrix(f)
joined_mtx[is.na(joined_mtx)] <- 0
joined_mtx=Matrix(joined_mtx, sparse = T)

#remove ERCCs


joined_mtx = joined_mtx[-grep("^ERCC-",rownames(joined_mtx)),]
rownames(joined_mtx)=sapply(strsplit(rownames(joined_mtx), "__"), "[", 1)

mtx_10x =joined_mtx[which(rownames(joined_mtx)%in%genes_10x$V1),]

rownames(mtx_10x)<-genes_10x$V2[match(rownames(mtx_10x),genes_10x$V1)]
rownames(mtx_10x) <- make.unique(rownames(mtx_10x))
dim(mtx_10x)

patient_ids =unique(Jan_nb$col.annot$PatientId)
cells_each_sample =sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})


cells_list = list("NB039_Tumour_organoids_in_20perc_HP_3D"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
                  "NB039_Tumour_organoids_in_OM_WNT_RSPO_TGFb_inh"=cells_each_sample$NB039$`Tumour organoids in OM_WNT_RSPO_TGFb inh`,
                  "NB039_Tumour_organoids_in_OM_WNT_RSPO"=cells_each_sample$NB039$`Tumour organoids in OM _WNT_RSPO`, 
                  "NB060_Fresh_biopsy"=as.character(cells_each_sample$NB060), 
                  "NB067_Tumour_organoids_in_OM"=as.character(cells_each_sample$NB067),
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
                  "NB152_Tumour_organoids_P1"=cells_each_sample$NB152$`Tumour organoids P1`, 
                  "000CGH_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK001", "KK002", "TM230", "TM231"))],
                  "000FQM_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK051", "KK052", "KK053", "KK054"))],
                  "000GGU_Fresh_biopsy"=colnames(mtx_10x)[which(sapply(strsplit(colnames(mtx_10x), "-"), "[", 1)%in%c("KK055", "KK056", "KK057", "KK058", "KK059"))])

meta_data=as.data.frame(colnames(mtx_10x))
colnames(meta_data)="CellId"
rownames(meta_data)=meta_data$CellId

meta_data$unique_sample <- "sample"
for (x in 1:length(cells_list)) {
  meta_data$unique_sample[which(meta_data$CellId%in%cells_list[[x]])] <-names(cells_list[x])
}

process_10x =function(mtx, meta_data){
  srat=CreateSeuratObject(mtx, meta.data = meta_data)
  srat[["percent.mt"]] = PercentageFeatureSet(srat, pattern = "^MT-")
  srat=subset(srat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)
  srat[["percent.hspGenes"]] = PercentageFeatureSet(srat, features=hspGenes[which(hspGenes%in%rownames(srat))])
  srat[["percent.riboGenes"]]=PercentageFeatureSet(srat, features=riboGenes[which(riboGenes%in%rownames(srat))])
  #exclude_genes=as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
  #keep_features=rownames(srat)[which(!rownames(srat)%in%exclude_genes)]
  #srat=subset(srat, features=keep_features)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features = rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


fresh_biopsies =unique(meta_data$unique_sample)[grep("Fresh",unique(meta_data$unique_sample))]
tumour_only_mtx =mtx_10x[,rownames(meta_data)[which(meta_data$unique_sample%in%fresh_biopsies)]]
tumour_only_meta =meta_data[which(meta_data$unique_sample%in%fresh_biopsies),]

srat_tumour =process_10x(tumour_only_mtx,tumour_only_meta)

DimPlot(srat_tumour, group.by = "unique_sample")

FeaturePlot(srat_tumour, c("PHOX2A", "PHOX2B", "SOX2", "SOX10"))


pdf("/home/jovyan/Neuroblastoma/dutch_standard.pdf",width=24,height=24)
#Clustering
plot(DimPlot(srat_tumour,label=TRUE, pt.size = 2, label.size = 6))
plot(DimPlot(srat_tumour,group.by='unique_sample', pt.size = 2))
plot(FeaturePlot(srat_tumour,c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.hspGenes', 'percent.riboGenes')))
plot(DimPlot(srat_tumour,group.by='Phase'))
plot(FeaturePlot(srat_tumour,c('PTPRC','PDGFRB','EPCAM','PECAM1','HBA1','HBB')))
plot(FeaturePlot(srat_tumour,c('PHOX2A','ISL1','CHRNA3','SOX2','SOX10','ERBB3','DBH','CHGB','TH')))
dev.off()


saveRDS(srat_tumour, "/lustre/scratch117/casm/team274/gk14/neuroblastoma_objects/srat_dutch.rds")
