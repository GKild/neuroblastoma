library(Matrix)
library(MASS)
library(dplyr)
library(Seurat)
#Processing Jan's data
#Loading data
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")
Jan_nb <- data
rm(data)
####split the data for separate processing, add metadata from relevant cells, save as rds objects.####
 
#get cols for each patient-sample. 
patient_ids <-unique(Jan_nb$col.annot$PatientId)
cells_each_sample <-sapply(patient_ids,function(x){
  sapply(unique(Jan_nb$col.annot[which(Jan_nb$col.annot$PatientId==x),]$SampleType), function(y){
    dplyr::filter(Jan_nb$col.annot, PatientId==x, SampleType==y)$CellId
  })
})

cells_list <- list("NB039_Tumour_organoids_in_20perc_HP_(3D)"=cells_each_sample$NB039$`Tumour organoids in 20perc HP (3D)`,
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


#####make a seurat object with metadata for each sample######
#get relevant cols for the mtx and meta, calc spike ins and add to meta, calc % reads in 10X genes and add to meta
sapply(1:length(cells_list), function(x){
  rel_mtx <- Jan_nb$data$counts[,which(colnames(Jan_nb$data$counts)%in%cells_list[[x]])]
  rel_meta <- Jan_nb$col.annot[which(Jan_nb$col.annot$CellId%in%cells_list[[x]]),]
  
})
x <- Jan_nb$data$counts[,which(colnames(Jan_nb$data$counts)%in%cells_list[[1]])]
obj <- CreateSeuratObject(x)

obj@assays 

#calculate fraction of reads from spike-ins


data$col.annot$spikeIn_fraction<- colSums(data$data$counts[grep("^ERCC-",rownames(data$data$counts)),])/colSums(data$data$counts)
data$col.annot$spikeIn_fraction

#remove spike-ins from the data

Jan_nb_no_spikeins <- data
Jan_nb_no_spikeins$data$counts <- Jan_nb_no_spikeins$data$counts[-grep("^ERCC-",rownames(Jan_nb_no_spikeins$data$counts)),]

#convert matrix rownames from ensembl_id to HGNC gene symbol (NOT ID)
rownames(Jan_nb_no_spikeins$data$counts)<-Jan_nb_no_spikeins$row.annot$gene[match(rownames(Jan_nb_no_spikeins$data$counts),Jan_nb_no_spikeins$row.annot$ensemblId)]
#make rownames unique
rownames(Jan_nb_no_spikeins$data$counts)<- make.unique(rownames(Jan_nb_no_spikeins$data$counts))
