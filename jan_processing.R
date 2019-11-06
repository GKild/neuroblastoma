#Processing Jan's data
#Loading data
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")

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

####split the data for separate processing, add metadata from relevant cells, save as rds objects.####
 
#get cols for each patient-sample. 

sapply(unique(Jan_nb_no_spikeins$))

