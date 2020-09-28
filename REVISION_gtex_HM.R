gen <- read.delim( "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm (1).gct", skip = 2, header = TRUE )

## First two columns corresponds to "Transcript name" and "Gene name"
gen_ann <- gen[ , c(1,2) ]
gen_clean <- as.matrix( gen[ , -seq( 2 ) ] )
rm(gen)
# Load map file
sample_map <- read.delim( "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" )
sample_map$SAMPID_new <- gsub( "-", "\\.", sample_map$SAMPID )
rownames( sample_map ) <- sample_map$SAMPID

sample_map=sample_map[which(sample_map$SAMPID_new%in%colnames(gen_clean)),]
srat_dutch=subset(srat_dutch, idents=15, invert=T)
srat_dutch@meta.data$new_idents="x"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(13,24,10,17))]="tumour"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(6))]="Schwannian stroma"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(7,29))]="Endothelium"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(28,18,23,5,26,0,8,20,9,3,22,25,4,19,21,27))]="Leukocytes"
srat_dutch@meta.data$new_idents[which(srat_dutch@meta.data$seurat_clusters%in%c(1,2,11,12,14,16))]="Mesenchyme"

dutch_tum_marks=quickMarkers(srat_dutch@assays$RNA@counts, srat_dutch@meta.data$new_idents, N = Inf, FDR=0.01)
inh_tum_marks=quickMarkers(srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$new_idents, N=Inf, FDR = 0.01)
adr_all_marks=quickMarkers(adr_all@assays$RNA@counts, adr_all@meta.data$new_clust, N=Inf, FDR = 0.01)

dutch_marks_filt=filter(dutch_tum_marks, cluster=="tumour", geneFrequency>0.5, geneFrequencySecondBest<0.2)
inh_marks_filt=filter(inh_tum_marks, cluster=="tumour", geneFrequency>0.5, geneFrequencySecondBest<0.2)


tum_genes_high=intersect(inh_marks_filt$gene, dutch_marks_filt$gene)

med_gene_frac=apply(just_med@assays$RNA@counts, 1, function(x){sum(x>0)/6451*100})
med_genes_high=names(med_gene_frac[which(med_gene_frac>10)])
med_and_tum=intersect(tum_genes_high, med_genes_high)

genes_of_i=c(med_and_tum,"MAGEA4", "CCDC144NL-AS1" , "PRAME","CRYBB2", "SIX3","LSMEM1") 
genes_of_i=intersect(genes_of_i, rownames(target_tpm))

length(which(genes_of_i%in%gen_ann$Description))
rownames(gen_clean)=gen_ann$Description

lr_target_tab=target_tpm[genes_of_i,target_clinical_data$TARGET_SHORT[which(target_clinical_data$COG.Risk.Group=="Low Risk")]]
lr_target_means=rowMeans(lr_target_tab)
hr_target_tab=target_tpm[genes_of_i,target_clinical_data$TARGET_SHORT[which(target_clinical_data$COG.Risk.Group=="High Risk")]]
hr_target_means=rowMeans(hr_target_tab)

gtex_goi=gen_clean[genes_of_i,]

bla=sapply(as.character(unique(sample_map$SMTSD)), function(x){
  rowMeans(gtex_goi[,sample_map$SAMPID_new[which(sample_map$SMTSD==x)]])
})

bla=as.data.frame(bla)

bla=cbind(lr_target_means, hr_target_means, bla)
bla_log2=log2(bla+1)
col_fun=colorRamp2(c(0, 5, 10), c("#edf8b1", "#7fcdbb", "#2c7fb8"))
Heatmap(as.matrix(bla_log2),col=col_fun, cluster_rows = T, cluster_columns = F, column_order = col_ord, row_split = row_ord, cluster_row_slices = F, show_row_dend = F )

col_ord=c("lr_target_means", "hr_target_means",colnames(bla_log2)[grep("Brain|Pituitary", colnames(bla_log2))],
          colnames(bla_log2)[-grep("Brain|target|Pituitary", colnames(bla_log2))])
row_ord=c(rep("tum-med shared", 91), rep("tum-med different", 6))

write.table(bla_log2[,col_ord], "gtex_table.txt", row.names = T, col.names = T, quote = F, sep = "\t")

which(genes_of_i%in%c("SIX3", "PAX6", "CCND1", "CCND2", "PDGFRA"))
