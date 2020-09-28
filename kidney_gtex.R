gen <- read.delim( "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm (1).gct", skip = 2, header = TRUE )

## First two columns corresponds to "Transcript name" and "Gene name"
gen_ann <- gen[ , c(1,2) ]
gen_clean <- as.matrix( gen[ , -seq( 2 ) ] )
rm(gen)
# Load map file
sample_map <- read.delim( "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" )
sample_map$SAMPID_new <- gsub( "-", "\\.", sample_map$SAMPID )
rownames( sample_map ) <- sample_map$SAMPID

adrenal_gland_ids=sample_map$SAMPID_new[which(sample_map$SMTS=="Kidney")]
adrenal_overlap=adrenal_gland_ids[which(adrenal_gland_ids%in%colnames(gen_clean))]
gtex_adrenal_gland=gen_clean[,adrenal_overlap]


length(unique(gen_ann$Description))
gen_ann$Description=make.unique(as.character(gen_ann$Description))

rownames(gtex_adrenal_gland)=gen_ann$Description


adrenal_log2tpm=log2(gtex_adrenal_gland+1)


plot(density(rowMeans(adrenal_log2tpm)))
hist(rowMeans(adrenal_log2tpm),30)

plot(x$breaks, x$counts, log='x', type='h')
adr_percent_above_threshold=apply(adrenal_log2tpm, 1, function(x){
  sum(x>3.5)/length(x)*100}) 
adr_mapk_cancer=data.frame(adr_percent_above_threshold[cancer_mapk_genes])
adr_mapk_general=data.frame(adr_percent_above_threshold[mapk_sc[which(mapk_sc%in%names(adr_percent_above_threshold))]])


colnames(adr_mapk_cancer)="value"
colnames(adr_mapk_general)="value"


adr_mapk_cancer$set="cancer_mapk"
adr_mapk_general$set="general_mapk"

adr_mapk=rbind(adr_mapk_general)





ggplot(data = adr_mapk, aes(x=set,y=value)) +
  geom_pointrange(mapping = aes(x = set, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  geom_quasirandom(shape=1) + ggtitle("% samples MAPK genes expressed at >4 log2(tpm) in kidney (GTEx)")

adr_mapk_general$gene=rownames(adr_mapk_general)

colnames(adr_mapk_general)=c("value", "gene")


kidney_mapk=data.frame(gene=adr_mapk_general$gene, risk="MAPK in kidney", value=adr_mapk_general$value, marker="cancer_mapk", cohort="GTEx kidney")



