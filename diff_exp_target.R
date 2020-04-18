#prep the target data
target_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv",
                         header = T, sep= "\t")
target_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t",
                        header = T)
target_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")
target_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv",
                                sep = "\t", header = T)
#sort out clinical data NAs
names(which(table(sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 3))==2))

df_bla=target_clinical_data[grep("PARBAJ|PARHAM|PASGAP|PASNPG|PATYIL|PAUDDK",target_clinical_data$TARGET.USI),]

target_clinical_data[39,2:ncol(target_clinical_data)]=target_clinical_data[40,2:ncol(target_clinical_data)]
target_clinical_data[42,2:ncol(target_clinical_data)]=target_clinical_data[78,2:ncol(target_clinical_data)]
target_clinical_data[76,2:ncol(target_clinical_data)]=target_clinical_data[125,2:ncol(target_clinical_data)]
target_clinical_data[96,2:ncol(target_clinical_data)]=target_clinical_data[41,2:ncol(target_clinical_data)]
target_clinical_data[102,2:ncol(target_clinical_data)]=target_clinical_data[80,2:ncol(target_clinical_data)]
target_clinical_data[142,2:ncol(target_clinical_data)]=target_clinical_data[1,2:ncol(target_clinical_data)]

#unify the IDs in some capacity
colnames(target_counts)=paste0(sapply(strsplit(colnames(target_counts), "\\."), "[", 4), "-", sapply(strsplit(colnames(target_counts),
                                                                                                              "\\."), "[", 5))
target_clinical_data$TARGET_SHORT=gsub("NA", "01A",paste0(sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 3),
                                                          "-", sapply(strsplit(as.character(target_clinical_data$TARGET.USI), "-"), "[", 4)))


#define risk groups
low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk_relapse"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"
all=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
extremes=rbind(low_risk_4s, high_risk_relapse)

#get tabs all in the same order 

count_tab=target_counts[,all$TARGET_SHORT]

extreme_count_tab=target_counts[,extremes$TARGET_SHORT]

#make an EdgeR object 


edg_file=DGEList(counts=count_tab, samples = all, genes = target_genes)

o <- order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file <- edg_file[o,]
d <- duplicated(edg_file$genes$external_gene_name)
edg_file <- edg_file[!d,]

keep <- filterByExpr(edg_file)
table(keep)

edg_file <- edg_file[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(edg_file)
design <- model.matrix(~edg_file$samples$risk)
y <- estimateDisp(y,design)
plotBCV(y)
edg_file <- calcNormFactors(edg_file)

edg_file$samples

edg_file <- estimateDisp(edg_file)

plotMD(cpm(edg_file, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)


et <- exactTest(edg_file, pair=c("high_risk_relapse","low_risk_4s"))

et$table
w=topTags(et,n = 5000) 
et$table[which(et$genes$external_gene_name=="ERBB3"),]

edg_file$samples

design <- model.matrix(~0+group, data=edg_file$samples)
colnames(design) <- levels(edg_file$samples$group)

quantile(edg_file$samples$norm.factors)

et <- exactTest(edg_file, dispersion=0.4^2, pair = c("high_risk_relapse","low_risk_4s"))

et$table[which(et$genes$external_gene_name%in%c("ERBB3", "SOX2", "MPZ", "PLP1")),]


#edgeR try with just relevant columns



edg_file=DGEList(extreme_count_tab, samples = extremes, genes = target_genes, group = extremes$risk)

o <- order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file <- edg_file[o,]
d <- duplicated(edg_file$genes$external_gene_name)
edg_file <- edg_file[!d,]
nrow(edg_file)


keep <- filterByExpr(edg_file)
edg_file <- edg_file[keep, , keep.lib.sizes=FALSE]
edg_file <- calcNormFactors(edg_file)
edg_file <- estimateDisp(edg_file)
plotBCV(edg_file)

et <- exactTest(edg_file, dispersion = 0.1^2, pair = c("high_risk_relapse","low_risk_4s"))
et <- exactTest(edg_file, pair = c("high_risk_relapse","low_risk_4s"))
et$table[which(et$genes$external_gene_name%in%c("ERBB3", "SOX2", "MPZ", "PLP1")),]

target_genes[which(target_genes$external_gene_name=="PLP1"),]

