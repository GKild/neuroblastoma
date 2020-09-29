#array data DE 

german_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_series_matrix.txt", header = T, sep= "\t", skip = 63, nrows =44708)

german_genes=read.delim("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_family_soft.txt", sep = "\t", header = T, skip=471, nrows = 44708, na.strings = "")


german_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/bulk_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)

german_sample_ids=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/German_bulkNB/GSE120572_series_matrix.txt", header = F, skip=26, nrows = 1)

x=as.data.frame(t(german_sample_ids))

colnames(german_counts)[2:395]=as.character(x$V1[2:395])


no_na_gs=german_genes[!is.na(german_genes$GeneSymbol),]
no_na_gs$ID
counts_gs=german_counts[which(german_counts$ID_REF%in%no_na_gs$ID),] 

counts_gs$gene_symbol=no_na_gs$GeneSymbol 
german_clinical_data=german_clinical_data[which(german_clinical_data$Tumor.ID%in%colnames(counts_gs)),]
german_clinical_data$Age.at.Diagnosis..d.=as.numeric(as.character(german_clinical_data$Age.at.Diagnosis..d.))

low_risk=dplyr::filter(german_clinical_data,
                       german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &!german_clinical_data$Status%in%c("relapse/progression","death of disease"))
high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified", german_clinical_data$Status%in%c("relapse/progression","death of disease"))
risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")
pt1=rbind(low_risk, high_risk, risk_4s)
pt2=german_clinical_data[-which(german_clinical_data$Tumor.ID%in%pt1$Tumor.ID),]
inter_death_or_relapse=dplyr::filter(pt2,Status%in%c("relapse/progression","death of disease"))
inter_no_event=pt2[-which(pt2$Tumor.ID%in%inter_death_or_relapse$Tumor.ID),]
low_risk$risk_level="low_risk"
inter_no_event$risk_level="inter_no_event"
inter_death_or_relapse$risk_level="inter_death_or_relapse"
high_risk$risk_level="high_risk"
risk_4s$risk_level="4s"

extremes=rbind(low_risk, high_risk)
all=rbind(low_risk, risk_4s, inter_no_event, inter_death_or_relapse, high_risk)

dim(extremes)

counts_gs$gene_symbol=make.unique(as.character(counts_gs$gene_symbol))
rownames(counts_gs)=counts_gs$gene_symbol
counts_gs=counts_gs[,2:(ncol(counts_gs)-1)]

extreme_counts_german=counts_gs[,extremes$Tumor.ID]


bla=apply(counts_gs, 1, function(x){unlist(x)})

bla_t=t(bla)

bla_t=bla_t[,extremes$Tumor.ID]


german_counts_log2=log2(bla_t)

plot(density(rowMeans(german_counts_log2)))
ncol(counts_log2)
hist(rowMeans(counts_log2),100)
threshold_10 <-apply(counts_log2, 1,
                     function(x){
                       sum(x >5)>=90 })

all_genes_lm=apply(counts_log2,1,function(x){
  exp_mod=lm(x~extremes$risk_level+extremes$Age.at.Diagnosis..d.)
  coefs=c(summary(exp_mod)$coefficients[1,1],summary(exp_mod)$coefficients[2,1],summary(exp_mod)$coefficients[2,4])
  return(coefs)
})

all_genes_lm=t(all_genes_lm)
colnames(all_genes_lm)=c("intercept", "slope", "pval")
all_genes_lm=as.data.frame(all_genes_lm)
all_genes_lm$qval=p.adjust(all_genes_lm$pval, method = "BH")
german_all_mark=all_genes_lm[all_marks[which(all_marks%in%rownames(all_genes_lm))],]

german_all_mark=as.data.frame(german_all_mark)
german_all_mark$qval=p.adjust(german_all_mark$pval, method = "BH")

german_all_mark[german_all_mark$qval<0.01,]

x1=german_all_mark[scp_marks[which(scp_marks%in%rownames(german_all_mark))],]

x1[x1$qval<0.01,]

dim(all_genes_lm[all_genes_lm$qval<0.01,])

hist(x1$slope,30)

table(table(german_genes$GeneSymbol))
hist(german_all_mark[which(!rownames(german_all_mark)%in%scp_marks),]$slope,30)




german_genes_above_threshold=names(rowMeans(german_counts_log2)[rowMeans(german_counts_log2)>10])
german_percent_above_threshold=apply(german_counts_log2, 1, function(x){
  sum(x>10)/length(x)*100
})


