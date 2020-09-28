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


target_length_norm =target_counts/(target_gene_length/1000)

scale_factor=colSums(target_length_norm)/1000000
target_tpm=target_length_norm/scale_factor
target_log2tpm=log2(target_tpm+1)

rownames(target_log2tpm)=make.unique(as.character(target_genes$external_gene_name))
rownames(target_tpm)=make.unique(as.character(target_genes$external_gene_name))
genes_above_threshold=names(rowMeans(target_log2tpm)[rowMeans(target_log2tpm)>5])
percent_above_threshold=apply(target_log2tpm, 1, function(x){
  sum(x>5)/length(x)*100
})
#define risk groups
low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"
all=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
extremes=rbind(low_risk_4s, high_risk, high_risk_relapse)

#get tabs all in the same order 

count_tab=target_counts[,all$TARGET_SHORT]

extreme_count_tab=target_counts[,extremes$TARGET_SHORT]

#make an EdgeR object 
dim(intermediate_risk)

edg_file=DGEList(counts=extreme_count_tab, samples = extremes, genes = target_genes, group=extremes$risk)

o = order(rowSums(edg_file$counts), decreasing=TRUE)
edg_file = edg_file[o,]
d = duplicated(edg_file$genes$external_gene_name)
edg_file = edg_file[!d,]

rownames(edg_file)=edg_file$genes$external_gene_name

keep = filterByExpr(edg_file)
table(keep)

edg_file = edg_file[keep, , keep.lib.sizes=FALSE]

edg_file=calcNormFactors(edg_file)

#et=exactTest(edg_file, dispersion = 0.4^2)
#et$table[which(et$genes$external_gene_name%in%c("ERBB3", "SOX2", "MPZ", "PLP1")),]

edg_file$samples$group=factor(edg_file$samples$group, levels = c("low_risk_4s", "high_risk"))

edg_file$samples$MYCN.status=factor(edg_file$samples$MYCN.status,levels=c("Not Amplified", "Amplified"))

design1=model.matrix(~group, edg_file$samples)

fit1 = glmQLFit(edg_file, design1, robust=TRUE, dispersion = 0.4^2)
qlf1 = glmQLFTest(fit1,coef=2)

qlf1$table[c("ERBB3", "SOX2", "MPZ", "PLP1"),]

qlf1$table$qval=p.adjust(qlf1$table$PValue, method = "BH")
length(which(qlf1$table$qval<0.01))


design2=model.matrix(~group+Age.at.Diagnosis.in.Days,edg_file$samples)

fit2 = glmQLFit(edg_file, design2, robust=TRUE, dispersion = 0.4^2)
qlf2=glmQLFTest(fit2,coef=2)

qlf2$table[c("ERBB3", "SOX2", "MPZ", "PLP1"),]

qlf2$table$qval=p.adjust(qlf2$table$PValue, method = "BH")
length(which(qlf2$table$qval<0.01))

design3=model.matrix(~group+Age.at.Diagnosis.in.Days+MYCN.status,edg_file$samples)
fit3=glmQLFit(edg_file, design3, robust=TRUE, dispersion = 0.4^2)
qlf3=glmQLFTest(fit3,coef=2)

qlf3$table[c("ERBB3", "SOX2", "MPZ", "PLP1"),]

qlf3_mycn=glmQLFTest(fit3,coef=4)
topTags(qlf3_mycn)
qlf3$table$qval=p.adjust(qlf3$table$PValue, method = "BH")
length(which(qlf3$table$qval<0.01))


tb=qlf3$table[which(qlf3$table$qval<0.01),]


length(which(hkGenes%in%rownames(tb)))


qlf3$table[scp_marks[which(scp_marks%in%rownames(qlf3))],]
tfidf_2=filter(tr2mark_filtered, tfidf>2)
scp_marks=as.character(tr2mark_filtered[which(tr2mark_filtered$cluster=="SCPs"),]$gene)
bridge_marks=as.character(tr2mark_filtered[which(tr2mark_filtered$cluster=="bridge"),]$gene)
right_marks=as.character(tr2mark_filtered[which(tr2mark_filtered$cluster=="right_clust"),]$gene)
left_marks=as.character(tr2mark_filtered[which(tr2mark_filtered$cluster=="left_clust"),]$gene)
all_marks=as.character(tr2mark_filtered$gene)

scp_marks=as.character(tfidf_2[which(tfidf_2$cluster=="SCPs"),]$gene)
bridge_marks=as.character(tfidf_2[which(tfidf_2$cluster=="bridge"),]$gene)
right_marks=as.character(tfidf_2[which(tfidf_2$cluster=="right_clust"),]$gene)
left_marks=as.character(tfidf_2[which(tfidf_2$cluster=="left_clust"),]$gene)
all_marks=as.character(tfidf_2$gene)

first_scp_mark=qlf3$table[scp_marks[which(scp_marks%in%rownames(qlf3))],]
first_bridge_mark=qlf3$table[bridge_marks[which(bridge_marks%in%rownames(qlf3))],]
first_right_mark=qlf3$table[right_marks[which(right_marks%in%rownames(qlf3))],]
first_left_mark=qlf3$table[left_marks[which(left_marks%in%rownames(qlf3))],]

all_first_mark=qlf3$table[all_marks[which(all_marks%in%rownames(qlf3))],]
first_scp_mark[first_scp_mark$qval<0.05,]
first_bridge_mark[first_bridge_mark$qval<0.05,]
first_right_mark[first_right_mark$qval<0.05,]
first_left_mark[first_left_mark$qval<0.05,]

all_first_mark$qval=p.adjust(all_first_mark$PValue, method = "BH")

dim(all_first_mark[all_first_mark$qval<0.01,])

get_target_marker_plots(target_tpm, scp_marks, right_marks , left_marks,"global_adr_markers")


length(which(podo_marks_filt$gene%in%genes_above_threshold))/length(podo_marks_filt$gene)

length(unique(all_marks))

low_risk_log2=target_log2tpm[,low_risk_4s$TARGET_SHORT]

target_lr_genes_above_threshold=names(rowMeans(low_risk_log2)[rowMeans(low_risk_log2)>5])
length(which(scp_marks%in%target_lr_genes_above_threshold))/length(scp_ma)

genes_above_threshold_3=names(rowMeans(target_log2tpm)[rowMeans(target_log2tpm)>3])

length(which(podo_marks_filt$gene%in%genes_above_threshold_3))


target_log2_10x_genes=target_log2tpm[genes_10x$V2,]

plot(density(rowMeans(target_log2_10x_genes)))
mean(rowMeans(target_log2_10x_genes))

t_scp=as.data.frame(percent_above_threshold[scp_marks[which(scp_marks%in%names(percent_above_threshold))]])
t_bridge=as.data.frame(percent_above_threshold[bridge_marks[which(bridge_marks%in%names(percent_above_threshold))]])
t_left=as.data.frame(percent_above_threshold[left_marks[which(left_marks%in%names(percent_above_threshold))]])
t_right=as.data.frame(percent_above_threshold[right_marks[which(right_marks%in%names(percent_above_threshold))]])
t_pods=as.data.frame(percent_above_threshold[as.character(podo_marks_filt$gene)[which(as.character(podo_marks_filt$gene)%in%names(percent_above_threshold))]])

t_scp$marker="scp"
t_bridge$marker="bridge"
t_left$marker="left"
t_right$marker="right"
t_pods$marker="pods"

t_scp$gene=rownames(t_scp)
t_bridge$gene=rownames(t_bridge)
t_left$gene=rownames(t_left)
t_right$gene=rownames(t_right)
t_pods$gene=rownames(t_pods)

rownames(t_scp)=1:nrow(t_scp)
rownames(t_bridge)=1:nrow(t_bridge)
rownames(t_left)=1:nrow(t_left)
rownames(t_right)=1:nrow(t_right)
rownames(t_pods)=1:nrow(t_pods)

colnames(t_scp)=c("value", "marker", "gene")
colnames(t_bridge)=c("value", "marker", "gene")
colnames(t_left)=c("value", "marker", "gene")
colnames(t_right)=c("value", "marker","gene")
colnames(t_pods)=c("value", "marker","gene")

t_all_m=rbind(t_scp, t_bridge,t_left,t_right, t_pods)


ggplot(t_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >5 log2(TPM) in TARGET cohort (N=498)")

which(left_marks%in%bridge_marks)



marks_de=qlf3$table[as.character(de_in_tg$gene)[which(as.character(de_in_tg$gene)%in%rownames(qlf3$table))],]

marks_de_ordered=marks_de[order(abs(marks_de$logFC), decreasing = T),]
