seqc_counts=read.table("seqc/GSE49711_SEQC_NB_TUC_G_log2.txt", sep = "\t", header = T, row.names = T)
seqc_meta=read.csv("seqc/seqc_meta.csv")

rownames(seqc_counts)=seqc_counts$X00gene_id

seqc_counts=seqc_counts[,2:ncol(seqc_counts)]
colnames(seqc_counts)=sapply(strsplit(colnames(seqc_counts), "_"), "[", 2)


seqc_unlist=apply(seqc_counts, 1, function(x){unlist(x)})

seqc_unlist_t=t(seqc_unlist)



untransf=(2^seqc_unlist_t)-1

threshold_10 <-apply(untransf, 1,
                     function(x){
                       sum(x >0)>=50 })
filtered_seqc=untransf[threshold_10,]

seqc_tpm=apply(filtered_seqc,2,function(x){
  (x/sum(x))*1000000
})


seqc_tpm_log2=log2(seqc_tpm+1)

seqc_lr=dplyr::filter(seqc_meta, Age<540 & MYCN.status=="1", INSS.Stage!="4S" & Class.label=="favorable")
seqc_hr=dplyr::filter(seqc_meta, Age>540 & MYCN.status=="Amp" &HighRisk=="HR")

seqc_lr$risk="low_risk"
seqc_hr$risk="high_risk"

seqc_extremes=rbind(seqc_lr, seqc_hr)

seqc_counts_extreme=seqc_tpm_log2[,seqc_extremes$NB.ID]


all_genes_lm_seqc=apply(seqc_counts_extreme,1,function(x){
  exp_mod=lm(x~seqc_extremes$risk+seqc_extremes$Age)
  coefs=c(summary(exp_mod)$coefficients[1,1],summary(exp_mod)$coefficients[2,1],summary(exp_mod)$coefficients[2,4])
  return(coefs)
})

all_genes_lm_seqc=t(all_genes_lm_seqc)
colnames(all_genes_lm_seqc)=c("intercept", "slope", "pval")
all_genes_lm_seqc=as.data.frame(all_genes_lm_seqc)
all_genes_lm_seqc$qval=p.adjust(all_genes_lm_seqc$pval, method = "BH")
seqc_all_mark=all_genes_lm_seqc[all_marks[which(all_marks%in%rownames(all_genes_lm_seqc))],]

plot(density(rowMeans(seqc_tpm_log2)))
seqc_genes_above_threshold=names(rowMeans(seqc_tpm_log2)[rowMeans(seqc_tpm_log2)>5])
length(which(podo_marks_filt$gene%in%seqc_genes_above_threshold))


seqc_low_risk_log2=seqc_tpm_log2[,as.character(seqc_lr$NB.ID)]
seqc_lr_genes_above_threshold=names(rowMeans(seqc_low_risk_log2)[rowMeans(seqc_low_risk_log2)>4])

length(which(scp_marks%in%seqc_lr_genes_above_threshold))/length(scp_marks)
rowMeans(seqc_low_risk_log2)["S100B"]


seqc_percent_above_threshold=apply(seqc_tpm_log2, 1, function(x){
  sum(x>5)/length(x)*100
})

x=as.data.frame(seqc_percent_above_threshold[scp_marks[which(scp_marks%in%names(seqc_percent_above_threshold))]])


seqc_scp=as.data.frame(seqc_percent_above_threshold[scp_marks[which(scp_marks%in%names(seqc_percent_above_threshold))]])
seqc_bridge=as.data.frame(seqc_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(seqc_percent_above_threshold))]])
seqc_left=as.data.frame(seqc_percent_above_threshold[left_marks[which(left_marks%in%names(seqc_percent_above_threshold))]])
seqc_right=as.data.frame(seqc_percent_above_threshold[right_marks[which(right_marks%in%names(seqc_percent_above_threshold))]])
seqc_pods=as.data.frame(seqc_percent_above_threshold[as.character(podo_marks_filt$gene)[which(as.character(podo_marks_filt$gene)%in%names(seqc_percent_above_threshold))]])

seqc_scp$marker="scp"
seqc_bridge$marker="bridge"
seqc_left$marker="left"
seqc_right$marker="right"
seqc_pods$marker="pods"

seqc_scp$gene=rownames(seqc_scp)
seqc_bridge$gene=rownames(seqc_bridge)
seqc_left$gene=rownames(seqc_left)
seqc_right$gene=rownames(seqc_right)
seqc_pods$gene=rownames(seqc_pods)

rownames(seqc_scp)=1:nrow(seqc_scp)
rownames(seqc_bridge)=1:nrow(seqc_bridge)
rownames(seqc_left)=1:nrow(seqc_left)
rownames(seqc_right)=1:nrow(seqc_right)
rownames(seqc_pods)=1:nrow(seqc_pods)

colnames(seqc_scp)=c("value", "marker", "gene")
colnames(seqc_bridge)=c("value", "marker", "gene")
colnames(seqc_left)=c("value", "marker", "gene")
colnames(seqc_right)=c("value", "marker","gene")
colnames(seqc_pods)=c("value", "marker","gene")

seqc_all_m=rbind(seqc_scp, seqc_bridge,seqc_left,seqc_right, seqc_pods)


ggplot(seqc_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() + ggtitle("% samples a marker is expressed at >5 log2(TPM) in SEQC cohort (N=498)")

