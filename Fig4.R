#prep the target data
target_counts=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkCountsFrag.tsv",
                         header = T, sep= "\t")
target_genes=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGenes.tsv", sep = "\t",
                        header = T)
target_meta=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkMetadata.tsv", sep = "\t")
target_clinical_data=read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET_NBL_ClinicalData_withClusters.tsv",
                                sep = "\t", header = T)
target_gene_length= read.table("/lustre/scratch119/casm/team274sb/ek12/Neuroblastoma/bulk_NB/TARGET-NBL_BulkGeneLengths.tsv", sep = "\t")
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

edg_file$samples$group=factor(edg_file$samples$group, levels = c("low_risk_4s", "high_risk"))
which(edg_file$samples$MYCN.status=="Unknown")
edg_file$samples$MYCN.status[80]="Not Amplified"
edg_file$samples$MYCN.status=factor(edg_file$samples$MYCN.status,levels=c("Not Amplified", "Amplified"))

design1=model.matrix(~group, edg_file$samples) 

fit1 = glmQLFit(edg_file, design1, robust=TRUE, dispersion = 0.4^2)
qlf1 = glmQLFTest(fit1,coef=2)

qlf1$table[c("ERBB3", "SOX2", "MPZ", "PLP1"),]

design3=model.matrix(~group+Age.at.Diagnosis.in.Days+MYCN.status,edg_file$samples)
fit3=glmQLFit(edg_file, design3, robust=TRUE, dispersion = 0.4^2)
qlf3=glmQLFTest(fit3,coef=2)

de_of_de=qlf3$table[as.character(final_de_sorted$gene)[which(as.character(final_de_sorted$gene)%in%rownames(qlf3$table))],]
de_of_de_sorted=de_of_de[order(abs(de_of_de$logFC), decreasing = T),]


de_of_de_pat=as.data.frame(percent_above_threshold[rownames(de_of_de_sorted)])
de_of_de_pat$gene=rownames(de_of_de_pat)

colnames(de_of_de_pat)=c("value", "gene")
de_of_de_pat=de_of_de_pat[rownames(de_of_de_sorted)[which(!rownames(de_of_de_sorted)%in%immune_sc)],]
de_of_de_pat=de_of_de_pat[order(de_of_de_pat$value, decreasing = T),]
de_of_de_pat$gene=factor(de_of_de_pat$gene, levels = de_of_de_pat$gene)


ggplot(de_of_de_pat, aes(x=reorder(gene, value), y=value))+geom_bar(stat="identity", width = 0.6) +coord_flip() +
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", face = "italic", size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 12)) + labs( y ="% of samples above expression threshold", x = "Gene")
above_log$qval=p.adjust(above)
above_log=de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),][which(abs(de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),]$logFC)>1),]
above_log=above_log[which(above_log$qval<0.05),]
genes_for_kap=rownames(above_log)


de_of_de_facet_tab=target_log2tpm[genes_for_kap, all$TARGET_SHORT]

de_tab=as.data.frame(t(de_of_de_facet_tab))  
de_tab$risk=all$risk  
de_tab$patient=rownames(de_tab)

de_tab_melted=melt(de_tab, id.vars=c("patient", "risk"))
de_tab_melted$risk=factor(de_tab_melted$risk, levels = c("low_risk_4s", "intermediate_risk","high_risk"))

de_tab_melted$variable=factor(de_tab_melted$variable, levels = c(rownames(de_of_de_sorted)))


ggplot(de_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("DE genes")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = de_tab_melted, aes(x=risk,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+facet_wrap(~variable,scales="free_y")+
  geom_pointrange(mapping = aes(x = risk, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +stat_compare_means()


x1=data.frame(group1="low_risk_4s", group2="high_risk", p.adj=above_log$qvalue, gene=rownames(above_log))
x1$p.adj=round(x1$p.adj, 4)
bp <- list()
for (i in genes_for_kap) {
  
  bp[[i]] <- ggboxplot(de_tab_melted[which(de_tab_melted$variable==i),], x = "risk", y = "value", width = 0.4, outlier.shape=NA) +
    stat_pvalue_manual(
      as.data.frame(x1[which(x1$gene==i),]), label = "p = {p.adj}",
      y.position = max(de_tab_melted[which(de_tab_melted$variable==i),]$value)-1) +
    ggtitle(i) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y="Expression (log2(TPM))")
}
do.call(grid.arrange, bp)






mapk1=read.table("mapk1.tsv", sep = "\t", header = T)
mapk2=read.table("mapk2.tsv", sep = "\t", header = T)


mapk_cancer_1=read.table("mapk_cancer1.tsv", sep = "\t", header = T)
mapk_cancer_2=read.table("mapk_cancer2.tsv", sep = "\t", header = T)
mapk_cancer_3=read.table("mapk_cancer3.tsv", sep = "\t", header = T)

nb_mapk=read.table("mapk_in_nb.txt", sep = "\t", header = F)

nb_mapk_genes=unique(as.character(nb_mapk$V1))

general_mapk_genes=unique(c(as.character(mapk1$Gene.Name), as.character(mapk2$Gene.Name)))

cancer_mapk_genes=unique(c(as.character(mapk_cancer_1$Gene.Name), as.character(mapk_cancer_2$Gene.Name),as.character(mapk_cancer_3$Gene.Name)))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk_relapse"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

target_risks=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
target_risks$patient=target_risks$TARGET_SHORT
get_mapk_per_risk=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
    })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  #scp_markers_in_rg=values_per_risk_group[general_mapk_genes[which(general_mapk_genes%in%rownames(values_per_risk_group))],]
  bridge_markers_in_rg=values_per_risk_group[mapk_sc[which(mapk_sc%in%rownames(values_per_risk_group))],]
  
  
  #scp_markers_in_rg=melt(scp_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  bridge_markers_in_rg=melt(bridge_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  
  #scp_markers_in_rg$marker="general_mapk"
  bridge_markers_in_rg$marker="cancer_mapk"
  
  #all_in_groups=rbind(scp_markers_in_rg, bridge_markers_in_rg)
  return(bridge_markers_in_rg)
  
}

target_mapk_per_risk=get_mapk_per_risk(target_log2tpm, target_risks,5)
german_mapk_per_risk=get_mapk_per_risk(german_counts_log2, german_risks, 10)
seqc_mapk_per_risk=get_mapk_per_risk(seqc_tpm_log2, seqc_risks, 5)

target_mapk_per_risk$cohort="target"
german_mapk_per_risk$cohort="german"
seqc_mapk_per_risk$cohort="seqc"


combined_mapk_per_risk=rbind(target_mapk_per_risk, german_mapk_per_risk, seqc_mapk_per_risk)

colnames(combined_mapk_per_risk)=c("gene", "risk", "value", "marker", "cohort")

combined_mapk_per_risk$risk=factor(combined_mapk_per_risk$risk, levels=c("low_risk","low_risk_4s", "intermediate_risk", "high_risk", "high_risk_relapse"))


combined_mapk_per_risk$cohort=factor(combined_mapk_per_risk$cohort, levels = c("target", "seqc", "german"))
ggplot(data = combined_mapk_per_risk, aes(x=risk,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3) +facet_wrap(vars(cohort)) + theme( axis.line = element_line(colour = "black"))+
  geom_pointrange(mapping = aes(x = risk, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))

mapk_tg_lr=target_mapk_per_risk$value[which(target_mapk_per_risk$variable=="low_risk_4s")]
mapk_tg_hr=target_mapk_per_risk$value[which(target_mapk_per_risk$variable=="high_risk")]
mapk_tg_hrl=target_mapk_per_risk$value[which(target_mapk_per_risk$variable=="high_risk_relapse")]

mapk_seqc_lr=seqc_mapk_per_risk$value[which(seqc_mapk_per_risk$variable=="low_risk")]
mapk_seqc_hr=seqc_mapk_per_risk$value[which(seqc_mapk_per_risk$variable=="high_risk")]

mapk_german_lr=german_mapk_per_risk$value[which(german_mapk_per_risk$variable=="low_risk")]
mapk_german_hr=german_mapk_per_risk$value[which(german_mapk_per_risk$variable=="high_risk")]
wilcox.test(mapk_tg_lr, mapk_tg_hr)$p.val
wilcox.test(mapk_tg_lr, mapk_tg_hrl)$p.val
wilcox.test(mapk_tg_hr, mapk_tg_hrl)$p.val
wilcox.test(mapk_seqc_lr, mapk_seqc_hr)$p.val
wilcox.test(mapk_german_lr, mapk_german_hr)$p.val

mapk_wilcox=data.frame(mod_1=c("mapk_tg_lr", "mapk_tg_lr", "mapk_tg_hr","mapk_seqc_lr", "mapk_german_lr"), 
                       mod2=c("mapk_tg_hr", "mapk_tg_hrl","mapk_tg_hrl","mapk_seqc_hr", "mapk_german_hr"),
                       p.val=c(3.050762e-05,0.0009250993,0.4321568,0.03545608,0.8379161))

wilcox.test(kidney_mapk$value, mapk_tg_lr)
wilcox.test(kidney_mapk$value, mapk_tg_hr)
wilcox.test(kidney_mapk$value, mapk_tg_hrl)

wilcox.test(kidney_mapk$value, mapk_seqc_lr)
wilcox.test(kidney_mapk$value, mapk_seqc_hr)

wilcox.test(kidney_mapk$value, mapk_german_lr)
wilcox.test(kidney_mapk$value, mapk_german_hr)

combined_mapk_per_risk

combined_mapk_per_risk


mapk_with_kidney=rbind(combined_mapk_per_risk, kidney_mapk)
mapk_with_kidney


ggplot(data = mapk_with_kidney, aes(x=risk,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3) +facet_wrap(vars(cohort), scales = "free_x") + theme( axis.line = element_line(colour = "black"))+
  geom_pointrange(mapping = aes(x = risk, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))


de_of_de_facet_tab=target_log2tpm[genes_for_kap, all$TARGET_SHORT]

de_tab=as.data.frame(t(de_of_de_facet_tab))  
de_tab$risk=all$risk  
de_tab$patient=rownames(de_tab)

de_tab_melted=melt(de_tab, id.vars=c("patient", "risk"))
de_tab_melted$risk=factor(de_tab_melted$risk, levels = c("low_risk_4s", "intermediate_risk","high_risk"))

de_tab_melted$variable=factor(de_tab_melted$variable, levels = c(rownames(de_of_de_sorted)))


ggplot(de_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("DE genes")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = de_tab_melted, aes(x=risk,y=value)) +
  geom_quasirandom(shape=1, alpha=0.3)+facet_wrap(~variable,scales="free_y")+
  geom_pointrange(mapping = aes(x = risk, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) +stat_compare_means() 
