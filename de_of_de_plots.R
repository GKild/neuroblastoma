all_genes=as.character(tum_n_de_final$gene)

t_de=as.data.frame(percent_above_threshold[all_genes[which(all_genes%in%names(percent_above_threshold))]])
seqc_de=as.data.frame(seqc_percent_above_threshold[all_genes[which(all_genes%in%names(seqc_percent_above_threshold))]])
german_de=as.data.frame(german_percent_above_threshold[all_genes[which(all_genes%in%names(german_percent_above_threshold))]])


t_de$gene=rownames(t_de)
seqc_de$gene=rownames(seqc_de)
german_de$gene=rownames(german_de)

rownames(t_de)=1:nrow(t_de)
rownames(seqc_de)=1:nrow(seqc_de)
rownames(german_de)=1:nrow(german_de)

colnames(t_de)=c("value", "gene")
colnames(seqc_de)=c("value", "gene")
colnames(german_de)=c("value", "gene")

t_de$dataset="TARGET"
seqc_de$dataset="SEQC"
german_de$dataset="GERMAN"


de_all=rbind(t_de, seqc_de, german_de)

nb_genes=c("MYCN", "ALK", "PHOX2B", "ATRX")
t_nb=as.data.frame(percent_above_threshold[nb_genes[which(nb_genes%in%names(percent_above_threshold))]])
seqc_nb=as.data.frame(seqc_percent_above_threshold[nb_genes[which(nb_genes%in%names(seqc_percent_above_threshold))]])
german_nb=as.data.frame(german_percent_above_threshold[nb_genes[which(nb_genes%in%names(german_percent_above_threshold))]])


t_nb$gene=rownames(t_nb)
seqc_nb$gene=rownames(seqc_nb)
german_nb$gene=rownames(german_nb)

rownames(t_nb)=1:nrow(t_nb)
rownames(seqc_nb)=1:nrow(seqc_nb)
rownames(german_nb)=1:nrow(german_nb)

colnames(t_nb)=c("value", "gene")
colnames(seqc_nb)=c("value", "gene")
colnames(german_nb)=c("value", "gene")

t_nb$dataset="TARGET"
seqc_nb$dataset="SEQC"
german_nb$dataset="GERMAN"


nb_all=rbind(t_nb, seqc_nb, german_nb)
nb_all$marker="NB_genes"

de_all$marker="singleCell_de_genes"

de_and_nb=rbind(de_all, nb_all)
ggplot(de_all, aes(x=dataset, y=value)) + 
  geom_boxplot(outlier.shape = NA) +geom_jitter(color="black", size=0.4, alpha=0.9)+ggtitle("% samples DE gene is expressed at above threshold in all cohorts")

ggplot(de_all, aes(x=gene, y=value))

de_of_de=qlf3$table[all_genes[which(all_genes%in%rownames(qlf3))],]

de_of_de$qval=p.adjust(de_of_de$PValue, method="fdr")


ggplot(de_and_nb, aes(x=marker, y=value, fill=dataset)) + 
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge(jitter.width = 0.2), size=0.75) +ggtitle("% samples a gene is expressed above a threshold across 3 cohorts (N=1050)")


ggplot(de_and_nb, aes(x = marker, y = value, fill = dataset)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.2))


de_of_marks=qlf3$table[as.character(podo_marks_filt$gene)[which(as.character(podo_marks_filt$gene)%in%rownames(qlf3))],]
de_of_marks$qval=p.adjust(de_of_marks$PValue, method = "BH")

table(podo_marks_filt[which(podo_marks_filt$gene%in%rownames(de_of_marks[de_of_marks$qval<0.01,])),]$cluster)




ra_genes=read.csv("ra_de.csv")


mapk_genes=c("ALK", "PTPN11", "FGFR1", "NF1", "NRAS", "HRAS", "KRAS", "BRAF")

german_percent_above_threshold[mapk_genes]


adr_all_med=subset(adr_all_new, idents=c(25,24,13,19,20))
adr_all_med@active.ident

new_idents=adr_all_med@meta.data$new_clust 
names(new_idents)=names(adr_all_med@active.ident)
new_idents=factor(new_idents, levels = c("SCPs", "bridge", "left_clust","right_clust" ))

new_med=adr_all_med

new_med@active.ident=new_idents

new_med@active.ident
new_avg=AverageExpression(new_med, features = rownames(new_med), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)

ra_genes_rel=ra_genes[which(ra_genes$External.Gene.ID%in%rownames(new_med)),]

ra_genes_rel_sorted=ra_genes_rel[order(abs(ra_genes_rel$LogFC), decreasing = T),]
col_fun2=colorRamp2(c(-2,-1,0,1,2), rev(brewer.pal(n = 5, name = "RdYlBu")))
col_log=colorRamp2(c(-4,0,4), c("blue", "white", "red"))
col_bh=colorRamp2(c(0.02,0), c("white", "red"))

ha_germ = rowAnnotation(log_fc=ra_genes_rel$LogFC, bh=ra_genes_rel$Benjamini.Hochberg,
                        col = list(log_fc = col_log,
                                   bh= col_bh))

Heatmap(new_avg@assays$RNA@scale.data[as.character(ra_genes_rel$External.Gene.ID),], 
        col=col_fun2, cluster_rows = T, cluster_columns = F, show_column_names = T, left_annotation = ha_germ, show_row_dend = F, show_row_names = F)
ha_germ = rowAnnotation(log_fc=ra_genes_rel_sorted$LogFC, bh=ra_genes_rel_sorted$Benjamini.Hochberg,
                        col = list(log_fc = col_log,
                                   bh= col_bh))
Heatmap(new_avg@assays$RNA@scale.data[as.character(ra_genes_rel_sorted$External.Gene.ID),], 
        col=col_fun2, cluster_rows = F, cluster_columns = F, show_column_names = T, left_annotation = ha_germ, show_row_dend = F, show_row_names = F)
