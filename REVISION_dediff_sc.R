merged_cells_gs=merged_cells
rownames(merged_cells_gs)=genes_10x$V2
tum_adr_comm_genes=intersect(rownames(tum_together), rownames(merged_cells_gs))

train_with_podos_sc=trainModel(merged_cells_gs[tum_adr_comm_genes,], cell_labs, workers=NULL)

ps_train_podos_tums=predictSimilarity(train_with_podos_sc, tum_together)

similarityHeatmap(ps_train_podos_tums, row_split=kmeans_2)
high_risk_symps=ps_train_podos_tums[rownames(ps_train_podos_tums)[which(tum_risk_labs=="high risk")],7]
high_risk_only=ps_train_podos_tums[rownames(ps_train_podos_tums)[which(tum_risk_labs=="high risk")],]
length(grep("STDY",names(high_risk_symps[which(high_risk_symps<1.5)])))
plot(density(high_risk_symps))
abline(v=1.5)
Dim


kmeans_2=kmeans(high_risk_symps,centers = 2)
kmeans_2$cluster


similarityHeatmap(high_risk_only, row_split=kmeans_2$cluster)

qc_df=data.frame(cells=names(kmeans_2$cluster), kmeans=kmeans_2$cluster)

qc_df$symp_cat="high_symp"
qc_df$symp_cat[which(qc_df$kmeans==2)]="low_symp"
qc_df$cohort="PMC"
qc_df$cohort[grep("STDY", qc_df$cells)]="GOSH"

qc_df$symp_and_cohort=paste0(qc_df$symp_cat, "_", qc_df$cohort)

qc_df$ngenes=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$nFeature_RNA,
               srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$nFeature_RNA)
qc_df$mito=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$percent.mt,
               srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$percent.mt)
qc_df$reads=c(srat_dutch@meta.data[rownames(qc_df)[which(qc_df$cohort=="PMC")],]$nCount_RNA,
              srat_inhouse@meta.data[rownames(qc_df)[which(qc_df$cohort=="GOSH")],]$nCount_RNA)

pdf("figs_for_revision/dediff_qc.pdf", height = 3, width = 5)
plot(ggplot(qc_df, aes(x=symp_cat, y=log10(ngenes))) +geom_boxplot())
plot(ggplot(qc_df, aes(x=symp_cat, y=log10(reads))) +geom_boxplot())
plot(ggplot(qc_df, aes(x=symp_cat, y=mito)) +geom_boxplot())
dev.off()

pdf("figs_for_revision/dediff_qc_cohort.pdf", height = 3, width = 5)
plot(ggplot(qc_df, aes(x=symp_and_cohort, y=log10(ngenes))) +geom_boxplot())
plot(ggplot(qc_df, aes(x=symp_and_cohort, y=log10(reads))) +geom_boxplot())
plot(ggplot(qc_df, aes(x=symp_and_cohort, y=mito)) +geom_boxplot())
dev.off()

just_high_mtx=tum_together[,rownames(qc_df)]
excludeGenes=c(riboGenes, mtGenes, riboGenes)
keep_genes=which(!rownames(just_high_mtx)%in%excludeGenes)
just_high_srat=CreateSeuratObject(just_high_mtx)
symp_sg=qc_df$symp_cat
names(symp_sg)=names(just_high_srat@active.ident)
just_high_srat@active.ident=factor(symp_sg)
just_high_srat=subset(just_high_srat, features=keep_genes)
just_high_srat=NormalizeData(just_high_srat)
just_high_srat=ScaleData(just_high_srat)
just_high_srat=FindVariableFeatures(just_high_srat)
just_high_srat=RunPCA(just_high_srat, npcs = 50)
just_high_srat = FindNeighbors(just_high_srat, dims = 1:50)
just_high_srat = FindClusters(just_high_srat, resolution = 1)
just_high_srat=RunUMAP(just_high_srat,dims=1:50, min.dist = 0.5, n.neighbors = 50)
high_low_symp_de=FindMarkers(just_high_srat, ident.1 = "high_symp", ident.2 = "low_symp", )
high_low_symp_de=high_low_symp_de[1:100,]
high_low_symp_de$gene=rownames(high_low_symp_de)



write.table(high_low_symp_de, "high_low_symp_de.txt", row.names = F, col.names = T, quote = F, sep="\t")

just_high_srat@meta.data$symp_and_cohort=qc_df$symp_and_cohort
qc_df$symp_signal=high_risk_symps
ggplot(qc_df, aes(x=symp_cat, y=symp_signal)) +geom_violin()


     