m13=read.delim("GSE99933_E13.5_counts.txt.gz")
m13=as.matrix(m13)

m13
m13_srat=CreateSeuratObject(m13)
m13_srat = NormalizeData(m13_srat)
m13_srat =FindVariableFeatures(m13_srat, selection.method = "vst", nfeatures = 2000)
m13_srat = ScaleData(m13_srat, features = rownames(m13_srat))
m13_srat = RunPCA(m13_srat, npcs = 25)
m13_srat = FindNeighbors(m13_srat, dims=1:25)
m13_srat = FindClusters(m13_srat, resolution = 0.6)
m13_srat = RunUMAP(m13_srat, dims=1:25, min.dist = 0.5, n.neighbors = 30)

DimPlot(m13_srat, label=T)
FeaturePlot(m13_srat, c("Sox10", "Mpz", "Dll3","Tlx2", "Th", "Cartpt", "Pnmt","Gap43", "Mapt"))

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(m13) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
dim(m13_srat)
dim(genesV2)
mouse_proper=convertMouseGeneList(rownames(m13_srat))

genesV2$MGI.symbol[grep("Mpz",genesV2$MGI.symbol)]


m13_srat@meta.data$new_clust="x"

m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==0)]="sympathoblast"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==1)]="chromaffin"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==2)]="scp"
m13_srat@meta.data$new_clust[which(m13_srat@meta.data$RNA_snn_res.0.6==3)]="bridge"
DimPlot(m13_srat, group.by = "new_clust", label=T, pt.size = 0.5) +theme(text = element_text(size=20), axis.text.x = element_blank(),axis.text.y = element_blank(),
                                                                           axis.ticks = element_blank() ,legend.position = "none")+ylab("UMAP 2") +xlab("UMAP 1") 

FeaturePlot(m13_srat, c("Sox10", "Erbb3", "Foxd3","Htr3a","Dll3","Thbd", "Th", "Chgb","Foxq1", "Cartpt", "Slc18a3"))


chromaffin_for_m=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="chromaffin"),]$gene)[1:20]
sympathoblastic_form_m=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="sympathoblastic"),]$gene)[1:20]

right_as_mouse=genesV2$MGI.symbol[which(genesV2$HGNC.symbol%in%chromaffin_for_m)]
left_as_mouse=genesV2$MGI.symbol[which(genesV2$HGNC.symbol%in%sympathoblastic_form_m)]


FeaturePlot(m13_srat, right_as_mouse)
FeaturePlot(m13_srat, left_as_mouse)

length(which(genesV2$MGI.symbol%in%rownames(m13_srat)))

t.first <- genesV2[match(unique(as.character(genesV2$HGNC.symbol)), as.character(genesV2$HGNC.symbol)),]
dim(t.first)

t.second=t.first[match(unique(as.character(t.first$MGI.symbol)), as.character(t.first$MGI.symbol)),]

final_overlap=t.second[which(t.second$HGNC.symbol%in%rownames(adr_all)),]

mouse_subs=subset(m13_srat, features=as.character(final_overlap$MGI.symbol))
adr_subs=subset(adr_all, features=as.character(final_overlap$HGNC.symbol))


adr_mtx=adr_subs@assays$RNA@counts


adr_right_order=adr_mtx[as.character(final_overlap$HGNC.symbol),]

rownames(adr_right_order)=as.character(final_overlap$MGI.symbol)


mouse_mtx=m13_srat@assays$RNA@counts
mouse_right_order=mouse_mtx[as.character(final_overlap$MGI.symbol),]


identical(rownames(mouse_right_order), rownames(adr_right_order))


fit_for_mouse = trainModel(adr_right_order, adr_all@meta.data$new_clust)

ps_mouse = predictSimilarity(fit_for_mouse, mouse_right_order,m13_srat@meta.data$new_clust)
similarityHeatmap(ps_mouse, row_order=c("scp", "bridge", "sympathoblast", "chromaffin"),
                  column_order=c("SCPs","Bridge","Sympathoblastic","Chromaffin","Mesenchyme","Cortex","Leukocytes", "Vascular endothelium", "Erythroblasts","Other"))


FeaturePlot(adr_med, c("SOX10", "ERBB3", "FOXD3","HTR3A","DLL3","THBD", "TH", "CHGB","FOXQ1", "CARTPT", "SLC18A3"), ncol=3) +
  theme(text = element_text(size=12), axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank(), legend.position = "none") +ylab("UMAP 2") +xlab("UMAP 1")


FeaturePlot(adr_med, c("CARTPT", "SOX10", "MPZ", "PLP1", "ERBB3", "DLL3", "TLX2", "TH","DBH", "PHOX2B", "CHGB", "GAP43", "BCL2", "NPY", "PNMT"), ncol=5)
FeaturePlot(m13_srat, c( "Sox2", "Sox10", "Mpz", "Plp1", "Erbb3", "Dll3", "Tlx2", "Th","Dbh", "Phox2b", "Chgb", "Gap43", "Bcl2", "Npy", "Pnmt"), ncol = 5)
