library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
judes_top_1000 <- read.table("/home/jovyan/Neuroblastoma/useful_files/judes_top_1000.txt", sep = "\t", header = F) 
colnames(judes_top_1000) <- c("gene","nMuts", "nSamples", "missense", "nonsense","splice_region", "frameshift", "splice", "proteindel", "proteinins", "CN_loss", "CN_gain", "structural_var")
cosmic=read.csv("useful_files/cancer_gene_census.csv")
judes_and_cosmic=judes_top_1000[which(judes_top_1000$gene%in%cosmic$Gene.Symbol),]
#genes in NB figure
x=as.character(judes_and_cosmic$gene[1:20])
germline_genes <-c("PHOX2B","ALK", "KIF1B", "HRAS", "CDKN1C", "RP53", "SDHB", "BRCA2", "PALB2", "BRIP1",
                   "PTPN11", "SOS1", "RAF1", "KRAS", "NF1", "H19", "IGF2", "KNBQOT1", "EZH2", "SDHAF2",
                   "SDHC", "SDHD", "FANCA", "FANCC", "FANCG", "BRCA1")
somatic_and_germline=c("ALK","PTPN11","RAF1","KRAS","NF1")
new_som=x[which(!x%in%c(germline_genes, somatic_and_germline))]

new_germ=germline_genes[which(!germline_genes%in%c(new_som, somatic_and_germline))]
neurocrist_all=c("AAAS", "IKBKAP", "RET", "GDNF","EDN3", "TRK", "PHOX2A", "PHOX2B", "GFRA1", "BMP2", "ECE1", "VHL", "RET", "SDHB",
                   "SDHD")
new_neuro=neurocrist_all[which(!neurocrist_all%in%c(somatic_and_germline, new_germ, new_som))]

adr_subs=subset(adr_all, ident=c(25,24,13,19,20, 4,7,12, 30, 23,18,11,5,22,28,0,1,2,10,8,9,6,21))
dim(adr_subs)

adr_subs <- ScaleData(adr_subs, features = rownames(adr_subs))
levels(adr_subs)=c(25,24,13,19,20, 4,7,12, 30, 23,18,11,5,22,28,0,1,2,10,8,9,6,21)

avg <- AverageExpression(adr_subs, features = rownames(adr_subs), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)
adr_scaled <- avg@assays$RNA@scale.data


col_types <-factor(c("SCP trajectory", "SCP trajectory","SCP trajectory","SCP trajectory", "SCP trajectory","Blood","Blood",
                     "Cortex","Cortex","Cortex","Cortex","Cortex", "Cortex"),
                   levels=c("SCP trajectory", "Blood","Cortex"))
all_cols=c("SCPs", "Bridge", "Neuroblasts", "Ganglia", "Ganglia", "Blood","Blood","Blood","Blood", "Mesenchyme",
  "Mesenchyme","Mesenchyme",rep("Cortex",11))
col_types=factor(all_cols, levels = c("SCPs", "Bridge", "Neuroblasts", "Ganglia","Mesenchyme", 
                 "Cortex", "Blood"))
germline_gene_tab=data.frame(gene=new_germ, type="germline")
somatic_gene_tab=data.frame(gene=new_som, type="somatic")
som_and_germ_tab=data.frame(gene=somatic_and_germline, type="somatic and germline")
neurocr_tab=data.frame(gene=new_neuro, type="neurocristopathy")

all=rbind(somatic_gene_tab, germline_gene_tab, som_and_germ_tab, neurocr_tab)
all$gene=as.character(all$gene)
all=all[which(all$gene%in%rownames(adr_scaled)),]
Heatmap(adr_scaled[all$gene,], name="scaled expression", col=col_fun, 
        column_order = all$gene, row_names_side = "left", column_names_rot = 0,
        column_names_side = "top", 
        row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
        column_split = col_types, column_title_gp = gpar(fontsize=12))
col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), rev(brewer.pal(n = 7, name = "RdBu")))
col_fun = colorRamp2(c(-3,0,3), c("blue", "black", "red"))
Heatmap(adr_scaled[all$gene,], col=col_fun,row_order = all$gene, column_order = colnames(adr_scaled[all$gene,]), 
        column_split = col_types, show_column_names = F, row_split = all$type, row_gap = unit(4, "mm"), column_gap =unit(7, "mm"),
        heatmap_legend_param = list(direction="horizontal",heatmap_legend_side = "bottom"))
