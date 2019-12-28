library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
judes_top_1000 <- read.table("/home/jovyan/Neuroblastoma/judes_top_1000.txt", sep = "\t", header = F) 
colnames(judes_top_1000) <- c("gene","nMuts", "nSamples", "missense", "nonsense","splice_region", "frameshift", "splice", "proteindel", "proteinins", "CN_loss", "CN_gain", "structural_var")
judes_top_100 <- judes_top_1000[1:100,]
top_mut <- data.frame(gene=judes_top_100$gene,top_mut=names(sapply(c(1:100), function(x){which.max(judes_top_100[x,4:13])})))
extra_rows <- data.frame(gene=c("TIAM1", "ARID1A", "ARID1B", "OR5T1", "PDE6G", "VANGL1", "NRAS"),
                         top_mut=c("missense","missense","structural_var", "missense","missense","missense","missense"))
top_mut <- rbind(top_mut, extra_rows)

adr <- readRDS("/home/jovyan/Neuroblastoma/babyAdrenalRef/babyAdrenal2.RDS")
adr <- NormalizeData(adr)
adr <- ScaleData(adr, features = rownames(adr))
levels(adr) <- c("14", "16", "10", "7","0","9","11","3","1","2", "6", "5","15", "12", "17","4","8","13")
avg <- AverageExpression(adr, features = rownames(adr), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)

adr_scaled <- avg@assays$RNA@scale.data

somatic_nb_genes <-top_mut[which(top_mut$gene%in%rownames(adr)),]

germline_genes <-c("PHOX2B","ALK", "KIF1B", "HRAS", "CDKN1C", "RP53", "SDHB", "BRCA2", "PALB2", "BRIP1",
                   "PTPN11", "SOS1", "RAF1", "KRAS", "NF1", "H19", "IGF2", "KNBQOT1", "EZH2", "SDHAF2",
                   "SDHC", "SDHD", "FANCA", "FANCC", "FANCG", "BRCA1")
germline_genes <- germline_genes[which(germline_genes%in%rownames(adr))]
gwas_genes <-c("TP53","CASC15", "NBAT1", "BARD1", "LMO1", "DUSP12", "DDX4", "IL31RA", "HSD17B12", "HACE1",
               "LIN28B", "TFCP2", "MLF1", "NEFL", "CDKN1B", "KIF15", "SPAG16")
gwas_genes <- gwas_genes[which(gwas_genes%in%rownames(adr))]
col_types <-factor(c("SCP trajectory", "SCP trajectory","SCP trajectory","SCP trajectory", "Proliferative cortex",
                     "Proliferative cortex", "Cortex","Cortex","Cortex","Cortex","Cortex", "Blood","Blood",
                     "V.E.", "CD45+", "Mesenchyme", "Mesenchyme","Mesenchyme") ,
                   levels=c("SCP trajectory", "Proliferative cortex","Cortex", "Blood","V.E.", "CD45+", "Mesenchyme"))

col_fun = colorRamp2(c(-4, 0, 4), c("purple", "black", "yellow"))
col_fun(seq(-3, 3))
avg_exp_col=colorRamp2(c(0, 1), c("white", "blue")) 

ha = rowAnnotation(avg_exp=anno_simple(rowMeans(adr@assays$RNA@data[somatic_nb_genes$gene,]), col=avg_exp_col))
ht =ComplexHeatmap::Heatmap(adr_scaled[somatic_nb_genes$gene,], name="mat", col=col_fun, 
                             column_order = colnames(adr_scaled[somatic_nb_genes$gene,]), row_names_side = "left", column_names_rot = 0,
                             column_names_side = "top", 
                             row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                             column_split = col_types, column_title_gp = gpar(fontsize=12), left_annotation = ha, row_split = somatic_nb_genes$top_mut,
                            row_names_max_width = max_text_width(
                              rownames(adr@assays$RNA@data[somatic_nb_genes$gene,]), 
                              gp = gpar(fontsize = 10)))

exp_legend = Legend(title = "avg exp", col = avg_exp_col, at = c(0,0.5, 1), 
                    labels = c("0", "0.5", "1"))
draw(ht, annotation_legend_list = exp_legend)
ha_ra = rowAnnotation(avg_exp=anno_simple(rowMeans(adr@assays$RNA@data[c("RARA", "RARB","RARG"),]), col=avg_exp_col))
hm_ra =ComplexHeatmap::Heatmap(adr_scaled[c("RARA", "RARB","RARG"),], name="scaled expression", col=col_fun, 
                        column_order = colnames(adr_scaled[c("RARA", "RARB","RARG"),]), row_names_side = "left", column_names_rot = 0,
                        column_names_side = "top", 
                        row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                        column_split = col_types, column_title_gp = gpar(fontsize=12), left_annotation = ha_ra,
                        row_names_max_width = max_text_width(
                          rownames(adr@assays$RNA@data[somatic_nb_genes$gene,]), 
                          gp = gpar(fontsize = 10)))
draw(hm_ra, annotation_legend_list=exp_legend)

ha_germ = rowAnnotation(avg_exp=anno_simple(rowMeans(adr@assays$RNA@data[germline_genes,]), col=avg_exp_col))
hm_germ =ComplexHeatmap::Heatmap(adr_scaled[germline_genes,], name="scaled expression", col=col_fun, 
                               column_order = colnames(adr_scaled[germline_genes,]), row_names_side = "left", column_names_rot = 0,
                               column_names_side = "top", 
                               row_names_gp = gpar(fontsize = 10), show_row_dend = F, 
                               column_split = col_types, column_title_gp = gpar(fontsize=12), left_annotation = ha_germ,
                               row_names_max_width = max_text_width(
                                 rownames(adr@assays$RNA@data[somatic_nb_genes$gene,]), 
                                 gp = gpar(fontsize = 10)))
pdf("exp_in_ref.pdf", height = 14, width = 14)
draw(ht, annotation_legend_list = exp_legend)
draw(hm_germ, annotation_legend_list=exp_legend)
draw(hm_ra, annotation_legend_list=exp_legend)
dev.off()


