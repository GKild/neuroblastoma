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

somatic_nb_genes <-top_mut[which(top_mut$gene%in%rownames(adr)),]

germline_genes <-c("PHOX2B","ALK", "KIF1B", "HRAS", "CDKN1C", "RP53", "SDHB", "BRCA2", "PALB2", "BRIP1",
                   "PTPN11", "SOS1", "RAF1", "KRAS", "NF1", "H19", "IGF2", "KNBQOT1", "EZH2", "SDHAF2",
                   "SDHC", "SDHD", "FANCA", "FANCC", "FANCG", "BRCA1")
gwas_genes <-c("TP53","CASC15", "NBAT1", "BARD1", "LMO1", "DUSP12", "DDX4", "IL31RA", "HSD17B12", "HACE1",
               "LIN28B", "TFCP2", "MLF1", "NEFL", "CDKN1B", "KIF15", "SPAG16")

adr_scaled <- avg@assays$RNA@scale.data

adr_scaled[somatic_nb_genes$gene,]
scaled_avg_nb_germline <-avg@assays$RNA@scale.data[which(rownames(avg@assays$RNA@scale.data)%in%c("PHOX2B","ALK", "KIF1B", "HRAS", "CDKN1C", "RP53", "SDHB", "BRCA2", "PALB2", "BRIP1", "PTPN11", "SOS1", "RAF1", "KRAS", "NF1", "H19", "IGF2", "KNBQOT1", "EZH2", "SDHAF2", "SDHC", "SDHD", "FANCA", "FANCC", "FANCG", "BRCA1")),]
gene_germ_list <-hclust(dist(scaled_avg_nb_germline))
gene_order_germline <-rownames(scaled_avg_nb_germline)[gene_germ_list$order]

scaled_avg_nb_gwas <-avg@assays$RNA@scale.data[which(rownames(avg@assays$RNA@scale.data)%in%c("TP53","CASC15", "NBAT1", "BARD1", "LMO1", "DUSP12", "DDX4", "IL31RA", "HSD17B12", "HACE1", "LIN28B", "TFCP2", "MLF1", "NEFL", "CDKN1B", "KIF15", "SPAG16")),]
gene_gwas_list <-hclust(dist(scaled_avg_nb_gwas))
gene_order_gwas <-rownames(scaled_avg_nb_gwas)[gene_gwas_list$order]
pdf("/home/jovyan/Neuroblastoma/heatmap_nb_genes.pdf", height = 15, width = 15)
gg=DoHeatmap(avg, features=gene_order_som, slot = "scale.data", draw.lines = F, raster = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 14)) + ggtitle("Average Expresion of somatically altered genes in NB in baby adrenal clusters")
plot(gg)
gg1=DoHeatmap(avg, features=gene_order_germline, slot = "scale.data", draw.lines = F, raster = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 14)) + ggtitle("Average Expresion of NB germline genes in baby adrenal clusters")
gg2=DoHeatmap(avg, features=gene_order_gwas, slot = "scale.data", draw.lines = F, raster = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 14)) + ggtitle("Average Expresion of genes w/ NB GWAS hits in baby adrenal clusters")
plot(plot_grid(gg1, gg2, ncol = 1, nrow = 2))
dev.off()


top100_mtx <-as.matrix(scaled_avg_nb_som)
judes_top_1000 <- read.table("/home/jovyan/Neuroblastoma/judes_top_1000.txt", sep = "\t", header = F) 
judes_top_100 <- judes_top_1000[1:100,]
colnames(judes_top_1000) <- c("gene","nMuts", "nSamples", "missense", "nonsense","splice_region", "frameshift", "splice", "proteindel", "proteinins", "CN_loss", "CN_gain", "structural_var")
top_mut <- data.frame(gene=judes_top_100$gene,top_mut=names(sapply(c(1:100), function(x){which.max(judes_top_100[x,4:13])})))
chr_1q <- read.table("/home/jovyan/Neuroblastoma/chr1q.txt", header = F, sep = "\t")
chr_2p <-read.table("/home/jovyan/Neuroblastoma/chr2p.txt", header = F, sep = "\t")
chr_7q <-read.table("/home/jovyan/Neuroblastoma/chr7q.txt", header = F, sep = "\t")
judes_top_1000[which(judes_top_1000$gene%in%chr_2p$V2),]

sub_mtx <- adr@assays$RNA@data[rownames(top100_mtx),]
adr@assays$RNA@

col_types <-factor(c("SCP trajectory", "SCP trajectory","SCP trajectory","SCP trajectory", "Proliferative cortex",
              "Proliferative cortex", "Cortex","Cortex","Cortex","Cortex","Cortex", "Blood","Blood",
              "V.E.", "CD45+", "Mesenchyme", "Mesenchyme","Mesenchyme") ,
              levels=c("SCP trajectory", "Proliferative cortex","Cortex", "Blood","V.E.", "CD45+", "Mesenchyme"))
col_fun = colorRamp2(c(-4, 0, 4), c("purple", "black", "yellow"))
col_fun(seq(-3, 3))
avg_exp_col=colorRamp2(c(0, 1), c("white", "blue")) 

ha = rowAnnotation(avg_exp=anno_simple(row_means, col=avg_exp_col))
ht <-ComplexHeatmap::Heatmap(top100_mtx, name="mat", col=col_fun, 
                        column_order = colnames(top100_mtx), row_names_side = "left", column_names_rot = 0,
                        column_names_side = "top", 
                        row_names_gp = gpar(fontsize = 6), show_row_dend = F, 
                        column_split = col_types, column_title_gp = gpar(fontsize=12), left_annotation = ha)

exp_legend = Legend(title = "avg exp", col = avg_exp_col, at = c(0,0.5, 1, 1.5, 2), 
                    labels = c("0", "0.5", "1", "1.5", "2"))
draw(ht, annotation_legend_list = exp_legend)

germline_ha=rowAnnotation(avg_exp=anno_simple(rowMeans(adr@assays$RNA@data[rownames(scaled_avg_nb_germline),]), col=avg_exp_col))
germline_ht =Heatmap(scaled_avg_nb_germline, column_order=colnames(scaled_avg_nb_germline), 
        col=col_fun,row_names_side = "left", column_names_rot = 0,
        column_names_side = "top", column_split = col_types, column_title_gp = gpar(fontsize=12),
        show_row_dend = F,left_annotation = germline_ha)
pdf("tr1.pdf", height = 15, width = 15)
draw(germline_ht, annotation_legend_list=exp_legend)
draw(ht,annotation_legend_list=exp_legend)
dev.off()


