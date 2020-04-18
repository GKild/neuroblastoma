
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
judes_top_1000 <- read.table("/home/jovyan/Neuroblastoma/useful_files/judes_top_1000.txt", sep = "\t", header = F) 
colnames(judes_top_1000) <- c("gene","nMuts", "nSamples", "missense", "nonsense","splice_region", "frameshift", "splice", "proteindel", "proteinins", "CN_loss", "CN_gain", "structural_var")
judes_top_100 <- judes_top_1000[1:100,]
top_mut <- data.frame(gene=judes_top_100$gene,top_mut=names(sapply(c(1:100), function(x){which.max(judes_top_100[x,4:13])})))
extra_rows <- data.frame(gene=c("TIAM1", "ARID1A", "ARID1B", "OR5T1", "PDE6G", "VANGL1", "NRAS"),
                         top_mut=c("missense","missense","structural_var", "missense","missense","missense","missense"))
top_mut <- rbind(top_mut, extra_rows)
germline_genes <-c("PHOX2B","ALK", "KIF1B", "HRAS", "CDKN1C", "RP53", "SDHB", "BRCA2", "PALB2", "BRIP1",
                   "PTPN11", "SOS1", "RAF1", "KRAS", "NF1", "H19", "IGF2", "KNBQOT1", "EZH2", "SDHAF2",
                   "SDHC", "SDHD", "FANCA", "FANCC", "FANCG", "BRCA1")
somatic_and_germline=c("ALK","PTPN11","RAF1","KRAS","NF1")
neurocrist_all=c("AAAS", "IKBKAP", "RET", "GDNF","EDN3", "TRK", "PHOX2A", "PHOX2B", "GFRA1", "BMP2", "ECE1", "VHL", "RET", "SDHB",
                 "SDHD")
new_som=top_mut[which(!top_mut$gene%in%c(germline_genes, somatic_and_germline)),]

new_germ=germline_genes[which(!germline_genes%in%c(new_som$gene, somatic_and_germline))]

new_neuro=neurocrist_all[which(!neurocrist_all%in%c(somatic_and_germline, new_germ, new_som$gene))]

germline_gene_tab=data.frame(gene=new_germ, type="germline")
som_and_germ_tab=data.frame(gene=somatic_and_germline, type="somatic and germline")
neurocr_tab=data.frame(gene=new_neuro, type="neurocristopathy")

all=rbind(germline_gene_tab, som_and_germ_tab, neurocr_tab)
all$gene=as.character(all$gene)
all=all[which(all$gene%in%rownames(adr_all)),]
all_som=new_som[which(new_som$gene%in%rownames(adr_all)),]

all_som$gene=as.character(all_som$gene)

avg <- AverageExpression(adr_all, features = rownames(adr_all), add.ident = NULL, return.seurat = TRUE, verbose = TRUE)
adr_scaled <- avg@assays$RNA@scale.data
new_levels=c("25","24","19","20","13", "11", "18", "23", 
             "5","22","28","0","1", "2", "10", "8", "6", "9", "21",
             "27","30","7","4","12", "31", "29", "26", "16", "3", "15", "14", "17")
ComplexHeatmap::Heatmap(adr_scaled[all$gene,], name="mat", col=col_fun, 
                        column_order = new_levels, row_names_side = "left", column_names_rot = 0,
                        column_names_side = "top", 
                        row_names_gp = gpar(fontsize = 10), show_row_dend = F, , column_title_gp = gpar(fontsize=12), row_split = all$type)

colnames(all_som)=c("gene", "type")

ComplexHeatmap::Heatmap(adr_scaled[all_som$gene,], name="mat", col=col_fun, 
                        column_order = new_levels, row_names_side = "left", column_names_rot = 0,
                        column_names_side = "top", 
                        row_names_gp = gpar(fontsize = 10), show_row_dend = F, , column_title_gp = gpar(fontsize=12), row_split = all_som$type)



blah=rbind(all_som, all)


rel_tab=adr_scaled[blah$gene, new_levels]






x=sapply(new_levels, function(x){
clust=WhichCells(adr_all, idents = x)
apply(adr_all@assays$RNA@counts[,clust],1,function(y){
  sum(y>0)/length(y)*100
})
})

avg_exp=sapply(new_levels, function(x){
  clust=WhichCells(adr_all, idents = x)
  apply(adr_all@assays$RNA@counts[,clust],1,function(y){
    mean(y)
  })
})
blah1=adr_scaled[,new_levels]
write.table(x, "percent_cells_exp_per_clust_adr_universe.txt", sep = "\t", quote=F, col.names = T, row.names = T)
write.table(blah1, "all_genes_scaled_avg_exp_adr_universe.txt", sep = "\t", quote=F, col.names = T, row.names = T)
write.table(avg_exp, "all_genes_avg_exp_adr_universe.txt", sep = "\t", quote=F, col.names = T, row.names = T)
write.table(blah, "gene_annot.txt", col.names = T, row.names = F, sep = "\t", quote = F)

