somatic_alts_nb <- read.table("somatic_alts_nb.txt", sep = "\t", header = F)          
heatmap_genes <- as.character(somatic_alts_nb$V1)
avg <- AverageExpression(adr, features = heatmap_genes_2, add.ident = NULL, return.seurat = TRUE, verbose = TRUE)
pdf("avg_exp_som_nb.pdf", width = 14, height = 21)
gg=DoHeatmap(avg, features=heatmap_genes[1:10], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[11:20], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[21:30], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[31:40], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[41:50], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[51:60], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[61:70], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[71:80], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[81:90], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
gg=DoHeatmap(avg, features=heatmap_genes[91:100], slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
plot(gg)
dev.off()


gg=DoHeatmap(adr, features=heatmap_genes_2, slot = "scale.data", draw.lines = F) +
  scale_fill_viridis(option = "viridis") +
  theme(text = element_text(size = 20))
adr=ScaleData(adr, features = rownames(adr@assays$RNA@counts))

which(adr$seurat_clusters==0)
scale_rowm <-rowMeans(adr@assays$RNA@scale.data[,which(adr$seurat_clusters==0)])
which(names(scale_rowm)=="MYCN")
scale_rowm[3304]
avg@assays$RNA@scale.data

avg_2 <- AverageExpression(adr, add.ident = NULL, return.seurat = TRUE, verbose = TRUE)
avg_2@assays$RNA@scale.data["MYCN",]
rm(adr)
adr <- ScaleData(adr, )
