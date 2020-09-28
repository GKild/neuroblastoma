DimPlot(adr_all)

library(umap)
iris.umap = umap(iris.data)

adr_umap=umap(Embeddings(adr_all, reduction="pca"))

Embeddings(adr_all, reduction="pca")

adr_umap$layout=Embeddings(adr_all, reduction="umap")

plot(x=adr_umap$layout[,1], y=adr_umap$layout[,2])

same_features=intersect(rownames(adr_all@reductions$pca@feature.loadings), rownames(srat_inhouse))

inh_scaled_var_ft=srat_inhouse@assays$RNA@scale.data[same_features,]
adr_all_loadings=adr_all@reductions$pca@feature.loadings[same_features,]
transf_inh=t(inh_scaled_var_ft)%*%adr_all_loadings
inh_predict = predict(adr_umap, transf_inh)
inh_predict_df=as.data.frame(inh_predict)
inh_predict_df$lab=paste0("tum_",srat_inhouse@meta.data$idents_for_plot)

adr_orig=as.data.frame(adr_umap$layout)
adr_orig$lab=paste0("ref_",adr_all@meta.data$new_clust)
identical(rownames(srat_inhouse@meta.data), rownames(inh_predict))


orig_and_proj=rbind(adr_orig, inh_predict_df)


ggplot(data = orig_and_proj, aes(UMAP_1, UMAP_2, color=lab)) +geom_point()

############ Seurat integration thing ###############

nb.query <-srat_inhouse
nb.anchors <- FindTransferAnchors(reference = adr_all, query = nb.query, 
                                  dims = 1:30)
predictions <- TransferData(anchorset = nb.anchors, refdata = adr_all$new_clust, 
                            dims = 1:30)
nb.query <- AddMetaData(nb.query, metadata = predictions)

grep("Tumour",nb.query@meta.data$idents_for_plot)
table(nb.query@meta.data$predicted.id[grep("Tumour",nb.query@meta.data$idents_for_plot)])


nb.query_dutch <-srat_dutch
nb.anchors <- FindTransferAnchors(reference = adr_all, query = nb.query_dutch, 
                                  dims = 1:30)
predictions <- TransferData(anchorset = nb.anchors, refdata = adr_all$new_clust, 
                            dims = 1:30)
nb.query_dutch <- AddMetaData(nb.query_dutch, metadata = predictions)

grep("Tumour",nb.query_dutch@meta.data$idents_for_plot)
table(nb.query_dutch@meta.data$predicted.id[grep("Tumour",nb.query_dutch@meta.data$idents_for_plot)])

DimPlot(nb.query_dutch, group.by = "predicted.id")
nb.query_dutch@meta.data$prediction.score.Sympathoblastic
