#filter the MAPK genes based on single cell expression
inhouse_tum_cells=rownames(srat_inhouse@meta.data)[which(srat_inhouse@meta.data$new_idents=="tumour")]

inhouse_tum_cellexp=apply(srat_inhouse@assays$RNA@counts[,inhouse_tum_cells],1,function(y){
  sum(y>0)/length(y)*100
})

dutch_tum_cells=rownames(srat_tumour@meta.data)[which(srat_tumour@meta.data$new_idents=="tumour")]
dutch_tum_cellexp=apply(srat_tumour@assays$RNA@counts[,dutch_tum_cells],1,function(y){
  sum(y>0)/length(y)*100
})

inhouse_immune_cells=rownames(srat_inhouse@meta.data)[which(srat_inhouse@meta.data$new_idents=="immune")]
inhouse_immune_cellexp=apply(srat_inhouse@assays$RNA@counts[,inhouse_immune_cells],1,function(y){
  sum(y>0)/length(y)*100
})

dutch_immune_cells=rownames(srat_tumour@meta.data)[which(srat_tumour@meta.data$new_idents=="immune")]
dutch_immune_cellexp=apply(srat_tumour@assays$RNA@counts[,dutch_immune_cells],1,function(y){
  sum(y>0)/length(y)*100
})
mapk_inhouse_t50=general_mapk_genes[which(general_mapk_genes%in%names(inhouse_tum_cellexp[inhouse_tum_cellexp>50]))]

mapk_dutch_t50=general_mapk_genes[which(general_mapk_genes%in%names(dutch_tum_cellexp[dutch_tum_cellexp>50]))]

tum_de_genes=as.character(tum_n_de_final_sorted$gene)

immune_inhouse_t50=tum_de_genes[which(tum_de_genes%in%names(inhouse_immune_cellexp[inhouse_immune_cellexp>25]))]
dutch_immune_t50=tum_de_genes[which(tum_de_genes%in%names(dutch_immune_cellexp[dutch_immune_cellexp>25]))]

immune_sc=unique(c(immune_inhouse_t50, dutch_immune_t50))



length(which(mapk_dutch_t50%in%mapk_inhouse_t50))

length(mapk_inhouse_t50)
length(mapk_dutch_t50)

mapk_sc=unique(c(mapk_inhouse_t50, mapk_dutch_t50))
length(mapk_sc)



de_of_de_pat=as.data.frame(percent_above_threshold[rownames(de_of_de_sorted)])
de_of_de_pat$gene=rownames(de_of_de_pat)

colnames(de_of_de_pat)=c("value", "gene")
de_of_de_pat=de_of_de_pat[order(de_of_de_pat$value, decreasing = T),]

de_of_de_pat$gene=factor(de_of_de_pat$gene, levels = de_of_de_pat$gene)
a=ifelse(de_of_de_pat$gene %in% immune_sc, "grey", "red")

ggplot(de_of_de_pat, aes(x=reorder(gene, value), y=value))+geom_bar(stat="identity") +coord_flip() +
  theme(axis.text.y = element_text(colour = rev(a))) +scale_x_discrete(limits = rev(levels(de_of_de_pat$gene)))

above_log=de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),][which(abs(de_of_de_sorted[which(!rownames(de_of_de_sorted)%in%immune_sc),]$logFC)>1),]

genes_for_kap=rownames(above_log)


######kaplan-meier curves#########

start_mtx=target_tpm[,all$TARGET_SHORT]

kap_meta=data.frame(patient=all$TARGET_SHORT, event_free_surv=all$Event.Free.Survival.Time.in.Days, vital_status=all$Vital.Status)


kap_meta$surv_num=0
kap_meta$surv_num[which(kap_meta$vital_status=="Alive")]=1

gene_df=as.data.frame(t(start_mtx[genes_for_kap,]))

kap_final=cbind(kap_meta, gene_df)
colnames(kap_final)[8]="CCDC144NLAS1"
sap_g=colnames(kap_final)[5:24]
sapply(sap_g, function(x){
BRCA_HNSC.surv_rnaseq.cut <- surv_cutpoint(
  kap_final,
  time = "event_free_surv",
  event = "surv_num",
  variables = c(x)
)
summary(BRCA_HNSC.surv_rnaseq.cut)
pdf(paste0("kaplan_meier/",x, "_dist.pdf"), height = 5, width = 7, onefile = F)
print(plot(BRCA_HNSC.surv_rnaseq.cut, x, palette = "npg"))
dev.off()
BRCA_HNSC.surv_rnaseq.cat <- surv_categorize(BRCA_HNSC.surv_rnaseq.cut) 

pdf(paste0("kaplan_meier/",x, "_km.pdf"), height = 5, width = 7, onefile = F)
print(RTCGA::kmTCGA(BRCA_HNSC.surv_rnaseq.cat, times="event_free_surv", status="surv_num", explanatory.names = x, pval = T))
dev.off()
})


length(no_org_final_sorted$no_org_lj.gene[which(!no_org_final_sorted$no_org_lj.gene%in%immune_sc)])
