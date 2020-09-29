de_of_de_facet_tab=target_log2tpm[genes_for_kap, all$TARGET_SHORT]

de_tab=as.data.frame(t(de_of_de_facet_tab))  
de_tab$risk=all$risk  
de_tab$patient=rownames(de_tab)

de_tab_melted=melt(de_tab, id.vars=c("patient", "risk"))
de_tab_melted$risk=factor(de_tab_melted$risk, levels = c("low_risk_4s", "intermediate_risk","high_risk"))

de_tab_melted$variable=factor(de_tab_melted$variable, levels = c(rownames(de_of_de_sorted)))

de_tab_melted


unique(seqc_risks$risk)

seqc_risk2=seqc_risks

seqc_risk2$risk[which(seqc_risk2$risk=="low_risk")]="low_risk_4s"
unique(seqc_risk2$risk)


genes_for_seqc=genes_for_kap[which(genes_for_kap%in%rownames(seqc_tpm_log2))]


seqc_facet_tab=seqc_tpm_log2[genes_for_seqc, seqc_risk2$patient]
seqc_facet_tab=as.data.frame(t(seqc_facet_tab))

seqc_facet_tab$risk=seqc_risk2$risk
seqc_facet_tab$patient=rownames(seqc_facet_tab)

seqc_facet_tab_melted=melt(seqc_facet_tab, id.vars=c("patient", "risk"))

seqc_facet_tab_melted$risk=factor(seqc_facet_tab_melted$risk, levels = c("low_risk_4s", "intermediate_risk","high_risk"))

de_tab_melted$dataset="TARGET"
seqc_facet_tab_melted$dataset="SEQC"
de_tab_melted
seqc_facet_tab_melted


de_of_de_two_datasets=rbind(de_tab_melted, seqc_facet_tab_melted)


ggplot(data = de_of_de_two_datasets, aes(x=risk,y=value, group=dataset)) +
  geom_quasirandom(shape=1, alpha=0.3)+facet_wrap(~variable,scales="free_y")+
  geom_pointrange(mapping = aes(x = risk, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black",
                  position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))

de_of_de_two_datasets$risk
de_of_de_two_datasets$variable=as.character(de_of_de_two_datasets$variable)
de_of_de_two_datasets$risk=as.character(de_of_de_two_datasets$risk)


de_of_de_two_datasets$risk=factor(de_of_de_two_datasets$risk, levels=c("low_risk_4s", "intermediate_risk", "high_risk"))
de_of_de_two_datasets$variable=factor(de_of_de_two_datasets$variable, levels=genes_for_kap)
de_of_de_two_datasets$dataset=factor(de_of_de_two_datasets$dataset, levels=c("TARGET", "SEQC"))

which(de_of_de_two_datasets$variable%in%c("ALDH1A2", "PHF24", "CCDC144NL-AS1", "NCAN", "PRAME", "CRYBB2", "SIX3", "CPEB1", "TMEM163", "LSMEM1", "EPB41L4B", "RAP1GAP2", "SLC9A7"))

de_of_de_two_datasets_cons=de_of_de_two_datasets[which(de_of_de_two_datasets$variable%in%c("MAGEA4", "PRAME", "CRYBB2", "SIX3")),]
ggplot(data=de_of_de_two_datasets_cons, aes(x=risk, y=value, group=dataset)) +
 geom_quasirandom(shape=19, alpha=0.2, cex=0.5,dodge.width=.8)+facet_wrap(~variable,scales="free_y")+
  geom_pointrange(mapping = aes(x = risk, y = value, shape=dataset, fill=dataset),
                  stat = "summary", colour="#1c2e4a",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black",
                  position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))+ 
  labs(y="Expression (log2(TPM))")

de_of_de_two_datasets$da

de_of_de_two_datasets_no_inter=de_of_de_two_datasets[which(de_of_de_two_datasets$risk!="intermediate_risk"),]
de_of_de_two_datasets_no_inter=de_of_de_two_datasets_cons[which(de_of_de_two_datasets_cons$risk!="intermediate_risk"),]

de_of_de_two_datasets_cons=de_of_de_two_datasets[which(de_of_de_two_datasets$variable%in%c("ALDH1A2", "PHF24", "CCDC144NL-AS1", "NCAN", "PRAME", "CRYBB2", "SIX3", "CPEB1", "TMEM163", "LSMEM1", "EPB41L4B", "RAP1GAP2", "SLC9A7")),]
ggplot(data=de_of_de_two_datasets_no_inter, aes(x=risk, y=value, group=dataset)) +
  geom_quasirandom(shape=19, alpha=0.2, cex=0.5,dodge.width=.8)+facet_wrap(~variable,scales="free_y")+
  geom_pointrange(mapping = aes(x = risk, y = value, shape=dataset, fill=dataset),
                  stat = "summary", colour="#1c2e4a",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black",
                  position=position_dodge(width=0.8)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1))+ 
  labs(y="Expression (log2(TPM))")


FindClusters(srat_inhouse@assays$RNA@counts, srat_inhouse@meta.data$idents_for_plot)



unique(srat_inhouse@meta.data$new_idents)
