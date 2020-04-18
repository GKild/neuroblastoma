#stratify by risk groups - target

low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))

low_risk_4s$risk="low_risk_4s"
high_risk_relapse$risk="high_risk"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

target_risks=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)

unique(target_risks$risk)

#stratify by risk group - german 

low_risk=dplyr::filter(german_clinical_data,
                       german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &!german_clinical_data$Status%in%c("relapse/progression","death of disease"))
high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified", german_clinical_data$Status%in%c("relapse/progression","death of disease"))
risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")
pt1=rbind(low_risk, high_risk, risk_4s)
intermediate_risk=german_clinical_data[-which(german_clinical_data$Tumor.ID%in%pt1$Tumor.ID),]

low_risk$risk="low_risk"
risk_4s$risk="low_risk_4s"
intermediate_risk$risk="intermediate_risk"
high_risk$risk="high_risk"

german_risks=rbind(low_risk, risk_4s, intermediate_risk, high_risk)

# stratify by risk group - seqc 


seqc_lr=dplyr::filter(seqc_meta, Age<540 & MYCN.status=="1", INSS.Stage!="4S")
seqc_4s=dplyr::filter(seqc_meta, INSS.Stage=="4S")
seqc_hr=dplyr::filter(seqc_meta, Age>540 & MYCN.status=="Amp" &INSS.Stage!="4S")
seqpt1=rbind(seqc_lr, seqc_4s, seqc_hr)
seqc_intermediate=seqc_meta[-which(seqc_meta$NB.ID%in%seqpt1$NB.ID),]

seqc_lr$risk="low_risk"
seqc_4s$risk="low_risk_4s"
seqc_hr$risk="high_risk"
seqc_intermediate$risk="intermediate_risk"

seqc_risks=rbind(seqc_lr, seqc_4s, seqc_hr, seqc_intermediate)
dim(german_risks)


german_risks$patient=german_risks$Tumor.ID
seqc_risks$patient=seqc_risks$NB.ID
target_risks$patient=target_risks$TARGET_SHORT


get_per_risk=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
      })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  scp_markers_in_rg=values_per_risk_group[scp_marks[which(scp_marks%in%rownames(values_per_risk_group))],]
  bridge_markers_in_rg=values_per_risk_group[bridge_marks[which(bridge_marks%in%rownames(values_per_risk_group))],]
  sympathoblastic_markers_in_rg=values_per_risk_group[left_marks[which(left_marks%in%rownames(values_per_risk_group))],]
  chromaffin_markers_in_rg=values_per_risk_group[right_marks[which(right_marks%in%rownames(values_per_risk_group))],]
  podo_markers_in_rg=values_per_risk_group[pod_marks[which(pod_marks%in%rownames(values_per_risk_group))],]
  
  scp_markers_in_rg=melt(scp_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  bridge_markers_in_rg=melt(bridge_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  sympathoblastic_markers_in_rg=melt(sympathoblastic_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  chromaffin_markers_in_rg=melt(chromaffin_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  podo_markers_in_rg=melt(podo_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  
  scp_markers_in_rg$marker="scp"
  bridge_markers_in_rg$marker="bridge"
  sympathoblastic_markers_in_rg$marker="sympathoblastic"
  chromaffin_markers_in_rg$marker="chromaffin"
  podo_markers_in_rg$marker="podocyte"
  
  all_in_groups=rbind(scp_markers_in_rg, bridge_markers_in_rg, sympathoblastic_markers_in_rg, chromaffin_markers_in_rg, podo_markers_in_rg)
  
  return(all_in_groups)
  
  }

target_per_risk=get_per_risk(target_log2tpm, target_risks, 5)
german_per_risk=get_per_risk(german_counts_log2, german_risks, 10)
seqc_per_risk=get_per_risk(seqc_tpm_log2, seqc_risks, 5)

target_per_risk$cohort="target"
german_per_risk$cohort="german"
seqc_per_risk$cohort="seqc"


combined_per_risk=rbind(target_per_risk, german_per_risk, seqc_per_risk)


colnames(combined_per_risk)=c("gene", "risk", "value", "marker", "cohort")

combined_per_risk$risk=factor(combined_per_risk$risk, levels=c("low_risk","low_risk_4s", "intermediate_risk", "high_risk", "high_risk_relapse"))

combined_per_risk$marker=factor(combined_per_risk$marker, levels=c("podocyte" , "scp","bridge","sympathoblastic","chromaffin"))

combined_per_risk$cohort=factor(combined_per_risk$cohort, levels = c("target", "seqc", "german"))



ggplot(combined_per_risk, aes(cohort, value)) + geom_boxplot(outlier.shape=NA) +
  geom_point() +facet_grid(vars(risk), vars(marker))

leftover_ids=as.character(target_clinical_data$ICDO.Description)[-grep("adrenal|abdominal|kidney|abdomen|unknown|retroperitoneum|other", tolower(as.character(target_clinical_data$ICDO.Description)))]
leftover_ids=leftover_ids[c(1:13,16:23)]

leftover_samples=target_clinical_data$TARGET_SHORT[which(target_clinical_data$ICDO.Description%in%leftover_ids)]



left_percent_above_threshold=apply(target_log2tpm[,leftover_samples], 1, function(y){
  sum(y>5)/length(y)*100
})





left_scp=as.data.frame(left_percent_above_threshold[scp_marks[which(scp_marks%in%names(left_percent_above_threshold))]])
left_bridge=as.data.frame(left_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(left_percent_above_threshold))]])
left_left=as.data.frame(left_percent_above_threshold[left_marks[which(left_marks%in%names(left_percent_above_threshold))]])
left_right=as.data.frame(left_percent_above_threshold[right_marks[which(right_marks%in%names(left_percent_above_threshold))]])
left_pods=as.data.frame(left_percent_above_threshold[pod_marks[which(pod_marks%in%names(left_percent_above_threshold))]])

left_scp$marker="scp"
left_bridge$marker="bridge"
left_left$marker="left"
left_right$marker="right"
left_pods$marker="pods"

left_scp$gene=rownames(left_scp)
left_bridge$gene=rownames(left_bridge)
left_left$gene=rownames(left_left)
left_right$gene=rownames(left_right)
left_pods$gene=rownames(left_pods)

rownames(left_scp)=1:nrow(left_scp)
rownames(left_bridge)=1:nrow(left_bridge)
rownames(left_left)=1:nrow(left_left)
rownames(left_right)=1:nrow(left_right)
rownames(left_pods)=1:nrow(left_pods)

colnames(left_scp)=c("value", "marker", "gene")
colnames(left_bridge)=c("value", "marker", "gene")
colnames(left_left)=c("value", "marker", "gene")
colnames(left_right)=c("value", "marker","gene")
colnames(left_pods)=c("value", "marker","gene")

left_all_m=rbind(left_scp, left_bridge,left_left,left_right, left_pods)

ggplot(left_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >5 log2(tpm) in TARGET non-medullary samples")


get_risk_de_of_de=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
    })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  de_markers_in_rg=values_per_risk_group[de_marks[which(de_marks%in%rownames(values_per_risk_group))],]
  
  de_markers_in_rg=melt(de_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
 
  
  
  return(de_markers_in_rg)
  
}


rsk1=get_risk_de_of_de(target_log2tpm, target_risks, 5)


ggplot(rsk1, aes(x=variable, y=value)) + 
  geom_boxplot() +ggtitle("% samples a DE gene is expressed at >5 log2(tpm) in TARGET")

unique(rsk1[rsk1$value>50,]$gene)


de_marks


de_of_de_sorted=de_of_de[order(abs(de_of_de$logFC), decreasing = T),]

which(rownames(de_of_de_sorted)%in%unique(rsk1[rsk1$value>50,]$gene))


de_of_de_facet_tab=target_log2tpm[rownames(de_of_de_sorted), all$TARGET_SHORT]

de_tab=as.data.frame(t(de_of_de_facet_tab))  
de_tab$risk=all$risk  
de_tab$patient=rownames(de_tab)

de_tab_melted=melt(de_tab, id.vars=c("patient", "risk"))
de_tab_melted$risk=factor(de_tab_melted$risk, levels = c("low_risk_4s", "intermediate_risk","high_risk"))

de_tab_melted$variable=factor(de_tab_melted$variable, levels = c(rownames(de_of_de_sorted)))


ggplot(de_tab_melted, aes(x=risk, y=value)) + geom_boxplot()+ ggtitle("DE markers")+facet_wrap(~variable,scales="free_y") + theme(axis.text.x = element_text(angle = 90, hjust = 1))




mapk1=read.table("mapk1.tsv", sep = "\t", header = T)
mapk2=read.table("mapk2.tsv", sep = "\t", header = T)


mapk_cancer_1=read.table("mapk_cancer1.tsv", sep = "\t", header = T)
mapk_cancer_2=read.table("mapk_cancer2.tsv", sep = "\t", header = T)
mapk_cancer_3=read.table("mapk_cancer3.tsv", sep = "\t", header = T)

nb_mapk=read.table("mapk_in_nb.txt", sep = "\t", header = F)

nb_mapk_genes=unique(as.character(nb_mapk$V1))

general_mapk_genes=unique(c(as.character(mapk1$Gene.Name), as.character(mapk2$Gene.Name)))

cancer_mapk_genes=unique(c(as.character(mapk_cancer_1$Gene.Name), as.character(mapk_cancer_2$Gene.Name),as.character(mapk_cancer_3$Gene.Name)))

get_mapk_per_risk=function(count_tab, risk_tab,th){
  values_per_risk_group=sapply(unique(risk_tab$risk), function(x){
    pts=risk_tab$patient[which(risk_tab$risk==x)]
    tab=count_tab[,pts]
    apply(tab, 1, function(y){
      sum(y>th)/length(y)*100
    })
  })
  
  values_per_risk_group=data.frame(values_per_risk_group)
  values_per_risk_group$gene=rownames(values_per_risk_group)
  
  #scp_markers_in_rg=values_per_risk_group[general_mapk_genes[which(general_mapk_genes%in%rownames(values_per_risk_group))],]
  bridge_markers_in_rg=values_per_risk_group[mapk_sc[which(mapk_sc%in%rownames(values_per_risk_group))],]
  
  
  #scp_markers_in_rg=melt(scp_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  bridge_markers_in_rg=melt(bridge_markers_in_rg, id.vars = "gene", measure.vars = as.character(unique(risk_tab$risk)))
  
  #scp_markers_in_rg$marker="general_mapk"
  bridge_markers_in_rg$marker="cancer_mapk"
  
  #all_in_groups=rbind(scp_markers_in_rg, bridge_markers_in_rg)
  return(bridge_markers_in_rg)
  
}

target_mapk_per_risk=get_mapk_per_risk(target_log2tpm, target_mapk_risks,5)
german_mapk_per_risk=get_mapk_per_risk(german_counts_log2, german_risks, 10)
seqc_mapk_per_risk=get_mapk_per_risk(seqc_tpm_log2, seqc_risks, 5)

target_mapk_per_risk$cohort="target"
german_mapk_per_risk$cohort="german"
seqc_mapk_per_risk$cohort="seqc"


combined_mapk_per_risk=rbind(target_mapk_per_risk, german_mapk_per_risk, seqc_mapk_per_risk)

colnames(combined_mapk_per_risk)=c("gene", "risk", "value", "marker", "cohort")

combined_mapk_per_risk$risk=factor(combined_mapk_per_risk$risk, levels=c("low_risk","low_risk_4s", "intermediate_risk", "high_risk", "high_risk_relapse"))


combined_mapk_per_risk$cohort=factor(combined_mapk_per_risk$cohort, levels = c("target", "seqc", "german"))

ggplot(data = combined_mapk_per_risk, aes(x=marker,y=value)) +
  geom_quasirandom(shape=1) +facet_grid(vars(risk), vars(cohort)) + theme( axis.line = element_line(colour = "black"))+
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black")









ggplot(combined_per_risk, aes(cohort, value)) + geom_boxplot(outlier.shape=NA) +
  geom_point() +facet_grid(vars(risk), vars(marker))

ggplot(data = combined_per_risk, aes(x=marker,y=value)) +
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +
  geom_quasirandom(shape=1)+facet_grid(vars(risk), vars(cohort)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
ggplot(joined_dat, aes(x=marker, y=value, fill=dataset)) + 
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge(jitter.width = 0.2),aes(group=dataset), size=0.75) +ggtitle("% samples a marker is expressed above a threshold across 3 cohorts (N=1050)")


ggplot(data = joined_dat, aes(x=marker,y=value, group=dataset)) +
  geom_quasirandom(mapping = aes(x=marker,y=value, group=dataset), dodge.width=.8, shape=1, alpha=0.3) + theme( axis.line = element_line(colour = "black"))+
  geom_pointrange(mapping = aes(x = marker, y = value, shape=dataset),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, fill="black", 
                  position=position_dodge(width=0.8)) +  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

ggplot(data = combined_mapk_per_risk, aes(x=marker,y=value)) +
  geom_quasirandom(shape=1) +facet_grid(vars(risk), vars(cohort)) + theme( axis.line = element_line(colour = "black"))+
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black")



ggplot(data = left_all_m, aes(x=marker,y=value)) +geom_quasirandom(shape=1) +
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22, fill="black") +ggtitle("% samples a marker is expressed at >5 log2(tpm) in TARGET non-medullary samples")





