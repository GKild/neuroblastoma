scp_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="SCPs"),]$gene)
bridge_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="bridge"),]$gene)
right_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="right_clust"),]$gene)
left_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="left_clust"),]$gene)
pod_marks=as.character(podo_marks_filt[which(podo_marks_filt$cluster=="pods"),]$gene)

t_scp=as.data.frame(percent_above_threshold[scp_marks[which(scp_marks%in%names(percent_above_threshold))]])
t_bridge=as.data.frame(percent_above_threshold[bridge_marks[which(bridge_marks%in%names(percent_above_threshold))]])
t_left=as.data.frame(percent_above_threshold[left_marks[which(left_marks%in%names(percent_above_threshold))]])
t_right=as.data.frame(percent_above_threshold[right_marks[which(right_marks%in%names(percent_above_threshold))]])
t_pods=as.data.frame(percent_above_threshold[pod_marks[which(pod_marks%in%names(percent_above_threshold))]])

t_scp$marker="scp"
t_bridge$marker="bridge"
t_left$marker="left"
t_right$marker="right"
t_pods$marker="pods"

t_scp$gene=rownames(t_scp)
t_bridge$gene=rownames(t_bridge)
t_left$gene=rownames(t_left)
t_right$gene=rownames(t_right)
t_pods$gene=rownames(t_pods)

rownames(t_scp)=1:nrow(t_scp)
rownames(t_bridge)=1:nrow(t_bridge)
rownames(t_left)=1:nrow(t_left)
rownames(t_right)=1:nrow(t_right)
rownames(t_pods)=1:nrow(t_pods)

colnames(t_scp)=c("value", "marker", "gene")
colnames(t_bridge)=c("value", "marker", "gene")
colnames(t_left)=c("value", "marker", "gene")
colnames(t_right)=c("value", "marker","gene")
colnames(t_pods)=c("value", "marker","gene")

t_all_m=rbind(t_scp, t_bridge,t_left,t_right, t_pods)


ggplot(t_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >5 log2(TPM) in TARGET cohort (N=154)")


seqc_scp=as.data.frame(seqc_percent_above_threshold[scp_marks[which(scp_marks%in%names(seqc_percent_above_threshold))]])
seqc_bridge=as.data.frame(seqc_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(seqc_percent_above_threshold))]])
seqc_left=as.data.frame(seqc_percent_above_threshold[left_marks[which(left_marks%in%names(seqc_percent_above_threshold))]])
seqc_right=as.data.frame(seqc_percent_above_threshold[right_marks[which(right_marks%in%names(seqc_percent_above_threshold))]])
seqc_pods=as.data.frame(seqc_percent_above_threshold[pod_marks[which(pod_marks%in%names(seqc_percent_above_threshold))]])

seqc_scp$marker="scp"
seqc_bridge$marker="bridge"
seqc_left$marker="left"
seqc_right$marker="right"
seqc_pods$marker="pods"

seqc_scp$gene=rownames(seqc_scp)
seqc_bridge$gene=rownames(seqc_bridge)
seqc_left$gene=rownames(seqc_left)
seqc_right$gene=rownames(seqc_right)
seqc_pods$gene=rownames(seqc_pods)

rownames(seqc_scp)=1:nrow(seqc_scp)
rownames(seqc_bridge)=1:nrow(seqc_bridge)
rownames(seqc_left)=1:nrow(seqc_left)
rownames(seqc_right)=1:nrow(seqc_right)
rownames(seqc_pods)=1:nrow(seqc_pods)

colnames(seqc_scp)=c("value", "marker", "gene")
colnames(seqc_bridge)=c("value", "marker", "gene")
colnames(seqc_left)=c("value", "marker", "gene")
colnames(seqc_right)=c("value", "marker","gene")
colnames(seqc_pods)=c("value", "marker","gene")

seqc_all_m=rbind(seqc_scp, seqc_bridge,seqc_left,seqc_right, seqc_pods)


ggplot(seqc_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() + ggtitle("% samples a marker is expressed at >5 log2(TPM) in SEQC cohort (N=498)")




german_scp=as.data.frame(german_percent_above_threshold[scp_marks[which(scp_marks%in%names(german_percent_above_threshold))]])
german_bridge=as.data.frame(german_percent_above_threshold[bridge_marks[which(bridge_marks%in%names(german_percent_above_threshold))]])
german_left=as.data.frame(german_percent_above_threshold[left_marks[which(left_marks%in%names(german_percent_above_threshold))]])
german_right=as.data.frame(german_percent_above_threshold[right_marks[which(right_marks%in%names(german_percent_above_threshold))]])
german_pods=as.data.frame(german_percent_above_threshold[pod_marks[which(pod_marks%in%names(german_percent_above_threshold))]])

german_scp$marker="scp"
german_bridge$marker="bridge"
german_left$marker="left"
german_right$marker="right"
german_pods$marker="pods"

german_scp$gene=rownames(german_scp)
german_bridge$gene=rownames(german_bridge)
german_left$gene=rownames(german_left)
german_right$gene=rownames(german_right)
german_pods$gene=rownames(german_pods)

rownames(german_scp)=1:nrow(german_scp)
rownames(german_bridge)=1:nrow(german_bridge)
rownames(german_left)=1:nrow(german_left)
rownames(german_right)=1:nrow(german_right)
rownames(german_pods)=1:nrow(german_pods)

colnames(german_scp)=c("value", "marker", "gene")
colnames(german_bridge)=c("value", "marker", "gene")
colnames(german_left)=c("value", "marker", "gene")
colnames(german_right)=c("value", "marker","gene")
colnames(german_pods)=c("value", "marker","gene")

german_all_m=rbind(german_scp, german_bridge,german_left,german_right, german_pods)


ggplot(german_all_m, aes(x=marker, y=value)) + 
  geom_boxplot() +ggtitle("% samples a marker is expressed at >10 log2(intensity) in GERMAN cohort (N=394)")

which(rownames(german_counts_log2)%in%"ERBB3.4")


t_all_m$dataset="TARGET"
german_all_m$dataset="GERMAN"
seqc_all_m$dataset="SEQC"


joined_dat=rbind(t_all_m, german_all_m, seqc_all_m)


joined_dat$marker=factor(joined_dat$marker, levels = c("scp", "bridge", "left", "right", "pods"))
joined_dat$dataset=factor(joined_dat$dataset, levels = c("TARGET", "SEQC", "GERMAN"))



ggplot(joined_dat, aes(x=marker, y=value, fill=dataset)) + 
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge(jitter.width = 0.2),aes(group=dataset), size=0.75) +ggtitle("% samples a marker is expressed above a threshold across 3 cohorts (N=1050)")



leftover_ids=as.character(target_clinical_data$ICDO.Description)[-grep("adrenal|abdominal|kidney|mediastinum|pelvis|abdomen", tolower(as.character(target_clinical_data$ICDO.Description)))]


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



write.table(podo_marks_filt, "medullary_and_podocyte_cluster_markers.txt", sep = "\t", quote = F, row.names = F)


ggplot(data = left_all_m, aes(x=marker,y=value)) +
  geom_pointrange(mapping = aes(x = marker, y = value),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, shape=22,fatten = 10, size = 0.5, fill="black") +
  geom_jitter(width = 0.1, shape=1)

