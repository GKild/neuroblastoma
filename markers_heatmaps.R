#get_corr

get_target_corr=function(exp_dat, scp_markers, neurobl_markers, mature_markers){
  #relevant ensembl IDs for the genes of interest
  scp_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%scp_markers),]$ensembl_gene_id)
  neurobl_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%neurobl_markers),]$ensembl_gene_id)
  mature_ensembl=as.character(target_genes[which(target_genes$external_gene_name%in%mature_markers),]$ensembl_gene_id)
  
  low_risk_4s=dplyr::filter(target_clinical_data, target_clinical_data$COG.Risk.Group=="Low Risk")
  high_risk_relapse=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & First.Event%in%c("Relapse", "Death"))
  intermediate_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="Intermediate Risk")
  high_risk=dplyr::filter(target_clinical_data, COG.Risk.Group=="High Risk" & !First.Event%in%c("Relapse", "Death"))
  
  low_risk_4s$risk="low_risk_4s"
  high_risk_relapse$risk="high_risk_relapse"
  intermediate_risk$risk="intermediate_risk"
  high_risk$risk="high_risk"
  
  extremes=rbind(low_risk_4s, high_risk_relapse)
  all=rbind(low_risk_4s,intermediate_risk, high_risk, high_risk_relapse)
  
  scp_tab=exp_dat[scp_ensembl, all$TARGET.USI]
  neurobl_tab=exp_dat[neurobl_ensembl, all$TARGET.USI]
  mature_tab=exp_dat[mature_ensembl, all$TARGET.USI]
  
  
  rownames(scp_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(scp_tab))]
  rownames(neurobl_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(neurobl_tab))]
  rownames(mature_tab)=target_genes$external_gene_name[which(target_genes$ensembl_gene_id%in%rownames(mature_tab))]
  
  scp_tab=as.data.frame(t(scp_tab))
  neurobl_tab=as.data.frame(t(neurobl_tab))
  mature_tab=as.data.frame(t(mature_tab))
  
  scp_tab$risk=all$risk
  neurobl_tab$risk=all$risk
  mature_tab$risk=all$risk
  
  scp_tab$patient=rownames(scp_tab)
  neurobl_tab$patient=rownames(neurobl_tab)
  mature_tab$patient=rownames(mature_tab)
  
  scp_tab$age=all$Age.at.Diagnosis.in.Days
  neurobl_tab$age=all$Age.at.Diagnosis.in.Days
  mature_tab$age=all$Age.at.Diagnosis.in.Days
  
  scp_tab[,1:(ncol(scp_tab)-3)]=log10(scp_tab[,1:(ncol(scp_tab)-3)]+1)
  neurobl_tab[,1:(ncol(neurobl_tab)-3)]=log10(neurobl_tab[,1:(ncol(neurobl_tab)-3)]+1)
  mature_tab[,1:(ncol(mature_tab)-3)]=log10(mature_tab[,1:(ncol(mature_tab)-3)]+1)
  
  risk_list=list(low_risk_4s="low_risk_4s",intermediate_risk="intermediate_risk", high_risk="high_risk",
                 high_risk_relapse="high_risk_relapse",
                 all=c("high_risk", "intermediate_risk", "low_risk_4s", "high_risk_relapse"))
  scp_cor=sapply(risk_list, function(x){sapply(colnames(scp_tab)[1:(ncol(scp_tab)-3)], function(y){
    cor(scp_tab[which(scp_tab$risk%in%x),y], scp_tab$age[which(scp_tab$risk%in%x)])
  })
  })
  
  neurobl_cor=sapply(risk_list, function(x){sapply(colnames(neurobl_tab)[1:(ncol(neurobl_tab)-3)], function(y){
    cor(neurobl_tab[which(neurobl_tab$risk%in%x),y], neurobl_tab$age[which(neurobl_tab$risk%in%x)])
  })
  })
  
  mature_cor=sapply(risk_list, function(x){sapply(colnames(mature_tab)[1:(ncol(mature_tab)-3)], function(y){
    cor(mature_tab[which(mature_tab$risk%in%x),y], mature_tab$age[which(mature_tab$risk%in%x)])
  })
  })
  
  scp_cor=as.matrix(scp_cor)
  neurobl_cor=as.matrix(neurobl_cor)
  mature_cor=as.matrix(mature_cor)
  
  
  col_fun=colorRamp2(c(-1, 0, 1), c("brown", "white", "green"))
  scp_hm=Heatmap(scp_cor, cluster_rows = F,col=col_fun, cluster_columns = F, 
                 column_names_side = "top", row_names_side = "left", row_order = scp_markers, column_title = "SCP")
  neurobl_hm=Heatmap(neurobl_cor, cluster_rows = F, col=col_fun,cluster_columns = F,
                     column_names_side = "top", row_names_side = "left", row_order = neurobl_mesen, column_title = "Neuroblast")
  mature_hm=Heatmap(mature_cor, cluster_rows = F,col=col_fun, cluster_columns = F,
                    column_names_side = "top", row_names_side = "left", column_title = "Mature")
  
  pdf("plts/new_plts/age_cor/target_age_cor.pdf", height = 8, width = 12)
  print(neurobl_hm)
  print(scp_hm)
  print(mature_hm)
  dev.off()
}
get_target_corr(target_tpm, scp_mesen, neurobl_mesen, mature_mesen)


get_german_corr=function(exp_dat, scp_markers, neurobl_markers, mature_markers){
  low_risk=dplyr::filter(german_clinical_data,
                         german_clinical_data$Age.at.Diagnosis..d.<540 & german_clinical_data$MYCN.Status=="not amplified" &german_clinical_data$Stage..INSS.!="4S" &!german_clinical_data$Status%in%c("relapse/progression","death of disease"))
  high_risk=dplyr::filter(german_clinical_data,german_clinical_data$Age.at.Diagnosis..d.>540&german_clinical_data$MYCN.Status=="amplified", german_clinical_data$Status%in%c("relapse/progression","death of disease"))
  risk_4s=dplyr::filter(german_clinical_data, german_clinical_data$Stage..INSS.=="4S")
  pt1=rbind(low_risk, high_risk, risk_4s)
  pt2=german_clinical_data[-which(german_clinical_data$Tumor.ID%in%pt1$Tumor.ID),]
  inter_death_or_relapse=dplyr::filter(pt2,Status%in%c("relapse/progression","death of disease"))
  inter_no_event=pt2[-which(pt2$Tumor.ID%in%inter_death_or_relapse$Tumor.ID),]
  low_risk$risk_level="low_risk"
  inter_no_event$risk_level="inter_no_event"
  inter_death_or_relapse$risk_level="inter_death_or_relapse"
  high_risk$risk_level="high_risk"
  risk_4s$risk_level="4s"
  
  extremes=rbind(low_risk, risk_4s, high_risk)
  all=rbind(low_risk, risk_4s, inter_no_event, inter_death_or_relapse, high_risk)
  
  scp_tab=exp_dat[which(counts_gs$gene_symbol%in%scp_markers),c(all$Tumor.ID, "gene_symbol")]
  neurobl_tab=exp_dat[which(counts_gs$gene_symbol%in%neurobl_markers),c(all$Tumor.ID, "gene_symbol")]
  mature_tab=exp_dat[which(counts_gs$gene_symbol%in%mature_markers),c(all$Tumor.ID, "gene_symbol")]
  
  scp_tab$gene_symbol = make.unique(as.character(scp_tab$gene_symbol))
  neurobl_tab$gene_symbol=make.unique(as.character(neurobl_tab$gene_symbol))
  mature_tab$gene_symbol=make.unique(as.character(mature_tab$gene_symbol))
  
  rownames(scp_tab)=scp_tab$gene_symbol
  rownames(neurobl_tab)=neurobl_tab$gene_symbol
  rownames(mature_tab)=mature_tab$gene_symbol
  scp_tab=scp_tab[,1:ncol(scp_tab)-1]
  neurobl_tab=neurobl_tab[,1:ncol(neurobl_tab)-1]
  mature_tab=mature_tab[,1:ncol(mature_tab)-1]
  
  scp_transposed=as.data.frame(t(scp_tab))
  neurobl_transposed=as.data.frame(t(neurobl_tab))
  mature_transposed=as.data.frame(t(mature_tab))
  
  scp_transposed$patient=rownames(scp_transposed)
  neurobl_transposed$patient=rownames(neurobl_transposed)
  mature_transposed$patient=rownames(mature_transposed)
  
  
  scp_transposed$risk= all$risk_level
  neurobl_transposed$risk=all$risk_level
  mature_transposed$risk=all$risk_level
  
  scp_transposed$age=all$Age.at.Diagnosis..d.
  neurobl_transposed$age=all$Age.at.Diagnosis..d.
  mature_transposed$age=all$Age.at.Diagnosis..d.
  
  scp_transposed[,1:(ncol(scp_transposed)-3)]=log10(scp_transposed[,1:(ncol(scp_transposed)-3)]+1)
  neurobl_transposed[,1:(ncol(neurobl_transposed)-3)]=log10(neurobl_transposed[,1:(ncol(neurobl_transposed)-3)]+1)
  mature_transposed[,1:(ncol(mature_transposed)-3)]=log10(mature_transposed[,1:(ncol(mature_transposed)-3)]+1)
  
  risk_list=list(low_risk="low_risk", risk_4s="4s",inter_no_event="inter_no_event",
                 inter_death_or_relapse="inter_death_or_relapse",high_risk="high_risk",
                 all=c("low_risk", "inter_no_event", "inter_death_or_relapse", "high_risk", "4s"))
  
  scp_cor=sapply(risk_list, function(x){sapply(colnames(scp_transposed)[1:(ncol(scp_transposed)-3)], function(y){
    cor(scp_transposed[which(scp_transposed$risk%in%x),y], scp_transposed$age[which(scp_transposed$risk%in%x)])
  })
  })
  
  neurobl_cor=sapply(risk_list, function(x){sapply(colnames(neurobl_transposed)[1:(ncol(neurobl_transposed)-3)], function(y){
    cor(neurobl_transposed[which(neurobl_transposed$risk%in%x),y], neurobl_transposed$age[which(neurobl_transposed$risk%in%x)])
  })
  })
  
  mature_cor=sapply(risk_list, function(x){sapply(colnames(mature_transposed)[1:(ncol(mature_transposed)-3)], function(y){
    cor(mature_transposed[which(mature_transposed$risk%in%x),y], mature_transposed$age[which(mature_transposed$risk%in%x)])
  })
  })
  
  scp_cor=as.matrix(scp_cor)
  neurobl_cor=as.matrix(neurobl_cor)
  mature_cor=as.matrix(mature_cor)
  
  scp_german_order=c("MPZ", "PLP1", "TMEM176A", "TMEM176B", "SCRG1", "CRABP1", "CRABP1.1", "CRABP1.2", "CRABP1.3",
                     "SEMA3B", "SEMA3B.1", "ERBB3", "ERBB3.1","ERBB3.2", "ERBB3.3", "SOX10", "SOX10.1", "SOX2", "SOX2.1", 
                     "OLFML2A", "OLFML2A.1", "CDH19","CDH19.1" ,"CMTM5", "HMGA2", "HMGA2.1" ,"HMGA2.2" )
  neurobl_german_order=c("SCG2", "SCG2.1", "CYB561", "CYB561.1", "CYB561.2", "CHRNA3", "SLC18A1", "FAIM2", "FAIM2.1", "FAIM2.2",
                         "DDC", "DDC.1","PHOX2A", "PHOX2A.1", "INSM1", "INSM1.1", "GATA2", "GATA2.1", "GATA2.2", "SHD", 
                         "C1QL1", "TAGLN3", "TAGLN3.1", "SYT5", "SYT5.1", "HAND2", "MIAT","MIAT.1","MIAT.2", "MIAT.3",
                         "PGAM2", "SYP", "KCNQ2", "KCNQ2.1", "KCNQ2.2","PCSK2", "SCG3", "SCG3.1", "GCH1", "GCH1.1",
                         "ICA1", "ICA1.1", "KRT19", "FAM57B", "SNAP25", "SNAP25.1", "ELAVL4", "ELAVL4.1","ELAVL4.2", "SST")
  
  col_fun=colorRamp2(c(-1, 0, 1), c("brown", "white", "green"))
  scp_hm=Heatmap(scp_cor, cluster_rows = F,col=col_fun, cluster_columns = F, 
                 column_names_side = "top", row_names_side = "left", column_title = "SCP", row_order = scp_german_order)
  neurobl_hm=Heatmap(neurobl_cor, cluster_rows = F, col=col_fun,cluster_columns = F,
                     column_names_side = "top", row_names_side = "left", column_title = "Neuroblast", row_order = neurobl_german_order)
  mature_hm=Heatmap(mature_cor, cluster_rows = F,col=col_fun, cluster_columns = F,
                    column_names_side = "top", row_names_side = "left", column_title = "Mature")
  
  pdf("plts/new_plts/age_cor/german_age_cor.pdf", height = 8, width = 12)
  print(neurobl_hm)
  print(scp_hm)
  print(mature_hm)
  dev.off()
}

get_german_corr(counts_gs, scp_mesen, neurobl_mesen, mature_mesen)

