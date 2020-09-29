adr_all <- readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/processedSeurat.RDS")
adr_all@meta.data$new_clust=as.character(adr_all@meta.data$seurat_clusters)
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("25"))]="SCPs"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("13"))]="Sympathoblastic"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("24"))]="Bridge"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("19", "20"))]="Chromaffin"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("3", "15", "14","17"))]="Endothelium"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("23", "11", "18"))]="Mesenchyme"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("22", "5", "28", "0", "1", "2", "10", "6", "8", "9","21"))]="Cortex"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("16", "26"))]="Leukocytes"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("30", "7", "4", "12"))]="Erythroblasts"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust%in%c("27", "31", "29"))]="Other"

#cell counts 

srat_inhouse@meta.data$GOSH_ID=as.character(srat_inhouse@meta.data$orig.ident)

srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7685340","STDY7685341","STDY7685342"))]="GOSH014"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7843576", "STDY7843577", "STDY7843578"))]="GOSH023"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7733084","STDY7733085", "STDY7733086"))]="GOSH019"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY8004894","STDY8004902", "STDY8004910"))]="GOSH025"
srat_inhouse@meta.data$GOSH_ID[which(srat_inhouse@meta.data$GOSH_ID%in%c("STDY7787237", "STDY7787238", "STDY7787239"))]="GOSH021"
DimPlot(srat_inhouse, group.by = "GOSH_ID")

DimPlot(adr_all, label=T)

table(srat_inhouse@meta.data$new_idents[which(srat_inhouse@meta.data$GOSH_ID=="GOSH-021")])

sapply(unique(srat_tumour@meta.data$unique_sample), function(x){
  table(srat_tumour@meta.data$new_idents[which(srat_tumour@meta.data$unique_sample==x)])
})

adr_cc=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Results/preProcess/fAdrenal/baseSeurat.RDS")


ba1=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal1_Seurat.RDS")
ba2=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
bilateral=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/bilateralAdrenals_Seurat.RDS")
bilat_tech=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/techComparison_Seurat.RDS")
adr3=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/seurat.RDS")

adr_cc = NormalizeData(adr_cc)
adr_cc =FindVariableFeatures(adr_cc, selection.method = "vst", nfeatures = 2000)
adr_cc = ScaleData(adr_cc, features = rownames(adr_cc))
adr_cc = RunPCA(adr_cc, npcs = 75)
adr_cc = FindNeighbors(adr_cc, dims=1:75)
adr_cc = FindClusters(adr_cc, resolution = 1)
adr_cc = RunUMAP(adr_cc, dims=1:75, min.dist = 0.5, n.neighbors = 50)

DimPlot(adr_cc, label=T)


adr_cc@meta.data$new_clust=as.character(adr_cc@meta.data$seurat_clusters)
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("19"))]="SCPs"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("28","25","33","11"))]="Sympathoblastic"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("30"))]="Bridge"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("40", "12"))]="Chromaffin"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("4", "5", "20", "22", "32", "38", "44", "46"))]="Vascular Endothelium"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("7", "21", "24", "27", "41"))]="Mesenchyme"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("0", "1", "3", "23", "18", "31", "6", "17",
                                                                 "14", "10","29","26", "9","35"))]="Cortex"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("16", "42", "36", "39"))]="Leukocytes"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("2", "8", "13", "15", "43"))]="Erythroblasts"
adr_cc@meta.data$new_clust[which(adr_cc@meta.data$new_clust%in%c("37", "34", "45", "47"))]="Other"

DimPlot(adr_cc, group.by = "new_clust")



w21_1=rownames(adr_cc@meta.data)[which(sapply(strsplit(names(adr_cc@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]
w21_2=rownames(adr_cc@meta.data)[which(sapply(strsplit(names(adr_cc@active.ident), "_"), "[", 6)%in%c("Adr8710634", "Adr8710635"))]


adr_cc@meta.data$sample_name="x"
adr_cc@meta.data$sample_name[which(sapply(strsplit(names(adr_cc@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]="w21_1"
adr_cc@meta.data$sample_name[which(sapply(strsplit(names(adr_cc@active.ident), "_"), "[", 6)%in%c("Adr8710634", "Adr8710635"))]="w21_2"
adr_cc@meta.data$sample_name[which(adr_cc@meta.data$orig.ident=="babyAdrenal1")]="w8"
adr_cc@meta.data$sample_name[which(adr_cc@meta.data$orig.ident=="babyAdrenal2")]="w8d6"
adr_cc@meta.data$sample_name[which(adr_cc@meta.data$orig.ident%in%c("5388STDY7717452","5388STDY7717453",
                                                                    "5388STDY7717454","5388STDY7717455"))]="w10d5_1"
adr_cc@meta.data$sample_name[which(adr_cc@meta.data$orig.ident%in%c("5388STDY7717456","5388STDY7717458",
                                                                    "5388STDY7717459"))]="w10d5_2"
adr_cc@meta.data$sample_name[which(adr_cc@meta.data$orig.ident%in%c("5698STDY7839907","5698STDY7839909",
                                                                    "5698STDY7839917"))]="w11"

as.data.frame(sapply(unique(adr_cc@meta.data$sample_name), function(x){
  table(adr_cc@meta.data$new_clust[which(adr_cc@meta.data$sample_name==x)])
}))

table(adr_cc@meta.data$sample_name)

WhichCells(adr_cc, idents=c(37,34,45,47))
