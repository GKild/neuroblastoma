#get all adrenal datasets
ba1=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal1_Seurat.RDS")
ba2=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/babyAdrenal2_Seurat.RDS")
bilateral=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/bilateralAdrenals_Seurat.RDS")
bilat_tech=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/preProcessSCP/output/techComparison_Seurat.RDS")
adr3=readRDS("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/SCP/Data/fAdrenal19wk/seurat.RDS")

DimPlot(adr3, label=T)
dim(ba1)
dim(ba2)
dim(bilateral)
dim(bilat_tech)


DimPlot(ba1, label=T)
DimPlot(ba2, label=T)
DimPlot(bilateral, label=T)
DimPlot(bilat_tech,label=T)

FeaturePlot(ba1, features=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB','MYCN','SOX10','SOX2','PHOX2A','CHGB','PHOX2B'))

FeaturePlot(ba2, features=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB','MYCN','SOX10','SOX2','PHOX2A','CHGB','PHOX2B'))
FeaturePlot(bilat_tech, features=c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB','MYCN','SOX10','SOX2','PHOX2A','CHGB','PHOX2B'))
length(WhichCells(ba2, idents = c(19,7,11,17,10)))
length(WhichCells(ba2, idents = c(0,1,2,4,6,8,12,13)))
length(WhichCells(ba2, idents = c(18)))
length(WhichCells(ba2, idents = c(14)))
length(WhichCells(ba2, idents = c(3,9,15)))
length(WhichCells(ba2, idents = c(5,16)))
FeaturePlot(bilateral, features=c("EPCAM"))

length(WhichCells(ba1, idents = c(14,19,20)))
length(WhichCells(ba1, idents = c(0,1,2,4,6,7,8,10,11,13)))
length(WhichCells(ba1, idents = c(5,17)))
length(WhichCells(ba1, idents = c(16,18)))
length(WhichCells(ba1, idents = c(3,9,15,12)))
            
length(WhichCells(bilateral, idents = c(12,15,16,20,22)))
length(WhichCells(bilateral, idents = c(0,2,4,5,9,10,23,24)))

length(WhichCells(bilat_tech, idents = c(0,2,3,4,9,10,16,11)))
length(WhichCells(bilat_tech, idents = c(6,7,14)))

FeaturePlot(bilat_tech, c("SOX2"))

length(WhichCells(adr3,idents=c(0,1,2,4,6,9,13,14)))
length(WhichCells(adr3,idents=c(23)))


length(WhichCells(bilateral, idents = c(22)))
length(WhichCells(ba2, idents = c(17,10)))
length(WhichCells(bilat_tech, idents = c(11)))

WhichCells(bilateral, orig.ident="5388STDY7717452")



left=rownames(bilateral@meta.data)[which(bilateral@meta.data$orig.ident%in%c("5388STDY7717452", "5388STDY7717453", "5388STDY7717454", "5388STDY7717455"))]
right=rownames(bilateral@meta.data)[which(!bilateral@meta.data$orig.ident%in%c("5388STDY7717452", "5388STDY7717453", "5388STDY7717454", "5388STDY7717455"))]

left_srat=subset(bilateral, cells=left)
right_srat=subset(bilateral, cells=right)



length(WhichCells(left_srat, idents = c(12,15,16,20,22)))
length(WhichCells(left_srat, idents = c(0,2,4,5,9,10,23,24)))
length(WhichCells(left_srat, idents = c(18,19)))
length(WhichCells(left_srat, idents = c(1,21,13)))
length(WhichCells(left_srat, idents = c(8,11)))
length(WhichCells(left_srat, idents = c(3,6,7,14,17)))

length(WhichCells(right_srat, idents = c(12,15,16,20,22)))
length(WhichCells(right_srat, idents = c(0,2,4,5,9,10,23,24)))
length(WhichCells(right_srat, idents = c(18,19)))
length(WhichCells(right_srat, idents = c(1,21,13)))
length(WhichCells(right_srat, idents = c(8,11)))
length(WhichCells(right_srat, idents = c(3,6,7,14,17)))

length(WhichCells(right_srat, idents = c(22)))
sapply(strsplit(names(adr3@active.ident), "_"), "[", 6)
old_left=rownames(adr3@meta.data)[which(sapply(strsplit(names(adr3@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]
old_right=rownames(adr3@meta.data)[which(!sapply(strsplit(names(adr3@active.ident), "_"), "[", 6)%in%c("Adr8710632", "Adr8710633"))]
old_left_srat=subset(adr3, cells=old_left)
old_right_srat=subset(adr3, cells=old_right)
length(WhichCells(old_left_srat,idents=c(0,1,2,4,6,9,13,14)))
length(WhichCells(old_left_srat,idents=c(23)))
length(WhichCells(old_left_srat,idents=c(12,17,27,19,20)))
length(WhichCells(old_left_srat,idents=c(3,5,7,24,21,25)))
length(WhichCells(old_left_srat,idents=c(18,15)))
length(WhichCells(old_left_srat,idents=c(22,8,11,10,16)))

length(WhichCells(old_right_srat,idents=c(0,1,2,4,6,9,13,14)))
length(WhichCells(old_right_srat,idents=c(23)))
length(WhichCells(old_right_srat,idents=c(12,17,27,19,20)))
length(WhichCells(old_right_srat,idents=c(3,5,7,24,21,25)))
length(WhichCells(old_right_srat,idents=c(18,15)))
length(WhichCells(old_right_srat,idents=c(22,8,11,10,16)))

length(WhichCells(bilat_tech,idents=c(1,5,15)))
length(WhichCells(bilat_tech,idents=c(17,13,18)))
length(WhichCells(bilat_tech,idents=c(19,8,12)))
length(WhichCells(bilat_tech,idents=c(20,21)))
DimPlot(ba1, label=T)

tab=read.csv("plts/Book3.csv", header = T)

tab_melt=melt(tab,id.vars = "Sample")

ggplot(tab_melt, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="fill", stat="identity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499"))

