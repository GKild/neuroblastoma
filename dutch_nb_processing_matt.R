#############
# Libraries #
#############
library(glmnet)
library(Matrix.utils)
library(ComplexHeatmap)
packrat::init('~/SeuratV3')
packrat::set_opts(local.repos = "~/R/packages")
library(Seurat)
library(SoupX)
library(ggplot2)
library(cowplot)
library(Matrix)
source('https://raw.githubusercontent.com/constantAmateur/MCVR/master/code.R')
source('~/trueHome/scratch/src/logisticRegression.R')
import('~/trueHome/Projects/Common/geneSets.R',as='geneSets')
geneSets$S=cc.genes.updated.2019$s.genes
geneSets$G2M=cc.genes.updated.2019$g2m.genes

#############
# Functions #
#############
#' Run standard analysis
standard10X = function(srat,nPCs=50,res=1.0){
  srat = NormalizeData(srat)
  srat = ScaleData(srat)
  srat = FindVariableFeatures(srat)
  srat = RunPCA(srat,approx=FALSE)
  srat = RunUMAP(srat,dims=seq(nPCs))
  srat = FindNeighbors(srat,dims=seq(nPCs))
  srat = FindClusters(srat,res=res)
  return(srat)
}



##############
# Parameters #
##############
srcFile = '~/trueHome/scratch/dutchNB/scRNAseq_NB_PMC_Full.rdata'
outDir = '~/trueHome/scratch/dutchNB/'
nPCs=50
mtMax=0.2
nGeneMin=300
nCntsMin=500
minCells=100
#Normal ref
normRef = '~/trueHome/scratch/fAdrenals/babyAdrenal2.RDS'
#Gene conversion table for adrenals
gnsFile = file.path('~/trueHome/scratch/HCA/CMN1//4602STDY7685338/outs','filtered_gene_bc_matrices','GRCh38','genes.tsv')





####################
# Basic processing #
####################
##################
# Process NB data
load(srcFile)
gns = data$row.annot
#Split into analysis units
mDat = data$col.annot
mDat = split(mDat,mDat$PatientId)
mDat = lapply(mDat,function(e) split(e,e$SampleType))
#Get data for each
cnts=list()
for(pid in names(mDat)){
  cnts[[pid]]=list()
  for(st in names(mDat[[pid]])){
    cnts[[pid]][[st]] = data$data$counts[,rownames(mDat[[pid]][[st]])]
  }
}
##################################
# Train logisitc regression model
gnsConv = read.table(gnsFile,sep='\t',header=FALSE)
adr = readRDS(normRef)
#Keep only the genes in common
dat = adr@assays$RNA@counts
dat = dat[rownames(dat) %in% gns$gene,]
fit = trainModel(dat,paste0('fAd',as.character(adr@active.ident)),maxCells=5000)
srats=list()
qc=list()
#Now process Seurat for each
for(pid in names(mDat)){
  message(pid)
  srats[[pid]]=list()
  qc[[pid]]=list()
  for(st in names(mDat[[pid]])){
    message(st)
    #Get the counts
    cc = cnts[[pid]][[st]]
    #Split out the spike-ins
    spikeCnts = cc[grep('ERCC',rownames(cc)),]
    geneCnts = cc[grep('ERCC',rownames(cc),invert=TRUE),]
    #Collapse into gene symbol counts
    s = gns[rownames(geneCnts),'gene']
    s[is.na(s)] = rownames(geneCnts)[is.na(s)]
    geneCnts = aggregate(geneCnts,s)
    #Get things that we will use to decide what to drop
    spikeFrac = colSums(spikeCnts)/colSums(cc)
    mtFrac = colSums(geneCnts[grep('^MT-',rownames(geneCnts)),])/colSums(geneCnts)
    mtFrac[is.na(mtFrac)]=0
    nUMIs = colSums(geneCnts)
    nGenes = colSums(geneCnts>0)
    #Decide which ones to keep
    w = which(mtFrac < mtMax & nGenes > nGeneMin & nUMIs > nCntsMin)
    mm = mDat[[pid]][[st]]
    mm$spikeFrac = spikeFrac
    mm$mtFrac = mtFrac
    passed = rep(FALSE,length(mtFrac))
    passed[w]=TRUE
    qc[[pid]][[st]] = cbind(mm,nUMIs,nGenes,pid,st,passed)
    #Create Seurat Object
    message(sprintf("After filtering, keeping %d cells of %d",length(w),ncol(cc)))
    if(length(w)<minCells){
      message(sprintf("Sample has only %d cells passing QC.  Too shit to continue.",length(w)))
      next
    }
    srat = CreateSeuratObject(geneCnts[,w],meta.data=mm[w,])
    srats[[pid]][[st]] = standard10X(srat,nPCs=nPCs)
    id = paste0(pid,'___',st)
    id = gsub(' ','_',id)
    #Run logistic regression
    ps = predictSimilarity(fit,srat@assays$RNA@counts,classes=as.character(srats[[pid]][[st]]@active.ident))
    pdf(file.path(outDir,paste0(id,'.pdf')),width=14,height=14)
    plot(DimPlot(srats[[pid]][[st]],label=TRUE)+guides(colour=FALSE))
    plot(FeaturePlot(srats[[pid]][[st]],c('HBB','HBA1','PECAM1','PTPRC','EPCAM','PDGFRB')))
    plot(FeaturePlot(srats[[pid]][[st]],c('MYCN','SOX10','SOX2','PHOX2A','CHGB','PHOX2B')))
    print(similarityHeatmap(ps,
                            row_order = rownames(ps)[order(as.numeric(gsub('[A-Za-z]','',rownames(ps))))],
                            column_order = colnames(ps)[order(as.numeric(gsub('[A-Za-z]','',colnames(ps))))]
    ))
    dev.off()
  }
}
qq = do.call(rbind,lapply(qc,function(e) do.call(rbind,e)))
qq$id = paste0(qq$pid,'\n',qq$st)
#Plot metrics
pdf(file.path(outDir,'QC.pdf'),width=21,height=14)
gg = ggplot(qq,aes(x=spikeFrac)) +
  geom_density() + 
  ggtitle('spikeFrac') +
  facet_wrap(~id)
plot(gg)
gg = ggplot(qq,aes(x=mtFrac)) +
  geom_density() + 
  ggtitle('mtFrac') +
  facet_wrap(~id)
plot(gg)
gg = ggplot(qq,aes(x=log10(nUMIs+1))) +
  geom_density() + 
  ggtitle('log10(nUMIs)') +
  facet_wrap(~id)
plot(gg)
gg = ggplot(qq,aes(x=log10(nGenes+1))) +
  geom_density() + 
  ggtitle('log10(nGenes)') +
  facet_wrap(~id)
plot(gg)
dev.off()
#Number of cells retained for each
filtCnts = table(qq$id,ifelse(qq$passed,'PASS','FAIL'))
write.table(filtCnts,file.path(outDir,'QC_cnts.tsv'),sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)











