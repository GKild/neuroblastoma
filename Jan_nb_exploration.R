#Running logistic regression on Neuroblastoma1 using baby adrenal as reference
load("Jan_NB/scRNAseq_NB_PMC_Full.RData")

unique(data$col.annot$BiopsyType)

 sapply(colnames(data$col.annot), function(x){
   unique(data$col.annot[,x])
 })

patient_nb039 <- data$col.annot[which(data$col.annot$PatientId=="NB039"),]
sapply(colnames(patient_nb039), function(x){
  unique(patient_nb039[,x])
})

CreateSeuratObject(data$data$counts, min.cells = 3, min.features = 200)

x <- data$data$counts
plot(density(log10(colSums(x))))

genes_10x <- read.table("Neuroblastoma1/4602STDY7685340/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv", sep = "\t", header = F)
genes_10x$V1 <- as.character(genes_10x$V1)
genes_10x$V2 <- as.character(genes_10x$V2)
x2 <-x[which(rownames(x)%in%genes_10x$V1),]
plot(density((colSums(x)-colSums(x2))/colSums(x), na.rm = T))

mean((colSums(x)-colSums(x2))/colSums(x), na.rm=T)
hist((colSums(x)-colSums(x2))/colSums(x), na.rm=T)
dim(genes_10x)
32254/33694
plot(table(match(genes_10x$V1,rownames(x)))/nrow(genes_10x))
a <-as.data.frame(table(match(genes_10x$V1, rownames(x))))
mtIDs = genes_10x$V1[grep("^MT-",genes_10x$V2)] 
length(which(a$Freq==1))
match(genes_10x$V1, rownames(x))
x2 <-x[which(rownames(x)%in%mtIDs),]
dim(x2)
plot(density(colSums(x2)/colSums(x), na.rm = T))
median(colSums(x2)/colSums(x), na.rm = T)



mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(x)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "external_gene_name"),values=genes,mart= mart)
listAttributes(mart)
is.data.frame(G_list)
dim(x)
dim(G_list)
length(which(genes_10x$V1%in%G_list$ensembl_gene_id))
length(unique(data$row.annot$gene))
dim(x)
length(unique(G_list$external_gene_name))

write.table(G_list, "ensemb_to_gene_id.tsv", col.names = T, quote = F, row.names = F, sep = "\t")
data$row.annot
x2 <-x[which(!rownames(x)%in%genes_10x$V1),]
dim(x2)
a <-rowSums(x2)
sort(a, decreasing = T)
mtIDs = data$row.annot$ensemblId[grep("^MT-",data$row.annot$gene)] 
x3 <- x[mtIDs,]
plot(density(colSums(x3)/colSums(x), na.rm = T))
