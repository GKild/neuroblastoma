import = function(src,as=NULL,searchDirs=.modPaths,all=FALSE,reattach=TRUE){
  #Force it to be a character
  src = as.character(substitute(src))
  #Get the calling environment
  #penv = parent.env(environment())
  #Make a new environment for the src file
  module = new.env()
  #Get the thing
  if(!file.exists(src)){
    #Try and find the thing in each of the search paths
    found=FALSE
    for(searchDir in searchDirs){
      if(file.exists(file.path(searchDir,sprintf('%s.r',src)))){
        src = file.path(searchDir,sprintf('%s.r',src))
        found=TRUE
        break
      }else if(file.exists(file.path(searchDir,sprintf('%s.R',src)))){
        src = file.path(searchDir,sprintf('%s.R',src))
        found=TRUE
        break
      }
    }
    if(!found)
      stop(sprintf("Could not find module %s on search path.",src))
  }
  sys.source(src,envir=module,keep.source=interactive())
  #Work out what to call it
  if(is.null(as))
    as=gsub('\\.[Rr]$','',basename(src))
  #Now do the assignment
  if(all){
    nom = sprintf("module:%s",as)
    if(nom %in% search()){
      if(reattach){
        detach(pos=match(nom,search()))
      }else{
        stop(sprintf("A module named %s is already attached.",as))
      }
    }
    attach(module,pos=2,name=nom)
  }else{
    penv = parent.frame()
    assign(as,module,envir=penv)
  }
  invisible(NULL)
}

##Core utilities and quality of life improvement functions
#importCore = function(src='~/Projects/Common/Code/INITR/core.R'){
#  source(src)
#}
#
##Functions for working with the Sanger infrastructure
#importSanger = function(src='~/Projects/Common/Code/INITR/sanger.R'){
#  source(src)
#}
#
##Function for working with genomic objects
#importGenome = function(src='~/Projects/Common/Code/INITR/genome.R'){
#  source(src)
#}
#
##Some customs plotting functions that are commonly useful
#importPlot = function(src='~/Projects/Common/Code/INITR/plot.R'){
#  source(src)
#}


quickCluster = function(dat,numPCs=75,clusteringRes=10.0){
  ss = CreateSeuratObject(dat)
  ss = NormalizeData(ss,verbose=FALSE)
  ss = FindVariableFeatures(ss,verbose=FALSE)
  ss = ScaleData(ss,verbose=FALSE)
  ss = RunPCA(ss,npcs=numPCs,approx=FALSE,verbose=FALSE)
  ss = FindNeighbors(ss,dims=seq(numPCs),verbose=FALSE)
  ss = FindClusters(ss,res=clusteringRes,verbose=FALSE)
  ss = RunUMAP(ss,dims=seq(numPCs),verbose=FALSE)
  return(ss)
}


runScrublet = function(dat10X,nPCs){
  fNom = tempfile()
  cmd = sprintf(
    "import scrublet as scr
import scipy.io
import numpy as np
import os
import gzip
os.chdir('%s')
input_dir = os.path.expanduser('%s')
output_file = '%s'
try:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
  genes = np.array(scr.load_genes(input_dir + '/genes.tsv', delimiter='\\t', column=1))
except FileNotFoundError:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
  #Copy relevant part of scr.load_genes
  gene_list = []
  gene_dict = {}
  with gzip.open(input_dir + '/features.tsv.gz','rt') as f:
    for l in f:
      gene = l.strip('\\n').split('\\t')[1]
      if gene in gene_dict:
        gene_dict[gene] += 1
        gene_list.append(gene + '__' + str(gene_dict[gene]))
        if gene_dict[gene] == 2:
          i = gene_list.index(gene)
          gene_list[i] = gene + '__1'
      else: 
        gene_dict[gene] = 1
        gene_list.append(gene)
  genes = np.array(gene_list)
#Convert to integers.  Will do nothing most of the time.
counts_matrix.data = np.round(counts_matrix.data).astype('int32')
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                        min_cells=3, 
                                                        min_gene_variability_pctl=85, 
                                                        n_prin_comps=%d)
with open(output_file,'w') as f:
  #Header row
  f.write('score\\tcall\\tsimScore1\\tsimScore2\\n')
  #As nSim = 2*nObs, we can split into two columns 
  offset = len(doublet_scores)
  f.write('\\n'.join([str(x)+'\\t'+str(predicted_doublets[i]).upper()+'\\t'+str(scrub.doublet_scores_sim_[i])+'\\t'+str(scrub.doublet_scores_sim_[offset+i]) for i,x in enumerate(doublet_scores)]))
  f.close()",getwd(),dat10X,fNom,nPCs)
  system(paste0('python -c "',cmd,'"'))
  out = read.table(fNom,header=TRUE,sep='\t')
  if(file.exists(file.path(dat10X,'barcodes.tsv'))){
    tmp = read.table(file.path(dat10X,'barcodes.tsv'),header=FALSE)
  }else{
    tmp = read.table(file.path(dat10X,'barcodes.tsv.gz'),header=FALSE)
  }
  out$barcode = tmp[,1]
  return(out)
}
