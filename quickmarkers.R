quickMarkers = function(toc,clusters,N=10,FDR=0.01,expressCut=0){
  #Convert to the more manipulable format
  toc = as(toc,'dgTMatrix')
  w = which(toc@x>expressCut)
  #Get the counts in each cluster
  clCnts = table(clusters)
  nObs = split(factor(rownames(toc))[toc@i[w]+1],clusters[toc@j[w]+1])
  nObs = sapply(nObs,table)
  #Calculate the observed and total frequency
  nTot = rowSums(nObs)
  tf = t(t(nObs)/as.integer(clCnts[colnames(nObs)]))
  ntf = t(t(nTot - nObs)/as.integer(ncol(toc)-clCnts[colnames(nObs)]))
  idf = log(ncol(toc)/nTot)
  score = tf*idf
  #Calculate p-values
  qvals = lapply(seq_len(ncol(nObs)),function(e)
    p.adjust(phyper(nObs[,e]-1,nTot,ncol(toc)-nTot,clCnts[colnames(nObs)[e]],lower.tail=FALSE),method='BH'))
  qvals = do.call(cbind,qvals)
  colnames(qvals) = colnames(nObs)
  #Now get the top N for each group
  w = lapply(seq_len(ncol(nObs)),function(e){
    o = order(score[,e],decreasing=TRUE)
    if(sum(qvals[,e]<FDR)>=N){
      o[seq(N)]
    }else{
      o[qvals[o,e]<FDR]
    }
  })
  #Now construct the data.frame with everything
  ww = cbind(unlist(w,use.names=FALSE),rep(seq_len(ncol(nObs)),lengths(w)))
  out = data.frame(gene = rownames(nObs)[ww[,1]],
                   cluster = colnames(nObs)[ww[,2]],
                   geneFrequency = tf[ww],
                   geneFrequencyOutsideCluster = ntf[ww],
                   geneFrequencyGlobal = nTot[ww[,1]]/ncol(toc),
                   tfidf = score[ww],
                   idf = idf[ww[,1]],
                   qval = qvals[ww])
  return(out)
}
