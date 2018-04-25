# Simulation functions
#
# Author: Xuran Wang
##################################################################################

#' Simulate Single cell read counts
#'
#' Simulate expected library sizes from a log-normal distribution
#'
#' @param N integer, number of subjects in total.
#' @param n.bulk integer, number of bulk subjects.
#' @param G integer, number of genes.
#' @param mu matrix, 2 by G, the mean expression of genes in each cell types.
#' @param sub.var matrix, 2 by G, scale of log-normal distribution for each subjects at genes. Part of them are zeros.
#' @param sort.bulk logical, default is FALSE. Sort the cell type proportion.
#'
#' @details
#' Generate simulation datasets. n.sc or n.bulk is zero, then we only generate
#' one dataset and cell type proportions.
#'
#' @return a list with elements
#'    * cell type proportions of generated bulk data
#'    * cell type proportions of generated single cell data
#'    * ExpressionSet of bulk data
#'    * ExpressionSet of single cell data
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats rlnorm
#'
#' @export
Twocelltype.Generator = function(n.sc = 10, n.bulk = 100, G = 100, mu, sub.var, sort.bulk = FALSE, sigma = NULL, ...){
  K = 2; #2 cell types
  N = n.sc + n.bulk

  #generate cell number and cell type proportion
  Ni = sample(700:2000, N, replace = T)
  p = rdirichlet(N, rep(4, K))
  Nk = cbind(round(p[,1]*Ni), round(p[,2]*Ni))
  p = Nk/Ni

  muA = mu[1, ]   #c(rep(0,10), runif(90, 20, 70))
  muB = mu[2, ]   #c(rep(50, 10), rep(0, 90)) + muA
  #sf = sample(1:10, 5); muA[sf] = 50; muB[sf] = 0;

  # Generate subject level variation
  sub.Fac = list(sapply(1:G, function(i){rlnorm(N, 0, sub.var[1, i])}),
                 sapply(1:G, function(i){rlnorm(N, 0, sub.var[2, i])}) )
  GN = NULL;
  for(g in 1:G){
    GN = c(GN, paste('gene', g, sep = ''))
  }
  if(G>9){
    GN[1:9] = c('gene01', 'gene02', 'gene03', 'gene04', 'gene05', 'gene06', 'gene07', 'gene08', 'gene09');
  }else{
    GN[1:G] = paste0('gene0', 1:G)
  }

  Bulk.null = FALSE; SC.null = FALSE;
  if(n.bulk == 0) Bulk.null = TRUE;
  if(n.sc == 0) SC.null = TRUE;

  if(is.null(sigma)){
    if(!Bulk.null){
      Xjg = NULL    # generate bulk tissue data
      for(i in 1:n.bulk){
        Xjg = cbind( Xjg, colSums( sapply(muA*sub.Fac[[1]][i, ], rpois, n = Nk[i,1]) )
                     + colSums( sapply(muB*sub.Fac[[2]][i, ], rpois, n = Nk[i,2]) ) )
      }
    }
    if(!SC.null){
      nk = NULL; #S = NULL; thetaA = NULL; thetaB = NULL;   # generate single cell statistics
      Xjgc = NULL; #Xsc.cov = NULL# generate single cell data
      for(i in (1+n.bulk):N){
        nk1 = max(rbinom(1, Nk[i, 1], 0.1), 2)
        Xjgc1 = t( sapply(muA*sub.Fac[[1]][i, ], rpois, n = nk1) )
        #thetaA = cbind(thetaA, rowSums(Xjgc1)/sum(Xjgc1))

        nk2 = max(rbinom(1, Nk[i, 2], 0.1), 2)
        Xjgc2 =t( sapply(muB*sub.Fac[[2]][i, ], rpois, n = nk2) )
        #thetaB = cbind(thetaB, rowSums(Xjgc2)/sum(Xjgc2))

        #S =rbind(S, c( mean( colSums(Xjgc1)), mean( colSums(Xjgc2))) )

        nk = rbind(nk, c(nk1, nk2));
        Xjgc = cbind(Xjgc, Xjgc1, Xjgc2);
        #Xsc.cov = cbind(Xsc.cov, rbind(rep(paste('indv', i-100, sep = '')), c(rep('A', nk1), rep('B', nk2) )) )
      }
    }
  }else{
    lambdaA = sapply(sigma/muA, rgamma, n = N, shape = sigma)
    lambdaB = sapply(sigma/muB, rgamma, n = N, shape = sigma)
    if(!Bulk.null){
      Xjg = NULL    # generate bulk tissue data
      for(i in 1:n.bulk){
        Xjg = cbind( Xjg, colSums( sapply(lambdaA[i, ] * sub.Fac[[1]][i, ], rpois, n = Nk[i,1]) )
                     + colSums( sapply(lambdaB[i, ] * sub.Fac[[2]][i, ], rpois, n = Nk[i,2]) ) )
      }
    }
    if(!SC.null){
      nk = NULL; #S = NULL; thetaA = NULL; thetaB = NULL;   # generate single cell statistics
      Xjgc = NULL; #Xsc.cov = NULL# generate single cell data
      for(i in (1+n.bulk):N){
        nk1 = max(rbinom(1, Nk[i, 1], 0.1), 2)
        Xjgc1 = t( sapply(lambdaA[i, ] * sub.Fac[[1]][i, ], rpois, n = nk1) )
        #thetaA = cbind(thetaA, rowSums(Xjgc1)/sum(Xjgc1))

        nk2 = max(rbinom(1, Nk[i, 2], 0.1), 2)
        Xjgc2 =t( sapply(lambdaB[i, ]*sub.Fac[[2]][i, ], rpois, n = nk2) )
        #thetaB = cbind(thetaB, rowSums(Xjgc2)/sum(Xjgc2))

        #S =rbind(S, c( mean( colSums(Xjgc1)), mean( colSums(Xjgc2))) )

        nk = rbind(nk, c(nk1, nk2));
        Xjgc = cbind(Xjgc, Xjgc1, Xjgc2);
        #Xsc.cov = cbind(Xsc.cov, rbind(rep(paste('indv', i-100, sep = '')), c(rep('A', nk1), rep('B', nk2) )) )
      }
    }
  }
  if(!SC.null){
    stop.idx = cumsum( rowSums(nk) )
    start.idx = c(1, stop.idx[-(N-n.bulk)]+1)
    idx.level = unlist(sapply(1:(N-n.bulk), function(i){rep(i, rowSums(nk)[i])}))
    ct.level = unlist(sapply(1:(N-n.bulk), function(i){c(rep(1, nk[i,1]), rep(2, nk[i, 2]))}))
  }
  if(!Bulk.null){
    colnames(p) = c('A', 'B');
    p.bulk = p[1:n.bulk, ];
    if(sort.bulk){
      op.bulk = order(p.bulk[,1])
      p.bulk = p.bulk[op.bulk, ]; rownames(p.bulk) = paste0('Bulk', 1:n.bulk)
      Xjg = Xjg[, op.bulk]; colnames(Xjg) <- paste0('Bulk', 1:n.bulk)
    }else{
      rownames(p.bulk) = paste0('Bulk', 1:n.bulk)
      colnames(Xjg) <- paste0('Bulk', 1:n.bulk)
    }
    rownames(Xjg) <- GN;

    Bulk.pData = data.frame(sampleID = 1:n.bulk, SubjectName = factor(colnames(Xjg), levels = colnames(Xjg)) )
    rownames(Bulk.pData) = colnames(Xjg);
    metadata <- data.frame(labelDescription = c("Sample ID", "Subject Name"), row.names = c("sampleID", "SubjectName"))
    Bulk.eset <- ExpressionSet(assayData = data.matrix(Xjg), phenoData = new("AnnotatedDataFrame", data = Bulk.pData, varMetadata = metadata))
  }

  if(!SC.null){
    p.sc = p[(n.bulk+1):N, ]; rownames(p.sc) = paste0("Sub", (n.bulk+1):N);
    rownames(Xjgc) <- GN
    colnames(Xjgc) <- paste0('cell', 1:ncol(Xjgc))
    SC.pData = data.frame(sampleID = ((n.bulk+1):N)[idx.level], SubjectName = factor(paste0("Sub", (n.bulk+1):N)[idx.level], levels = paste0("Sub", (n.bulk+1):N)),
                          cellTypeID = ct.level, cellType = c('A', 'B')[ct.level])

    rownames(SC.pData) = colnames(Xjgc)

    metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"),
                           row.names=c("sampleID", "SubjectName","cellTypeID", "cellType"))

    SC.eset = ExpressionSet(assayData = data.matrix(Xjgc), phenoData = new("AnnotatedDataFrame", data = SC.pData, varMetadata = metadata))
  }
  if(SC.null){
    p.sc = NULL; SC.eset = NULL;
  }
  if(Bulk.null){
    p.bulk = NULL; Bulk.eset = NULL;
  }
  out = list(p.bulk = p.bulk, p.sc = p.sc, Bulk.eset = Bulk.eset, SC.eset = SC.eset)
  return(out)
}


#' Estimate cell type proportions by different methods and evaluate estimations
#'
#'
#' @param bulk.eset
#' @param sc.eset
#' @param p.bulk
#'
#' @details
#'
#' @return a matrix of evaluations
#'
Est_sim = function(bulk.eset, sc.eset, p.bulk){
  simMarkers.m = list(A = rownames(sc.eset)[( rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'A', drop = F]) ) - rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'B', drop = F]) ) )/ ( rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'A', drop = F]) ) + rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'B', drop = F]) ) ) > 0.2],
                      B = rownames(sc.eset)[( rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'A', drop = F]) ) - rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'B', drop = F]) ) )/ ( rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'A', drop = F]) ) + rowMeans( exprs(sc.eset[ , sc.eset$cellType == 'B', drop = F]) ) ) < -0.2])
  cat('length of proportion selected genes: ', length(unlist(simMarkers.m)), '\n')

  Theta = w_nnls_Theta(x = sc.eset, clusters = 'cellType', samples = 'sampleID')
  simMarkers.s = list(A = rownames(sc.eset)[sapply(1:G, function(g){t.test(Theta[G*(0:9)+g, 1], Theta[G*(0:9)+g, 2], alternative = 'greater')$p.value})< 0.1],
                      B = rownames(sc.eset)[sapply(1:G, function(g){t.test(Theta[G*(0:9)+g, 1], Theta[G*(0:9)+g, 2], alternative = 'less')$p.value})< 0.1])
  cat('length of t.test selected genes: ', length(unlist(simMarkers.s)), '\t')
  # Try to filter the estimated
  sim.out1.est = w_nnls_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, samples = 'SubjectName',
                             clusters = 'cellType', verbose = F)
  sim.out1.marker.est = w_nnls_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, samples = 'SubjectName',
                                    clusters = 'cellType', markers = unlist(simMarkers), verbose = F)
  sim.out1.marker.w.est = w_nnls_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, samples = 'SubjectName',
                                      clusters = 'cellType', markers = unlist(simMarkers.w), verbose = F)
  sim.out1.marker.m.est = w_nnls_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, samples = 'SubjectName',
                                      clusters = 'cellType', markers = unlist(simMarkers.m), verbose = F)
  sim.out1.marker.s.est = w_nnls_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, samples = 'SubjectName',
                                      clusters = 'cellType', markers = unlist(simMarkers.s), verbose = F)
  B.sim1 = bseqsc_basis(sc.eset, simMarkers, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)
  fit.sim1 <- bseqsc_proportions(exprs(bulk.eset), B.sim1, verbose = F)
  Est.prop.bseq.sim1 = t(coef(fit.sim1))
  B.sim1.w = bseqsc_basis(sc.eset, simMarkers.w, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)
  fit.sim1 <- bseqsc_proportions(exprs(bulk.eset), B.sim1.w, verbose = F)
  Est.prop.bseq.sim1.w = t(coef(fit.sim1))
  B.sim1.m = bseqsc_basis(sc.eset, simMarkers.m, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)
  fit.sim1 <- bseqsc_proportions(exprs(bulk.eset), B.sim1.m, verbose = F)
  Est.prop.bseq.sim1.m = t(coef(fit.sim1))
  B.sim1.s = bseqsc_basis(sc.eset, simMarkers.s, clusters = 'cellType', samples = 'sampleID', ct.scale = TRUE)
  fit.sim1 <- bseqsc_proportions(exprs(bulk.eset), B.sim1.s, verbose = F)
  Est.prop.bseq.sim1.s = t(coef(fit.sim1))

  Eval.sim1 = Eval_multi(prop.real = p.bulk,
                         prop.est = list(sim.out1.est$Est.prop.allgene, sim.out1.est$Est.prop.weighted, Est.prop.bseq.sim1,
                                         sim.out1.marker.est$Est.prop.allgene, sim.out1.marker.est$Est.prop.weighted,
                                         sim.out1.marker.w.est$Est.prop.allgene, sim.out1.marker.w.est$Est.prop.weighted,
                                         Est.prop.bseq.sim1.w, sim.out1.marker.m.est$Est.prop.allgene, sim.out1.marker.m.est$Est.prop.weighted,
                                         Est.prop.bseq.sim1.m, sim.out1.marker.s.est$Est.prop.allgene, sim.out1.marker.s.est$Est.prop.weighted,
                                         Est.prop.bseq.sim1.s),
                         method.name = c('NNLS', 'W-NNLS', 'BSEQ-sc',
                                         'NNLS.\n marker', 'W-NNLS.\n marker', 'NNLS.\n wmarker',
                                         'W-NNLS.\n wmarker', 'BSEQ.w', 'NNLS.\n pmarker', 'W-NNLS.\n pmarker', 'BSEQ.p',
                                         'NNLS.\n smarker', 'W-NNLS\n.smarker', 'BSEQ.s'))
  return(Eval.sim1)
}

Eval_sim_boxplot = function(list, title.sup = ''){
  par(mfrow = c(3, 1))
  boxplot(t(list$R.sim), main = paste0('Boxplot of Pearson Correlation, ', title.sup))
  lines(c(-10, 100), rep(max(apply(list$R.sim, 1 , median)), 2), col = 2 )

  boxplot(t(list$RMSD.sim), main = paste0('Boxplot of RMSD, ', title.sup))
  lines(c(-10, 100), rep(min(apply(list$RMSD.sim, 1 , median)), 2), col = 2 )

  boxplot(t(list$mAD.sim), main = paste0('Boxplot of mAD, ', title.sup))
  lines(c(-10, 100), rep(min(apply(list$mAD.sim, 1 , median)), 2), col = 2 )
}
