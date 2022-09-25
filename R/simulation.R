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
#' \itemize{
#'    \item {p.bulk: cell type proportions of generated bulk data;}
#'    \item {p.sc: cell type proportions of generated single cell data;}
#'    \item {bulk.mtx: Expression Matrix of bulk data;}
#'    \item {sc.sce: SingleCellExperiment of single cell data;}
#'}
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
    
    bulk.mtx = data.matrix(Xjg)
  }
  
  if(!SC.null){
    p.sc = p[(n.bulk+1):N, ]; rownames(p.sc) = paste0("Sub", (n.bulk+1):N);
    rownames(Xjgc) <- GN
    colnames(Xjgc) <- paste0('cell', 1:ncol(Xjgc))
    SC.pData = data.frame(sampleID = ((n.bulk+1):N)[idx.level], SubjectName = factor(paste0("Sub", (n.bulk+1):N)[idx.level], levels = paste0("Sub", (n.bulk+1):N)),
                          cellTypeID = ct.level, cellType = c('A', 'B')[ct.level])
    
    rownames(SC.pData) = colnames(Xjgc)
    sc.sce = SingleCellExperiment(list(counts = data.matrix(Xjgc)), colData = SC.pData)
  }
  if(SC.null){
    p.sc = NULL; sc.sce = NULL;
  }
  if(Bulk.null){
    p.bulk = NULL; bulk.mtx = NULL;
  }
  out = list(p.bulk = p.bulk, p.sc = p.sc, bulk.mtx = bulk.mtx, sc.sce = sc.sce)
  return(out)
}

#' Evaluate simulation results by boxplots
#' 
Eval_sim_boxplot = function(list, title.sup = ''){
  par(mfrow = c(3, 1))
  boxplot(t(list$R.sim), main = paste0('Boxplot of Pearson Correlation, ', title.sup))
  lines(c(-10, 100), rep(max(apply(list$R.sim, 1 , median)), 2), col = 2 )
  
  boxplot(t(list$RMSD.sim), main = paste0('Boxplot of RMSD, ', title.sup))
  lines(c(-10, 100), rep(min(apply(list$RMSD.sim, 1 , median)), 2), col = 2 )
  
  boxplot(t(list$mAD.sim), main = paste0('Boxplot of mAD, ', title.sup))
  lines(c(-10, 100), rep(min(apply(list$mAD.sim, 1 , median)), 2), col = 2 )
}