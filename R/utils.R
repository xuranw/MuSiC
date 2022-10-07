# Utility functions
#
# Author: Xuran Wang
# Copyright Xuran Wang (2018)
##########################################################################


#' Calculate relative abundance
#' @param X non-negative matrix for calculate relative abundance
#' @param by.col logical, default as TRUE
#' @export
relative.ab = function(X, by.col = TRUE){
  if(sum(X < 0) > 0){
    stop('Negative entry appears!')
  }
  if(by.col == T){
    RX = sweep(X, 2, colSums(X), '/')
  }else{
    RX = sweep(X, 1, rowSums(X), '/')
  }
  return(RX)
}

#' Calculate square root of all values
#' @return  sqrt(x)   x >= 0
#'         -sqrt(-x)  x < 0
#' @param x real number
#' @export
my.sqrt = function(x){
  sign(x)*sqrt(abs(x))
}

#' Calculate fpkm to tpm
fpkmToTpm <- function(fpkm){
  apply(fpkm, 2, function(x){
    exp(log(x) - log(sum(x)) + log(1e6))
  })
}

#' log transformation with regulation
#'
#' @param x numeric
#' @param nu regulation
#'
#' @export
my.log = function(x, nu){
  sign(x)*( log(abs(x) + nu)- log(nu) )
}

#' Calculate Row Means with NA remove
#'
#' @param x matrix
#' @param na.rm default FALSE
#'
#' @return vector of row means
#'
#' @export
my.rowMeans = function (x, na.rm = FALSE, dims = 1L){
  if (is.data.frame(x))
    x <- as.matrix(x)
  if(length(dn <- dim(x)) < 2L){
    return(x)
  }
  if (!is.array(x) || length(dn <- dim(x)) < 2L )
    stop("'x' must be an array of at least two dimensions")
  if (dims < 1L || dims > length(dn) - 1L)
    stop("invalid 'dims'")
  p <- prod(dn[-(id <- seq_len(dims))])
  dn <- dn[id]
  z <- if (is.complex(x))
    .Internal(rowMeans(Re(x), prod(dn), p, na.rm)) + (0+1i) *
    .Internal(rowMeans(Im(x), prod(dn), p, na.rm))
  else .Internal(rowMeans(x, prod(dn), p, na.rm))
  if (length(dn) > 1L) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[id]
  }
  else names(z) <- dimnames(x)[[1L]]
  z
}

#' Get upper triangle matrix
#'
#' @param cormat square matrix
#' @return upper triangle matrix
#'
#' @export
get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#' MuSiC Deconvolution
#'
#' This function is to calculate the MuSiC deconvolution proportions
#'
#' @param bulk.mtx matrix of expression for bulk data
#' @param sc.sce SingleCellExperiment for single cell data
#' @param markers vector or list of gene names, default as NULL. If NULL, use all genes that provided by both bulk and single cell dataset.
#' @param clusters character, the colData of single cell dataset used as clusters;
#' @param samples character,the colData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param verbose logical, default as TRUE.
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking reciprocal
#' @param eps Threshold of convergence
#' @param centered logic, subtract avg of Y and D
#' @param normalize logic, divide Y and D by their standard deviation
#' @return a list with elements:
#' \itemize{
#'    \item {Estimates of MuSiC;}
#'    \item {Estimates of NNLS;}
#'    \item {Weight of MuSiC;}
#'    \item {r.squared of MuSiC;}
#'    \item {Variance of MuSiC estimates.}
#'    }
#' @seealso
#' \code{\link{music_basis}}
#' @export
music_prop = function(bulk.mtx, sc.sce, markers = NULL, clusters, samples, select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE,
                      iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE, ... ){
  bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx = bulk.mtx[bulk.gene, ]
  if(is.null(markers)){
    sc.markers = bulk.gene
  }else{
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = sc.markers, clusters = clusters, samples = samples, select.ct = select.ct, cell_size = cell_size, ct.cov = ct.cov, verbose = verbose)
  cm.gene = intersect( rownames(sc.basis$Disgn.mtx), bulk.gene )
  if(is.null(markers)){
    if(length(cm.gene)< 0.2*min(length(bulk.gene), nrow(sc.sce)) )
      stop("Too few common genes!")
  }else{
    if(length(cm.gene)< 0.2*length(unlist(markers)))
      stop("Too few common genes!")
  }
  if(verbose){message(paste('Used', length(cm.gene), 'common genes...'))}
  
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx)); m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]; 
  M.S = colMeans(sc.basis$S, na.rm = T);
  
  if(!is.null(cell_size)){
    if(!is.data.frame(cell_size)){
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if(sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))){
      stop("Cell type names in cell_size must match clusters")
    }else if (any(is.na(as.numeric(cell_size[, 2])))){
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  
  Yjg = relative.ab(bulk.mtx[m.bulk, ]); N.bulk = ncol(bulk.mtx);
  if(ct.cov){
    Sigma.ct = sc.basis$Sigma.ct[, m.sc];
    
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    
    for(i in 1:N.bulk){
      if(sum(Yjg[, i] == 0) > 0){
        D1.temp = D1[Yjg[, i]!=0, ];
        Yjg.temp = Yjg[Yjg[, i]!=0, i];
        Sigma.ct.temp = Sigma.ct[, Yjg[,i]!=0];
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
      }else{
        D1.temp = D1;
        Yjg.temp = Yjg[, i];
        Sigma.ct.temp = Sigma.ct;
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
      }
      
      lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S, Sigma.ct.temp, iter.max = iter.max,
                                     nu = nu, eps = eps, centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }else{
    Sigma = sc.basis$Sigma[m.sc, ];
    
    valid.ct = (colSums(is.na(Sigma)) == 0)&(colSums(is.na(D1)) == 0)&(!is.na(M.S))
    
    if(sum(valid.ct)<=1){
      stop("Not enough valid cell type!")
    }
    
    if(verbose){message(paste('Used', sum(valid.ct), 'cell types in deconvolution...' ))}
    
    D1 = D1[, valid.ct]; M.S = M.S[valid.ct]; Sigma = Sigma[, valid.ct];
    
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for(i in 1:N.bulk){
      if(sum(Yjg[, i] == 0) > 0){
        D1.temp = D1[Yjg[, i]!=0, ];
        Yjg.temp = Yjg[Yjg[, i]!=0, i];
        Sigma.temp = Sigma[Yjg[,i]!=0, ];
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
      }else{
        D1.temp = D1;
        Yjg.temp = Yjg[, i];
        Sigma.temp = Sigma;
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
      }
      
      lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, Sigma.temp, iter.max = iter.max,
                                  nu = nu, eps = eps, centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)
  
  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
              Weight.gene = Weight.gene, r.squared.full = r.squared.full, Var.prop = Var.prop))
}

#' MuSiC Deconvolution with Clusters
#'
#' This function is to calculate the MuSiC deconvolution proportions with clusters
#'
#' @param bulk.mtx Matrix of expression for bulk data
#' @param sc.sce SingleCellExperiment for single cell data
#' @param group.markers list of gene names. The list include differential expressed genes within groups. 
#'        List name must be the same as `clusters.type`.
#' @param groups character, the colData of single cell data used as groups;
#' @param clusters character, the colData of single cell dataset used as clusters;
#' @param samples character,the colData of single cell dataset used as samples;
#' @param clusters.type list of cell types. The list identify groups of similar cell types.
#' @param verbose logical, default as TRUE.
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking reciprocal
#' @param eps Threshold of convergence
#' @param centered logic, subtract avg of Y and D
#' @param normalize logic, divide Y and D by their standard deviation
#' @return matrix of estimated proportions by MuSiC with cluster information.
#' @seealso
#' \code{\link{music_basis}}; \code{\link{music_prop}}
#' @export
music_prop.cluster = function(bulk.mtx, sc.sce, group.markers, groups, clusters, samples, clusters.type,
                              verbose = TRUE, iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE, ... ){
  bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx = bulk.mtx[bulk.gene, ]
  select.ct = unlist(clusters.type)
  
  if(length(setdiff(names(group.markers), names(clusters.type))) > 0 || length(setdiff(names(clusters.type), names(group.markers))) > 0){
    stop("Cluster number is not matching!")
  }else{
    group.markers = group.markers[names(clusters.type)]
  }
  
  
  if(verbose){message('Start: cluster estimations...')}
  cluster.sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = NULL, clusters = groups, samples = samples, select.ct = names(clusters.type), verbose = verbose)
  if(verbose){message('Start: cell type estimations...')}
  sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = NULL, clusters = clusters, samples = samples, select.ct = select.ct, verbose = verbose)
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  
  if(length(cm.gene)< 0.2*min(length(bulk.gene), nrow(sc.sce)) ){
    stop("Too few common genes!")
  }
  
  if(verbose){message(paste('Used', length(cm.gene), 'common genes...'))}
  
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx)); m.bulk = match(cm.gene, bulk.gene)
  group.markers = lapply(group.markers, intersect, cm.gene)
  
  D1 = sc.basis$Disgn.mtx[m.sc,]; M.S = sc.basis$M.S; Sigma = sc.basis$Sigma[m.sc, ]
  
  cluster.select = setdiff(rownames(D1), unique(unlist(group.markers)))
  cluster.diff = unique(unlist(group.markers))
  
  D1.cluster = cluster.sc.basis$Disgn.mtx[cluster.select, ]; M.S.cluster = cluster.sc.basis$M.S;
  Yjg = relative.ab(bulk.mtx[m.bulk, ]); N.bulk = ncol(bulk.mtx);
  
  Sigma.cluster = cluster.sc.basis$Sigma[cluster.select, ];
  
  D1.sub = cluster.sc.basis$Disgn.mtx[cluster.diff, ]; Sigma.sub = cluster.sc.basis$Sigma[cluster.diff, ];
  
  Est.prop.weighted.cluster = NULL
  for(i in 1:N.bulk){
    if(sum(Yjg[, i] == 0) > 0){
      name.temp = rownames(Yjg)[Yjg[, i] != 0]
      D1.cluster.temp = D1.cluster[rownames(D1.cluster)%in%name.temp, ];
      D1.sub.temp = D1.sub[rownames(D1.sub)%in%name.temp, ];
      Yjg.temp = Yjg[Yjg[, i]!=0, i];
      Sigma.cluster.temp = Sigma.cluster[rownames(Sigma.cluster)%in%name.temp, ];
      Sigma.sub.temp = Sigma.sub[rownames(Sigma.sub)%in%name.temp, ];
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
    }else{
      D1.cluster.temp = D1.cluster;
      D1.sub.temp = D1.sub;
      Yjg.temp = Yjg[, i];
      Sigma.cluster.temp = Sigma.cluster;
      Sigma.sub.temp = Sigma.sub;
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
    }
    lm.D1.cluster = music.iter(Yjg.temp, D1.cluster.temp, M.S.cluster, Sigma.cluster.temp, iter.max = iter.max,
                               nu = nu, eps = eps, centered = centered, normalize = normalize)
    p.weight = NULL
    p.cluster.weight = lm.D1.cluster$p.weight
    for(j in 1:length(clusters.type)){
      if(length(clusters.type[[j]]) == 1){
        p.weight = c(p.weight, p.cluster.weight[j])
      }else{
        if(p.cluster.weight[j] == 0){
          p.weight = c(p.weight, rep(0, length(clusters.type[[j]])))
        }else{
          c.marker = intersect(group.markers[[j]], names(Yjg.temp))
          Y.sub = D1.sub.temp[c.marker, j]*p.cluster.weight[j] +
            (Yjg.temp[c.marker] - D1.sub.temp[c.marker, ]%*% p.cluster.weight) * p.cluster.weight[j]
          names(Y.sub) = c.marker
          Y.sub = Y.sub[Y.sub > 0]
          lm.D1.sub = music.iter(Y.sub, D1[c.marker, clusters.type[[j]]], M.S[clusters.type[[j]]],
                                 Sigma[c.marker, clusters.type[[j]]])
          p.weight = c(p.weight, p.cluster.weight[j] * lm.D1.sub$p.weight)
        }
      }
    }
    Est.prop.weighted.cluster = rbind(Est.prop.weighted.cluster, p.weight)
  }
  
  colnames(Est.prop.weighted.cluster) = unlist(clusters.type)
  rownames(Est.prop.weighted.cluster) = colnames(Yjg)
  
  return(list(Est.prop.weighted.cluster = Est.prop.weighted.cluster))
}

################### ANOVA select Marker Genes ######################
#' ANOVA Test for each gene
#'
#' This function is to test the two-way significance of cell type specific expression cross samples
#'
#' @param sce SingleCellExperiment single cell dataset
#' @param non.zero logical, default as TRUE. If TRUE, we only use gene that have non-zero expression
#' @param markers vector of characters, default as NULL. If NULL, use all genes in \code{sce}
#' @param clusters character, the name of colData used as clusters
#' @param samples character, the name of colData used as samples
#' @param select.ct vector of cell types included, default as \code{NULL}. If \code{NULL}, include all cell types in \code{x}
#' @param num.info numeric, number of selected gene for each cell type. Default at 25
#'
#' @return a list of
#'   \itemize{
#'        \item {F statistics of \code{samples};}
#'        \item {F statistics of \code{clusters};}
#'        \item {F statistics of two-way anova;}
#'        \item {selected informative genes: high F statictis for \code{clusters} compare to \code{samples};}
#'}
#' @importFrom stats aov
#' @importFrom Matrix rowSums
#' @export
Anova_info = function(sce, non.zero = TRUE, markers = NULL, clusters, samples, select.ct = NULL, num.info = 25, ... ){
  if(!is.null(select.ct)){
    sce = sce[, sce@colData[, clusters] %in% select.ct]
    message(paste('Selected', length(select.ct), 'cell type(s) ...' ))
  }
  if(non.zero){  ## eliminate non expressed genes
    sce <- sce[Matrix::rowSums(counts(sce))>0, ]
    message(paste('Eliminating non expressed gene(s) ...' ))
  }
  if(!is.null(markers)){
    sce <- sce[unlist(markers), ]
    message(paste('Selected', length(unlist(markers)), 'marker gene(s) ...' ))
  }
  sce@assays$logcounts = log(relative.ab(counts(sce))+10^{-10}) #change counts to relative abundance
  message('Transform counts to log scaled relative abundance ...')
  clusters <- as.character(colData(sce)[ , clusters])
  samples <- as.character(colData(sce)[ , samples])
  cell.type.select = unique(as.character(colData(sce)[ , clusters]))
  ngene = nrow(sce);
  j = 1
  #message(paste0('Selected ', length(cell.type.select), ' cell type(s) ...'))
  F.indv = NULL; F.cell.type = NULL; F.inter = NULL;
  for(i in 1:nrow(sce)){
    f.indv = NULL; f.cell.type = NULL; f.inter = NULL;
    marker.data = data.frame(log.RA = sce@assays$logcounts[i, ], indv = samples, cell.type = clusters)
    s.aov = unlist(summary(aov(log.RA ~ indv*cell.type, data = marker.data)))[13:15]
    f.indv = c(f.indv, s.aov[1]); f.cell.type = c(f.cell.type, s.aov[2]); f.inter = c(f.inter, s.aov[3])
    
    for(k in 1:length(cell.type.select)){
      marker.data.1ct = marker.data;
      cell.type.1ct = rep(cell.type.select[k], ncol(sce));
      cell.type.1ct[(marker.data$cell.type != cell.type.select[k])] = paste('non', cell.type.select[k]);
      marker.data.1ct$cell.type = cell.type.1ct;
      ct.aov = unlist(summary(aov(log.RA ~ indv*cell.type, data = marker.data.1ct)) )[13:15];
      f.indv = c(f.indv, ct.aov[1]); f.cell.type = c(f.cell.type, ct.aov[2]); f.inter = c(f.inter, ct.aov[3]);
    }
    F.indv = rbind(F.indv, f.indv); F.cell.type = rbind(F.cell.type, f.cell.type);
    F.inter = rbind(F.inter, f.inter);
    if(ngene >=  10000){
      if(i == round(ngene*j*0.1)){
        message(paste(j*10, '% is done ...'));
        j = j+ 1;
      }
    }
  }
  #cat(dim(F.indv), dim(F.cell.type), dim(F.inter))
  rownames(F.indv) <- rownames(F.cell.type) <- rownames(F.inter) <- rownames(sce);
  colnames(F.indv) <- colnames(F.cell.type) <- colnames(F.inter) <- c('full', cell.type.select)
  marker.score = log(F.cell.type/F.indv)
  info.gene = NULL; num.info = min(num.info, nrow(F.indv));
  for(k in 1:ncol(F.indv)){
    info.gene = cbind(info.gene, rownames(marker.score)[order(marker.score[, k], decreasing = T)[1:num.info]])
  }
  colnames(info.gene) = colnames(marker.score)
  return( list(F.indv = F.indv, F.cell.type = F.cell.type, F.inter = F.inter, informative.gene = info.gene) )
}

