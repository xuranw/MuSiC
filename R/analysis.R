# Functions to perform the deconvolution analysis
#
# Author: Xuran Wang
#####################################################################################################

#' Estimate cell type proportion with MuSiC and NNLS
#'
#' @param Y vector of bulk tissue expression
#' @param X matrix, Signature matrix
#' @param S vector of Avg. Library size
#' @param Sigma matrix of Subject level variation (gene * cell type)
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking recipical
#' @param eps Thredshold of convergence
#'
#' @return a list with elements:
#'   * p.nnls: vector of cell type proportions estimated by nnls (add up to 1)
#'   * q.nnls: vector of original estimation from nnls
#'   * fit.nnls: fitted value of Y from nnls
#'   * resid.nnls: residual value from nnls
#'   * p.weight: vector of cell type proportions estimated by weighted-nnls (add up to 1)
#'   * q.weight: vector of original estimation from weight-nnls
#'   * fit.weight: fitted value of Y from weighted-nnls
#'   * resid.weight: residual value from weighted-nnls
#'   * weight.gene: weight calculated from weighted-nnls
#'   * converge: 'Reach Maxiter', 'Converge at m'
#'   * rsd: residual value from weighted-nnls transfromed by weight
#'   * R.squared: R square of weighted-nnls
#'   * var.p: variance of weighted-nnls estimated cell type proportions
#'
#' @export
#' @importFrom nnls nnls
music.basic = function(Y, X, S, Sigma, iter.max, nu, eps){
  k = ncol(X)

  lm.D = nnls(X, Y)
  r = resid(lm.D);
  weight.gene = 1/(nu + r^2 + colSums( (lm.D$x*S)^2*t(Sigma)  ))
  Y.weight = Y*sqrt(weight.gene)
  D.weight = sweep(X, 1, sqrt(weight.gene), '*')
  lm.D.weight = nnls(D.weight, Y.weight)
  p.weight = lm.D.weight$x/sum(lm.D.weight$x)
  p.weight.iter = p.weight
  r = resid(lm.D.weight)
  for(iter in 1:iter.max){
    weight.gene = 1/(nu + r^2 + colSums( (lm.D.weight$x*S)^2*t(Sigma)  ))
    Y.weight = Y*sqrt(weight.gene)
    D.weight = X * as.matrix(sqrt(weight.gene))[,rep(1,k)]
    lm.D.weight = nnls(D.weight, Y.weight )
    p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
    r.new = resid(lm.D.weight)
    if(sum(abs(p.weight.new - p.weight)) < eps){
      p.weight = p.weight.new;
      r = r.new
      R.squared = 1 - var(Y - X%*%as.matrix(lm.D.weight$x))/var(Y)
      fitted = X%*%as.matrix(lm.D.weight$x)
      var.p = diag(solve(t(D.weight)%*%D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
      return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
                  p.weight = p.weight, q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight =  Y - X%*%as.matrix(lm.D.weight$x),
                  weight.gene = weight.gene, converge = paste0('Converge at ', iter),
                  rsd = r, R.squared = R.squared, var.p = var.p));
    }
    p.weight = p.weight.new;
    r = r.new;
  }
  fitted = X%*%as.matrix(lm.D.weight$x)
  R.squared = 1 - var(Y - X%*%as.matrix(lm.D.weight$x))/var(Y)
  var.p = diag(solve(t(D.weight)%*%D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
  return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
              p.weight = p.weight, q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight =  Y - X%*%as.matrix(lm.D.weight$x),
              weight.gene = weight.gene, converge = 'Reach Maxiter', rsd = r,
              R.squared = R.squared, var.p = var.p))
}


#' Scaling bulk data and signature matrix and estimate cell type proportion
#'
#' @inheritParams music.basic
#' @param Y vector of bulk tissue expression
#' @param D matrix, Signature matrix
#' @param S vector of Avg. Library size
#' @param Sigma matrix of subject level variation (gene * cell type)
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, when take recipical
#' @param eps thredshold of converagence
#' @param centered logic, substract avg of Y and D
#' @param normalize logic, divide Y and D by their standard deviation
#'
#'
#' @return a list same as nnls.weight.basic
#' @export
music.iter = function(Y, D, S, Sigma, iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE){
  if(length(S)!=ncol(D)){
    common.cell.type = intersect(colnames(D), names(S))
    if(length(common.cell.type)<=1){
      stop('Not enough cell types!')
    }
    D = D[,match(common.cell.type, colnames(D))]
    S = S[match(common.cell.type, names(S))]
  }
  if(ncol(Sigma) != ncol(D)){
    common.cell.type = intersect(colnames(D), colnames(Sigma))
    if(length(common.cell.type)<=1){
      stop('Not enough cell type!')
    }
    D = D[, match(common.cell.type, colnames(D))]
    Sigma = Sigma[, match(common.cell.type, colnames(Sigma))]
    S = S[match(common.cell.type, names(S))]
  }
  k = ncol(D); # number of cell types

  common.gene = intersect(names(Y), rownames(D))
  if(length(common.gene)< 0.1*min(length(Y), nrow(D))){
    stop('Not enough common genes!')
  }
  Y = Y[match(common.gene, names(Y))];
  D = D[match(common.gene, rownames(D)), ]
  Sigma = Sigma[match(common.gene, rownames(Sigma)), ]

  X = D
  ## First, no intercept and no normalization
  if(centered){
    X = X - mean(X)
    Y = Y - mean(Y)
  }
  if(normalize){
    X = X/sd(as.vector(X));
    S = S*sd(as.vector(X));
    Y = Y/sd(Y)
  }else{
    Y = Y*100
  }
  lm.D = music.basic(Y, X, S, Sigma, iter.max = iter.max, nu = nu, eps = eps)
  return(lm.D)
}

#' Calculate weight with cross cell type covariance
#'
#' @param Sp vector, library size of each cell type
#' @param Sigma.ct matrix, cell type^2 by genes. Each columns store the cross cell-type covariance matrix of a gene.
#'
#' @return vector of weight
#'
#' @export
weight.cal.ct = function(Sp, Sigma.ct){
  nGenes = ncol(Sigma.ct); n.ct = length(Sp);
  weight = sapply(1:nGenes, function(g){
    sum(Sp%*%t(Sp)*matrix(Sigma.ct[,g], n.ct))
  })
  return(weight)
}

#' Calculate weight with cross-subject variance for each cell types
#'
#' @param p vector, cell type proportions
#' @param M.S vector of average library size for each cell type
#' @param Sigma gene by cell type matrix. Each entry is the cross-subject variance.
#'
#' @return vector of weight
#'
#' @export
Weight_cal <- function(p, M.S, Sigma){
  Sigma%*%(M.S * t(p))^2
}

#' Estimate cell type proportion with MuSiC and NNLS
#'
#' weight is estimated with cell type covariance
#'
#' @param Y vector of bulk tissue expression
#' @param X matrix, Signature matrix
#' @param S vector of Avg. Library size
#' @param Sigma.ct matrix of Subject level variation with
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking recipical
#' @param eps Thredshold of convergence
#'
#'
#' @return a list with elements:
#'   * p.nnls: vector of cell type proportions estimated by nnls (add up to 1)
#'   * q.nnls: vector of original estimation from nnls
#'   * fit.nnls: fitted value of Y from nnls
#'   * resid.nnls: residual value from nnls
#'   * p.weight: vector of cell type proportions estimated by weighted-nnls (add up to 1)
#'   * q.weight: vector of original estimation from weight-nnls
#'   * fit.weight: fitted value of Y from weighted-nnls
#'   * resid.weight: residual value from weighted-nnls
#'   * weight.gene: weight calculated from weighted-nnls
#'   * converge: 'Reach Maxiter', 'Converge at m'
#'   * rsd: residual value from weighted-nnls transfromed by weight
#'   * R.squared: R square of weighted-nnls
#'   * var.p: variance of weighted-nnls estimated cell type proportions
#'
#' @export
#' @importFrom nnls nnls
music.basic.ct = function(Y, X, S, Sigma.ct, iter.max, nu, eps){
  k = ncol(X)

  lm.D = nnls(X, Y)
  r = resid(lm.D);
  weight.gene = 1/(nu + r^2 + weight.cal.ct(lm.D$x*S, Sigma.ct))
  Y.weight = Y*sqrt(weight.gene)
  D.weight = sweep(X, 1, sqrt(weight.gene), '*')
  lm.D.weight = nnls(D.weight, Y.weight)
  p.weight = lm.D.weight$x/sum(lm.D.weight$x)
  p.weight.iter = p.weight
  r = resid(lm.D.weight)
  for(iter in 1:iter.max){
    weight.gene = 1/(nu + r^2 + weight.cal.ct(lm.D$x*S, Sigma.ct))
    Y.weight = Y*sqrt(weight.gene)
    D.weight = X * as.matrix(sqrt(weight.gene))[,rep(1,k)]
    lm.D.weight = nnls(D.weight, Y.weight )
    p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
    r.new = resid(lm.D.weight)
    if(sum(abs(p.weight.new - p.weight)) < eps){
      p.weight = p.weight.new;
      r = r.new
      R.squared = 1 - var(Y - X%*%as.matrix(lm.D.weight$x))/var(Y)
      fitted = X%*%as.matrix(lm.D.weight$x)
      var.p = diag(solve(t(D.weight)%*%D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
      return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
                  p.weight = p.weight, q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight =  Y - X%*%as.matrix(lm.D.weight$x),
                  weight.gene = weight.gene, converge = paste0('Converge at ', iter),
                  rsd = r, R.squared = R.squared, var.p = var.p));
      break;
    }
    p.weight = p.weight.new;
    r = r.new;
  }
  fitted = X%*%as.matrix(lm.D.weight$x)
  R.squared = 1 - var(Y - X%*%as.matrix(lm.D.weight$x))/var(Y)
  var.p = diag(solve(t(D.weight)%*%D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
  return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
              p.weight = p.weight, q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight =  Y - X%*%as.matrix(lm.D.weight$x),
              weight.gene = weight.gene, converge = 'Reach Maxiter', rsd = r,
              R.squared = R.squared, var.p = var.p))
}



#' Scaling bulk data and signature matrix and estimate cell type proportion
#'
#' @inheritParams music.basic.ct
#'
#' @param Y vector of bulk tissue expression
#' @param X matrix, Signature matrix
#' @param S vector of Avg. Library size
#' @param Sigma.ct matrix of Subject level variation with cross-cell-type covariance.
#' @param iter.max numeric, maximum iteration number. Default is 1000
#' @param nu regulation parameter, take care of weight when taking recipical
#' by default nu = 0.0001.
#' @param eps Thredshold of convergence. Default is 0.01.
#' @param centered logic, substract avg of Y and D.
#' Default is FALSE.
#' @param normalize logic, divide Y and D by their standard deviation.
#' Default is FALSE
#'
#' @return a list same as nnls.weight.basic
#' @export
music.iter.ct = function(Y, D, S, Sigma.ct, iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE){
  if(length(S)!=ncol(D)){
    common.cell.type = intersect(colnames(D), names(S))
    if(length(common.cell.type)<=1){
      stop('Not enough cell types!')
    }
    D = D[,match(common.cell.type, colnames(D))]
    S = S[match(common.cell.type, names(S))]
  }

  k = ncol(D); # number of cell types

  common.gene = intersect(names(Y), rownames(D))
  common.gene = intersect(common.gene, colnames(Sigma.ct))
  if(length(common.gene)< 0.1*min(length(Y), nrow(D), ncol(Sigma.ct))){
    stop('Not enough common genes!')
  }
  Y = Y[match(common.gene, names(Y))];
  D = D[match(common.gene, rownames(D)), ]
  Sigma.ct = Sigma.ct[, match(common.gene, colnames(Sigma))]

  X = D
  ## First, no intercept and no normalization
  if(centered){
    X = X - mean(X)
    Y = Y - mean(Y)
  }
  if(normalize){
    X = X/sd(as.vector(X));
    S = S*sd(as.vector(X));
    Y = Y/sd(Y)
  }else{
    Y = Y*100
  }
  lm.D = music.basic.ct(Y, X, S, Sigma.ct, iter.max = iter.max, nu = nu, eps = eps)
  return(lm.D)
}
