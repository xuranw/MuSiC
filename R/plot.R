# Plotting functions
#
# Author: Xuran Wang
############################################################################################################

#' Convert list of real and estimated cell type proportions to data frame
#'
#' This is a function for converting real and estimated cell type proportions to data frame
#'
#' @param prop.real, real cell type proportions
#' @param prop.est, estimated cell type proportions
#' @param method.name vector of estimation methods, default as NULL, name will be provided as 'Est.Method1', 'Est.Method2', ...
#' @return data.frame of real and estimated cell type proportion, cell type, subject name.
#'
#' @export
Prop_convert = function(prop.real, prop.est, method.name = NULL, ...){
  ct.real = colnames(prop.real); sub.real = rownames(prop.real);
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est)
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub.est = Reduce(intersect, l.sub.est)
  ct.est = Reduce(intersect, l.ct.est)
  celltype = intersect(ct.real, ct.est); sub = intersect(sub.real, sub.est);
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
  m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, each = N), levels = celltype),
                           Sub = factor(rep(sub,K), levels = sub), Method = rep('Real', N*K) )
  m.prop.est = NULL;
  for(l in 1:L){
    m.prop.temp = m.prop.real
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp$Prop = c( cm.prop.est )
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    m.prop.est = rbind(m.prop.est, m.prop.temp)
  }
  m.prop = rbind(m.prop.real, m.prop.est);
  return(m.prop)
}


#' Plot heatmap of real and estimated cell type proportions
#'
#' Generate heatmaps of real and estimated cell type proportions side by side. Pearson correlation can be calculated and printed.
#'
#' @param prop.real a matrix of real cell type proportions
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estimation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the title of the plot
#' @param eval logical, default as TRUE. If TRUE, pearson correlation evaluation will be printed on the figure.
#'
#' @return a 'ggplot' object with [ggplot2::geom_tile]
#'
#' @seealso
#' \code{\link{Abs_diff_multi}}, \code{\link{Eval_multi}},
#' \code{\link{Scatter_multi}}
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' # generate random data
#'
Prop_comp_multi = function(prop.real, prop.est, method.name = NULL, title = NULL, eval = TRUE, ... ){
  ct.real = colnames(prop.real); sub.real = rownames(prop.real);
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est)
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub.est = Reduce(intersect, l.sub.est)
  ct.est = Reduce(intersect, l.ct.est)
  celltype = intersect(ct.real, ct.est); sub = intersect(sub.real, sub.est);
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
  m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, each = N), levels = celltype),
                           Sub = factor(rep(sub,K), levels = sub), Method = rep('Real', N*K) )
  m.prop.est = NULL;
  ann = data.frame(metric='',Method='Real')
  for(l in 1:L){
    m.prop.temp = m.prop.real
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp$Prop = c( cm.prop.est )
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    # ann = c(ann, paste0('R = ', round(cor(c(cm.prop.real), c(cm.prop.est)), digits = 4)))
    ann <- rbind(ann,data.frame(metric=paste0('R = ', round(cor(c(cm.prop.real), c(cm.prop.est)), digits = 4)),Method=method.name[l]))
    m.prop.est = rbind(m.prop.est, m.prop.temp)
  }
  m.prop = rbind(m.prop.real, m.prop.est);
  if(is.null(title)){
    title = 'Heatmap of Estimated and Real Proportion'
  }
  if(eval){
    ggplot(m.prop, aes(CellType, Sub)) + geom_tile(aes(fill = Prop), colour = 'white') + scale_fill_gradient2(
      low = 'steelblue', high = "red", mid = 'white', midpoint = 0.5, limit = c(0, 1), name = 'Est Prop\n')  + facet_wrap(~Method, nrow = 1) + 
      theme_minimal() + theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0) ) + ggtitle(title) + 
      geom_text(data = ann,aes(x = round(4*K/5), y = N,label=metric),  size = 2.5, colour = "black")
  }else{
    ggplot(m.prop, aes(CellType, Sub)) + geom_tile(aes(fill = Prop), colour = 'white') + scale_fill_gradient2(
      low = 'steelblue', high = "red", mid = 'white', midpoint = 0.5, limit = c(0, 1), name = 'Est Prop\n')  + facet_wrap(~Method, nrow = 1) + 
      theme_minimal() + theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0) ) + 
      ggtitle(title)
  }
}


#' Plot heatmap of Absolute difference between estimated and real cell type proportions
#'
#' Generate a series of heatmap of absolute difference between estimated and real cell type proportion with option of printing RMSD and mAD between real and estimated cell type proportions
#'
#' @param prop.real a matrix of real cell type proportions
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the tile of the plot
#' @param eval logical, default as TRUE. If TRUE, RMSD and mAD will be printed with figures
#'
#' @return a 'ggplot' object with [ggplot2::geom_tile]
#'
#' @seealso
#' \code{\link{Prop_comp_multi}}, \code{\link{Eval_multi}},
#' \code{\link{Scatter_multi}}
#'
#' @import ggplot2
#' @export
#'
Abs_diff_multi = function(prop.real, prop.est, method.name = NULL, title = NULL, eval = TRUE, ... ){
  ct.real = colnames(prop.real); sub.real = rownames(prop.real);
  if(is.list(prop.est)){
    L = length(prop.est)
    if(is.null(method.name)){
      method.name = paste0('Est.Method', 1:L)
    }else{
      if(length(method.name) < L){
        method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
      }else{
        method.name = method.name[1:L]
      }
    }
    l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
    sub.est = Reduce(intersect, l.sub.est)
    ct.est = Reduce(intersect, l.ct.est)
    celltype = intersect(ct.real, ct.est); sub = intersect(sub.real, sub.est);
    N = length(sub); K = length(celltype);
    if(N<1){
      stop("No common Subjects! Check rowname!")
    }
    if(K <=1 ){
      stop("Not enough cell types!")
    }
    cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
    m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, each = N), levels = celltype),
                             Sub = factor(rep(sub,K), levels = sub), Method = rep('Real', N*K) )
    
    abs.diff = NULL
    ann <- data.frame(metric = character(0), Method = character(0))
    for(l in 1:L){
      abs.diff.temp = m.prop.real; colnames(abs.diff.temp)[1] = 'Abs.Diff';
      cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
      abs.diff.temp$Abs.Diff = abs( c( cm.prop.est ) - m.prop.real$Prop)
      abs.diff.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
      ann <- rbind(ann, data.frame(metric = paste0("RMSD = ", round(sqrt(mean((cm.prop.real - cm.prop.est)^2)), digits = 5), "\n mAD = ",
                                                   round(mean(abs.diff.temp$Abs.Diff), digits = 5)), Method = method.name[l]))
      abs.diff = rbind(abs.diff, abs.diff.temp)
    }
  }else{
    ct.est = colnames(prop.est); sub.est = rownames(prop.est);
    celltype = intersect(ct.real, ct.est); sub = intersect(sub.real, sub.est);
    N = length(sub); K = length(celltype);
    if(N<1){
      stop("No common Subjects! Check rowname!")
    }
    if(K <=1 ){
      stop("Not enough cell types!")
    }
    cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
    cm.prop.est = prop.est[match(sub, sub.real), match(celltype, ct.real)]
    m.prop.real = data.frame(Prop = c(cm.prop.real), CellType = factor(rep(celltype, each = N), levels = celltype), Sub = factor(rep(sub,K), levels = sub), Method = rep('Real', N*K) )
    abs.diff = m.prop.real; colnames(abs.diff)[1] = 'Abs.Diff';
    abs.diff$Abs.Diff = abs( c(cm.prop.est-cm.prop.real) );
    abs.diff$Method = rep(method.name, K*N);
    ann =paste0('RMSD = ', round( sqrt(mean((cm.prop.real-cm.prop.est)^2)), digits = 5 ),
                '\n mAD = ', round( mean( abs.diff$Abs.Diff), digits = 5 ));
  }
  if(is.null(title)){
    title = 'Heatmap of Absolute Difference |p - Est.p|'
  }
  if(eval){
    ggplot(abs.diff, aes(CellType, Sub)) + geom_tile(aes(fill = Abs.Diff), colour = 'white')+ scale_fill_gradient(
      low = 'white', high = 'steelblue', name = 'Abs.Diff\n') + facet_wrap( ~ Method, nrow = 1) + theme_minimal() + 
      theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0) ) + #geom_text(aes(label = round(Abs.Diff, 2))) +
      ggtitle(title) + 
      geom_text(data = ann, aes(x = round(4 * K / 5), y = N - 0.5, label = metric), size = 2.5, colour = "black")
  }else{
    ggplot(abs.diff, aes(CellType, Sub)) + geom_tile(aes(fill = Abs.Diff), colour = 'white')+ scale_fill_gradient(
      low = 'white', high = 'steelblue', name = 'Abs.Diff\n') + facet_wrap( ~ Method, nrow = 1) + theme_minimal() +
      theme(axis.text.x = element_text(angle = -90, size = 10, vjust = 0) ) + #geom_text(aes(label = round(Abs.Diff, 2))) +
      ggtitle(title)
  }
}

#' Evaluate estimation methods
#'
#' Calculate Pearson correlation, RMSE and mAD between real and estimated cell type proportions
#'
#' @param prop.real a matrix of real cell type proportions
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param by.subject logical, default is FALSE. If TRUE, Pearson correlation is estimated subject by subject.
#' @param select.ct a vector of cell types selected for evaluation. Default is NULL. If select.ct is NULL, we evaluate
#' real and estimated cell type proportion with all common cell types. Otherwise use only selected cell types.
#'
#' @return a matrix of evaluation
#'
#' @export
Eval_multi = function(prop.real, prop.est, method.name = NULL, by.subject = FALSE, select.ct = NULL, ... ){
  ct.real = colnames(prop.real); sub.real = rownames(prop.real);
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est)
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub.est = Reduce(intersect, l.sub.est)
  ct.est = Reduce(intersect, l.ct.est)
  celltype = intersect(ct.real, ct.est);
  if(!is.null(select.ct)){
    celltype = intersect(celltype, select.ct)
  }
  sub = intersect(sub.real, sub.est);
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
  cm.prop.real = relative.ab(cm.prop.real, by.col = F)
  Eval.matrix = NULL;
  for(l in 1:L){
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    cm.prop.est = relative.ab(cm.prop.est, by.col = F);
    RMSD = round( sqrt(mean((cm.prop.real - cm.prop.est)^2)), digits = 5);
    mAD = round( mean( abs(cm.prop.real - cm.prop.est) ) , digits = 5);
    Pearson = round(cor(c(cm.prop.real), c(cm.prop.est)), digits = 4);
    
    if(by.subject){
      Pearson.by.subject = round(mean( sapply(1:N, function(i){
        cor(cm.prop.real[i, ], cm.prop.est[i, ])}) ), digits = 4)
      Eval.matrix = rbind(Eval.matrix, c(RMSD,mAD, Pearson, Pearson.by.subject))
    }else{
      Eval.matrix = rbind(Eval.matrix, c(RMSD,mAD, Pearson))
    }
  }
  rownames(Eval.matrix) = method.name
  if(by.subject){
    colnames(Eval.matrix) = c('RMSD', 'mAD', 'R', 'R.by.sub')
  }else{
    colnames(Eval.matrix) = c('RMSD', 'mAD', 'R')
  }
  return(Eval.matrix)
}



#' Plot Scatter plot of real and estimated cell type proportions
#'
#' Generate a series of scatter plot for each cell type between real and estimated cell type
#' proportions
#'
#' @param prop.real a matrix of real cell type proportions
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the tile of the plot
#' @param oneline logical, default as FALSE. If TRUE, set \code{nrow = 1} in \code{facet_wrap()}, all plots are in one row.
#'
#' @return a 'ggplot' object with [ggplot2::geom_point]
#'
#' @seealso
#' \code{\link{Prop_comp_multi}}, \code{\link{Abs_diff_multi}},
#' \code{\link{Eval_multi}}
#'
#' @import ggplot2
#' @export
#'
Scatter_multi = function(prop.real, prop.est, method.name = NULL, title = NULL, oneline = FALSE, ... ){
  ct.real = colnames(prop.real); sub.real = rownames(prop.real);
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est);
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub.est = Reduce(intersect, l.sub.est)
  ct.est = Reduce(intersect, l.ct.est)
  celltype = intersect(ct.real, ct.est); sub = intersect(sub.real, sub.est);
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  cm.prop.real = prop.real[match(sub, sub.real), match(celltype, ct.real)]
  m.prop.real = data.frame(Real.prop = c(cm.prop.real), CellType = factor(rep(celltype, each = N), levels = celltype),
                           Sub = factor(rep(sub,K), levels = sub))
  m.prop.est = NULL;
  for(l in 1:L){
    m.prop.temp = m.prop.real
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp$Est.prop = c( cm.prop.est )
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    m.prop.est = rbind(m.prop.est, m.prop.temp)
  }
  if(is.null(title)){
    title = 'Scatter plot of Real and Est. Proportion'
  }
  if(oneline){
    ggplot(m.prop.est, aes(Real.prop, Est.prop)) + geom_point(aes(color = Method, shape = Method)) + 
      facet_wrap(~CellType, scale = 'free', nrow = 1) + geom_abline(slope = 1, color = 'gray', size = 1) + 
      theme_minimal() + ggtitle(title)
  }else{
    ggplot(m.prop.est, aes(Real.prop, Est.prop)) + geom_point(aes(color = Method, shape = Method)) + 
      facet_wrap(~CellType, scale = 'free') + geom_abline(slope = 1, color = 'gray', size = 1) + 
      theme_minimal() + ggtitle(title)
  }
  #return(m.prop.est)
}

#' Boxplot of estimated cell type proportions
#'
#' Generate boxplot of estimated cell type proportions cell type by cell type
#'
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the tile of the plot
#'
#' @return a 'ggplot' object with [ggplot2::geom_boxplot]
#'
#' @seealso
#' \code{\link{Jitter_Est}}, \code{\link{Prop_heat_Est}}
#'
#' @import ggplot2
#' @export
Boxplot_Est = function(prop.est, method.name = NULL, title = NULL, ... ){
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est);
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub = Reduce(intersect, l.sub.est)
  celltype = Reduce(intersect, l.ct.est)
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  m.prop.est = NULL
  for(l in 1:L){
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp = data.frame(Prop = c(cm.prop.est),  CellType = factor(rep(celltype, each = N), levels = celltype),
                             Sub = factor(rep(sub,K), levels = sub));
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    m.prop.est = rbind(m.prop.est, m.prop.temp);
  }
  #return(m.prop.est)
  if(is.null(title)){
    title = 'Boxplot of Estimated Cell Type Proportions'
  }
  ggplot(m.prop.est, aes(Method, Prop))+ geom_boxplot(aes(color = Method)) + facet_wrap(~ CellType, scales = 'free') + 
    theme_minimal() + theme(axis.text.x=element_text(angle=30, size = 10, vjust=0.5)) + ggtitle(title)
}

#' Jitter plot of estimated cell type proportions
#'
#' Generate jitter plots of estimated cell type proportions cell type by cell type
#'
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the tile of the plot
#'
#' @return a 'ggplot' object with [ggplot2::geom_jitter]
#'
#' @seealso
#' \code{\link{Boxplot_Est}}, \code{\link{Prop_heat_Est}}
#'
#' @import ggplot2
#' @export
Jitter_Est = function(prop.est, method.name = NULL, title = NULL, ... ){
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est);
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub = Reduce(intersect, l.sub.est)
  celltype = Reduce(intersect, l.ct.est)
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  m.prop.est = NULL
  for(l in 1:L){
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp = data.frame(Prop = c(cm.prop.est),  CellType = factor(rep(celltype, each = N), levels = celltype),
                             Sub = factor(rep(sub,K), levels = sub));
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    m.prop.est = rbind(m.prop.est, m.prop.temp);
  }
  #return(m.prop.est)
  if(is.null(title)){
    title = 'Jitter plot of Estimated Cell Type Proportions'
  }
  ggplot(m.prop.est, aes(Method, Prop))+ geom_jitter(aes(color = Method), width = 0.2, height = 0) + facet_wrap(~ CellType, scales = 'free') + 
    theme_minimal() + theme(axis.text.x=element_text(angle=30, size = 10, vjust=0.5)) + ggtitle(title)
}

#' Heatmap of estimated cell type proportions
#'
#' Generate heatmap of estimated cell type proportions cell type by cell type
#'
#' @param prop.est a matrix or a list of matrices of estimated cell type proportions
#' @param method.name vector of the names of estmation methods. Default is NULL and the names will be
#' generated automatically as 'Est1', 'Est2', ...
#' @param title a string of the tile of the plot
#'
#' @return a 'ggplot' object with [ggplot2::geom_tile]
#'
#' @seealso
#' \code{\link{Boxplot_Est}}, \code{\link{Jitter_Est}}
#'
#' @import ggplot2
#' @export
Prop_heat_Est = function(prop.est, method.name = NULL, title = NULL, ... ){
  if(!is.list(prop.est)){
    prop.est = list(prop.est)
  }
  L = length(prop.est);
  if(is.null(method.name)){
    method.name = paste0('Est.Method', 1:L)
  }else{
    if(length(method.name) < L){
      method.name = c(method.name, paste0('Est.Method', 1:(L-length(method.name))))
    }else{
      method.name = method.name[1:L]
    }
  }
  l.sub.est = lapply(prop.est, rownames); l.ct.est = lapply(prop.est, colnames);
  sub = Reduce(intersect, l.sub.est)
  celltype = Reduce(intersect, l.ct.est)
  N = length(sub); K = length(celltype);
  if(N<1){
    stop("No common Subjects! Check rowname!")
  }
  if(K <=1 ){
    stop("Not enough cell types!")
  }
  m.prop.est = NULL
  for(l in 1:L){
    cm.prop.est = prop.est[[l]][match(sub, l.sub.est[[l]]), match(celltype, l.ct.est[[l]])]
    m.prop.temp = data.frame(Prop = c(cm.prop.est),  CellType = factor(rep(celltype, each = N), levels = celltype),
                             Sub = factor(rep(sub,K), levels = sub));
    m.prop.temp$Method = factor(rep(method.name[l], K*N), levels = method.name);
    m.prop.est = rbind(m.prop.est, m.prop.temp);
  }
  #return(m.prop.est)
  if(is.null(title)){
    title = 'Heatmap of Estimated Cell Type Proportions'
  }
  ggplot(m.prop.est, aes(CellType, Sub)) + geom_tile(aes(fill = Prop), colour = 'white') + scale_fill_gradient2(
    low = 'steelblue', high = "red", mid = 'white', midpoint = 0.1, limit = c(0, 1), name = 'Est Prop')  + facet_wrap(~Method, nrow = 1) + 
    theme_minimal() + theme(axis.text.x = element_text(angle = 50, size = 10, vjust = 0.5) ) + ggtitle(title)
}
############## Compare two datasets #######################
#' Compare cell type specific total expression (library size) between 2 dataset
#'
#'
CellTotal.df = function(sce, clusters, samples){
  df <- colData(sce)
  df$Total = Matrix::colSums(counts(sce))
  mdf <- t(data.matrix(sapply(unique(df[, clusters]), function(cl){
    msample <- sapply(unique(df[, samples]), function(sid){
      mean(df$Total[df[, clusters] == cl & df[, samples] == sid])
    })
    x = data.frame(avgtot = mean(msample), avgtot.sd = sd(msample))
    x$lb = x$avgtot - x$avgtot.sd
    x$ub = x$avgtot + x$avgtot.sd
    return(x)
  })))
  
  mdf = as.data.frame(mdf)
  mdf$cellType = unique(df[, clusters])
  rownames(mdf) = unique(df[, clusters])
  mdf = mdf[Matrix::rowSums(is.na(mdf)) == 0 ,]
  return(mdf)
}

#' Plot the cell type specific library size of 2 single cell datasets
#'
#' This function is to plot the cell type specific library size of the common cell types of 2 single cell datasets
#'
#' @param sce1 SingleCellExpression of the first single cell dataset
#' @param sce2 SingleCellExpression of the second single cell dataset
plotCellTotal.two = function(sce1, sce2, clusters = 'cellType', samples = 'sampleID', name1, name2){
  mdf1 = CellTotal.df(sce1, clusters, samples)
  mdf2 = CellTotal.df(sce2, clusters, samples)
  
  ids = intersect(mdf1$cellType, mdf2$cellType)
  mdf1 <- mdf1[match(ids, mdf1$cellType), ]
  mdf2 <- mdf2[match(ids, mdf2$cellType), ]
  mdf1 <- data.frame(avgtot = mdf1$avgtot, cellType = mdf1$cellType, lb = mdf1$lb, ub = mdf1$ub, Name = rep(name1, length(ids)))
  mdf2 <- data.frame(avgtot = mdf2$avgtot, cellType = mdf2$cellType, lb = mdf2$lb, ub = mdf2$ub, Name = rep(name2, length(ids)))
  mdf <- rbind(mdf1, mdf2)
  #mdf$Name = factor( rep(c(name1, name2), each = length(ids)), levels = c(name1, name2) )
  ggplot(mdf, aes(y = avgtot, x = cellType, fill = cellType)) +
    geom_bar(stat = "identity") + guides(fill = guide_legend("")) +
    geom_errorbar(aes_string(ymin = "lb", ymax = "ub"), color = "#555555",
                  width = 0.25) + facet_wrap( ~Name, scales = 'free') + theme_minimal() + 
    theme(axis.text.x = element_text(angle = 320, hjust = 0)) + ylab("Average total count") + xlab("")
}

#' Boxplot of relative abundance
#'
#' Generate boxplot of cell type specific relative abundance for each subjects
#'
#' @param sc.sce ExpressionSet for single cell data
#' @param gene.name character, for gene name
#' @param nu numeric regulator for log transformation
#' @param marker.id numeric indicator of cell type marker genes. Order in the select cell types.
#' Default is NULL. If NULL, do not print out the cell type.
#' @param log.trans logical, default as TRUE. If FALSE, do not take log transformation of relative abundance.
#' @param select.ct vector of characters. Default is NULL. If NULL use all cell types provided. Otherwise select
#' cell types provided.
#'
#' @return a 'ggplot' object with [ggplot2::geom_boxplot]
#' @importFrom Matrix rowSums
#' @import ggplot2
#' @export
#'
Relative_gene_boxplot = function(sc.sce, gene.name, nu = 10^{-10}, marker.id = NULL, log.trans = TRUE, select.ct = NULL, ... ){
  ## eliminate non expressed genes
  x <- sc.sce[Matrix::rowSums(counts(sc.sce))>0 , ]
  
  if(sum(rownames(x) %in% gene.name) != 1){
    stop(paste0('No such gene: ', gene.name, '!'))
  }
  
  if(is.null(select.ct)){
    sc.ra = relative.ab(counts(x))[gene.name, ]
    m.sc.ra = colData(sc.sce)
    m.sc.ra$Relative.Ab = sc.ra
  }else{
    sc.ra = relative.ab(counts(x))[gene.name, sc.sce$cellType %in% select.ct]
    m.sc.ra = colData(sc.sce)[sc.sce$cellType%in% select.ct, ]
    m.sc.ra$Relative.Ab = sc.ra
    m.sc.ra$cellType  = factor(m.sc.ra$cellType, levels = select.ct)
  }
  if(is.null(marker.id)){
    title.gg = paste0('Boxplot of relative abundance at ', gene.name)
  }else{
    cell.type.marker = levels(m.sc.ra$cellType)[marker.id]
    title.gg = paste0('Boxplot of relative abundance at ', gene.name, '\n marker for ', cell.type.marker, ' cell' )
  }
  if(log.trans){
    ggplot(m.sc.ra, aes(SubjectName, log(Relative.Ab + nu), color = cellType)) + theme_minimal() + 
      theme(axis.text.y=element_text(angle=20, size = 10, vjust=0.5)) + labs(y = paste0('log(Rel.Ab + ', nu, ')')) + geom_boxplot() + facet_wrap( ~ cellType, ncol = 1)+ coord_flip() + ggtitle(title.gg)
  }else{
    ggplot(m.sc.ra, aes(SubjectName, Relative.Ab, color = cellType)) + theme_minimal() +
      theme(axis.text.y=element_text(angle=20, size = 10, vjust=0.5)) + geom_boxplot() + facet_wrap( ~ cellType, ncol = 1)+ coord_flip() + ggtitle(title.gg)
  }
}
