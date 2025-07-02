# Utility Function
#
#
# Author: Jiaxin Fan, Xuran Wang
###################################################
#' @title MuSiC2_Deconvolution
#'
#' @description This function is used to deconvolve bulk RNA-seq data using single-cell reference generated under a different condition.
#' @param bulk.control.mtx Matrix of expression for bulk data, control group
#' @param bulk.case.mtx Matrix of expression for bulk data, case group
#' @param sc.sce SingleCellExperiment for single cell data
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character, the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by the single cell dataset;
#' @param expr_low numeric, cutoff for defining lowly expressed genes in bulk data. Genes with mean expression across samples in bulk data < expr_low will be excluded from cell-type-specific DE gene detection. Default is 20;
#' @param prop_r numeric, cutoff for defining rare cell types. Cell types with mean proportion across samples in bulk data < prop_r will be characterized as rare cell types. Otherwise, will be characterized as common cell types. Default is 0.1;
#' @param eps_c numeric, convergence cutoff for common cell types. The cell type proportion estimate is converged if absolute relative change of proportion estimates for the current iteration against the previous iteration < eps_c. Default is 0.05;
#' @param eps_r numeric, convergence cutoff for rare cell types. The cell type proportion estimate is converged if absolute change of proportion estimates for the current iteration against the previous iteration < eps_r. Default is 0.01;
#' @param n_resample numeric, number of resamples used for detecting cell-type-specific DE genes. Default is 20;
#' @param sample_prop numeric, proportion of samples to be randomly sampled without replacement under each condition at each resampling. Default is 0.5;
#' @param cutoff_expr numeric, cutoff for defining lowly expressed genes over resamples. Genes with average cell-type-specific expression calculated over all resamples in the lower cutoff_expr quantile are excluded from cell-type-specific DE gene detection. Default is 0.05;
#' @param cutoff_c numeric, cutoff for defining cell-type-specific DE genes for common cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_c quantile are considered as cell-type-specific DE genes. Default is 0.05;
#' @param cutoff_r numeric, cutoff for defining cell-type-specific DE genes for rare cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_r quantile are considered as cell-type-specific DE genes. Default is 0.01;
#' @param maxiter numeric, maximum number of iterations. Default is 200;
#' @param markers vector or list of gene names. Default as NULL, i.e., use all genes that provided by both bulk and single cell datasets;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param centered logic, subtract avg of Y and D;
#' @param normalize logic, divide Y and D by their standard deviation;
#' @return If MuSiC2 converges, return:
#' \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {n.iter: numeric, number of iterations.}
#'    \item {DE.genes: vector, cell-type-specific DE genes being removed.}
#'    }
#'  Or if MuSiC2 does not converge, return:
#'  \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {id.not.converge: vector, sample ids that failed to converge.}
#'    }
#' @seealso
#' \code{\link{music_prop}}
#' @export

music2_prop = function(bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, samples, select.ct, expr_low=20, 
                       prop_r=0.1, eps_c=0.05, eps_r=0.01, n_resample=20, sample_prop=0.5,cutoff_expr=0.05, 
                       cutoff_c=0.05, cutoff_r=0.01, maxiter = 200, markers = NULL, cell_size = NULL, ct.cov = FALSE, 
                       centered = FALSE, normalize = FALSE){
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if(length(gene.bulk) < 0.1*min(nrow(bulk.control.mtx), nrow(bulk.case.mtx))){
    stop('Not enough genes for bulk data! Please check gene annotations.')
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, ])
  
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if(length(gene_all) < 0.2*min(length(gene.bulk), nrow(sc.sce))){
    stop('Not enough genes between bulk and single-cell data! Please check gene annotations.')
  }
  
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  
  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr>=expr_low])
  
  # Analyse separately based on their disease status
  
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)];
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)];
  
  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  # tutorial for basic music see https://xuranw.github.io/MuSiC/articles/MuSiC.html
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce,
                            clusters = clusters, samples = samples, select.ct = select.ct,
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                            nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case.mtx, sc.sce = sc.sce,
                             clusters = clusters, samples = samples, select.ct = select.ct,
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                             nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control,prop_CASE)
  
  # start iteration
  iter=1
  ncell = length(select.ct)
  id_conv = NULL
  while(iter <= maxiter){
    # print(iter)
    # step 2: identify cell-type-specific DE genes
    # calculate mean/ sd of log fold change as an indicator of differential expression
    LOGFC = NULL
    MOD0=MOD1=matrix(0L, nrow = ncell, ncol = length(exp_genel))
    
    for(i in 1:n_resample){
      id_h = sample(colnames(bulk.control), round(ncol(bulk.control)*sample_prop))
      control_s = bulk.control[exp_genel, colnames(bulk.control) %in% id_h]
      prop_h = prop_control[colnames(control_s),]
      
      mod0 = apply(control_s, 1, function(x){
        mod = nnls(prop_h,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD0=MOD0+mod0
      
      id_d = sample(colnames(bulk.case),round(ncol(bulk.case)*sample_prop))
      case_s = bulk.case[exp_genel, colnames(bulk.case) %in% id_d]
      prop_d = prop_CASE[colnames(case_s),]
      
      mod1 = apply(case_s, 1, function(x){
        mod = nnls(prop_d,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD1=MOD1+mod1
      LOGFC = rbind(LOGFC,log1p(mod1)-log1p(mod0))
    }
    
    rcv=NULL
    for(i in 1:ncell){
      s = LOGFC[seq(from=i,to=nrow(LOGFC), by=ncell),]
      rcv = rbind(rcv,apply(s,2,function(x){ifelse(mean(x)==0, 0, mean(x)/sd(x))}))
    }
    abs_rcv_logfc = abs(rcv)
    MOD0 = MOD0/n_resample
    MOD1 = MOD1/n_resample
    rownames(MOD0)=rownames(MOD1)=rownames(abs_rcv_logfc)=select.ct
    
    # average cell type proportion
    mex = apply(prop_all,2,mean)
    lr = NULL
    for(celltype in select.ct){
      m = mex[celltype]
      rh = MOD0[celltype,]
      rd = MOD1[celltype,]
      # for genes with average expression within lower cutoff_expr for both conditions are removed from cell-type-specific DE genes detection
      llr = unique(intersect(names(rd[rd <= quantile(rd,prob=cutoff_expr)]),
                             names(rh[rh <= quantile(rh,prob=cutoff_expr)])))
      x = abs_rcv_logfc[celltype,]
      x = x[!names(x) %in% llr]
      # select genes with large mean/cv of log fc as DE
      if(m >= prop_r){
        lr = c(lr, names(x[x >= quantile(x,prob=1-cutoff_c)]))
      }else{
        lr = c(lr, names(x[x >= quantile(x,prob=1-cutoff_r)]))
      }
    }
    lr = unique(lr)
    
    # step 3: update sc gene list
    # remove identified DE genes from sc rna-seq data
    l = setdiff(gene_all,lr)
    sc.iter.sce = sc.sce[l,]
    
    # step 1: update cell type proportion based on new gene list
    if(length(id_conv)>0){
      case_sample = bulk.case[ , !colnames(bulk.case) %in% id_conv]
    }else{
      case_sample = bulk.case
    }
    
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce,
                           clusters=clusters, samples=samples, select.ct=select.ct,
                           markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                           nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize,verbose = F)$Est.prop.weighted
    
    prop_CASE = rbind(prop_case,prop_case_fix)
    
    if(length(id_conv)==1){
      rownames(prop_CASE) = c(rownames(prop_case),id_conv)
    }
    
    prop_all = rbind(prop_control,prop_CASE)
    
    # check convergence, by cell type
    prop_case=prop_case[rownames(prop_case_ini),]
    pc = abs(prop_case-prop_case_ini)
    conv = pc
    conv[,] = 1
    # use difference if rare cell type
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    # use percent change if common cell type
    pc[prop_case_ini > prop_r] = pc[prop_case_ini>prop_r]/prop_case_ini[prop_case_ini>prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > prop_r] < eps_c,0,1)
    convf = apply(conv,1,function(x){all(x==0)})
    
    # if an id converged, not updating anymore
    all_converge=FALSE
    id_conv = c(id_conv,names(convf[convf==TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% id_conv,]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv,]
    
    # if all converged or if only one subjects not converging--> music2 converged
    if(is.vector(prop_case_ini)){
      all_converge=TRUE
      break
    }else if(nrow(prop_case_ini)==0){
      all_converge=TRUE
      break
    }
    iter=iter+1
  }
  # return estimated proportion
  if(all_converge){
    return(list('Est.prop' = prop_all,'convergence'=TRUE,'n.iter'=iter,'DE.genes'=lr))}
  else{
    return(list('Est.prop' = prop_all,'convergence'=FALSE,'id.not.converge'=rownames(prop_case_ini)))}
}




#' MuSiC2 Deconvolution with T Statistics
#' @description
#' This function is used to deconvolve bulk RNA-seq data using single-cell reference generated under a different clinical condition. Cell-type-specific differentially expressed (DE) genes are defined using a T statistics, calculated as the ratio of the mean to the standard deviation of the log fold change of cell-type-specific expression between conditions over all resampling iterations. 
#' 
#' @param bulk.control.mtx Matrix of expression for bulk data, control group
#' @param bulk.case.mtx Matrix of expression for bulk data, case group
#' @param sc.sce SingleCellExperiment for single cell data
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character, the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by the single cell dataset;
#' @param expr_low numeric, cutoff on gene expression of the bulk data. Genes with mean expression across samples in bulk data < expr_low will be excluded from cell-type-specific DE gene detection. Default is 20;
#' @param prop_r numeric, cutoff on cell type proportions for defining rare cell types. Cell types with mean proportion across samples in bulk data < prop_r will be characterized as rare cell types. Otherwise, will be characterized as common cell types. Default is 0.1;
#' @param eps_c numeric, convergence cutoff for common cell types. The cell type proportion estimate is converged if absolute relative change of proportion estimates for the current iteration against the previous iteration < eps_c. Default is 0.05;
#' @param eps_r numeric, convergence cutoff for rare cell types. The cell type proportion estimate is converged if absolute change of proportion estimates for the current iteration against the previous iteration < eps_r. Default is 0.01;
#' @param n_resample numeric, number of resamples used for detecting cell-type-specific DE genes. Default is 20;
#' @param sample_prop numeric, proportion of samples to be randomly sampled without replacement under each condition at each resampling iteration. Default is 0.5;
#' @param cutoff_expr numeric, cutoff on gene expression over resamples. Genes with average cell-type-specific expression calculated over all resamples in the lower cutoff_expr quantile are excluded from cell-type-specific DE gene detection. Default is 0.05;
#' @param cutoff_fc numeric, cutoff on fold change over resamples. Genes with absolute value of the mean fold change calculated over all resamples < cutoff_fc are excluded from cell-type-specific DE gene detection. Default is 2;
#' @param cutoff_c numeric, cutoff on T statistics for defining cell-type-specific DE genes for common cell types. Genes with the value of T statistic in the upper cutoff_c quantile are considered as cell-type-specific DE genes. Default is 0.05;
#' @param cutoff_r numeric, cutoff on T statistics for defining cell-type-specific DE genes for rare cell types. Genes with the value of T statistic in the upper cutoff_r quantile are considered as cell-type-specific DE genes. Default is 0.01;
#' @param maxiter numeric, maximum number of iterations. Default is 200;
#' @param markers vector or list of gene names. Default as NULL, i.e., use all genes that provided by both bulk and single cell datasets;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param centered logic, subtract avg of Y and D;
#' @param normalize logic, divide Y and D by their standard deviation;
#' 
#' @return If MuSiC2 converges, return:
#' \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {n.iter: numeric, number of iterations.}
#'    \item {DE.genes: vector, cell-type-specific DE genes being removed.}
#'    }
#'  Or if MuSiC2 does not converge, return:
#'  \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {id.not.converge: vector, sample ids that failed to converge.}
#'    }
#'    
#' @seealso
#' \code{\link{music_prop}}
#' 
#' @export

music2_prop_t_statistics = function(bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, samples, select.ct, expr_low=20, 
                                    prop_r=0.1, eps_c=0.05, eps_r=0.01, n_resample=20, sample_prop=0.5,cutoff_expr=0.05, 
                                    cutoff_fc=2, cutoff_c=0.05, cutoff_r=0.01, maxiter = 200, markers = NULL, 
                                    cell_size = NULL, ct.cov = FALSE, centered = FALSE, normalize = FALSE){
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if(length(gene.bulk) < 0.1*min(nrow(bulk.control.mtx), nrow(bulk.case.mtx))){
    stop('Not enough genes for bulk data! Please check gene annotations.')
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, ])
  
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if(length(gene_all) < 0.2*min(length(gene.bulk), nrow(sc.sce))){
    stop('Not enough genes between bulk and single-cell data! Please check gene annotations.')
  }
  
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  
  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr>=expr_low])
  
  # Analyse separately based on their disease status
  
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)];
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)];
  
  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce,
                            clusters = clusters, samples = samples, select.ct = select.ct,
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                            nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx =bulk.case, sc.sce = sc.sce,
                             clusters = clusters, samples = samples, select.ct = select.ct,
                             markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                             nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control, prop_CASE)
  
  # start iteration
  iter=1
  ncell = length(select.ct)
  id_conv = NULL
  while(iter <= maxiter){
    # print(iter)
    # step 2: identify cell-type-specific DE genes
    # calculate mean/ sd of log fold change as an indicator of differential expression
    LOGFC = NULL
    MOD0=MOD1=matrix(0L, nrow = ncell, ncol = length(exp_genel))
    
    for(i in 1:n_resample){
      id_h = sample(colnames(bulk.control),round(ncol(bulk.control)*sample_prop))
      control_s = bulk.control[exp_genel, colnames(bulk.control) %in% id_h]
      prop_h = prop_control[colnames(control_s),]
      
      mod0 = apply(control_s, 1, function(x){
        mod = nnls(prop_h,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD0=MOD0+mod0
      
      id_d = sample(colnames(bulk.case),round(ncol(bulk.case)*sample_prop))
      case_s = bulk.case[exp_genel, colnames(bulk.case) %in% id_d]
      prop_d = prop_CASE[colnames(case_s),]
      
      mod1 = apply(case_s, 1, function(x){
        mod = nnls(prop_d,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD1=MOD1+mod1
      LOGFC = rbind(LOGFC,log1p(mod1)-log1p(mod0))
    }
    
    rcv=NULL
    mfc=NULL
    for(i in 1:ncell){
      s = LOGFC[seq(from=i,to=nrow(LOGFC), by=ncell),]
      rcv = rbind(rcv,apply(s,2,function(x){ifelse(mean(x)==0, 0, mean(x)/sd(x))}))
      mfc = rbind(mfc,apply(s,2,function(x){mean(x)}))
    }
    abs_rcv_logfc = abs(rcv)
    MOD0 = MOD0/n_resample
    MOD1 = MOD1/n_resample
    rownames(MOD0)=rownames(MOD1)=rownames(abs_rcv_logfc)=rownames(mfc)=select.ct
    
    # average cell type proportion
    mex = apply(prop_all,2,mean)
    lr = NULL
    for(celltype in select.ct){
      m = mex[celltype]
      rh = MOD0[celltype,]
      rd = MOD1[celltype,]
      fc = mfc[celltype,]
      # for genes with average expression within lower cutoff_expr for both conditions are removed from cell-type-specific DE genes detection
      llr = unique(intersect(names(rd[rd <= quantile(rd,prob=cutoff_expr)]),
                             names(rh[rh <= quantile(rh,prob=cutoff_expr)])))
      x = abs_rcv_logfc[celltype,]
      x = x[!names(x) %in% llr]
      # select genes with large mean/SD of log fc and with at least moderate log fc as DE
      if(m >= prop_r){
        lr = c(lr, intersect(names(x[x >= quantile(x,prob=1-cutoff_c)]), names(fc[abs(fc)>=log1p(cutoff_fc)])))
      }else{
        lr = c(lr, intersect(names(x[x >= quantile(x,prob=1-cutoff_r)]), names(fc[abs(fc)>=log1p(cutoff_fc)])))
      }
    }
    lr = unique(lr)
    
    # step 3: update sc gene list
    # remove identified DE genes from sc rna-seq data
    l = setdiff(gene_all,lr)
    sc.iter.sce = sc.sce[l,]
    
    # step 1: update cell type proportion based on new gene list
    if(length(id_conv)>0){
      case_sample = bulk.case[ , !colnames(bulk.case) %in% id_conv]
    }else{
      case_sample = bulk.case
    }
    
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce,
                           clusters=clusters, samples=samples, select.ct=select.ct,
                           markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                           nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize,verbose = F)$Est.prop.weighted
    
    prop_CASE = rbind(prop_case, prop_case_fix)
    
    if(length(id_conv)==1){
      rownames(prop_CASE) = c(rownames(prop_case), id_conv)
    }
    
    prop_all = rbind(prop_control,prop_CASE)
    
    # check convergence, by cell type
    prop_case=prop_case[rownames(prop_case_ini),]
    pc = abs(prop_case-prop_case_ini)
    conv = pc
    conv[,] = 1
    # use difference if rare cell type
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    # use percent change if common cell type
    pc[prop_case_ini > prop_r] = pc[prop_case_ini>prop_r]/prop_case_ini[prop_case_ini>prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > prop_r] < eps_c,0,1)
    convf = apply(conv,1,function(x){all(x==0)})
    
    # if an id converged, not updating anymore
    all_converge=FALSE
    id_conv = c(id_conv,names(convf[convf==TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% id_conv,]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv,]
    
    # if all converged or if only one subjects not converging--> music2 converged
    if(is.vector(prop_case_ini)){
      all_converge=TRUE
      break
    }else if(nrow(prop_case_ini)==0){
      all_converge=TRUE
      break
    }
    iter=iter+1
  }
  # return estimated proportion
  if(all_converge){
    return(list('Est.prop' = prop_all,'convergence'=TRUE,'n.iter'=iter,'DE.genes'=lr))}
  else{
    return(list('Est.prop' = prop_all,'convergence'=FALSE,'id.not.converge'=rownames(prop_case_ini)))}
}




#' MuSiC2 Deconvolution with TOAST
#' @description This function is used to deconvolve bulk RNA-seq data using single-cell reference generated under 
#' a different clinical condition. Cell-type-specific DE genes are defined using TOAST (Li et al. 2019 
#'  \href{https://academic.oup.com/bioinformatics/article/35/20/3898/5418952}{Bioinformatics}) 
#'  based on cell-type-specific FDR adjusted p-values.
#' 
#' @param bulk.control.mtx Matrix of expression for bulk data, control group
#' @param bulk.case.mtx Matrix of expression for bulk data, case group
#' @param sc.sce SingleCellExperiment for single cell data
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character, the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by the single cell dataset;
#' @param expr_low numeric, cutoff on gene expression of the bulk data. Genes with mean expression across samples in bulk data < expr_low will be excluded from cell-type-specific DE gene detection. Default is 20;
#' @param prop_r numeric, cutoff on cell type proportions for defining rare cell types. Cell types with mean proportion across samples in bulk data < prop_r will be characterized as rare cell types. Otherwise, will be characterized as common cell types. Default is 0.1;
#' @param eps_c numeric, convergence cutoff for common cell types. The cell type proportion estimate is converged if absolute relative change of proportion estimates for the current iteration against the previous iteration < eps_c. Default is 0.05;
#' @param eps_r numeric, convergence cutoff for rare cell types. The cell type proportion estimate is converged if absolute change of proportion estimates for the current iteration against the previous iteration < eps_r. Default is 0.01;
#' @param cutoff_c numeric, cutoff on FDR adjusted p-values for defining cell-type-specific DE genes for common cell types. Genes with FDR adjusted p-value <= cutoff_c are considered as cell-type-specific DE genes. Default is 10^(-3);
#' @param cutoff_r numeric, cutoff on FDR adjusted p-values for for defining cell-type-specific DE genes for rare cell types. Genes with FDR adjusted p-value <= cutoff_r are considered as cell-type-specific DE genes. Default is 10^(-3);
#' @param cap numeric, cutoff on maximum number of genes removed for each cell type. For each cell type, at the maximum, genes with FDR adjusted p-value within the lower cap quantile are removed. Default is 0.3;
#' @param maxiter numeric, maximum number of iterations. Default is 200;
#' @param markers vector or list of gene names. Default as NULL, i.e., use all genes that provided by both bulk and single cell datasets;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param centered logic, substract avg of Y and D;
#' @param normalize logic, divide Y and D by their standard deviation;
#' 
#' @return If MuSiC2 converges, return:
#' \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {n.iter: numeric, number of iterations.}
#'    \item {DE.genes: vector, cell-type-specific DE genes being removed.}
#'    }
#'  Or if MuSiC2 does not converge, return:
#'  \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {id.not.converge: vector, sample ids that failed to converge.}
#'    }
#'    
#' @seealso
#' \code{\link{music_prop}}
#' 
#' @export
#' @import TOAST

music2_prop_toast = function(bulk.control.mtx, bulk.case.mtx, sc.sce, clusters, samples, select.ct, expr_low=20, prop_r=0.1, 
                        eps_c=0.05, eps_r=0.01, cutoff_c=10^(-3), cutoff_r=10^(-3), cap=0.3, maxiter = 200){
  gene.bulk = intersect(rownames(bulk.control.mtx), rownames(bulk.case.mtx))
  if(length(gene.bulk) < 0.1*min(nrow(bulk.control.mtx), nrow(bulk.case.mtx))){
    stop('Not enough genes for bulk data! Please check gene annotations.')
  }
  bulk.mtx = cbind(bulk.control.mtx[gene.bulk, ], bulk.case.mtx[gene.bulk, ])
  Pheno = data.frame(condition = factor(c(rep('control', ncol(bulk.control.mtx)), rep('case', ncol(bulk.case.mtx))), 
                                        levels = c('control', 'case')))
  rownames(Pheno) = colnames(bulk.mtx)
  
  gene_all = intersect(gene.bulk, rownames(sc.sce))
  if(length(gene_all) < 0.2*min(length(gene.bulk), nrow(sc.sce))){
    stop('Not enough genes between bulk and single-cell data! Please check gene annotations.')
  }
  
  bulk.mtx = bulk.mtx[gene_all, ]
  sc.iter.sce = sc.sce[gene_all, ]
  
  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  expr = apply(bulk.mtx, 1, mean)
  exp_genel = names(expr[expr>=expr_low])
  
  # Analyse separately based on their disease status
  
  bulk.control = bulk.mtx[, colnames(bulk.control.mtx)];
  bulk.case = bulk.mtx[, colnames(bulk.case.mtx)];
  
  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  prop_control = music_prop(bulk.mtx = bulk.control, sc.sce = sc.sce, clusters = clusters, samples = samples, 
                            select.ct = select.ct, markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                            iter.max = 1000, nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, 
                            verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.mtx = bulk.case, sc.sce = sc.sce, clusters = clusters, samples = samples, 
                             select.ct = select.ct, markers = markers, cell_size = cell_size, ct.cov = ct.cov, 
                             iter.max = 1000, nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, 
                             verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control,prop_CASE)
  
  # start iteration
  iter=1
  ncell = length(select.ct)
  id_conv = NULL
  while(iter <= maxiter){
    # step 2: identify cell-type-specific DE genes using TOAST
    # first, log transformed the bulk expression
    Y_raw = log1p(bulk.mtx)
    #Pheno = pData(bulk.eset)
    design = Pheno
    Prop <- prop_all[rownames(Pheno),]
    Design_out <- makeDesign(design, Prop)
    fitted_model <- fitModel(Design_out, Y_raw)
    # run TOAST to detect DE between conditions
    res_table<- csTest(fitted_model, coef = "group", verbose = F)
    
    # average cell type proportion
    mex = apply(prop_all, 2, mean)
    lr = NULL
    for(celltype in select.ct){
      m = mex[celltype]
      DE = res_table[[celltype]]
      pval = DE$fdr
      names(pval)=rownames(DE)
      pval = pval[names(pval) %in% exp_genel]
      # select DE genes 
      if(m >= prop_r){
        lr = c(lr, names(pval[pval <= cutoff_c & pval <= quantile(pval, prob=cap)]))
      }else{
        lr = c(lr, names(pval[pval <= cutoff_r & pval <= quantile(pval, prob=cap)]))
      }
    }
    lr = unique(lr)
    
    # step 3: update sc gene list
    # remove identified DE genes from sc rna-seq data
    l = setdiff(gene_all,lr)
    sc.iter.sce = sc.sce[l,]
    
    # step 1: update cell type proportion based on new gene list
    if(length(id_conv)>0){
      case_sample = bulk.case[ , !colnames(bulk.case) %in% id_conv]
    }else{
      case_sample = bulk.case
    }
    
    prop_case = music_prop(bulk.mtx = case_sample, sc.sce = sc.iter.sce,
                           clusters=clusters, samples=samples, select.ct=select.ct,verbose = F)$Est.prop.weighted
    
    prop_CASE = rbind(prop_case, prop_case_fix)
    
    if(length(id_conv)==1){
      rownames(prop_CASE) = c(rownames(prop_case),id_conv)
    }
    
    prop_all = rbind(prop_control, prop_CASE)
    
    # check convergence, by cell type
    prop_case=prop_case[rownames(prop_case_ini),]
    pc = abs(prop_case - prop_case_ini)
    conv = pc
    conv[,] = 1
    # use difference if rare cell type
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    # use percent change if common cell type
    pc[prop_case_ini > prop_r] = pc[prop_case_ini>prop_r]/prop_case_ini[prop_case_ini>prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > prop_r] < eps_c,0,1)
    convf = apply(conv,1,function(x){all(x==0)})
    
    # if an id converged, not updating anymore
    all_converge=FALSE
    id_conv = c(id_conv,names(convf[convf==TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% id_conv,]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv,]
    
    # if all converged or if only one subjects not converging--> music2 converged
    if(is.vector(prop_case_ini)){
      all_converge=TRUE
      break
    }else if(nrow(prop_case_ini)==0){
      all_converge=TRUE
      break
    }
    iter=iter+1
  }
  # return estimated proportion
  if(all_converge){
    return(list('Est.prop' = prop_all,'convergence'=TRUE,'n.iter'=iter,'DE.genes'=lr))}
  else{
    return(list('Est.prop' = prop_all,'convergence'=FALSE,'id.not.converge'=rownames(prop_case_ini)))}
}
