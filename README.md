Multi-sample Single Cell deconvolution (MuSiC)
=============================================

`MuSiC` is a deconvolution method that utilizes cross-subject scRNA-seq to estimate cell type proportions in bulk RNA-seq data.
![MuSiC\_pipeline](FigureMethod.jpg)

How to cite `MuSiC`
-------------------
This work is not published yet, please see [bioRxiv](https://www.biorxiv.org/content/early/2018/06/26/354944).

Installation
------------

``` r
# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)
```

More Information
-----------------
Please see [vignette](https://github.com/xuranw/MuSiC/blob/master/vignettes/vignette.Rmd).
