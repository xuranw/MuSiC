MUlti-sample SIgle Cell deconvolution (MuSiC)
=============================================

`MuSiC` is a deconvolution method that utilizes cross-subject scRNA-seq to estimate cell type proportions in bulk RNA-seq data.
![MuSiC\_pipeline](image/FigureMethod.jpg)

How to cite `MuSiC`
-------------------
This work is not published yet, please see bioRxiv.

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
Please see vignette.