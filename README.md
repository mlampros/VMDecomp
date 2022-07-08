
[![tic](https://github.com/mlampros/VMDecomp/workflows/tic/badge.svg?branch=master)](https://github.com/mlampros/VMDecomp/actions)
[![docs: passing](https://img.shields.io/badge/docs-passing-success.svg)](https://mlampros.github.io/VMDecomp/reference/index.html)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/VMDecomp)](http://cran.r-project.org/package=VMDecomp)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/VMDecomp?color=blue)](http://www.r-pkg.org/pkg/VMDecomp)
[![](https://img.shields.io/docker/automated/mlampros/vmdecomp.svg)](https://hub.docker.com/r/mlampros/vmdecomp)
[![Dependencies](https://tinyverse.netlify.com/badge/VMDecomp)](https://cran.r-project.org/package=VMDecomp)
[![codecov.io](https://codecov.io/github/mlampros/VMDecomp/coverage.svg?branch=master)](https://codecov.io/github/mlampros/VMDecomp?branch=master)


## VMDecomp

<br>

The **VMDecomp** R package is the *RcppArmadillo* implementation of the ["Variational Mode Decomposition" Matlab code](https://math.montana.edu/dzosso/code/) in R. More details on the functionality of VMDecomp can be found in the package [Documentation](https://mlampros.github.io/VMDecomp/reference/index.html), [Vignette](https://mlampros.github.io/VMDecomp/articles/variatonal_mode_decomposition.html) and [blog post](http://mlampros.github.io/2022/06/11/variatonal_mode_decomposition/).

<br>

### Examples:

<br>

**Variational Mode Decomposition** (including residuals)

<br>

```R

require(VMDecomp)
data(arrhythmia)

alpha = 2000       # moderate bandwidth constraint
tau = 0            # noise-tolerance (no strict fidelity enforcement)
K = 9              # 9 modes
DC = FALSE         # no DC part imposed
init = 1           # initialize omegas uniformly
tol = 1e-6

vec_arrhythmia = arrhythmia[['MLII']]

set.seed(1)
arr_vmd = vmd(data = vec_arrhythmia, 
              alpha = alpha,
              tau = tau, 
              K = K, 
              DC = DC,
              init = init, 
              tol = tol,
              verbose = TRUE)
              
imfs = data.table::data.table(arr_vmd$u)
colnames(imfs) = glue::glue("IMF_{1:ncol(imfs)}")
imfs$residual = rowSums(imfs) - vec_arrhythmia

round(imfs, digits = 5)

#           IMF_1    IMF_2    IMF_3   IMF_4   IMF_5    IMF_6    IMF_7    IMF_8    IMF_9 residual
#     1: -0.06947  0.07831  0.13355 0.14031 0.10371  0.05622 -0.00143 -0.09686 -0.00629 -0.01194
#     2: -0.06971  0.07765  0.13199 0.13698 0.09920  0.05226 -0.00062 -0.07683 -0.01055 -0.00963
#     3: -0.07016  0.07639  0.12896 0.13047 0.09046  0.04475  0.00097 -0.04075 -0.01487 -0.00377
#     4: -0.07068  0.07468  0.12466 0.12114 0.07810  0.03447  0.00329  0.00427 -0.01331  0.00662
#     5: -0.07108  0.07273  0.11937 0.10947 0.06297  0.02246  0.00622  0.04929 -0.00211  0.01932
#    ---                                                                                        
#  9996: -0.07001 -0.13154 -0.24738 0.18826 0.03381 -0.07354  0.00076  0.03773  0.01426  0.01234
#  9997: -0.06980 -0.13333 -0.25498 0.21256 0.03400 -0.11214  0.05833 -0.00481  0.01450  0.00432
#  9998: -0.06951 -0.13452 -0.26056 0.23154 0.03414 -0.14316  0.10925 -0.04580  0.00791 -0.00570
#  9999: -0.06934 -0.13534 -0.26432 0.24447 0.03413 -0.16496  0.14706 -0.07813 -0.00137 -0.00780
# 10000: -0.06932 -0.13581 -0.26626 0.25095 0.03402 -0.17629  0.16708 -0.09601 -0.00797  0.00040

```

<br>

**Estimation of the K-modes Parameter** (*correlation threshold* of 0.1 and a *minimum K* of 2)

<br>

```R

require(VMDecomp)
data(arrhythmia)

default_vmd_params = list(alpha = 2000,
                          tau = 0,
                          DC = FALSE,
                          init = 1,
                          tol = 1e-6)

res_k = estimate_k_modes(signal_1d = arrhythmia[['MLII']],
                         cor_thresh = 0.1,
                         default_vmd_params = default_vmd_params,
                         min_K = 2,
                         seed = 1,
                         verbose = TRUE)
                         
# VMD based on a K of '2' will be computed ... 
# VMD based on a K of '3' will be computed ... 
# VMD based on a K of '4' will be computed ... 
# VMD based on a K of '5' will be computed ... 
# VMD based on a K of '6' will be computed ... 
# VMD based on a K of '7' will be computed ... 
# VMD based on a K of '8' will be computed ... 
# VMD based on a K of '9' will be computed ... 
# Optimal K parameter: '8'  Pre-specified correlation coefficient threshold: '0.1'
# Elapsed time: 0 hours and 1 minutes and 19 seconds.

res_k
# [1] 8

```

<br>

### Installation:

<br>

To install the package from CRAN use, 

```R
install.packages("VMDecomp")

```
<br>

and to download the latest version of the package from Github,

```R
remotes::install_github('mlampros/VMDecomp')

```

<br>

#### **Docker Image**

<br>

**Docker images** of the *VMDecomp* package are available to download from my [dockerhub](https://hub.docker.com/r/mlampros/vmdecomp) account. The images come with *Rstudio* and the *R-development* version (latest) installed. The whole process was tested on Ubuntu 18.04. To **pull** & **run** the image do the following,

<br>

```R

docker pull mlampros/vmdecomp:rstudiodev

docker run -d --name rstudio_dev -e USER=rstudio -e PASSWORD=give_here_your_password --rm -p 8787:8787 mlampros/vmdecomp:rstudiodev

```

<br>

The user can also **bind** a home directory / folder to the image to use its files by specifying the **-v** command,

<br>

```R

docker run -d --name rstudio_dev -e USER=rstudio -e PASSWORD=give_here_your_password --rm -p 8787:8787 -v /home/YOUR_DIR:/home/rstudio/YOUR_DIR mlampros/vmdecomp:rstudiodev


```

<br>

The **USER** defaults to *rstudio* but you have to give your **PASSWORD** of preference (see [www.rocker-project.org](https://www.rocker-project.org/) for more information).

<br>

Open your web-browser and depending where the docker image was *build / run* give, 

<br>

**1st. Option** on your personal computer,

<br>

```R
http://0.0.0.0:8787 

```

<br>

**2nd. Option** on a cloud instance, 

<br>

```R
http://Public DNS:8787

```

<br>

to access the Rstudio console in order to give your username and password.

<br>

### **Similar Projects:**

* https://github.com/vrcarva/vmdpy   (*Variational Mode Decomposition in Python*)
* https://github.com/helske/Rlibeemd  (*ensemble empirical mode decomposition (EEMD) and its complete variant (CEEMDAN)*) 

<br>

### **Citation:**

If you use the **VMDecomp** R package in your paper or research please cite both **VMDecomp** and the **original articles / software** `https://CRAN.R-project.org/package=VMDecomp`:

<br>

```R
@Manual{,
  title = {{VMDecomp}: Variational Mode Decomposition using R},
  author = {Lampros Mouselimis},
  year = {2022},
  note = {R package version 1.0.1},
  url = {https://CRAN.R-project.org/package=VMDecomp},
}
```

<br>
