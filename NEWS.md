
## VMDecomp 1.0.3

* I removed the `SystemRequirements` from the DESCRIPTION file.


## VMDecomp 1.0.2

* I updated the `Makevars` and `Makevars.win` files by adding `-DARMA_USE_CURRENT` (see issue: https://github.com/RcppCore/RcppArmadillo/issues/476)
* I removed the "CXX_STD = CXX11" from the "Makevars" files, and the "[[Rcpp::plugins(cpp11)]]" from the "vmd.cpp" file due to the following NOTE from CRAN, "NOTE Specified C++11: please drop specification unless essential" (see also: https://www.tidyverse.org/blog/2023/03/cran-checks-compiled-code/#note-regarding-systemrequirements-c11)


## VMDecomp 1.0.1

* I've added the *estimate_K()* function (Correlation Coefficient Method) which is based on the paper https://doi.org/10.1155/2020/8304903


## VMDecomp 1.0.0
