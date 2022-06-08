FROM rocker/rstudio:devel

LABEL maintainer='Lampros Mouselimis'

RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update && \
 apt-get install -y zlib1g-dev pandoc pandoc-citeproc make libfftw3-dev libicu-dev libpng-dev && \
 apt-get install -y sudo && \
 apt-get install -y libarmadillo-dev && \
 apt-get install -y libblas-dev && \
 apt-get install -y liblapack-dev && \
 apt-get install -y libarpack++2-dev && \
 apt-get install -y gfortran  && \
 apt-get install -y libxml2-dev && \
 apt-get install -y libssh2-1-dev && \
 apt-get install -y zlib1g-dev && \
 R -e "install.packages('devtools', dependencies = TRUE, repos = 'https://cloud.r-project.org/')" && \
 R -e "install.packages(c( 'Rcpp', 'RcppArmadillo', 'OpenImageR', 'R.matlab', 'glue', 'testthat', 'rmarkdown', 'knitr', 'remotes' ), repos =  'https://cloud.r-project.org/' )" && \
 R -e "remotes::install_github('mlampros/VMDecomp', upgrade = 'never', dependencies = FALSE, repos = 'https://cloud.r-project.org/')" && \
 apt-get autoremove -y && \
 apt-get clean

ENV USER rstudio


