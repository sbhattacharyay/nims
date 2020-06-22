#!/bin/bash
#SBATCH --job-name=R_grid_search
#SBATCH --time=72:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#
#---------------------------------------------------------------------
# SLURM job script to run serial R
#---------------------------------------------------------------------

ml stack
ml r/3.6.1
ml r-assertthat
ml r-callr
ml r-cli
ml r-rlang
ml r-crayon
ml r-rcpp
ml r-usethis
ml r-fs
ml r-devtools
ml r-stringr
ml r-purrr
ml r-curl
ml r-magrittr
ml r-glue
ml r-desc
ml r-digest
ml r-git2r
ml r-remotes
ml r-httr
ml r-jsonlite
ml r-memoise
ml r-pkgbuild
ml r-pkgload
ml r-rcmdcheck
ml r-roxygen2
ml r-rstudioapi
ml r-sessioninfo
ml r-testthat
ml r-withr
ml r-yaml
ml r-sys
ml r-r6
ml r-askpass
ml r-backports
ml r-brew
ml r-clipr
ml r-clisymbols
ml r-commonmark
ml r-evaluate
ml r-gh
ml r-ini
ml r-mime
ml r-openssl
ml r-praise
ml r-prettyunits
ml r-processx
ml r-ps
ml r-rprojroot
ml r-stringi 
ml r-whisker
ml r-xml2
ml r-xopen
ml r-rgdal
ml r-class
ml r-classint
ml r-colorspace
ml r-dbi
ml r-e1071
ml r-ellipsis
ml r-fansi
ml r-ggplot2
ml r-gtable
ml r-igraph
ml r-kernsmooth
ml r-labeling
ml r-lattice
ml r-lazyeval
ml r-labeling
ml r-mass
ml r-matrix
ml r-mgcv
ml r-munsell
ml r-nlme
ml r-pillar
ml r-pkgconfig
ml r-rcolorbrewer
ml r-reshape2
ml r-scales
ml r-sf
ml r-sp
ml r-tibble
ml r-units
ml r-utf8
ml r-vctrs
ml r-viridislite
ml r-zeallot
ml # confirm modules used
Rscript ~/work/nims/scripts/master10_3_ML_grid_search.R