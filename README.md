# sarscov2Rutils
Scripts to support preparation and analysis of SARS CoV 2 phylogenetic analyses.   

*Archived version of sarscov2Rutils package as used for the production of skyline plots and importation dates for the paper 'Saudi Arabian SARS-CoV-2 genomes implicate a mutant Nucleoplasmid protein in modulating host interactions and increased viral load in COVID-19 patients'*


[![DOI](https://zenodo.org/badge/335401763.svg)](https://zenodo.org/badge/latestdoi/335401763)

**Version:** 0.1.4  
**Date:** 2020-06-17  
**Code Contributors:** Erik M Volz, Manon Ragonnet, David Jorgensen, Lily Geidelberg, Olivia Boyd, Robert Johnson, Igor Siveroni  
**Maintainer:** Erik M Volz <erik.volz@gmail.com>  
**License:** LGPL

### System requirements
**Required R packages:** ape, treedater(>=0.5.0), lubridate, ggtree, ggplot2, treeio, knitr, coda, phangorn, Hmisc, yaml, glue  
**Suggested R packages:** pika, dplyr, skygrowth, seqinr, treestructure, Biostrings  
**Tested R version:** 3.6.3  
**Tested OS:** Windows 10.0, MacOS 10.14, Linux CentOS 7  
**Tested dependencies:** ape 5.4, treedater 0.5.0, lubridate 1.7.9, ggtree 1.16.6, ggplot2 3.3.2, treeio 1.8.2, knitr 1.29, coda 0.19-3, phangorn 2.5.5, Hmisc 4.4-0, yaml 2.2.1, glue 1.4.1  

Certain functions in this R package also require [IQ-TREE](iqtree.org) (version 1.6.12-Windows tested) and [MAFFT](https://mafft.cbrc.jp/alignment/software/) (version 7.450-win64 tested) and interact with these programs via system calls from R. Adding these programs to PATH is reccomended, or modifications can be made to these functions.

### Installation

Download [tar file](https://github.com/emvolz-phylodynamics/sarscov2Rutils/blob/sarscov2Rutils/sarscov2_0.1.4.tar.gz) and install from package archive file.   
```install.packages(path_to_tar_file, repos = NULL, type="source")```  

Typical install time \< 1 minute (not incluing R installation and package dependencies).

