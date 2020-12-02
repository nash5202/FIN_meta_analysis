# Mixture inverse-normal (mixIN) method for meta analysis
Mixture inverse-normal method for meta-analysis of RNA-seq data.

This repository contains an implementation of mixture inverse-normal (mixIN) method for meta-analysis of multiple but related RNA sequencing (RNA-seq) studies. The current code demonstrates the method for combination of p-values from two different real glioblastoma RNA-seq studies.

## Published article
xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## Citation
xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## Prerequisites:
* Results from per-study differential analysis for RNA-seq studies. Popular methods such as DESeq and edgeR can be used for this step. These results should at least contain the raw p-values and logFC from the individual differential analaysis.  
* Assessment of the underlying assumption that p-values for all genes obtained from per-study differential analysis are uniformly distributed under the null hypothesis. Usually this assumption is not satisfied in case of RNA-seq data but filtering of weakly expressed genes using the method described in Rau et al. (2014)[1] or using the counts per million criteria described in Raithel et al. (2016)[2] circumvents this difficulty to a significant extent.
* R version 3.6.0 or above.

## References
1. Rau A, Marot G, Jaffrézic F. Differential meta-analysis of RNA-seq data from multiple studies. BMC bioinformatics. 2014 Dec 1;15(1):91.
2. Raithel S, Johnson L, Galliart M, Brown S, Shelton J, Herndon N, Bello NM. Inferential considerations for low-count RNA-seq transcripts: a case study on the dominant prairie grass Andropogon gerardii. BMC genomics. 2016 Dec 1;17(1):140.

