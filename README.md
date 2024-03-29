# Fused inverse-normal (FIN) method for meta analysis
Fused inverse-normal method for meta-analysis of RNA-seq data.

This repository contains an implementation of fused inverse-normal (FIN) method for meta-analysis of multiple but related RNA sequencing (RNA-seq) studies. The current code demonstrates the method for combination of p-values from three different real glioblastoma RNA-seq studies and can be adapted to any number of RNA-seq gene expression studies.

Although designed for meta-analysis in the context of RNA-seq data, the method can be applied in other contexts when the variables of interest assume two distinct states across multiple studies. 

## Published article and citation
Prasad, B., Li, X. Fused inverse-normal method for integrated differential expression analysis of RNA-seq data. BMC Bioinformatics 23, 320 (2022). https://doi.org/10.1186/s12859-022-04859-9

## Prerequisites:
* Results from per-study differential analysis for RNA-seq studies. Popular methods such as DESeq and edgeR can be used for this step. These results should at least contain the raw p-values and log_2(FC) (logFC) from the individual differential analaysis.  
* Assessment of the underlying assumption that p-values for all genes obtained from per-study differential analysis are uniformly distributed under the null hypothesis. Usually this assumption is not satisfied in case of RNA-seq data but filtering of weakly expressed genes using the method described in Rau et al. (2014)[1] or using the counts per million criteria described in Raithel et al. (2016)[2] circumvents this difficulty to a significant extent.
* R version 3.6.0 or above.

## References
1. Rau A, Marot G, Jaffrézic F. Differential meta-analysis of RNA-seq data from multiple studies. BMC bioinformatics. 2014 Dec 1;15(1):91. [DOI](https://doi.org/10.1186/1471-2105-15-91)
2. Raithel S, Johnson L, Galliart M, Brown S, Shelton J, Herndon N, Bello NM. Inferential considerations for low-count RNA-seq transcripts: a case study on the dominant prairie grass Andropogon gerardii. BMC genomics. 2016 Dec 1;17(1):140. [DOI](https://doi.org/10.1186/s12864-016-2442-7)

