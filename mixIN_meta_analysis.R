## combination of p-values for meta-analysis using the mixture inverse-normal method.
## @uthor: B Prasad
## created on : March 12, 2021
path <- "Specify the current working directory path"
setwd(path)
data_path <- "/Data" # Specify the path to the individual differential expression analysis results
## required packages
if(!require("pacman")) install.packages("pacman")
pacman::p_load(org.Hs.eg.db, annotate)

## RNA-seq datasets and sample size of each dataset.
datasets <- c("GSE123892", "GSE151352", "TCGA_GBM")          
n_samp <- c(7, 24, 165)
sig_cutoff <- 1e-16       # the lower limit of system cutoff for a number.
# Load individual study p-value data obtained from edgeR. 
# Other methods such as DESeq2 or others can be used for individual 
# differential expression analysis. However the results from individual
# analysis should be put in the input format as shown in the input datasets.
study_1 <- read.csv(file = sprintf("%s/GSE123892_edgeR_results.csv", data_path), header = TRUE, stringsAsFactors = FALSE)
study_2 <- read.csv(file = sprintf("%s/GSE151352_edgeR_results.csv", data_path), header = TRUE, stringsAsFactors = FALSE)
study_3 <- read.csv(file = sprintf("%s/TCGA_GBM_edgeR_results.csv", data_path), header = TRUE, stringsAsFactors = FALSE)

## mix_in function for mixture inverse-normal method
mix_in <- function(study_1, study_2, study_3, n_samp, datasets){
  ## Assess the direction of expression of each gene in each study
  study_1$dir <- sign(study_1$logFC)
  study_2$dir <- sign(study_2$logFC)
  study_3$dir <- sign(study_3$logFC)
  # get all unique genes among the studies  
  unique_genes <- Reduce(union, list(study_1$entrez_id, study_2$entrez_id, study_3$entrez_id))
  # for all the unique genes get the direction status of expression if they are conflicting or not.
  sign <- matrix(data = 0, nrow = length(unique_genes), ncol = length(n_samp)+1)
  row.names(sign) <- unique_genes
  sign[, 1] <- study_1$dir[match(unique_genes, study_1$entrez_id)]
  sign[, 2] <- study_2$dir[match(unique_genes, study_2$entrez_id)]
  sign[, 3] <- study_3$dir[match(unique_genes, study_3$entrez_id)]
  for (l in 1:length(unique_genes))
  {
    if (1 %in% sign[l, c(1:length(n_samp))] & -1 %in% sign[l, c(1:length(n_samp))])
    {
      sign[l, (length(n_samp)+1)] <- 1
    }
  }
  ## computation of N_g
  # 1. Estimation of weights.
  weights <- matrix(0, nrow = length(unique_genes), ncol = length(n_samp))
  weights[which(unique_genes %in% study_1$entrez_id == TRUE), 1] <- n_samp[1]
  weights[which(unique_genes %in% study_2$entrez_id == TRUE), 2] <- n_samp[2]
  weights[which(unique_genes %in% study_3$entrez_id == TRUE), 3] <- n_samp[3]
  denom <- apply(weights, 1, sum)
  weights <- weights/denom
  weights <- sqrt(weights)
  row.names(weights) <- as.character(unique_genes)
  colnames(weights) <- c("study_1", "study_2", "study_3")
  weights <- as.data.frame(weights, stringsAsFactors = FALSE)
  #write.csv(weights, file = sprintf("%s/weights.csv", path))          # to save weights as a csv file
  # 2. calculation of each term of Ng. 
  ng_terms <- matrix(0, nrow = nrow(weights), ncol = ncol(weights))
  for(j in 1:nrow(ng_terms))
  {
    if (sign[j, ncol(sign)] == 0)    # for the non-conflicting genes
    {
      # dataset 1
      if (unique_genes[j] %in% study_1$entrez_id)
      {
        k = which(study_1$entrez_id == unique_genes[j])
        p_val=min(max(study_1$PValue[k],1e-16),1-1e-16)
        ng_terms[j, 1] <- weights$study_1[j] * qnorm((1-p_val), mean = 0, sd = 1)
      }
      # dataset 2
      if (unique_genes[j] %in% study_2$entrez_id)
      {
        k = which(study_2$entrez_id == unique_genes[j])
        p_val=min(max(study_2$PValue[k],sig_cutoff),1-sig_cutoff)
        ng_terms[j, 2] <- weights$study_2[j] * qnorm((1-p_val), mean = 0, sd = 1)
      }
      # dataset 3
      if (unique_genes[j] %in% study_3$entrez_id)
      {
        k = which(study_3$entrez_id == unique_genes[j])
        p_val=min(max(study_3$PValue[k],sig_cutoff),1-sig_cutoff)
        ng_terms[j, 3] <- weights$study_3[j] * qnorm((1-p_val), mean = 0, sd = 1)
      }
    }
    
    if (sign[j, ncol(sign)] == 1) # for the conflicting genes
    {
      # dataset 1
      if (unique_genes[j] %in% study_1$entrez_id)
      {
        k = which(study_1$entrez_id == unique_genes[j])
        p_val=min(max(study_1$PValue[k],sig_cutoff),1-sig_cutoff)
        ind_sign <- sign(study_1$logFC[k])
        ng_terms[j, 1] <- weights$study_1[j]* ind_sign * abs(qnorm((1-p_val), mean = 0, sd = 1))
      }
      #dataset 2
      if (unique_genes[j] %in% study_2$entrez_id)
      {
        k = which(study_2$entrez_id == unique_genes[j])
        p_val=min(max(study_2$PValue[k],sig_cutoff),1-sig_cutoff)
        ind_sign <- sign(study_2$logFC[k])
        ng_terms[j, 2] <- weights$study_2[j]* ind_sign *abs(qnorm((1-p_val), mean = 0, sd = 1))
      } 
      # dataset 3
      if (unique_genes[j] %in% study_3$entrez_id)
      {
        k = which(study_3$entrez_id == unique_genes[j])
        p_val=min(max(study_3$PValue[k],sig_cutoff),1-sig_cutoff)
        ind_sign <- sign(study_3$logFC[k])
        ng_terms[j, 3] <- weights$study_3[j]* ind_sign *abs(qnorm((1-p_val), mean = 0, sd = 1))
      } 
    }
  }
  
  colnames(ng_terms) <- datasets
  ng_terms <- as.data.frame(ng_terms, stringsAsFactors = FALSE)
  row.names(ng_terms) <- row.names(weights)
  # sum all ng_terms row-wise.
  ng <- as.data.frame(rowSums(ng_terms))
  colnames(ng) <- c("ng")
  row.names(ng) <- row.names(ng_terms)
  ## 3. hypothesis testing
  ng$mix_in_p_val <- 1-pnorm(ng$ng)                   # first do one-sided for all and then replace with two sided for conflicting direction genes
  conf_ind <- which(sign[, ncol(sign)] == 1)          # index of conflicting direction genes
  ng$mix_in_p_val[conf_ind] <- 2 * (1-pnorm(abs(ng$ng[conf_ind])))
  ng$BH_p_value <- p.adjust(ng$mix_in_p_val, method = "BH", n = length(ng$mix_in_p_val))
  # 4. Annotations
  ng$entrez_id <- as.numeric(as.character(row.names(ng)))
  # now use the entrez ids to get the symbols.
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  match_SY <- match(row.names(ng), as.character(egSYMBOL$gene_id))
  ng$symbol <- as.character(egSYMBOL$symbol[match_SY])
  return(ng)
}
