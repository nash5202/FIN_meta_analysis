## combination of p-values for meta-analysis using the mixture inverse-normal method.
## @uthor: B Prasad
## created on : December 2, 2020
path <- "path specified by the user"
setwd(path)
data_path <- "path to the results from per-study differential analysis"
## packages required.
if(!require("pacman")) install.packages("pacman")
pacman::p_load(org.Hs.eg.db, annotate)
## dataset information
datasets <- c("GSE123892", "GSE151352")
n_samp <- c(7, 24)
## load individual study p-value data obtained from per-study differential analysis.
study_1 <- read.csv(file = sprintf("%s/GSE151352_edgeR_result.csv", data_path), header = TRUE, stringsAsFactors = FALSE)
study_2 <- read.csv(file = sprintf("%s/GSE123892_edgeR_result.csv", data_path), header = TRUE, stringsAsFactors = FALSE)

## function for mixture inverse-normal (mixIN) method for meta-analysis of RNA-seq studies
mix_in <- function(study_1, study_2, n_samp, datasets){
  ## check 1: distribution of raw p-values for a study for uniform distribution
  ## under H_0 assumption.
  h <- hist(study_1$p_value)                     
  plot( h, col= "red", xlim=c(0,1), main = "Study_1", xlab = "p-value")  #histogram
  ### p-value combination using individual studies (mixture inverse_normal approach)
  ## Assess the direction of expression of each gene in each study.
  dir_exp <- function(study){
    for (c in 1: nrow(study))
    {
      if(study$logFC[c] > 0)
      {
        study$dir[c] <- 1
      }
      else
      {
        study$dir[c] <- -1
      }
    }
    return(study)
  }
  study_1 <- dir_exp(study_1)
  study_2 <- dir_exp(study_2)
  # get all unique genes among the studies  
  unique_genes <- Reduce(union, list(study_1$entrez_id, study_2$entrez_id))
  # for all the unique genes get the direction status if they are conflicting or not.
  sign <- matrix(data = 0, nrow = length(unique_genes), ncol = length(n_samp)+2)
  for (l in 1:length(unique_genes))
  {
    sign[l, 1] <- unique_genes[l]
    if(unique_genes[l] %in% study_1$entrez_id)
    {
      ind <- which(study_1$entrez_id == unique_genes[l])
      sign[l, 2] <- study_1$dir[ind]
    }
    if(unique_genes[l] %in% study_2$entrez_id)
    {
      ind <- which(study_2$entrez_id == unique_genes[l])
      sign[l, 3] <- study_2$dir[ind]
    }
    if (1 %in% sign[l, ] & -1 %in% sign[l, ])
    {
      sign[l, 4] <- 1
    }
  }
  ## computation of N_g
  # 1. Estimation of weights.
  weights <- matrix(0, nrow = length(unique_genes), ncol = length(n_samp))
  for (i in 1:nrow(weights))
  {
    if(unique_genes[i] %in% study_1$entrez_id)
    {
      weights[i, 1] <- n_samp[1]
    }
    if(unique_genes[i] %in% study_2$entrez_id)
    {
      weights[i, 2] <- n_samp[2]
    }
  }
  denom <- apply(weights, 1, sum)
  weights <- weights/denom
  weights <- sqrt(weights)
  row.names(weights) <- as.character(unique_genes)
  colnames(weights) <- c("study_1", "study_2")
  weights <- as.data.frame(weights, stringsAsFactors = FALSE)
  #write.csv(weights, file = sprintf("%s/weights.csv", path))          # to save weights as a csv file
  # 3. calculation of each term of Ng. 
  ng_terms <- matrix(0, nrow = nrow(weights), ncol = ncol(weights))
  for(j in 1:nrow(ng_terms))
  {
    if (sign[j, ncol(sign)] == 0)
    {
      if (unique_genes[j] %in% study_1$entrez_id)
      {
        k = which(study_1$entrez_id == unique_genes[j])
        p_val=min(max(study_1$p_value[k],1e-16),1-1e-16)
        ng_terms[j, 1] <- weights$study_1[j] * qnorm((1-p_val), mean = 0, sd = 1)
      }
      if (unique_genes[j] %in% study_2$entrez_id)
      {
        k = which(study_2$entrez_id == unique_genes[j])
        p_val=min(max(study_2$p_value[k],1e-16),1-1e-16)
        ng_terms[j, 2] <- weights$study_2[j] * qnorm((1-p_val), mean = 0, sd = 1)
      }
    }
    if (sign[j, ncol(sign)] == 1)
    {
      if (unique_genes[j] %in% study_1$entrez_id)
      {
        k = which(study_1$entrez_id == unique_genes[j])
        p_val=min(max(study_1$p_value[k],1e-16),1-1e-16)
        ind_sign <- sign(study_1$logFC[k])
        ng_terms[j, 1] <- (weights$study_1[j]* ind_sign * abs(((qnorm((1-p_val), mean = 0, sd = 1)))))
      }
      if (unique_genes[j] %in% study_2$entrez_id)
      {
        k = which(study_2$entrez_id == unique_genes[j])
        p_val=min(max(study_2$p_value[k],1e-16),1-1e-16)
        ind_sign <- sign(study_2$logFC[k])
        ng_terms[j, 2] <- (weights$study_2[j]* ind_sign *abs((qnorm((1-p_val), mean = 0, sd = 1))))
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
  ## hypothesis testing
  ng$mix_in_p_val <- 1-pnorm(ng$ng)                   # first do one-sided for all and then replace with two sided for conflicting direction genes
  conf_ind <- which(sign[, ncol(sign)] == 1)          # index of conflicting direction genes
  ng$mix_in_p_val[conf_ind] <- 2 * (1-pnorm(abs(ng$ng[conf_ind])))
  ng$BH_p_value <- p.adjust(ng$mix_in_p_val, method = "BH", n = length(ng$mix_in_p_val))
  # Annotations
  ng$entrez_id <- as.numeric(as.character(row.names(ng)))
  # now use the entrez ids to get the symbols.
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  match_SY <- match(row.names(ng), as.character(egSYMBOL$gene_id))
  ng$symbol <- as.character(egSYMBOL$symbol[match_SY])
  return(ng)
}











