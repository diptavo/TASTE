## 1. Overview
Increasing evidence suggests that related cancers share alterations of
common regulatory programs.  Trans-associations of cancer risk variants mediated via
molecular phenotypes, such as gene expression and protein levels, can help uncover
these downstream mechanisms.  In this paper we introduce TASTE (Trans Association using Shared factorization and
TEsting), a summary statistic-based framework to identify protein sets that are trans-
regulated by genetic variants associated with sets of biologically related cancers.
TASTE consists of three steps: (1) **TASTE-D**, a low-rank matrix factorization to estimate
shared and group-specific trans-association patterns across cancers; (2) **TASTE-S**, a
sparse singular value decomposition to identify proteins driving shared effects; (3)
**TASTE-T**, a competitive testing strategy for evaluating significance of trans-associations
captured by the identified protein-set.  This vignettes will demonstrate how to use **TASTE** package in R to extract the  set of proteins that are  associated with sets of biologically related cancers.

## 2. Input Data Format

```r

harmonize_matrices <- function(matrices) {
  
  com_cols <- Reduce(intersect, lapply(matrices, colnames))
  
  if (length(com_cols) == 0) {
    message("No common columns found.")
    return(NULL)
  }
 
  updated_matrices <- lapply(matrices, function(mat) {
    mat[, com_cols, drop = FALSE]
  })
  return(updated_matrices)
}

,,,

The inputs are list of matrices corresponding to a group of cancer( example- genitourinary cancers).  Each matrices of the list corresponds to a cancer type ((renal cell carcinoma [RCC] and bladder cancer [BLCA]).  Each column of the matrix should represent Protein ID and each row represents summary statistics (Z-values) from standard trans-pQTL analysis.  Different matrices may contain overlapping but not identical protein panels.  The function  **harmonize_matrices()** alings all matrices with the common set of Proteins.  If no such common proteins are founds, the function stops with the message **No common columns found**.

## 3.


