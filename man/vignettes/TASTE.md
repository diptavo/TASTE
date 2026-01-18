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

'''r
#' function to harmonize the matrices (Step1)
#' @param matrices A list of two or more matrices of z-statistics of corresponding to two or more groups.
harmonize_matrices <- function(matrices) {
  # Find the common column names
  com_cols <- Reduce(intersect, lapply(matrices, colnames))
  # If no common columns are found
  if (length(com_cols) == 0) {
    message("No common columns found.")
    return(NULL)
  }
  # Subset matrices to include only the common columns
  updated_matrices <- lapply(matrices, function(mat) {
    mat[, com_cols, drop = FALSE]
  })
  return(updated_matrices)
}

,,,


