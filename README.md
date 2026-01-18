## 1.  TASTE
Increasing evidence suggests that related cancers share alterations of common regulatory programs. Trans-associations of cancer risk variants mediated via molecular phenotypes, such as gene expression and protein levels, can help uncover these downstream mechanisms. In this paper we introduce TASTE (Trans Association using Shared factorization and TEsting), a summary statistic-based framework to identify protein sets that are trans- regulated by genetic variants associated with sets of biologically related cancers. TASTE consists of three steps: (1) TASTE-D, a low-rank matrix factorization to estimate shared and group-specific trans-association patterns across cancers; (2) TASTE-S, a sparse singular value decomposition to identify proteins driving shared effects; (3) TASTE-T, a competitive testing strategy for evaluating significance of trans-associations captured by the identified protein-set.  **TASTE** is a R package for extracting protein sets that are trans-regulated by genetic variants associated with sets of biologically related cancers.

## 2.  Installation
To install the package use the following code.
```r
install.packages("devtools")
devtools::install_github("diptavo /TASTE")

#Load package
library(TASTE)
```
## Data for Illustration:
**TASTE** comes with an illustrative data set that can be loaded directly after installing and loading the package.

```r
dat <- readRDS(system.file("extdata", "data_illustration.rds", package = "TASTE"))
```

A detailed description of the method, input format, and example analysis is provided in the package tutorial: **[TASTE Tutorial](TASTE.md)**
