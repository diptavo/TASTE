## 1. Overview:

Increasing evidence suggests that related cancers share alterations of
common regulatory programs. Trans-associations of cancer risk variants
mediated via molecular phenotypes, such as gene expression and protein
levels, can help uncover these downstream mechanisms. In this paper we
introduce TASTE (Trans Association using Shared factorization and
TEsting), a summary statistic-based framework to identify protein sets
that are trans- regulated by genetic variants associated with sets of
biologically related cancers. TASTE consists of three steps: (1)
**TASTE-D**, a low-rank matrix factorization to estimate shared and
group-specific trans-association patterns across cancers; (2)
**TASTE-S**, a sparse singular value decomposition to identify proteins
driving shared effects; (3) **TASTE-T**, a competitive testing strategy
for evaluating significance of trans-associations captured by the
identified protein-set. This vignettes will demonstrate how to use
**TASTE** package in R to extract the set of proteins that are
associated with sets of biologically related cancers.  


## 2. Input Data Format:

``` r
harmonize_matrices <- function(matrices) {
  # Find the common column names
  com_cols <- Reduce(intersect, lapply(matrices, colnames))
  # If no common columns are found
  if (length(com_cols) == 0) {
    message("No common columns found.")
    return(NULL)
  }
   # Subset matrices to include only the common columns
  um <- lapply(matrices, function(mat) {
    mat[, com_cols, drop = FALSE]
  })
  updated_matrices=lapply(um,na.omit)  
  return(updated_matrices)
}
```

The inputs are list of matrices corresponding to a group of cancer(
example- genitourinary cancers). Each matrices of the list corresponds
to a cancer type ((renal cell carcinoma \[RCC\] and bladder cancer
\[BLCA\]). Each column of the matrix should represent Protein ID and
each row should represent summary statistics (Z-values) from standard
trans-pQTL analysis. Different matrices may contain overlapping but not
identical protein panels. The function **harmonize_matrices()** alings
all matrices with the common set of Proteins. If no such common proteins
are founds, the function stops with the message **No common columns
found**. In addition, this function also removes the rows which have
missing entries.

## 3. Low Rank Decomposition for the Estimation of the Joint Structure:

``` r
TASTE.D=function(matrices,mi,status)
{
  library(r.jive)
  rr=jive(matrices,maxiter=mi,showProgress=status)
  jr=rr$rankJ
  jm=do.call(rbind, rr$joint)#rbind the list of matrices.
  result=list(jr,jm)
  return(result)
}
```

This function is used for estimating the joint structure via a low rank
decomposition using JIVE in r.jive package. It takes input as
(harmonized) list of matrices and parameters of JIVE (**maxiter**: The
maximum number of iterations for each instance of the JIVE algorithm
,**showProgress** :  
A boolean indicating whether or not to give output showing the progress
of the algorithm.) and returns the **rank of the joint effect matrices**
and the **joint matrix** as a list.

## 4. Protein Selection:

``` r
TASTE.S=function(pro_name,jefs,lower,upper,gamma)
{
  library(PMA)
  R=jefs[[2]] #rbind the list of matrices.
  
  pro.list <- vector("list", length = jefs[[1]])
  n1=pro_name
  c1=lower
  c2=upper
  for(i in 1:(jefs[[1]]))
  {
    pp=length(n1)
    l=c2+1
    z=1
    while (l<c1 || l>c2 && z< pp^0.25) {
      tt1=SPC(R,sumabsv=z, K=1)
      l=length(tt1$v[tt1$v!=0])
      z=z+gamma
    }
    ind=which(tt1$v!=0) # indices of protein previously selected
    pro.list[[i]]=n1[ind]
    
    x=(tt1$d)*(tt1$u)%*%t((tt1$v))
    R.star=R-x
    R=R.star[,-ind]
    n1=n1[-ind]
  }
  return(pro.list)
}
```

This function takes into input **pro_name**: common protein ID, **jefs**
: list that contains rank of the joint effect matrices and joint matrix
as its first and second elemnt, **lower**: lower limit of the $l_0$ norm
of the vector v a sparse loading vector, corresponding to columns of
joint matrix. , **upper**: upper limit of the $l_0$ norm of the vector v
a sparse loading vector, corresponding to columns of joint matrix. ,
**gamma**: step size. It extracts the proteins that drives the joint
structure by performing a variant of sparse principle component analysis
via PMA package.

## 5. Full Pipeline Function:

``` r
TASTE=function(matrices,mi=20,status="False",lower=30,upper=100,gamma=0.1)
{
  # Check if the input is a list of matrices
  if (!all(sapply(matrices, is.matrix))) {
    stop("All inputs must be matrices.")
  }
  
  hm=harmonize_matrices(matrices)
  
  if (is.null(hm)) {
    stop("No common columns found.")
  }
  #pro_name=colnames(hm[[1]])
  
  jx=TASTE.D(hm,mi,status)
  
  protein_list=TASTE.S(colnames(hm[[1]]),jx,lower,upper,gamma)
  return(protein_list)
}
```

The complete workflow is implemented in the function TASTE().

## 6. Data Illustration:

``` r
library(TASTE)

dat <- readRDS(system.file("extdata", "data_illustration.rds", package = "TASTE"))

res <- TASTE(dat)
```

    ## 12
    ## 1234
    ## 12345
    ## 12345
    ## 1234
    ## 1234
    ## 1234
    ## 12345
    ## 123456
    ## 123456
    ## 123456
    ## 12345
    ## 123456
    ## 123456
    ## 123456
    ## 123456
    ## 123456
    ## 1234567
    ## 1234567
    ## 1234567
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920
    ## 1234567891011121314151617181920

``` r
res
```

    ## [[1]]
    ##  [1] "BTN2A1"  "BTN3A2"  "CAPS"    "CCL19"   "CD22"    "CD72"    "CD74"   
    ##  [8] "CD8A"    "CPVL"    "CXCL1"   "CXCL3"   "FOLR2"   "GZMA"    "GZMH"   
    ## [15] "IGF1R"   "IL10"    "IL12A"   "IL12B"   "IL5RA"   "IL7R"    "KIR2DL2"
    ## [22] "KIR2DL3" "KLRB1"   "LAMB1"   "LILRB1"  "LILRB4"  "LRPAP1"  "RNASET2"
    ## [29] "SFTPD"   "SGSH"    "SIGLEC9"
    ## 
    ## [[2]]
    ##  [1] "AGR2"     "APOA4"    "CA4"      "CCL28"    "CD160"    "CD300E"  
    ##  [7] "CFP"      "CGREF1"   "CLEC1A"   "CLEC4C"   "CRNN"     "ENTPD6"  
    ## [13] "FAP"      "FCRL6"    "FOLR3"    "IGF2R"    "IL15RA"   "LILRB2"  
    ## [19] "LTBR"     "MANSC4"   "MMP3"     "OPTC"     "PCDH1"    "PDGFRA"  
    ## [25] "RNASE6"   "TF"       "TIMD4"    "TNFRSF1B" "TSHB"     "VCAM1"   
    ## [31] "WASL"     "XCL1"
