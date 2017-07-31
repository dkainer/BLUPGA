###############################################

#' Construct a weighted G matrix from genotypes
#'
#' This function ...
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @param weights vector of weights for each of the SNPs in the \code{genomat}. Use a vector of 1s to avoid weighting.
#' @return a weighted Genomic Relationship Matrix with cols and rows named according to the rownames in \code{genomat}. Uses the 'cpgen' package implementation of Van Radens GRM method
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' make_weighted_GRM()
make_weighted_GRM <- function(genomat, weights)
{
  G1 <- cpgen::cgrm(genomat, weights) #uses the fast, low mem cpgen implementation of Van Raden (2008).
  colnames(G1) <- rownames(genomat)
  rownames(G1) <- rownames(genomat)
  return(G1)
}

#' Construct the G matrix from genotypes
#'
#' This function ....
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @return additive Genomic Relationship Matrix with cols and rows named according to the rownames in \code{genomat}. Uses the 'cpgen' package implementation of Van Radens GRM method
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' make_GRM()
make_GRM <- function(genomat)
{
  G1 <- cpgen::cgrm.A(genomat) #uses the fast, low mem cpgen implementation of Van Raden (2008).
  colnames(G1) <- rownames(genomat)
  rownames(G1) <- rownames(genomat)
  return(G1)
}

#' estimate the additive effect of each SNP upon a trait (i.e. quick GWAS)
#'
#' estimate the additive effect of each SNP upon a trait (i.e. quick GWAS). SNPs are provided in a genotype matrix.
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices Defaults to NULL.
#' @return a vector of squared SNP effects, one per SNP in \code{genomat}, estimated with the EMMAX GWAS method from the 'cpgen' package.
#' @export
#' @examples
#' est_SNPeffects()
est_SNPeffects <- function(phenodata, genomat, valset, fixmat=NULL)
{
  cat("estimating marker effects with EMMAX\n")

  pheno.blup            <- phenodata
  pheno.blup$y[valset]  <- NA
  selected    <- phenodata[valset,]$ID
  trainset    <- which(!phenodata$ID %in% selected )
  eff     <- (cpgen::cGWAS.emmax(phenodata$y[trainset], genomat[trainset,], X=fixmat[trainset,])$beta)^2
  return(eff)
}
