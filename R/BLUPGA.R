#' Standard GBLUP
#'
#' This function runs the BLUP|GA method where certain SNPs are weighted in the a special GRM (S) prior to GBLUP.
#' @param G GRM constructed using all available SNPs and all samples Defaults to NULL. Use \code{make_GRM()} to get this.
#' @param S Weighted G-matrix constructed using only selected SNPs and all samples Defaults to NULL. Use \code{make_weighted_GRM()} to get this.
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices that defines which rows in \code{phenodata} will be set to NA and used for cross-validation Defaults to NULL.
#' @return a data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}. Each row contains that correlation at increasing values of the weighting value W.
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' blupga()
blupga <- function(G, S, phenodata, valset, verbose=T)
{
  if (verbose==T) cat("POS\n")
  pheno.blup              <- phenodata
  pheno.blup$y[ valset ]  <- NA

  out <- data.frame(W=numeric(), COR=numeric())

  for(wt in seq(0.00, 1, 0.10))
  {
    Tmat <- wt*S + (1-wt)*G
    colnames(Tmat) <- rownames(G)
    rownames(Tmat) <- rownames(G)
    GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=Tmat, PEV=F)
    cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
    row     <- data.frame(W=wt, COR=cg)
    out     <- rbind(out, row)
  }
  return(out)
}


#' Standard GBLUP
#'
#' This function runs straightforward GBLUP using the rrBLUP package. This is the equivalent of BLUP|GA where W=0 (i.e. no SNPs are weighted)
#' @param G GRM constructed using all available SNPs and all samples Defaults to NULL. Use \code{make_GRM()} to get this.
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices that defines which rows in \code{phenodata} will be set to NA and used for cross-validation Defaults to NULL.
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' blupga_GBLUP()
blupga_GBLUP <- function(G, phenodata, valset, verbose=T)
{
  if (verbose==T) cat("GBLUP\n")
  wt                  <- 0
  pheno.blup          <- phenodata
  pheno.blup$y[ valset ] <- NA

  GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=G, PEV=F)
  cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
  out     <- data.frame(W=wt, COR=cg)
  return(out)
}


#' BLUP|GA using SNPs in Candidate Genes
#'
#' This function runs the BLUP|GA method where SNPs in candidate genes are weighted in the GRM prior to GBLUP
#' @param G G-matrix constructed using all available SNPs and samples Defaults to NULL.
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices that defines which rows in \code{phenodata} will be set to NA and used for cross-validation Defaults to NULL.
#' @param GRMs List of G-matrices, each constructed from just the SNPs in one candidate gene/region. Defaults to NULL.
#' @return a data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}. Each row contains that correlation at increasing values of the weighting value W.
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' blupga_CAND()
blupga_CAND <- function(G=NULL, phenodata=NULL, valset=NULL, GRMs=NULL, verbose=T)
{
  if (verbose==T) cat("CAND\n")
  out <- data.frame(W=numeric(), COR=numeric())

  pheno.blup            <- phenodata
  pheno.blup$y[valset]  <- NA

  for(wt in seq(0.00, 1, 0.10))
  {
    wt_per_gene <- wt/length(GRMs$grm)
    tmp         <- lapply(GRMs$grm, function(x){ wt_per_gene*x })
    Tmat           <- Reduce("+", tmp) + (1-wt)*G
    colnames(Tmat) <- phenodata$ID
    rownames(Tmat) <- phenodata$ID
    GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=Tmat, PEV=F)
    cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
    row     <- data.frame(W=wt, COR=cg)
    out     <- rbind(out, row)
  }
  return(out)
}


# blupga_para_CAND <- function(G=NULL, phenodata=NULL, valset=NULL, GRMs=NULL, verbose=T)
# {
#   if (verbose==T) cat("CAND\n")
#   out <- data.frame(W=numeric(), COR=numeric())
#
#   pheno.blup            <- phenodata
#   pheno.blup$y[valset]  <- NA
#
#   cl<-makeCluster(7, type="SOCK", outfile="")
#   registerDoSNOW(cl)
#
#   out <- foreach(wt = seq(0.00, 1, 0.10), .packages=c("rrBLUP"),
#                         .verbose=TRUE, .combine=rbind, .inorder = FALSE) %dopar%
#   {
#     wt_per_gene <- wt/length(GRMs$grm)
#     tmp         <- lapply(GRMs$grm, function(x){ wt_per_gene*x })
#     Tmat           <- Reduce("+", tmp) + (1-wt)*G
#     colnames(Tmat) <- phenodata$ID
#     rownames(Tmat) <- phenodata$ID
#     GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=Tmat, PEV=F)
#     cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
#     row     <- data.frame(W=wt, COR=cg)
#     row
#   }
#
#   stopCluster(cl)
#   closeAllConnections()
#   return(out)
# }

#' BLUP|GA using top SNPS according to estimated SNP effect
#'
#' This function runs the BLUP|GA method where SNPs with the greatest squared effect size are weighted in the GRM prior to GBLUP.
#' @param G G-matrix constructed using all available SNPs and samples Defaults to NULL.
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices that defines which rows in \code{phenodata} will be set to NA and used for cross-validation Defaults to NULL.
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states).
#' @param bsq vector of squared marker effects the same length as the number of SNPs in genomat. Can be obtained with \code{est_SNPeffects()}.
#' @param perc proportion of SNPs to be weighted (between 0 and 1), where 0.05 means the top 5 percent of SNPs will be weighted. Defaults to NULL.
#' @param flank choose to include the immediate SNPs to the left and right of a top 'perc' SNP. Defaults to TRUE.
#' @param verbose just leave this for now!
#' @return a data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}. Each row contains that correlation at increasing values of the weighting value W.
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' blupga_EFF()
blupga_EFF <- function(G, phenodata, valset=NULL, genomat, bsq, perc, flank=TRUE, verbose=TRUE)
{
  if (verbose==T) cat("EFF\n")
  out <- data.frame(W=numeric(), COR=numeric())

  pheno.blup            <- phenodata
  pheno.blup$y[valset]  <- NA

  top     <- which( bsq > quantile(bsq, 1-perc) )

  if (flank == TRUE)
    top   <- sort(unique(c(top, top-1, top+1)))   # add flanking SNPs

  cat("Using ", length(top),"weighted SNPs\n")

  if(max(top) > ncol(genomat)) { top <- top[top<length(top)]  }
  if(min(top) < 0) { top <- top[top>0]  }

  top.bsq     <- bsq[top]
  top.genomat <- genomat[,top]

  # make the S GRM and remaining unweighted GRM for weighting scheme 1
  W   <- top.bsq * (length(top.bsq)) / sum(top.bsq)
  S   <- cpgen::cgrm(top.genomat, W)

  for(wt in seq(0.00, 1, 0.10))
  {
    Tmat <- wt*S + (1-wt)*G
    colnames(Tmat) <- rownames(genomat)
    rownames(Tmat) <- rownames(genomat)
    GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=Tmat, PEV=F)
    cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
    row     <- data.frame(W=wt, COR=cg)
    out     <- rbind(out, row)
  }
  return(out)
}


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

#' run a GWAS using a training set of samples to estimate marker effects
#'
#' This function runs the BLUP|GA method where SNPs with the greatest squared effect size are weighted in the GRM prior to GBLUP
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


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("BLUP|GA genomic Prediction package is now loaded")
}
