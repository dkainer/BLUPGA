devtools::use_package("cpgen")
devtools::use_package("rrBLUP")

#' Basic BLUP|GA
#'
#' This function runs the BLUP|GA method where certain SNPs are weighted in the a special GRM (S) prior to GBLUP.
#' @param G GRM constructed using all available SNPs and all samples Defaults to NULL. Use \code{\link{make_GRM}} to get this.
#' @param S Weighted G-matrix constructed using only selected SNPs and all samples Defaults to NULL. Use \code{\link{make_weighted_GRM}} to get this.
#' @param phenodata data frame with 2 columns. One col must be named 'ID' and contain sample IDs. Another column must be named 'y' and contain the phenotypes. Defaults to NULL.
#' @param valset vector of indices that defines which rows in \code{phenodata} will be set to NA and used for cross-validation Defaults to NULL.
#' @return A data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}.
#'
#'  Since BLUP|GA is run for each value of omega (W) from 0.0 to 1.0 in increments of 0.10, each row of the returned data frame
#'   contains the cross-validation correlation at one value of omega (W). This allows the user to find the value of W
#'   at which the predictive ability (COR) is maximised.
#'
#'   Columns:
#' \describe{
#'   \item{W}{omega weighting for selected SNPS in candidate genes (0.0--1.0)}
#'   \item{COR}{cross validation predictive ability (0.0--1.0)}
#'   }
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' # get example genotype matrix and phenotype data
#' data(M)
#' data(pheno)
#'
#' # select some 'special' SNPs from M to be weighted
#' GAsnps <- sample(1:ncol(M), 20)
#'
#' # generate random weights for the 'special' SNPs
#' wt <- runif(length(GAsnps), min = 1, max = 10)
#' # make a weighted GRM for the 'special' SNPs
#' S <- make_weighted_GRM(M[,GAsnps], wt)
#' # make a standard GRM for all SNPs
#' G <- make_GRM(M)
#'
#' # choose a validation set of 20 random individuals
#' val <- sample(1:nrow(pheno), 20)
#' results <- blupga(G, S, pheno, val)
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
#' @return a data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the valset.
#' \describe{
#'   \item{W}{omega weighting for selected SNPS in candidate genes (0.0--1.0)}
#'   \item{COR}{cross validation predictive ability (0.0--1.0)}
#'   }
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' # get example genotype matrix and phenotype data
#' data(M)
#' data(pheno)
#'
#' G <- make_GRM(M)
#'
#' # choose a validation set of 20 random individuals
#' val <- sample(1:nrow(pheno), 20)
#' results <- blupga_GBLUP(G, pheno, val)
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
#' @return A data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}.
#'
#'  Since BLUP|GA is run for each value of omega (W) from 0.0 to 1.0 in increments of 0.10, each row of the returned data frame
#'   contains the cross-validation correlation at one value of omega (W). This allows the user to find the value of W
#'   at which the predictive ability (COR) is maximised.
#' \describe{
#'   \item{W}{omega weighting for selected SNPS in candidate genes (0.0--1.0)}
#'   \item{COR}{cross validation predictive ability (0.0--1.0)}
#'   }
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
#' @return A data frame containing the correlation between the predicted phenotype and the true phenotype of the individuals in the \code{valset}.
#'
#'  Since BLUP|GA is run for each value of omega (W) from 0.0 to 1.0 in increments of 0.10, each row of the returned data frame
#'   contains the cross-validation correlation at one value of omega (W). This allows the user to find the value of W
#'   at which the predictive ability (COR) is maximised.
#' \describe{
#'   \item{W}{omega weighting for selected SNPS in candidate genes (0.0--1.0)}
#'   \item{COR}{cross validation predictive ability (0.0--1.0)}
#'   }
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' # get example genotype matrix and phenotype data
#' data(M)
#' data(pheno)
#'
#'
#' # choose a validation set of 20 random individuals
#' val <- sample(1:nrow(pheno), 20)
#'
#' # Run a GWAS to estimate squared SNP effects for all SNPs.
#' # By calling est_SNPeffects() this process will be done using the individuals NOT in your 'val' set.
#' # Of course you can use your own method to generate the full vector of SNP effects using a training set of individuals.
#' bsq <- est_SNPeffects(pheno, M, val)
#'
#' # make a standard GRM for all SNPs
#' G <- make_GRM(M)
#'
#' # run BLUPGA where the top 1 percent of SNPs (according to effect size) are weighted
#' results <- blupga_EFF(G, pheno, val, M, bsq, 0.01, flank=FALSE)
#'
#' # run BLUPGA where the top 0.1 percent of SNPs (according to effect size) are weighted and
#' # SNPs immediately upstream and downstream of those top SNPs are also weighted.
#' results <- blupga_EFF(G, pheno, val, M, bsq, 0.001, flank=TRUE)
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




.onAttach <- function(libname, pkgname) {
  packageStartupMessage("BLUP|GA genomic Prediction package is now loaded")
}
