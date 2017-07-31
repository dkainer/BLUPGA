#' Standard GBLUP
#'
#' This function runs the BLUP|GA method where SNPs in candidate genes are weighted in the GRM prior to GBLUP
#' @param valset vector of indices Defaults to NULL.
#' @param G G-matrix constructed using all available SNPs and all samples Defaults to NULL.
#' @param S Weighted G-matrix constructed using only selected SNPs and all samples Defaults to NULL.
#' @param phenodata data frame with 2 columns: ID and y Defaults to NULL.
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
    T <- wt*S + (1-wt)*G
    colnames(T) <- rownames(G)
    rownames(T) <- rownames(G)
    GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=T, PEV=F)
    cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
    row     <- data.frame(W=wt, COR=cg)
    out     <- rbind(out, row)
  }
  return(out)
}


#' Standard GBLUP
#'
#' This function runs the BLUP|GA method where SNPs in candidate genes are weighted in the GRM prior to GBLUP
#' @param valset vector of indices Defaults to NULL.
#' @param G G-matrix constructed using all available SNPs and samples Defaults to NULL.
#' @param phenodata data frame with 2 columns: ID and y Defaults to NULL.
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
#' @param valset vector of indices Defaults to NULL.
#' @param G G-matrix constructed using all available SNPs and samples Defaults to NULL.
#' @param phenodata data frame with 2 columns: ID and y Defaults to NULL.
#' @param GRMs List of G-matrices, each constructed from just the SNPs in one candidate gene/region Defaults to NULL.
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

blupga_para_CAND <- function(G=NULL, phenodata=NULL, valset=NULL, GRMs=NULL, verbose=T)
{
  if (verbose==T) cat("CAND\n")
  out <- data.frame(W=numeric(), COR=numeric())

  pheno.blup            <- phenodata
  pheno.blup$y[valset]  <- NA

  cl<-makeCluster(7, type="SOCK", outfile="")
  registerDoSNOW(cl)

  out <- foreach(wt = seq(0.00, 1, 0.10), .packages=c("rrBLUP"),
                        .verbose=TRUE, .combine=rbind, .inorder = FALSE) %dopar%
  {
    wt_per_gene <- wt/length(GRMs$grm)
    tmp         <- lapply(GRMs$grm, function(x){ wt_per_gene*x })
    Tmat           <- Reduce("+", tmp) + (1-wt)*G
    colnames(Tmat) <- phenodata$ID
    rownames(Tmat) <- phenodata$ID
    GBLUP   <- rrBLUP::kin.blup(data=pheno.blup, geno="ID", pheno="y", GAUSS=F, K=Tmat, PEV=F)
    cg      <- cor(phenodata$y[valset], GBLUP$pred[valset])
    row     <- data.frame(W=wt, COR=cg)
    row
  }

  stopCluster(cl)
  closeAllConnections()
  return(out)
}

#' BLUP|GA using top SNPS according to estimated SNP effect
#'
#' This function runs the BLUP|GA method where SNPs with the greatest squared effect size are weighted in the GRM prior to GBLUP
#' @param valset vector of indices Defaults to NULL.
#' @param G G-matrix constructed using all available SNPs and samples Defaults to NULL.
#' @param phenodata data frame with 2 columns: ID and y Defaults to NULL.
#' @param perc proportion of SNPs to be weighted (between 0 and 1), where 0.05 means the top 5% of SNPs will be weighted. Defaults to NULL.
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' blupga_EFF()
blupga_EFF <- function(G, phenodata, valset=NULL, genomat, bsq, perc, flank=TRUE, verbose=T)
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
  S   <- cgrm(top.genomat, W)

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
#' @param weights vector of weights for each of the SNPs in the geno matrix
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' make_weighted_GRM()
make_weighted_GRM <- function(genomat, weights)
{
  G1 <- cgrm(genomat, weights) #uses the fast, low mem cpgen implementation of Van Raden (2008).
  colnames(G1) <- rownames(genomat)
  rownames(G1) <- rownames(genomat)
  return(G1)
}

#' Construct the G matrix from genotypes
#'
#' This function ....
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' make_GRM()
make_GRM <- function(genomat)
{
  G1 <- cgrm.A(genomat) #uses the fast, low mem cpgen implementation of Van Raden (2008).
  colnames(G1) <- rownames(genomat)
  rownames(G1) <- rownames(genomat)
  return(G1)
}

#' run a GWAS using a training set of samples to estimate marker effects
#'
#' This function runs the BLUP|GA method where SNPs with the greatest squared effect size are weighted in the GRM prior to GBLUP
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @param phenodata data frame with 2 columns: ID and y Defaults to NULL.
#' @param valset vector of indices Defaults to NULL.
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
