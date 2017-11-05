###############################################

#' Construct a weighted G matrix from genotypes
#'
#' This function ...
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @param weights vector of weights for each of the SNPs in the \code{genomat}. Use a vector of 1s to avoid weighting.
#' @return a weighted Genomic Relationship Matrix with cols and rows named according to the rownames in \code{genomat}. Uses the \pkg{cpgen} package implementation of Van Radens GRM method
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' # get an example genotype matrix
#' data(M)
#' # generate weights for the SNPs
#' wt <- runif(ncol(M), min = 1, max = 10)
#'
#' make_weighted_GRM(M, wt)
make_weighted_GRM <- function(genomat, weights)
{
  G <- cpgen::cgrm(genomat, weights) #uses the fast, low mem cpgen implementation of Van Raden (2008).
  colnames(G) <- rownames(genomat)
  rownames(G) <- rownames(genomat)
  return(G)
}

#' Construct the G matrix from genotypes
#'
#' This function ....
#' @param genomat matrix of genotypes in -1,0,1 format (i.e. 0 is heterozygous, 1 and -1 are opposing homozygous states)
#' @return additive Genomic Relationship Matrix with cols and rows named according to the rownames in \code{genomat}.Uses the \pkg{cpgen} package implementation of Van Radens GRM method
#' @keywords GBLUP,BLUP|GA,SNP selection
#' @export
#' @examples
#' # get an example genotype matrix
#' data(M)
#' make_GRM(M)
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
#' @return a vector of squared SNP effects, one per SNP in \code{genomat}, estimated with the EMMAX GWAS method from the \pkg{cpgen} package.
#' @export
#' @examples
#' # get an example genotype matrix
#' data(M)
#' # get an example phenotype data frame
#' data(pheno)
#' # choose a validation set of 20 random individuals
#' val <- sample(1:nrow(pheno), 20)
#' est_SNPeffects(pheno, M, val)
est_SNPeffects <- function(phenodata, genomat, valset, fixmat=NULL, method="emmax", R2=0.5)
{
  cat("estimating marker effects using: ", method,"\n")

  selected    <- phenodata[valset,]$ID
  trainset    <- which(!phenodata$ID %in% selected )

  if(method=="BayesA")
  {
    set_num_threads(1)
    z_list <- list(genomat[trainset,])
    par_random <- list(list(method="BayesA",scale=var(pheno$y[trainset])/2, df=4, name="add"))
    fit <- cpgen::clmm(y = phenodata$y[trainset], X=fixmat[trainset,], Z = z_list, par_random=par_random, niter=2000, burnin=500)
    eff <- fit$add$posterior$estimates_mean^2
  }
  else if(method=="BayesC")
  {
    FE.mat <- model.matrix(ID ~ 1 + as.factor(FE), data=phenodata)
    hyp <- VIGoR::hyperpara(Geno=genomat[trainset,], Method = method, Mvar=0.5, f=0.1, Nu=c(4,8,12), Kappa=c(0.05,0.01,0.001))
    eff <- VIGoR::vigor(phenodata$y[trainset], genomat[trainset,], Method="BayesC",
                        Hyperparameters = hyp, Function = "tuning", Covariates = FE.mat[trainset,])$Beta^2
  }
  else if(method=="BayesB")
  {
    FE.mat <- model.matrix(ID ~ 1 + as.factor(FE), data=phenodata)
    #hyp <- VIGoR::hyperpara(Geno=genomat[trainset,], Method = method, Mvar=0.5, f=0.1, Nu=c(4,8,12), Kappa=c(0.05,0.01,0.001))
    hyp <- c(16,0.000533,0.01)
    eff <- VIGoR::vigor(phenodata$y[trainset], genomat[trainset,], Method="BayesB",
                        Hyperparameters = c(5,1,0.01), Function = "fitting", Covariates = FE.mat[trainset,])$Beta^2
  }
  else if(method=="BL")
  {
    FE.mat <- model.matrix(ID ~ 1 + as.factor(FE), data=phenodata)
    eff <- VIGoR::vigor(phenodata$y[trainset], genomat[trainset,], Method=method,
                        Hyperparameters = matrix(c(1,0.1, 1,0.01, 1,0.001), byrow = T, nrow = 2),
                        Function = "tuning", Covariates = FE.mat[trainset,])$Beta^2

  }
  else if(method=="emRR")
  {
    eff <- bWGR::emRR(y = phenodata$y[trainset], genomat[trainset,], R2)$b^2
  }
  else if(method=="emBL2")
  {
    eff <- bWGR::emDE(y = phenodata$y[trainset], genomat[trainset,])$b^2
  }
  else if(method=="Ridge")
  {
    eff <- rrBLUP::mixed.solve(y=phenodata$y[trainset], Z=genomat[trainset,], X=fixmat[trainset,])$u^2
    eff <- as.numeric(eff)
  }
  else if(method=="GWAS")
  {
    eff     <- (cpgen::cGWAS(phenodata$y[trainset], genomat[trainset,], X=fixmat[trainset,])$beta)^2
  }
  else if(method=="LMM")
  {
    phenodata %<>% dplyr::rename(scanID=ID)
    nullmod   <- fitNullMM( phenodata[trainset,], outcome="y", covars=c("Q1","Q2","Q3","Q4","Q5"), covMatList=as.matrix(G1[trainset,trainset]) )
    assoc <- assocTestMM(GenotypeData(whole.genesis), nullmod, impute.geno = F, snp.block.size = 20000, snp.include = pruned.id)
    eff <- -log10(assoc$Wald.pval)
  }
  else  # default is EMMAX
  {
    cat("estimating marker effects using: emmax \n")
    eff     <- (cpgen::cGWAS.emmax(phenodata$y[trainset], genomat[trainset,], X=fixmat[trainset,])$beta)^2
  }
  return(eff)
}

est_SNPpvals <- function(phenodata, genomat, valset, fixmat=NULL, method="emmax")
{
  cat("estimating marker effects using: ", method)

  #pheno.blup            <- phenodata
  #pheno.blup$y[valset]  <- NA
  selected    <- phenodata[valset,]$ID
  trainset    <- which(!phenodata$ID %in% selected )

  # default is EMMAX
  pval     <- -log10(cpgen::cGWAS.emmax(phenodata$y[trainset], genomat[trainset,], X=fixmat[trainset,])$p_value)

  return(pval)
}

simulate <- function(nind, nsnp, prop_qtl, h2, nval, seed=999, effmethod="emmax", model="EFF0.1")
{
  cpgen::rand_data(n=nind, p_marker = nsnp, h2=h2, prop_qtl = prop_qtl, seed = seed)
  pheno <- data.frame(ID=seq(1,nind), y=y)
  val <- sample(1:nrow(pheno), nval)
  G1 <- make_GRM(M)
  colnames(G1) <- pheno$ID
  rownames(G1) <- pheno$ID

  if(effmethod=="BayesA")
    eff <- est_SNPeffects(pheno, M, val, method = "BayesA")
  else
    eff <- est_SNPpvals(pheno, M, val, method="emmax")
  plot(eff)
  if(model == "EFF0.1")
    res <- blupga_EFF(G1, pheno, val, M, eff, perc=0.001, flank=TRUE)
  else
    res <- blupga_EFF(G1, pheno, val, M, eff, perc=0.01, flank=TRUE)

  print(res)
}

simulate_mixture <- function(nind, nsnp, prop_qtl, h2, nval, seed=999, effmethod="emmax", model="EFF0.1")
{

  # major effect QTLs
  cpgen::rand_data(n=nind, p_marker = 10, h2=h2, prop_qtl = 1, seed = seed)
  M.big <- M

  # small effect QTLs
  cpgen::rand_data(n=nind, p_marker = 100, h2=h2, prop_qtl = 1, seed = seed)
  M.small <- M

  # tiny effect QTLs
  cpgen::rand_data(n=nind, p_marker = 500, h2=h2, prop_qtl = 1, seed = seed)
  M.tiny <- M

}
