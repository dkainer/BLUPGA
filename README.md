# BLUP|GA
R implementation of BLUP|GA genomic prediction.

## installation
Can be installed directly from github with:
```R
# install.packages("devtools")
devtools::install_github("dkainer/BLUPGA") 
```

## what is BLUP|GA?
BLUP|GA means *Best Linear Unbiased Prediction given Genetic Architecture*. It is a genomic prediction model that is basically a weighted version of GBLUP. The original BLUP|GA publication can be found here:

Accuracy of Whole-Genome Prediction Using a Genetic Architecture-Enhanced Variance-Covariance Matrix
Zhang et al (2015) [https://doi.org/10.1534/g3.114.016261]

In a nutshell, GBLUP uses a genomic relationship matrix (GRM) to describe the covariance between individuals in the study population. All SNPs used in calculating the GRM are treated the same as each other. With BLUP|GA the calculation of the GRM is done with two separate sub-GRMs:

* **S-matrix:** this is a GRM made with selected SNPs that have an (assumed) additive effect on the trait. These are the *genetic architecture SNPs*. They are given weights according to a metric such as their squared effect size from a GWAS model.
* **G-matrix:** this is a GRM made with all SNPs. These are the *background* SNPs that probably don't affect the trait.

A further weighting factor **(w)** is applied to S so that the genetic architecture SNPs can be given more or less importance relative to the background SNPs in G. The S and G matrices are then recombined into the final GRM (called the **T-matrix**):

**T** = w(**S**) + (1-w)(**G**)

#### the effect of w
w ranges from 0 to 1. 
- When w=0 the selected SNPs in the S matrix are given no overall weight and the model is simply GBLUP
- When w=1 the selected SNPs in the S matrix are given full weight and the background SNPs are ignored.

## Using this R package
To run BLUP|GA you need:
- a genotype matrix (in -1,0,1} format. Both columns and rows should be named with your sample IDs.
- a phenotype data frame with the first column containing your sample IDs..
- a vector of weightings for SNPs.
